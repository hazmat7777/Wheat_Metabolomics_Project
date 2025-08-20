# HPC-Optimized Random Forest Script for Phyloseq GT Status Prediction
# Fixed version based on local testing - removes probability=TRUE and sparse matrix issues

# Set up personal library path
lib_path <- "~/R_libs"
.libPaths(c(lib_path, .libPaths()))

# Packages ----------------------------------
library(tidyverse)
library(phyloseq)
library(ranger)
library(parallel)
library(doSNOW)       # For PBS clusters  
library(foreach)
library(caret)
library(data.table)

# HPC Setup ----------------------------------
# Detect HPC environment and set up appropriate parallel backend
slurm_env <- Sys.getenv("SLURM_JOB_ID")
pbs_env <- Sys.getenv("PBS_JOBID")
mpi_env <- Sys.getenv("OMPI_COMM_WORLD_SIZE")

if(nzchar(pbs_env)) {
    # PBS environment
    cat("Detected PBS environment - Job ID:", pbs_env, "\n")
    
    # First try NCPUS (common PBS variable)
    ncpus <- Sys.getenv("NCPUS", "")
    if(nzchar(ncpus)) {
        total_cores <- as.numeric(ncpus)
        cat("Using NCPUS:", total_cores, "\n")
    } else {
        # Get PBS allocation details
        pbs_nodefile <- Sys.getenv("PBS_NODEFILE", "")
        if(nzchar(pbs_nodefile) && file.exists(pbs_nodefile)) {
            nodes <- readLines(pbs_nodefile)
            total_cores <- length(nodes)
            unique_nodes <- unique(nodes)
            cat("PBS allocated", total_cores, "cores across", length(unique_nodes), "nodes\n")
        } else {
            # Try PBS_NP or PBS_NCPUS
            total_cores <- as.numeric(Sys.getenv("PBS_NP", 
                                     Sys.getenv("PBS_NCPUS", 
                                     Sys.getenv("NCPUS", detectCores()))))
            cat("Using PBS_NP/PBS_NCPUS or detected cores:", total_cores, "\n")
        }
    }
    
    # Create SNOW cluster for PBS
    cl <- makeCluster(min(total_cores, 16), type = "SOCK")  # Cap at requested cores
    registerDoSNOW(cl)
    opts <- list(preschedule = FALSE)  # Remove chunkSize - it's not recognized
}

cat("Total cores available:", total_cores, "\n")

# Set up nested parallelization strategy
n_cv_folds <- 5
cores_per_fold <- max(1, floor(total_cores / n_cv_folds))
cat("Will use", cores_per_fold, "cores per CV fold\n")

# Memory management
options(warn = 1)
gc()

# Functions ----------------------------------
preparePhyloseqForRF_Simple <- function(phyloseq_obj, target_variable) {
    cat("Preparing data for Random Forest...\n")
    
    # Extract sample data
    sample_df <- data.frame(sample_data(phyloseq_obj))
    
    # Extract OTU/ESV table and transpose (samples as rows)
    otu_mat <- as.matrix(otu_table(phyloseq_obj))
    if(taxa_are_rows(phyloseq_obj)) {
        otu_df <- data.frame(t(otu_mat))
    } else {
        otu_df <- data.frame(otu_mat)
    }
    
    # Ensure sample names match
    common_samples <- intersect(rownames(sample_df), rownames(otu_df))
    sample_df <- sample_df[common_samples, ]
    otu_df <- otu_df[common_samples, ]
    
    # Combine sample data with abundances
    combined_df <- cbind(sample_df, otu_df)
    
    # Remove samples with missing target variable
    combined_df <- combined_df[!is.na(combined_df[[target_variable]]), ]
    
    # Convert target to factor
    combined_df[[target_variable]] <- as.factor(combined_df[[target_variable]])
    
    # Clean problematic columns that cause 0% accuracy
    combined_df$sample_name <- NULL
    combined_df$og_sample_names <- NULL
    if("tillage" %in% names(combined_df)) {
        combined_df$tillage <- as.factor(combined_df$tillage)
    }
    
    # Remove any remaining character columns (except target)
    char_cols <- sapply(combined_df, is.character)
    char_cols[target_variable] <- FALSE
    if(any(char_cols)) {
        cat("Removing character columns:", names(combined_df)[char_cols], "\n")
        combined_df <- combined_df[, !char_cols]
    }
    
    return(combined_df)
}

# Aggressive filtering for 316k taxa (keep rare taxa but remove impossible ones)
smartFilter <- function(data, target_col) {
    cat("Applying smart filtering for 316k taxa...\n")
    
    # Get ESV columns
    sample_data_cols <- which(names(data) == target_col)
    metadata_cols <- 1:max(sample_data_cols)
    esv_cols <- setdiff(1:ncol(data), metadata_cols)
    
    cat("Original ESVs:", length(esv_cols), "\n")
    
    # Strategy 1: Remove completely absent ESVs
    esv_sums <- colSums(data[, esv_cols, drop = FALSE])
    present_esvs <- esv_cols[esv_sums > 0]
    
    cat("ESVs with any presence:", length(present_esvs), "\n")
    
    # Strategy 2: Keep top abundant + prevalent taxa to manage memory
    if(length(present_esvs) > 50000) {
        cat("Too many ESVs, applying additional filtering...\n")
        
        # Get abundance and prevalence
        esv_data <- data[, present_esvs, drop = FALSE]
        esv_abundance <- colSums(esv_data)
        esv_prevalence <- colSums(esv_data > 0) / nrow(esv_data)
        
        # Keep ESVs that are:
        # - In top 20k most abundant, OR
        # - Present in ≥10% of samples, OR  
        # - In top 5% by abundance
        top_abundant <- names(sort(esv_abundance, decreasing = TRUE))[1:min(20000, length(esv_abundance))]
        prevalent <- names(esv_prevalence)[esv_prevalence >= 0.10]
        top_5pct <- names(esv_abundance)[esv_abundance >= quantile(esv_abundance, 0.95)]
        
        keep_esvs_names <- unique(c(top_abundant, prevalent, top_5pct))
        keep_esvs <- present_esvs[names(data)[present_esvs] %in% keep_esvs_names]
        
        cat("ESVs after smart filtering:", length(keep_esvs), "\n")
    } else {
        keep_esvs <- present_esvs
    }
    
    # Keep all metadata + filtered ESVs
    filtered_data <- data[, c(metadata_cols, keep_esvs)]
    
    return(filtered_data)
}

# Parallel cross-validation function (FIXED - no probability=TRUE)
parallelCV <- function(data, target_col, n_folds = 5, n_trees = 2000, 
                      mtry = 50, min_node_size = 1, sample_fraction = 0.632) {
    
    cat("Setting up parallel cross-validation...\n")
    
    # Create stratified folds
    folds <- createFolds(data[[target_col]], k = n_folds, list = TRUE, returnTrain = FALSE)
    
    
    cat("Using", cores_per_fold, "threads per Random Forest model\n")
    cat("Using mtry =", mtry, "\n")
    
    # Parallel CV loop with PBS-optimized settings
    cv_results <- foreach(i = 1:n_folds, 
                         .packages = c("ranger", "caret"),
                         .combine = rbind,
                         .errorhandling = "pass",
                         .options.snow = opts,
                         .export = c("cores_per_fold")) %dopar% {
        
        # Create train/test splits
        test_indices <- folds[[i]]
        train_data <- data[-test_indices, ]
        test_data <- data[test_indices, ]
        
        # Fit Random Forest - NO probability=TRUE (this was causing 0% accuracy)
        rf_model <- ranger(
            formula = as.formula(paste(target_col, "~ .")),
            data = train_data,
            importance = 'permutation',
            mtry = mtry,
            min.node.size = min_node_size,
            sample.fraction = sample_fraction,
            num.trees = n_trees,
            num.threads = cores_per_fold,
            seed = 123 + i,
            verbose = FALSE
        )
        
        # Predictions - now returns direct class predictions
        predictions <- predict(rf_model, test_data, num.threads = cores_per_fold)
        pred_classes <- predictions$predictions
        
        # Ensure factor levels match for proper comparison
        actual <- test_data[[target_col]]
        pred_classes <- factor(pred_classes, levels = levels(actual))
        
        # Debug- check the factor levels match now:
        if(!all(levels(pred_classes) == levels(actual))) {
            warning(paste("Factor levels mismatch in fold", i))
        }

        # Calculate metrics
        accuracy <- mean(pred_classes == actual, na.rm = TRUE)
        
        # Return results (no list columns to avoid combine issues)
        data.frame(
            fold = i,
            accuracy = accuracy,
            n_train = nrow(train_data),
            n_test = nrow(test_data),
            oob_error = rf_model$prediction.error
        )
    }
      
    return(cv_results)
}

# Feature importance function (simplified)
getFeatureImportance <- function(data, target_col, n_trees = 500, mtry = 50) {
    cat("Calculating feature importance on full dataset...\n")
    
    
    # Fit RF on full dataset for importance
    rf_full <- ranger(
        formula = as.formula(paste(target_col, "~ .")),
        data = data,
        importance = 'permutation',
        mtry = mtry,
        num.trees = n_trees,
        num.threads = total_cores,
        verbose = FALSE
    )
    
    # Get importance
    importance_vec <- rf_full$variable.importance
    importance_df <- data.frame(
        feature = names(importance_vec),
        importance = as.numeric(importance_vec)
    ) %>%
        arrange(desc(importance))
    
    return(importance_df)
}

# Main Analysis ----------------------------------
cat("=== HPC RANDOM FOREST ANALYSIS FOR GT PREDICTION ===\n")
cat("Fixed version - no probability=TRUE, proper factor handling\n")

# Load data
ps_highdiv_absolute <- readRDS("ps_highdiv_relative.rds")
cat("Phyloseq object loaded\n")
print(ps_highdiv_absolute)

# Check memory usage
cat("Initial memory usage:", format(object.size(ps_highdiv_absolute), units = "MB"), "\n")

# Prepare data 
rf_data <- preparePhyloseqForRF_Simple(ps_highdiv_absolute, "gt")

cat("Data dimensions after preparation:", dim(rf_data), "\n")
cat("GT status distribution:\n")
print(table(rf_data$gt))

# Apply smart filtering for 316k taxa
rf_data_filtered <- smartFilter(rf_data, "gt")

cat("Final data dimensions:", dim(rf_data_filtered), "\n")
cat("Memory usage of filtered data:", format(object.size(rf_data_filtered), units = "MB"), "\n")

# Remove original data to free memory
rm(rf_data)
gc()

# HPC-Optimized Parameters ----------------------------------
n_trees <- 1000  
n_cv_folds <- 5
# # Use more conservative mtry for high-dimensional data
# optimal_mtry <- max(50, floor(sqrt(ncol(rf_data_filtered) - 1)))
# optimal_mtry <- min(optimal_mtry, 500)  # Cap at 500 for memory

cat("RF Parameters:\n")
cat("- Trees:", n_trees, "\n")
#cat("- mtry:", optimal_mtry, "\n")
cat("- CV folds:", n_cv_folds, "\n")
cat("- Features:", ncol(rf_data_filtered) - 1, "\n")

# Run Parallel Cross-Validation ----------------------------------
start_time <- Sys.time()
cat("Starting parallel cross-validation at", format(start_time), "\n")

cv_results <- parallelCV(
    data = rf_data_filtered,
    target_col = "gt",
    n_folds = n_cv_folds,
    n_trees = n_trees,
    #mtry = optimal_mtry,
    min_node_size = 3,
    sample_fraction = 0.632
)

# debug- check for any NA accuracies
if(any(is.na(cv_results$accuracy))) {
    cat("WARNING: Some folds returned NA accuracy\n")
    print(cv_results)
}

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
cat("Cross-validation completed in", round(runtime, 2), "minutes\n")

# Results Analysis ----------------------------------
cat("\n=== CROSS-VALIDATION RESULTS ===\n")

if(is.data.frame(cv_results) && nrow(cv_results) > 0) {
    mean_accuracy <- mean(cv_results$accuracy, na.rm = TRUE)
    sd_accuracy <- sd(cv_results$accuracy, na.rm = TRUE)
    
    cat("Mean CV Accuracy:", round(mean_accuracy, 4), "±", round(sd_accuracy, 4), "\n")
    cat("Individual fold accuracies:", round(cv_results$accuracy, 4), "\n")
    
    # Feature Importance Analysis
    importance_summary <- getFeatureImportance(rf_data_filtered, "gt", n_trees)
    
    cat("\nTop 20 most important features:\n")
    print(head(importance_summary, 20))
    
    # Save results
    write_csv(importance_summary, "hpc_gt_prediction_importance.csv")
    write_csv(cv_results, "hpc_gt_prediction_cv_results.csv")
    
    # Taxonomic analysis of important features
    if(exists("ps_highdiv_absolute")) {
        cat("\n=== TAXONOMIC ANALYSIS ===\n")
        
        # Get taxonomy table
        tax_df <- data.frame(tax_table(ps_highdiv_absolute))
        tax_df$feature <- rownames(tax_df)
        
        # Match important features to taxonomy
        important_taxa <- merge(importance_summary, tax_df, 
                               by = "feature", all.x = TRUE) %>%
            filter(!is.na(Phylum)) %>%
            select(feature, importance, Phylum, Class, Order, Family, Genus) %>%
            head(50)
        
        if(nrow(important_taxa) > 0) {
            cat("Top 20 taxonomically-identified important ESVs:\n")
            print(head(important_taxa, 20))
            
            write_csv(important_taxa, "hpc_gt_prediction_important_taxa.csv")
        }
    }
    
} else {
    cat("ERROR: Cross-validation failed or returned no results\n")
    print(cv_results)
}

# Cleanup ----------------------------------
if(!is.null(cl)) {
    stopCluster(cl)
}
gc()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Total runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes\n")
cat("Files saved:\n")
cat("- hpc_gt_prediction_importance.csv\n")
cat("- hpc_gt_prediction_cv_results.csv\n") 
cat("- hpc_gt_prediction_important_taxa.csv\n")

# Final summary
cat("\nFINAL SUMMARY:\n")
cat("Dataset: 316k taxa reduced to", ncol(rf_data_filtered) - 1, "features\n")
if(exists("mean_accuracy")) {
    cat("Cross-validation accuracy:", round(mean_accuracy, 3), "\n")
    if(exists("importance_summary")) {
        cat("Top predictor:", importance_summary$feature[1], "\n")
    }
}
cat("Memory-optimized HPC analysis completed successfully!\n")