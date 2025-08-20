# HPC-Optimized Random Forest Script for Phyloseq GT Status Prediction
# FIXED VERSION - Addresses protection stack overflow and core detection
# HPC-Optimized Random Forest Script for Phyloseq GT Status Prediction

# ============================================
# CRITICAL FIXES - MUST BE FIRST!
# ============================================
# CRITICAL FIX 1: Prevent protection stack overflow
options(expressions = 500000)

# Set up personal library path
lib_path <- "~/R_libs"
.libPaths(c(lib_path, .libPaths()))

# ============================================
# NOW LOAD PACKAGES
# ============================================
library(tidyverse)
library(phyloseq)
library(ranger)
library(parallel)
library(doParallel)  # or doSNOW
library(foreach)
library(caret)
library(data.table)

# ============================================
# CRITICAL FIX 2: Force core detection (AFTER packages)
# ============================================
# Check multiple possible environment variables
ncpus_env <- Sys.getenv("NCPUS", "")
pbs_np <- Sys.getenv("PBS_NP", "")
omp_threads <- Sys.getenv("OMP_NUM_THREADS", "")

if(nzchar(ncpus_env)) {
    total_cores <- as.numeric(ncpus_env)
    cat("Using NCPUS:", total_cores, "\n")
} else if(nzchar(pbs_np)) {
    total_cores <- as.numeric(pbs_np)
    cat("Using PBS_NP:", total_cores, "\n")
} else if(nzchar(omp_threads)) {
    total_cores <- as.numeric(omp_threads)
    cat("Using OMP_NUM_THREADS:", total_cores, "\n")
} else {
    total_cores <- 16  # Hard default
    cat("No environment variables found, defaulting to:", total_cores, "\n")
}

# Sanity check
if(total_cores == 1) {
    cat("WARNING: Only 1 core detected, forcing to 16\n")
    total_cores <- 16
}

cat("Total cores available:", total_cores, "\n")

# Set up parallelization
registerDoParallel(cores = total_cores)

cores_per_fold <- max(1, floor(total_cores / 5))  # Divide cores among CV folds
cat("Cores per fold:", cores_per_fold, "\n")

# Continue with rest of script...
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

# CRITICAL FIX 3: More aggressive filtering to prevent memory issues
smartFilter <- function(data, target_col, max_features = 20000) {  # Reduced from 20000
    cat("Applying smart filtering for high-dimensional data...\n")
    
    # Get ESV columns
    sample_data_cols <- which(names(data) == target_col)
    metadata_cols <- 1:max(sample_data_cols)
    esv_cols <- setdiff(1:ncol(data), metadata_cols)
    
    cat("Original ESVs:", length(esv_cols), "\n")
    
    # Strategy 1: Remove completely absent ESVs
    esv_sums <- colSums(data[, esv_cols, drop = FALSE])
    present_esvs <- esv_cols[esv_sums > 0]
    
    cat("ESVs with any presence:", length(present_esvs), "\n")
    
    # Strategy 2: Aggressive filtering to prevent protection stack overflow
    if(length(present_esvs) > max_features) {
        cat(sprintf("Reducing to %d features to prevent memory issues...\n", max_features))
        
        # Get abundance and prevalence
        esv_data <- data[, present_esvs, drop = FALSE]
        esv_abundance <- colSums(esv_data)
        esv_prevalence <- colSums(esv_data > 0) / nrow(esv_data)
        
        # Create a combined score (abundance * prevalence)
        esv_scores <- esv_abundance * esv_prevalence
        
        # Keep top features by combined score
        top_features <- names(sort(esv_scores, decreasing = TRUE))[1:max_features]
        keep_esvs <- present_esvs[names(data)[present_esvs] %in% top_features]
        
        cat("ESVs after filtering:", length(keep_esvs), "\n")
    } else {
        keep_esvs <- present_esvs
    }
    
    # Keep all metadata + filtered ESVs
    filtered_data <- data[, c(metadata_cols, esv_cols)] # change esv_cols to keep_esvs if you want to filter.
    
    return(filtered_data)
}

# CRITICAL FIX 4: Simplified parallel CV without complex options
parallelCV <- function(data, target_col, n_folds = 5, n_trees = 2000,  # Reduced trees
                      mtry = NULL, min_node_size = 5, sample_fraction = 0.632) {
    
    cat("Setting up parallel cross-validation...\n")
    
    # Auto-calculate mtry if not provided
    if(is.null(mtry)) {
        n_features <- ncol(data) - 1
        mtry <- floor(sqrt(n_features))
        mtry <- min(mtry, 100)  # Cap at 100 for stability
    }
    
    cat("Using mtry =", mtry, "\n")
    cat("Using", n_trees, "trees per model\n")
    
    # Create stratified folds
    folds <- createFolds(data[[target_col]], k = n_folds, list = TRUE, returnTrain = FALSE)
    
    # Simplified parallel CV loop
    cv_results <- foreach(i = 1:n_folds, 
                         .packages = c("ranger"),
                         .combine = rbind,
                         .errorhandling = "remove") %dopar% {  # Changed to "remove" bad results
        
        tryCatch({
            # Create train/test splits
            test_indices <- folds[[i]]
            train_data <- data[-test_indices, ]
            test_data <- data[test_indices, ]
            
            # Fit Random Forest with memory-efficient settings
            rf_model <- ranger(
                formula = as.formula(paste(target_col, "~ .")),
                data = train_data,
                importance = 'none',  # Skip importance in CV to save memory
                mtry = mtry,
                min.node.size = min_node_size,
                sample.fraction = sample_fraction,
                num.trees = n_trees,
                num.threads = cores_per_fold,
                seed = 123 + i,
                verbose = FALSE,
                respect.unordered.factors = 'ignore'  # Faster processing
            )
            
            # Predictions
            predictions <- predict(rf_model, test_data, num.threads = cores_per_fold)
            pred_classes <- predictions$predictions
            
            # Ensure factor levels match
            actual <- test_data[[target_col]]
            pred_classes <- factor(pred_classes, levels = levels(actual))
            
            # Calculate metrics
            accuracy <- mean(pred_classes == actual, na.rm = TRUE)
            
            # Return results
            data.frame(
                fold = i,
                accuracy = accuracy,
                n_train = nrow(train_data),
                n_test = nrow(test_data),
                oob_error = rf_model$prediction.error,
                stringsAsFactors = FALSE
            )
        }, error = function(e) {
            cat("Error in fold", i, ":", e$message, "\n")
            NULL
        })
    }
    
    return(cv_results)
}

# Simplified feature importance function
getFeatureImportance <- function(data, target_col, n_trees = 500, mtry = NULL) {
    cat("Calculating feature importance on full dataset...\n")
    
    if(is.null(mtry)) {
        n_features <- ncol(data) - 1
        mtry <- floor(sqrt(n_features))
        mtry <- min(mtry, 100)
    }
    
    # Fit RF on full dataset for importance
    rf_full <- ranger(
        formula = as.formula(paste(target_col, "~ .")),
        data = data,
        importance = 'permutation',
        mtry = mtry,
        num.trees = n_trees,
        num.threads = total_cores,
        verbose = FALSE,
        respect.unordered.factors = 'ignore'
    )
    
    # Get importance
    importance_vec <- rf_full$variable.importance
    importance_df <- data.frame(
        feature = names(importance_vec),
        importance = as.numeric(importance_vec),
        stringsAsFactors = FALSE
    ) %>%
        arrange(desc(importance))
    
    return(importance_df)
}

# Main Analysis ----------------------------------
cat("=== HPC RANDOM FOREST ANALYSIS FOR GT PREDICTION ===\n")
cat("Fixed version with protection stack and core detection fixes\n")
cat("Date:", format(Sys.time()), "\n")

# Load data
ps_highdiv_relative <- readRDS("ps_highdiv_relative.rds")
cat("Phyloseq object loaded\n")
print(ps_highdiv_relative)

# Check memory usage
cat("Initial memory usage:", format(object.size(ps_highdiv_relative), units = "MB"), "\n")

# Prepare data 
rf_data <- preparePhyloseqForRF_Simple(ps_highdiv_relative, "gt")

cat("Data dimensions after preparation:", dim(rf_data), "\n")
cat("GT status distribution:\n")
print(table(rf_data$gt))

# Apply smart filtering - REDUCED FEATURES
rf_data_filtered <- smartFilter(rf_data, "gt", max_features = 5000)  # Reduced from 20000

cat("Final data dimensions:", dim(rf_data_filtered), "\n")
cat("Memory usage of filtered data:", format(object.size(rf_data_filtered), units = "MB"), "\n")

# Remove original data to free memory
rm(rf_data)
gc()

# HPC-Optimized Parameters ----------------------------------
n_trees <- 2000  # Reduced from 1000 for faster testing
n_cv_folds <- 5

cat("\nRF Parameters:\n")
cat("- Trees:", n_trees, "\n")
cat("- CV folds:", n_cv_folds, "\n")
cat("- Features:", ncol(rf_data_filtered) - 1, "\n")
cat("- Cores:", total_cores, "\n")

# Run Parallel Cross-Validation ----------------------------------
start_time <- Sys.time()
cat("\nStarting parallel cross-validation at", format(start_time), "\n")

cv_results <- parallelCV(
    data = rf_data_filtered,
    target_col = "gt",
    n_folds = n_cv_folds,
    n_trees = n_trees,
    mtry = 50,  # Will be auto-calculated
    min_node_size = 5,
    sample_fraction = 0.632
)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
cat("Cross-validation completed in", round(runtime, 2), "minutes\n")

# Results Analysis ----------------------------------
cat("\n=== CROSS-VALIDATION RESULTS ===\n")

if(!is.null(cv_results) && is.data.frame(cv_results) && nrow(cv_results) > 0) {
    cat("Successfully completed", nrow(cv_results), "out of", n_cv_folds, "folds\n")
    
    mean_accuracy <- mean(cv_results$accuracy, na.rm = TRUE)
    sd_accuracy <- sd(cv_results$accuracy, na.rm = TRUE)
    
    cat("Mean CV Accuracy:", round(mean_accuracy * 100, 2), "% ±", round(sd_accuracy * 100, 2), "%\n")
    cat("Individual fold accuracies:", round(cv_results$accuracy * 100, 2), "%\n")
    
    # Feature Importance Analysis
    cat("\nCalculating feature importance...\n")
    importance_summary <- getFeatureImportance(rf_data_filtered, "gt", n_trees = 500)
    
    cat("\nTop 20 most important features:\n")
    print(head(importance_summary, 20))
    
    # Save results
    write_csv(importance_summary, "hpc_gt_prediction_importance.csv")
    write_csv(cv_results, "hpc_gt_prediction_cv_results.csv")
    
    # Taxonomic analysis of important features
    if(exists("ps_highdiv_relative")) {
        cat("\n=== TAXONOMIC ANALYSIS ===\n")
        
        # Get taxonomy table
        tax_df <- data.frame(tax_table(ps_highdiv_relative))
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
    cat("ERROR: Cross-validation failed\n")
    if(!is.null(cv_results)) {
        print(cv_results)
    }
}

# Cleanup ----------------------------------
gc()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Total runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes\n")

# Final summary
cat("\nFINAL SUMMARY:\n")
cat("Dataset: 316k taxa reduced to", ncol(rf_data_filtered) - 1, "features\n")
if(exists("mean_accuracy")) {
    cat("Cross-validation accuracy:", round(mean_accuracy * 100, 2), "%\n")
    if(exists("importance_summary") && nrow(importance_summary) > 0) {
        cat("Top predictor:", importance_summary$feature[1], "\n")
    }
}
cat("Script completed successfully!\n")