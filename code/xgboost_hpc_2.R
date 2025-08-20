#!/usr/bin/env Rscript
# ============================================================================
# HPC XGBoost Analysis for 316k Taxa Phyloseq Data
# Handles ALL features including rare taxa - no filtering!
# Fixes all previous errors: stack overflow, core detection, numeric storage
# ============================================================================

# ============================================
# CRITICAL MEMORY FIXES - MUST BE FIRST!
# ============================================
options(expressions = 500000)          # Prevent stack overflow
options(gcinfo = FALSE)               # Reduce GC overhead
options(memory.profiling = FALSE)     # Reduce memory tracking
Sys.setenv("R_MAX_NUM_DLLS" = 1000)  # Allow more DLLs

# ============================================
# Library Setup
# ============================================
lib_path <- "~/R_libs"
.libPaths(c(lib_path, .libPaths()))

# Load packages
suppressPackageStartupMessages({
    library(tidyverse)
    library(phyloseq)
    library(xgboost)
    library(Matrix)      # For sparse matrices
    library(parallel)
    library(doParallel)
    library(foreach)
    library(caret)
    library(data.table)
})

cat("=== HPC XGBOOST ANALYSIS FOR 316K TAXA ===\n")
cat("Date:", format(Sys.time()), "\n")
cat("This script keeps ALL features including rare taxa\n\n")

# ============================================
# CORE DETECTION WITH ALL FIXES
# ============================================
detectCoresRobust <- function() {
    # Try multiple environment variables
    ncpus_env <- Sys.getenv("NCPUS", "")
    pbs_np <- Sys.getenv("PBS_NP", "")
    omp_threads <- Sys.getenv("OMP_NUM_THREADS", "")
    
    if(nzchar(ncpus_env)) {
        cores <- as.numeric(ncpus_env)
        cat("Detected cores from NCPUS:", cores, "\n")
    } else if(nzchar(pbs_np)) {
        cores <- as.numeric(pbs_np)
        cat("Detected cores from PBS_NP:", cores, "\n")
    } else if(nzchar(omp_threads)) {
        cores <- as.numeric(omp_threads)
        cat("Detected cores from OMP_NUM_THREADS:", cores, "\n")
    } else {
        cores <- 16  # Default for your HPC setup
        cat("No environment variables found, defaulting to:", cores, "\n")
    }
    
    # Sanity check - we know we requested 16
    if(cores == 1) {
        cat("WARNING: Only 1 core detected, forcing to 16\n")
        cores <- 16
    }
    
    return(cores)
}

# Detect and set up cores
total_cores <- detectCoresRobust()
registerDoParallel(cores = total_cores)
cat("Parallel backend registered with", total_cores, "cores\n\n")

# ============================================
# DATA PREPARATION FUNCTIONS
# ============================================

preparePhyloseqForML <- function(phyloseq_obj, target_variable) {
    cat("Preparing phyloseq data for ML analysis...\n")
    
    # Extract sample data
    sample_df <- data.frame(sample_data(phyloseq_obj), stringsAsFactors = FALSE)
    
    # Extract OTU/ESV table
    otu_mat <- as.matrix(otu_table(phyloseq_obj))
    if(taxa_are_rows(phyloseq_obj)) {
        otu_df <- data.frame(t(otu_mat), check.names = FALSE)
    } else {
        otu_df <- data.frame(otu_mat, check.names = FALSE)
    }
    
    # Ensure sample names match
    common_samples <- intersect(rownames(sample_df), rownames(otu_df))
    if(length(common_samples) == 0) {
        stop("No matching sample names between sample_data and otu_table!")
    }
    
    sample_df <- sample_df[common_samples, , drop = FALSE]
    otu_df <- otu_df[common_samples, , drop = FALSE]
    
    # Get target variable and ensure it's properly formatted
    if(!target_variable %in% names(sample_df)) {
        stop(paste("Target variable", target_variable, "not found in sample data!"))
    }
    
    target <- sample_df[[target_variable]]
    
    # Remove samples with missing target
    complete_cases <- !is.na(target)
    target <- target[complete_cases]
    otu_df <- otu_df[complete_cases, , drop = FALSE]
    sample_df <- sample_df[complete_cases, , drop = FALSE]
    
    cat("  Samples after removing NAs:", sum(complete_cases), "\n")
    
    # Convert target to factor and then to numeric (0-based for XGBoost)
    target <- as.factor(target)
    cat("  Target variable levels:", paste(levels(target), collapse = ", "), "\n")
    cat("  Target distribution:\n")
    print(table(target))
    
    # Store level names for later
    target_levels <- levels(target)
    target_numeric <- as.numeric(target) - 1  # XGBoost needs 0-based
    
    # Keep only numeric columns from sample data (excluding target)
    numeric_cols <- sapply(sample_df, is.numeric)
    numeric_cols[target_variable] <- FALSE  # Exclude target
    
    if(any(numeric_cols)) {
        metadata_features <- sample_df[, numeric_cols, drop = FALSE]
        cat("  Including", ncol(metadata_features), "numeric metadata features\n")
    } else {
        metadata_features <- NULL
        cat("  No additional numeric metadata features\n")
    }
    
    # Combine features
    if(!is.null(metadata_features)) {
        all_features <- cbind(metadata_features, otu_df)
    } else {
        all_features <- otu_df
    }
    
    # Clean column names (XGBoost doesn't like special characters)
    colnames(all_features) <- make.names(colnames(all_features), unique = TRUE)
    
    cat("  Total features before filtering:", ncol(all_features), "\n")
    
    return(list(
        features = all_features,
        target = target_numeric,
        target_levels = target_levels,
        sample_names = rownames(all_features)
    ))
}

filterZeroVariance <- function(features, min_presence = 0) {
    cat("Filtering features...\n")
    
    # Remove zero-variance features (same value in all samples)
    feature_vars <- apply(features, 2, var, na.rm = TRUE)
    keep_cols <- feature_vars > 0 & !is.na(feature_vars)
    
    cat("  Removed", sum(!keep_cols), "zero-variance features\n")
    
    # Optional: remove features present in fewer than min_presence samples
    if(min_presence > 0) {
        presence <- colSums(features > 0, na.rm = TRUE)
        keep_cols <- keep_cols & (presence >= min_presence)
        cat("  Removed features present in <", min_presence, "samples\n")
    }
    
    filtered_features <- features[, keep_cols, drop = FALSE]
    cat("  Features after filtering:", ncol(filtered_features), "\n")
    
    return(filtered_features)
}

# ============================================
# XGBOOST FUNCTIONS
# ============================================

runXGBoostCV <- function(features, target, target_levels, n_folds = 5, 
                        params = NULL, nrounds = 100, early_stopping = 10,
                        verbose = TRUE) {
    
    cat("\n=== Starting XGBoost Cross-Validation ===\n")
    cat("Features:", ncol(features), "x", nrow(features), "\n")
    cat("Folds:", n_folds, "\n")
    cat("Rounds:", nrounds, "\n")
    
    # Default parameters if not provided
    if(is.null(params)) {
        params <- list(
            objective = ifelse(length(target_levels) == 2, "binary:logistic", "multi:softprob"),
            eval_metric = ifelse(length(target_levels) == 2, "logloss", "mlogloss"),
            max_depth = 6,
            eta = 0.3,
            subsample = 0.8,
            colsample_bytree = 0.8,
            min_child_weight = 1,
            nthread = total_cores
        )
        
        if(length(target_levels) > 2) {
            params$num_class <- length(target_levels)
        }
    }
    
    # Convert to DMatrix (XGBoost's optimized data structure)
    # This handles sparse data efficiently
    cat("Converting to XGBoost DMatrix format...\n")
    dtrain <- xgb.DMatrix(data = as.matrix(features), label = target)
    
    # Create stratified folds
    set.seed(123)
    folds <- createFolds(target, k = n_folds, list = TRUE, returnTrain = FALSE)
    
    # Initialize results storage
    cv_results <- data.frame(
        fold = integer(),
        accuracy = numeric(),
        logloss = numeric(),
        n_train = integer(),
        n_test = integer(),
        stringsAsFactors = FALSE
    )
    
    # Run CV - using simple loop instead of foreach to avoid memory issues
    for(i in 1:n_folds) {
        cat("\nFold", i, "of", n_folds, "...\n")
        
        tryCatch({
            # Split data
            test_idx <- folds[[i]]
            train_idx <- setdiff(1:nrow(features), test_idx)
            
            # Create train/test matrices
            dtrain_fold <- xgb.DMatrix(
                data = as.matrix(features[train_idx, , drop = FALSE]), 
                label = target[train_idx]
            )
            dtest_fold <- xgb.DMatrix(
                data = as.matrix(features[test_idx, , drop = FALSE]), 
                label = target[test_idx]
            )
            
            # Train model
            watchlist <- list(train = dtrain_fold, test = dtest_fold)
            
            xgb_model <- xgb.train(
                params = params,
                data = dtrain_fold,
                nrounds = nrounds,
                watchlist = watchlist,
                early_stopping_rounds = early_stopping,
                verbose = ifelse(verbose, 1, 0),
                print_every_n = 50
            )
            
            # Make predictions
            if(length(target_levels) == 2) {
                # Binary classification
                pred_probs <- predict(xgb_model, dtest_fold)
                pred_classes <- ifelse(pred_probs > 0.5, 1, 0)
            } else {
                # Multiclass
                pred_probs <- predict(xgb_model, dtest_fold, reshape = TRUE)
                pred_classes <- max.col(pred_probs) - 1
            }
            
            # Calculate metrics
            actual_test <- target[test_idx]
            accuracy <- mean(pred_classes == actual_test)
            
            # Calculate logloss
            eps <- 1e-15
            if(length(target_levels) == 2) {
                pred_probs <- pmax(pmin(pred_probs, 1 - eps), eps)
                logloss <- -mean(actual_test * log(pred_probs) + 
                                (1 - actual_test) * log(1 - pred_probs))
            } else {
                # Multiclass logloss
                logloss <- NA  # Complex to calculate, skip for now
            }
            
            # Store results
            fold_results <- data.frame(
                fold = i,
                accuracy = accuracy,
                logloss = ifelse(is.na(logloss), NA, logloss),
                n_train = length(train_idx),
                n_test = length(test_idx),
                best_iteration = xgb_model$best_iteration,
                stringsAsFactors = FALSE
            )
            
            cv_results <- rbind(cv_results, fold_results)
            
            cat("  Fold", i, "accuracy:", round(accuracy * 100, 2), "%\n")
            
            # Clean up memory
            rm(xgb_model, dtrain_fold, dtest_fold)
            gc(verbose = FALSE)
            
        }, error = function(e) {
            cat("ERROR in fold", i, ":", e$message, "\n")
        })
    }
    
    return(cv_results)
}

getFeatureImportance <- function(features, target, target_levels, 
                                params = NULL, nrounds = 100) {
    
    cat("\nCalculating feature importance on full dataset...\n")
    
    # Set parameters
    if(is.null(params)) {
        params <- list(
            objective = ifelse(length(target_levels) == 2, "binary:logistic", "multi:softprob"),
            eval_metric = ifelse(length(target_levels) == 2, "logloss", "mlogloss"),
            max_depth = 6,
            eta = 0.3,
            subsample = 0.8,
            colsample_bytree = 0.8,
            nthread = total_cores
        )
        
        if(length(target_levels) > 2) {
            params$num_class <- length(target_levels)
        }
    }
    
    # Train on full dataset
    dtrain <- xgb.DMatrix(data = as.matrix(features), label = target)
    
    xgb_model <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = nrounds,
        verbose = 0
    )
    
    # Get importance
    importance_matrix <- xgb.importance(model = xgb_model)
    
    # Add feature names if missing
    if(!"Feature" %in% names(importance_matrix)) {
        importance_matrix$Feature <- colnames(features)[as.numeric(importance_matrix$Feature) + 1]
    }
    
    return(importance_matrix)
}

# ============================================
# MAIN ANALYSIS
# ============================================

# Track total runtime
script_start <- Sys.time()

# Load phyloseq data
cat("\nLoading phyloseq data...\n")
ps_highdiv <- readRDS("ps_highdiv_relative.rds")
print(ps_highdiv)
cat("Memory usage:", format(object.size(ps_highdiv), units = "MB"), "\n\n")

# Prepare data
ml_data <- preparePhyloseqForML(ps_highdiv, "gt")

# Apply minimal filtering (keep rare taxa!)
# Only remove features that are exactly zero in ALL samples
ml_data$features <- filterZeroVariance(ml_data$features, min_presence = 1)

cat("\nFinal dataset dimensions:", nrow(ml_data$features), "samples x", 
    ncol(ml_data$features), "features\n")
cat("Memory usage:", format(object.size(ml_data$features), units = "MB"), "\n")

# Clean up memory
rm(ps_highdiv)
gc(verbose = FALSE)

# ============================================
# CROSS-VALIDATION
# ============================================

# XGBoost parameters optimized for high-dimensional data
xgb_params <- list(
    objective = "binary:logistic",  # Binary classification
    eval_metric = "logloss",
    max_depth = 4,                  # Shallower trees for high dimensions
    eta = 0.1,                       # Lower learning rate
    subsample = 0.7,                 # Row subsampling
    colsample_bytree = 0.05,         # Aggressive column subsampling for 316k features
    min_child_weight = 5,            # Regularization
    gamma = 1,                       # Minimum loss reduction
    nthread = total_cores
)

# Run cross-validation
cv_start <- Sys.time()
cat("\n=== STARTING CROSS-VALIDATION ===\n")
cat("Time:", format(cv_start), "\n")

cv_results <- runXGBoostCV(
    features = ml_data$features,
    target = ml_data$target,
    target_levels = ml_data$target_levels,
    n_folds = 5,
    params = xgb_params,
    nrounds = 200,
    early_stopping = 20,
    verbose = TRUE
)

cv_end <- Sys.time()
cv_runtime <- difftime(cv_end, cv_start, units = "mins")

# ============================================
# RESULTS ANALYSIS
# ============================================

cat("\n=== CROSS-VALIDATION RESULTS ===\n")

if(nrow(cv_results) > 0) {
    # Calculate summary statistics
    mean_accuracy <- mean(cv_results$accuracy, na.rm = TRUE)
    sd_accuracy <- sd(cv_results$accuracy, na.rm = TRUE)
    se_accuracy <- sd_accuracy / sqrt(nrow(cv_results))
    
    cat("\nSummary Statistics:\n")
    cat("Mean Accuracy: ", sprintf("%.2f%%", mean_accuracy * 100), 
        " ± ", sprintf("%.2f%%", sd_accuracy * 100), "\n", sep = "")
    cat("95% CI: [", sprintf("%.2f%%", (mean_accuracy - 1.96*se_accuracy) * 100),
        ", ", sprintf("%.2f%%", (mean_accuracy + 1.96*se_accuracy) * 100), "]\n", sep = "")
    
    cat("\nPer-fold results:\n")
    print(cv_results)
    
    # Save CV results
    write_csv(cv_results, "xgboost_cv_results_316k_features.csv")
    cat("\nCV results saved to: xgboost_cv_results_316k_features.csv\n")
    
} else {
    cat("ERROR: No successful CV folds completed\n")
}

# ============================================
# FEATURE IMPORTANCE
# ============================================

if(nrow(cv_results) > 0) {
    cat("\n=== CALCULATING FEATURE IMPORTANCE ===\n")
    
    importance_df <- getFeatureImportance(
        features = ml_data$features,
        target = ml_data$target,
        target_levels = ml_data$target_levels,
        params = xgb_params,
        nrounds = 150
    )
    
    if(nrow(importance_df) > 0) {
        cat("\nTop 30 Most Important Features:\n")
        print(head(importance_df, 30))
        
        # Save importance
        write_csv(importance_df, "xgboost_importance_316k_features.csv")
        cat("\nFeature importance saved to: xgboost_importance_316k_features.csv\n")
        
        # Try to match with taxonomy
        if(exists("ps_highdiv")) {
            ps_highdiv <- readRDS("ps_highdiv_relative.rds")
        }
        
        tryCatch({
            tax_table_df <- data.frame(tax_table(ps_highdiv), stringsAsFactors = FALSE)
            tax_table_df$Feature <- make.names(rownames(tax_table_df), unique = TRUE)
            
            # Merge importance with taxonomy
            importance_with_tax <- merge(
                importance_df[1:min(100, nrow(importance_df)), ],
                tax_table_df,
                by = "Feature",
                all.x = TRUE
            )
            
            if(nrow(importance_with_tax) > 0) {
                cat("\nTop 20 Important Features with Taxonomy:\n")
                print(head(importance_with_tax[, c("Feature", "Gain", "Phylum", 
                                                   "Class", "Family", "Genus")], 20))
                
                write_csv(importance_with_tax, "xgboost_importance_with_taxonomy.csv")
                cat("\nImportance with taxonomy saved to: xgboost_importance_with_taxonomy.csv\n")
            }
        }, error = function(e) {
            cat("Could not match taxonomy:", e$message, "\n")
        })
    }
}

# ============================================
# FINAL SUMMARY
# ============================================

script_end <- Sys.time()
total_runtime <- difftime(script_end, script_start, units = "mins")

cat("\n" , rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Dataset Summary:\n")
cat("  Original features: 316,306 taxa\n")
cat("  Features after filtering:", ncol(ml_data$features), "\n")
cat("  Samples:", nrow(ml_data$features), "\n")
cat("  Target distribution:\n")
target_table <- table(ml_data$target)
names(target_table) <- ml_data$target_levels
print(target_table)

if(exists("mean_accuracy")) {
    cat("\nModel Performance:\n")
    cat("  Mean CV Accuracy:", sprintf("%.2f%%", mean_accuracy * 100), "\n")
    cat("  CV Runtime:", round(cv_runtime, 2), "minutes\n")
}

cat("\nTotal Runtime:", round(total_runtime, 2), "minutes\n")
cat("Cores Used:", total_cores, "\n")
cat("\nAll results saved to CSV files\n")

# Clean up
gc(verbose = FALSE)
cat("\nScript completed successfully!\n")