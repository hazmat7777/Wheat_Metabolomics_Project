# ============================================
# CRITICAL MEMORY FIXES - MUST BE FIRST!
# ============================================
options(expressions = 500000)
options(gcinfo = FALSE)
Sys.setenv("R_MAX_NUM_DLLS" = 1000)

# ============================================
# Library Setup
# ============================================
lib_path <- "~/R_libs"
.libPaths(c(lib_path, .libPaths()))

suppressPackageStartupMessages({
    library(tidyverse)
    library(phyloseq)
    library(xgboost)
    library(Matrix)
    library(parallel)
    library(doParallel)
    library(foreach)
    library(data.table)
})

cat("=== XGBOOST HYPERPARAMETER TUNING FOR 316K TAXA ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# ============================================
# CORE DETECTION
# ============================================
total_cores <- as.numeric(Sys.getenv("NCPUS", "16"))
if(total_cores == 1) total_cores <- 16
cat("Using", total_cores, "cores\n\n")
registerDoParallel(cores = total_cores)

# ============================================
# LOAD AND PREPARE DATA (Same as before)
# ============================================
prepareData <- function() {
    cat("Loading and preparing data...\n")
    
    # Load phyloseq
    ps <- readRDS("ps_highdiv_relative.rds")
    
    # Extract data
    sample_df <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
    otu_mat <- as.matrix(otu_table(ps))
    if(taxa_are_rows(ps)) {
        otu_df <- data.frame(t(otu_mat), check.names = FALSE)
    } else {
        otu_df <- data.frame(otu_mat, check.names = FALSE)
    }
    
    # Get target
    target <- sample_df[["gt"]]
    complete_cases <- !is.na(target)
    target <- as.factor(target[complete_cases])
    otu_df <- otu_df[complete_cases, ]
    
    # Convert to numeric
    target_numeric <- as.numeric(target) - 1
    
    # Remove zero-variance features
    feature_vars <- apply(otu_df, 2, var, na.rm = TRUE)
    keep_cols <- feature_vars > 0 & !is.na(feature_vars)
    features <- otu_df[, keep_cols]
    
    cat("  Features:", ncol(features), "\n")
    cat("  Samples:", nrow(features), "\n\n")
    
    return(list(
        features = as.matrix(features),
        target = target_numeric,
        target_levels = levels(target)
    ))
}

# ============================================
# HYPERPARAMETER GRID
# ============================================
createParameterGrid <- function() {
    # Define parameter ranges to test
    param_grid <- expand.grid(
        max_depth = c(3, 4, 6, 8),                    # Tree depth
        eta = c(0.01, 0.05, 0.1, 0.2),                # Learning rate
        subsample = c(0.5, 0.7, 0.9),                 # Row sampling
        colsample_bytree = c(0.01, 0.05, 0.1, 0.2),   # Column sampling (critical for 316k!)
        min_child_weight = c(1, 5, 10),               # Regularization
        gamma = c(0, 1, 5),                           # Min split loss
        stringsAsFactors = FALSE
    )
    
    # Add fixed parameters
    param_grid$nrounds <- 500  # Max rounds (early stopping will handle)
    param_grid$early_stopping_rounds <- 20
    
    cat("Testing", nrow(param_grid), "parameter combinations\n\n")
    return(param_grid)
}

# ============================================
# SINGLE PARAMETER EVALUATION
# ============================================
evaluateParameters <- function(params, features, target, n_folds = 3) {
    # Create XGBoost parameters
    xgb_params <- list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        max_depth = params$max_depth,
        eta = params$eta,
        subsample = params$subsample,
        colsample_bytree = params$colsample_bytree,
        min_child_weight = params$min_child_weight,
        gamma = params$gamma,
        nthread = 1  # Use 1 thread per model since we parallelize across parameters
    )
    
    # Create DMatrix
    dtrain <- xgb.DMatrix(data = features, label = target)
    
    # Run cross-validation
    cv_results <- xgb.cv(
        params = xgb_params,
        data = dtrain,
        nrounds = params$nrounds,
        nfold = n_folds,
        early_stopping_rounds = params$early_stopping_rounds,
        verbose = 0,
        stratified = TRUE,
        seed = 123
    )
    
    # Get best score
    best_iteration <- cv_results$best_iteration
    best_score <- cv_results$evaluation_log$test_logloss_mean[best_iteration]
    
    # Calculate accuracy from logloss (approximate)
    # For binary classification: accuracy ≈ 1 - logloss (rough approximation)
    approx_accuracy <- 1 / (1 + exp(best_score))
    
    return(list(
        logloss = best_score,
        accuracy = approx_accuracy,
        best_iteration = best_iteration,
        params = params
    ))
}

# ============================================
# PARALLEL GRID SEARCH
# ============================================
runGridSearch <- function(param_grid, features, target, sample_size = NULL) {
    cat("Starting parallel grid search...\n")
    
    # Optional: Use a subset for faster tuning
    if(!is.null(sample_size) && sample_size < nrow(param_grid)) {
        param_grid <- param_grid[sample(nrow(param_grid), sample_size), ]
        cat("Using random sample of", sample_size, "parameter combinations\n")
    }
    
    # Run parallel evaluation
    results <- foreach(i = 1:nrow(param_grid), 
                      .packages = c("xgboost"),
                      .combine = rbind,
                      .errorhandling = "remove") %dopar% {
        
        params <- param_grid[i, ]
        
        # Evaluate this parameter set
        result <- tryCatch({
            evaluateParameters(params, features, target, n_folds = 3)
        }, error = function(e) {
            list(logloss = NA, accuracy = NA, best_iteration = NA, params = params)
        })
        
        # Return as data frame row
        data.frame(
            max_depth = params$max_depth,
            eta = params$eta,
            subsample = params$subsample,
            colsample_bytree = params$colsample_bytree,
            min_child_weight = params$min_child_weight,
            gamma = params$gamma,
            logloss = result$logloss,
            accuracy = result$accuracy,
            best_iteration = result$best_iteration,
            stringsAsFactors = FALSE
        )
    }
    
    return(results)
}

# ============================================
# FULL EVALUATION WITH BEST PARAMETERS
# ============================================
evaluateBestModel <- function(best_params, features, target, target_levels) {
    cat("\n=== EVALUATING BEST MODEL ===\n")
    
    # Create parameter list
    xgb_params <- list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        max_depth = best_params$max_depth,
        eta = best_params$eta,
        subsample = best_params$subsample,
        colsample_bytree = best_params$colsample_bytree,
        min_child_weight = best_params$min_child_weight,
        gamma = best_params$gamma,
        nthread = total_cores
    )
    
    # 5-fold CV for final evaluation
    set.seed(123)
    n_folds <- 5
    folds <- cut(1:length(target), breaks = n_folds, labels = FALSE)
    folds <- sample(folds)  # Shuffle
    
    cv_accuracies <- numeric(n_folds)
    
    for(i in 1:n_folds) {
        train_idx <- which(folds != i)
        test_idx <- which(folds == i)
        
        dtrain <- xgb.DMatrix(features[train_idx, ], label = target[train_idx])
        dtest <- xgb.DMatrix(features[test_idx, ], label = target[test_idx])
        
        model <- xgb.train(
            params = xgb_params,
            data = dtrain,
            nrounds = best_params$best_iteration,
            verbose = 0
        )
        
        pred <- predict(model, dtest)
        pred_class <- ifelse(pred > 0.5, 1, 0)
        cv_accuracies[i] <- mean(pred_class == target[test_idx])
        
        cat("  Fold", i, "accuracy:", round(cv_accuracies[i] * 100, 2), "%\n")
    }
    
    cat("\nMean CV Accuracy:", round(mean(cv_accuracies) * 100, 2), "%\n")
    cat("SD:", round(sd(cv_accuracies) * 100, 2), "%\n")
    
    return(cv_accuracies)
}

# ============================================
# MAIN EXECUTION
# ============================================
main <- function() {
    start_time <- Sys.time()
    
    # Load data
    data <- prepareData()
    
    # Create parameter grid
    param_grid <- createParameterGrid()
    
    # Option 1: Test ALL combinations (comprehensive but slow)
    # results <- runGridSearch(param_grid, data$features, data$target)
    
    # Option 2: Random search - test subset (faster)
    results <- runGridSearch(param_grid, data$features, data$target, sample_size = 500)
    
    # Remove NA results
    results <- results[!is.na(results$logloss), ]
    
    # Sort by performance
    results <- results[order(results$logloss), ]
    
    # Display top 10 parameter sets
    cat("\n=== TOP 10 PARAMETER COMBINATIONS ===\n")
    print(head(results, 10))
    
    # Save all results
    write_csv(results, "xgboost_tuning_results.csv")
    cat("\nAll results saved to: xgboost_tuning_results.csv\n")
    
    # Get best parameters
    best_params <- results[1, ]
    cat("\n=== BEST PARAMETERS ===\n")
    print(best_params)
    
    # Full evaluation with best parameters
    cv_accuracies <- evaluateBestModel(best_params, data$features, data$target, data$target_levels)
    
    # Save best parameters
    best_model_info <- list(
        parameters = best_params,
        cv_accuracies = cv_accuracies,
        mean_accuracy = mean(cv_accuracies),
        sd_accuracy = sd(cv_accuracies)
    )
    
    saveRDS(best_model_info, "best_xgboost_params.rds")
    
    # Final summary
    end_time <- Sys.time()
    runtime <- difftime(end_time, start_time, units = "mins")
    
    cat("\n============================================================\n")
    cat("TUNING COMPLETE\n")
    cat("============================================================\n")
    cat("Total parameter combinations tested:", nrow(results), "\n")
    cat("Best accuracy achieved:", round(mean(cv_accuracies) * 100, 2), "%\n")
    cat("Total runtime:", round(runtime, 2), "minutes\n")
    cat("\nBest parameters saved to: best_xgboost_params.rds\n")
    cat("All results saved to: xgboost_tuning_results.csv\n")
}

# Run the analysis
main()