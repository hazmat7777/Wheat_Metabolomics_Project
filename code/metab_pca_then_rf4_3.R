# PCA + Random Forest Cross-Validation Script
# Runs PCA preprocessing followed by Random Forest 10 times
# POSITIVE CLASS = 1 (GT_present)

# Load required libraries
library(caret)
library(randomForest)
library(dplyr)

# Load data
boost_data <- readRDS("../data/fk_metabolomics_boost_data.RDS")

###########################################################################
# PCA PREPROCESSING FUNCTION (Same as before)
###########################################################################

preprocess_metabolomics_pca <- function(data, target_col, train_idx, 
                                          variance_threshold = 0.8,
                                          max_components = 10,
                                          min_variance_filter = 0.01) {
  
  # Separate features and target
  if(is.character(target_col)) {
    target <- data[[target_col]]
    features <- data[, !names(data) %in% target_col]
  } else {
    target <- data[, target_col]
    features <- data[, -target_col]
  }
  
  cat("=== PCA PREPROCESSING (TRAINING DATA ONLY) ===\n")
  cat("Original data dimensions:", dim(features), "\n")
  
  # STEP 1: Robust zero imputation using TRAINING statistics
  features_imputed <- features
  zero_replacement_values <- numeric(length(names(features)))
  names(zero_replacement_values) <- names(features)
  
  for(col in names(features)) {
    zeros_all <- features[[col]] == 0
    if(sum(zeros_all) > 0) {
      train_nonzero <- features[train_idx, col][features[train_idx, col] > 0]
      if(length(train_nonzero) > 0) {
        min_val <- quantile(train_nonzero, 0.1, na.rm = TRUE)
        replacement_val <- min_val / 2
        features_imputed[[col]][zeros_all] <- replacement_val
        zero_replacement_values[col] <- replacement_val
      }
    }
  }
  
  # STEP 2: Log transform and robust scaling using TRAINING statistics
  features_log <- log10(features_imputed + 1)
  
  train_log <- features_log[train_idx, ]
  train_medians <- sapply(train_log, median, na.rm = TRUE)
  train_mads <- sapply(train_log, mad, na.rm = TRUE)
  
  features_scaled <- features_log
  for(i in 1:ncol(features_log)) {
    if(train_mads[i] > 0) {
      features_scaled[, i] <- (features_log[, i] - train_medians[i]) / train_mads[i]
    }
  }
  
  # STEP 3: Basic variance filtering
  train_features_scaled <- features_scaled[train_idx, ]
  train_var <- sapply(train_features_scaled, var, na.rm = TRUE)
  
  valid_features <- train_var > min_variance_filter & is.finite(train_var)
  
  if(sum(valid_features) == 0) {
    features_for_pca <- features_scaled
    valid_feature_names <- names(features_scaled)
  } else {
    features_for_pca <- features_scaled[, valid_features]
    valid_feature_names <- names(features_scaled)[valid_features]
  }
  
  # STEP 4: Apply PCA using TRAINING data only
  train_features_for_pca <- features_for_pca[train_idx, ]
  
  pca_result <- tryCatch({
    prcomp(train_features_for_pca, center = FALSE, scale. = FALSE)
  }, error = function(e) {
    prcomp(train_features_for_pca, center = TRUE, scale. = FALSE)
  })
  
  # Calculate explained variance
  explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_var <- cumsum(explained_var)
  
  # Determine number of components
  if(variance_threshold < 1) {
    n_components <- min(which(cumulative_var >= variance_threshold))
    n_components <- min(n_components, max_components, ncol(train_features_for_pca))
  } else {
    n_components <- min(max_components, ncol(train_features_for_pca))
  }
  
  n_components <- max(2, n_components)
  
  cat("Selected", n_components, "principal components\n")
  cat("Explained variance:", round(cumulative_var[n_components], 3), "\n")
  
  # STEP 5: Transform ALL data using the PCA fitted on training data
  features_pca <- predict(pca_result, features_for_pca)[, 1:n_components, drop = FALSE]
  colnames(features_pca) <- paste0("PC", 1:n_components)
  
  # Create final dataset
  processed_data <- as.data.frame(features_pca)
  processed_data$gt <- target
  
  return(list(
    data = processed_data,
    n_components = n_components,
    explained_variance = cumulative_var[n_components]
  ))
}

###########################################################################
# SINGLE PCA + RF RUN FUNCTION
###########################################################################

run_pca_rf_single <- function(data, target_col = "gt", test_prop = 0.3, seed = 123,
                              variance_threshold = 0.8, max_components = 10) {
  
  set.seed(seed)
  
  # Create train/test split
  train_idx <- createDataPartition(data[[target_col]], p = 1 - test_prop, list = FALSE)
  test_idx <- setdiff(1:nrow(data), train_idx)
  
  # Apply PCA preprocessing
  preprocessed <- preprocess_metabolomics_pca(
    data, target_col, train_idx,
    variance_threshold = variance_threshold,
    max_components = max_components
  )
  
  processed_data <- preprocessed$data
  
  # Split processed data
  train_data <- processed_data[train_idx, ]
  test_data <- processed_data[test_idx, ]
  
  # Train Random Forest with hyperparameter tuning
  cat("Training Random Forest with hyperparameter tuning...\n")
  
  # Set up cross-validation for hyperparameter tuning
  train_control <- trainControl(
    method = "cv",
    number = 3,  # 3-fold CV for speed
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    verboseIter = FALSE
  )
  
  # Prepare data for caret (needs factor with proper levels)
  train_data_caret <- train_data
  train_data_caret[[target_col]] <- factor(train_data_caret[[target_col]], 
                                           levels = c(0, 1), 
                                           labels = c("Class0", "Class1"))
  
  # Create comprehensive tuning grid for RF
  n_features <- ncol(train_data) - 1
  
  # More comprehensive mtry search
  if(n_features <= 5) {
    # For small feature sets (like PCA), test all possibilities
    mtry_values <- 1:n_features
  } else if(n_features <= 10) {
    # For medium feature sets, test most values
    mtry_values <- c(1, 2, 3, floor(sqrt(n_features)), 
                    floor(n_features/3), floor(n_features/2), 
                    floor(2*n_features/3), n_features)
  } else {
    # For larger feature sets, use strategic sampling
    mtry_values <- c(1, 2, 3, floor(sqrt(n_features)), 
                    floor(n_features/4), floor(n_features/3), 
                    floor(n_features/2), floor(2*n_features/3))
  }
  
  # Remove duplicates and ensure valid range
  mtry_values <- unique(mtry_values[mtry_values <= n_features & mtry_values >= 1])
  
  # More comprehensive tuning grid
  tune_grid <- expand.grid(
    mtry = mtry_values,
    splitrule = c("gini", "extratrees"),  # Try both splitting rules
    min.node.size = c(1, 2, 3, 5, 10)    # More nodesize options
  )
  
  cat("  Comprehensive tuning grid size:", nrow(tune_grid), "combinations\n")
  cat("  mtry values being tested:", paste(mtry_values, collapse = ", "), "\n")
  
  # Train RF with hyperparameter tuning
  rf_tuned <- train(
    as.formula(paste(target_col, "~ .")),
    data = train_data_caret,
    method = "ranger",
    trControl = train_control,
    tuneGrid = tune_grid,
    num.trees = 200,
    importance = 'impurity',
    metric = "ROC",
    num.threads = 1
  )
  
  # Train final RF model with best parameters
  best_mtry <- rf_tuned$bestTune$mtry
  best_nodesize <- rf_tuned$bestTune$min.node.size
  
  rf_model <- randomForest(
    as.factor(gt) ~ ., 
    data = train_data,
    ntree = 200,
    mtry = best_mtry,
    nodesize = best_nodesize,
    importance = TRUE
  )
  
  # Make predictions on test set
  rf_pred_prob <- predict(rf_model, test_data, type = "prob")[,2]  # Probability of class 1
  rf_pred_class <- predict(rf_model, test_data, type = "class")
  
  # Convert predictions to factor with correct levels
  rf_pred_class <- factor(as.numeric(as.character(rf_pred_class)), levels = c(0, 1))
  y_test <- factor(test_data[[target_col]], levels = c(0, 1))
  
  # Calculate confusion matrix and metrics
  cm <- confusionMatrix(rf_pred_class, y_test, positive = "1")
  
  # Extract metrics
  metrics <- list(
    accuracy = as.numeric(cm$overall['Accuracy']),
    balanced_accuracy = as.numeric(cm$byClass['Balanced Accuracy']),
    sensitivity = as.numeric(cm$byClass['Sensitivity']),
    specificity = as.numeric(cm$byClass['Specificity']),
    kappa = as.numeric(cm$overall['Kappa']),
    precision = as.numeric(cm$byClass['Pos Pred Value']),
    n_components = preprocessed$n_components,
    explained_variance = preprocessed$explained_variance,
    best_mtry = best_mtry,
    best_nodesize = best_nodesize
  )
  
  return(metrics)
}

###########################################################################
# MAIN 10-FOLD CV FUNCTION
###########################################################################

run_pca_rf_cv <- function(data, target_col = "gt", n_runs = 10, test_prop = 0.3,
                          variance_threshold = 0.8, max_components = 10) {
  
  cat("=== RUNNING", n_runs, "RUNS OF PCA + RANDOM FOREST ===\n")
  
  # Convert gt to numeric if needed
  data <- data %>%
    mutate(gt = ifelse(gt == "GT_present", 1,
                ifelse(gt == "GT_absent", 0, gt))) %>%
    mutate(gt = as.numeric(gt))
  
  # Verify class distribution
  cat("=== CLASS DISTRIBUTION ===\n")
  print(table(data$gt))
  cat("0 = GT_absent (negative class)\n")
  cat("1 = GT_present (positive class)\n\n")
  
  # Initialize storage
  all_metrics <- list()
  
  # Run multiple iterations
  for(i in 1:n_runs) {
    cat("Run", i, "...")
    
    # Use different seed for each run
    metrics <- run_pca_rf_single(
      data, 
      target_col = target_col, 
      test_prop = test_prop,
      seed = 123 + i * 10,
      variance_threshold = variance_threshold,
      max_components = max_components
    )
    
    all_metrics[[i]] <- metrics
    
    cat(" Balanced Acc:", round(metrics$balanced_accuracy, 3),
        " PCs:", metrics$n_components,
        " mtry:", metrics$best_mtry, "\n")
  }
  
  # Extract results into vectors
  accuracies <- sapply(all_metrics, function(x) x$accuracy)
  balanced_accuracies <- sapply(all_metrics, function(x) x$balanced_accuracy)
  sensitivities <- sapply(all_metrics, function(x) x$sensitivity)
  specificities <- sapply(all_metrics, function(x) x$specificity)
  kappas <- sapply(all_metrics, function(x) x$kappa)
  precisions <- sapply(all_metrics, function(x) x$precision)
  n_components_vec <- sapply(all_metrics, function(x) x$n_components)
  explained_variances <- sapply(all_metrics, function(x) x$explained_variance)
  best_mtrys <- sapply(all_metrics, function(x) x$best_mtry)
  best_nodesizes <- sapply(all_metrics, function(x) x$best_nodesize)
  
  # Print individual results
  cat("\n=== INDIVIDUAL RUN RESULTS ===\n")
  for(i in 1:n_runs) {
    cat(sprintf("Run %2d: Acc = %.3f, Bal.Acc = %.3f, Kappa = %.3f, PCs = %d, mtry = %d\n",
                i, accuracies[i], balanced_accuracies[i], kappas[i], 
                n_components_vec[i], best_mtrys[i]))
  }
  
  # Calculate summary statistics
  cat("\n=== SUMMARY STATISTICS ===\n")
  cat(sprintf("Accuracy: %.3f ± %.3f (range: %.3f - %.3f)\n", 
              mean(accuracies), sd(accuracies), min(accuracies), max(accuracies)))
  cat(sprintf("Balanced Accuracy: %.3f ± %.3f (range: %.3f - %.3f)\n", 
              mean(balanced_accuracies), sd(balanced_accuracies), 
              min(balanced_accuracies), max(balanced_accuracies)))
  cat(sprintf("Sensitivity: %.3f ± %.3f\n", mean(sensitivities), sd(sensitivities)))
  cat(sprintf("Specificity: %.3f ± %.3f\n", mean(specificities), sd(specificities)))
  cat(sprintf("Kappa: %.3f ± %.3f\n", mean(kappas), sd(kappas)))
  cat(sprintf("Precision: %.3f ± %.3f\n", mean(precisions, na.rm = TRUE), sd(precisions, na.rm = TRUE)))
  
  # PCA component analysis
  cat("\n=== PCA COMPONENT ANALYSIS ===\n")
  cat(sprintf("Number of PCs: %.1f ± %.1f (range: %d - %d)\n", 
              mean(n_components_vec), sd(n_components_vec), 
              min(n_components_vec), max(n_components_vec)))
  cat(sprintf("Explained Variance: %.3f ± %.3f\n", 
              mean(explained_variances), sd(explained_variances)))
  
  cat("PC distribution:\n")
  print(table(n_components_vec))
  
  # RF hyperparameter analysis
  cat("\n=== RF HYPERPARAMETER ANALYSIS ===\n")
  cat("mtry distribution:\n")
  print(table(best_mtrys))
  cat("nodesize distribution:\n")
  print(table(best_nodesizes))
  
  # Compare to No Information Rate
  nir <- max(table(data[[target_col]])) / nrow(data)
  cat(sprintf("\nNo Information Rate: %.3f\n", nir))
  cat(sprintf("Mean improvement over NIR: %.3f (%.1f percentage points)\n", 
              mean(balanced_accuracies) - nir,
              (mean(balanced_accuracies) - nir) * 100))
  
  # Statistical significance test
  if(length(balanced_accuracies) > 1) {
    t_test <- t.test(balanced_accuracies, mu = nir, alternative = "greater")
    cat(sprintf("T-test vs NIR: p-value = %.4f\n", t_test$p.value))
    if(t_test$p.value < 0.05) {
      cat("Performance is significantly better than NIR (p < 0.05)\n")
    } else {
      cat("Performance is NOT significantly better than NIR (p >= 0.05)\n")
    }
  }
  
  # Create results dataframe
  results_df <- data.frame(
    Run = 1:n_runs,
    Accuracy = round(accuracies, 3),
    Balanced_Accuracy = round(balanced_accuracies, 3),
    Sensitivity = round(sensitivities, 3),
    Specificity = round(specificities, 3),
    Kappa = round(kappas, 3),
    Precision = round(precisions, 3),
    N_Components = n_components_vec,
    Explained_Variance = round(explained_variances, 3),
    Best_mtry = best_mtrys,
    Best_nodesize = best_nodesizes,
    stringsAsFactors = FALSE
  )
  
  # Summary statistics dataframe
  summary_df <- data.frame(
    Metric = c("Accuracy", "Balanced_Accuracy", "Sensitivity", "Specificity", "Kappa", "Precision"),
    Mean = round(c(mean(accuracies), mean(balanced_accuracies), mean(sensitivities), 
                   mean(specificities), mean(kappas), mean(precisions, na.rm = TRUE)), 3),
    SD = round(c(sd(accuracies), sd(balanced_accuracies), sd(sensitivities), 
                 sd(specificities), sd(kappas), sd(precisions, na.rm = TRUE)), 3),
    Min = round(c(min(accuracies), min(balanced_accuracies), min(sensitivities), 
                  min(specificities), min(kappas), min(precisions, na.rm = TRUE)), 3),
    Max = round(c(max(accuracies), max(balanced_accuracies), max(sensitivities), 
                  max(specificities), max(kappas), max(precisions, na.rm = TRUE)), 3),
    stringsAsFactors = FALSE
  )
  
  return(list(
    detailed_results = results_df,
    summary_stats = summary_df,
    nir = nir,
    t_test_p_value = if(length(balanced_accuracies) > 1) t_test$p.value else NA,
    all_metrics = all_metrics
  ))
}

###########################################################################
# USAGE
###########################################################################

# Run PCA + RF cross-validation
cv_results <- run_pca_rf_cv(
  boost_data, 
  target_col = "gt", 
  n_runs = 10, 
  test_prop = 0.3,
  variance_threshold = 0.8,
  max_components = 10
)

# Save results
write.csv(cv_results$detailed_results, "../results/metab_rf/pca_rf_cv_detailed.csv", row.names = FALSE)
write.csv(cv_results$summary_stats, "../results/metab_rf/pca_rf_cv_summary.csv", row.names = FALSE)

cat("\n=== FINAL RESULTS SUMMARY ===\n")
print(cv_results$summary_stats)

cat("\nResults saved to CSV files\n")