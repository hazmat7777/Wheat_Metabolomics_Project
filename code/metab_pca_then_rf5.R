# this script runs machine learning models on metabolomic data using PCA
# POSITIVE CLASS = 1 (GT_present)
# UPDATED: Now includes hyperparameter tuning for Random Forest

# Load required libraries
library(gbm)
library(caret)
library(glmnet)
library(randomForest)
library(dplyr)
library(doParallel)

# Set up parallel processing for caret
registerDoParallel(cores = detectCores() - 1)

# load data
boost_data <- readRDS("../data/fk_metabolomics_boost_data.RDS")

###########################################################################
# PCA-BASED PREPROCESSING FUNCTION (UNCHANGED)
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
  cat("Number of zeros:", sum(features == 0), "\n")
  
  # STEP 1: Robust zero imputation using TRAINING statistics
  cat("Performing robust zero imputation on training data...\n")
  
  features_imputed <- features
  zero_replacement_values <- numeric(length(names(features)))
  names(zero_replacement_values) <- names(features)
  
  # Impute zeros using training data statistics
  for(col in names(features)) {
    zeros_all <- features[[col]] == 0
    if(sum(zeros_all) > 0) {
      # Calculate replacement value from TRAINING data only
      train_nonzero <- features[train_idx, col][features[train_idx, col] > 0]
      if(length(train_nonzero) > 0) {
        # Use median of smallest 10% instead of minimum for robustness
        min_val <- quantile(train_nonzero, 0.1, na.rm = TRUE)
        replacement_val <- min_val / 2
        features_imputed[[col]][zeros_all] <- replacement_val
        zero_replacement_values[col] <- replacement_val
      }
    }
  }
  
  # STEP 2: Log transform and robust scaling using TRAINING statistics
  cat("Performing log transform and robust scaling...\n")
  
  features_log <- log10(features_imputed + 1)
  
  # Calculate robust scaling parameters from training data only
  train_log <- features_log[train_idx, ]
  train_medians <- sapply(train_log, median, na.rm = TRUE)
  train_mads <- sapply(train_log, mad, na.rm = TRUE)
  
  # Apply scaling to ALL data using training parameters
  features_scaled <- features_log
  for(i in 1:ncol(features_log)) {
    if(train_mads[i] > 0) {
      features_scaled[, i] <- (features_log[, i] - train_medians[i]) / train_mads[i]
    }
  }
  
  # STEP 3: Basic variance filtering (optional - remove truly invariant features)
  cat("Performing basic variance filtering...\n")
  
  train_features_scaled <- features_scaled[train_idx, ]
  train_var <- sapply(train_features_scaled, var, na.rm = TRUE)
  
  # Keep features with some variance
  valid_features <- train_var > min_variance_filter & is.finite(train_var)
  
  if(sum(valid_features) == 0) {
    cat("Warning: No features with sufficient variance. Using all features.\n")
    features_for_pca <- features_scaled
    valid_feature_names <- names(features_scaled)
  } else {
    features_for_pca <- features_scaled[, valid_features]
    valid_feature_names <- names(features_scaled)[valid_features]
    cat("Removed", sum(!valid_features), "features with very low variance\n")
  }
  
  cat("Features going into PCA:", ncol(features_for_pca), "\n")
  
  # STEP 4: Apply PCA using TRAINING data only
  cat("Applying PCA to training data...\n")
  
  train_features_for_pca <- features_for_pca[train_idx, ]
  
  # Perform PCA on training data
  pca_result <- tryCatch({
    prcomp(train_features_for_pca, center = FALSE, scale. = FALSE)  # Already scaled
  }, error = function(e) {
    cat("PCA failed with error:", e$message, "\n")
    cat("Trying with centering...\n")
    prcomp(train_features_for_pca, center = TRUE, scale. = FALSE)
  })
  
  # Calculate explained variance
  explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_var <- cumsum(explained_var)
  
  # Determine number of components
  if(variance_threshold < 1) {
    # Select components explaining desired variance threshold
    n_components <- min(which(cumulative_var >= variance_threshold))
    n_components <- min(n_components, max_components, ncol(train_features_for_pca))
  } else {
    # Use fixed number of components
    n_components <- min(max_components, ncol(train_features_for_pca))
  }

  n_components <- max(2, n_components)
  
  cat("Selected", n_components, "principal components\n")
  cat("Explained variance:", round(cumulative_var[n_components], 3), "\n")
  cat("Individual PC variances:", round(explained_var[1:n_components], 3), "\n")
  
  # STEP 5: Transform ALL data using the PCA fitted on training data
  cat("Transforming all data using training-fitted PCA...\n")
  
  # Transform all data (both training and test)
  features_pca <- predict(pca_result, features_for_pca)[, 1:n_components, drop = FALSE]
  
  # Create column names for PCs
  colnames(features_pca) <- paste0("PC", 1:n_components)
  
  # Create final dataset
  processed_data <- as.data.frame(features_pca)
  processed_data$gt <- target
  
  # Store preprocessing parameters
  preprocessing_params <- list(
    pca_model = pca_result,
    n_components = n_components,
    explained_variance = explained_var[1:n_components],
    cumulative_variance = cumulative_var[n_components],
    train_medians = train_medians,
    train_mads = train_mads,
    zero_replacement_values = zero_replacement_values,
    valid_features = valid_feature_names
  )
  
  # Feature loadings for interpretation
  loadings <- pca_result$rotation[, 1:n_components, drop = FALSE]
  
  cat("Final processed data dimensions:", dim(processed_data), "\n")
  cat("Principal components:", colnames(features_pca), "\n")
  
  return(list(
    data = processed_data,
    params = preprocessing_params,
    loadings = loadings,
    selected_features = paste0("PC", 1:n_components)  # For compatibility
  ))
}

###########################################################################
# UPDATED MODEL TRAINING FUNCTION WITH RF TUNING
###########################################################################

train_multiple_models <- function(train_data, target_col = "gt") {
  
  cat("=== TRAINING MODELS ===\n")
  
  # Prepare data
  x_train <- as.matrix(train_data[, !names(train_data) %in% target_col])
  y_train <- train_data[[target_col]]
  
  # Model 1: GBM with CV for optimal trees
  cat("Training GBM...\n")
  gbm_model <- gbm(
    gt ~ ., 
    data = train_data,
    distribution = "bernoulli",
    n.trees = 200,
    interaction.depth = 2,
    shrinkage = 0.1,
    n.minobsinnode = 3,
    bag.fraction = 0.8,
    cv.folds = 5,
    verbose = FALSE
  )
  optimal_trees <- gbm.perf(gbm_model, method = "cv", plot.it = FALSE)
  
  # Model 2: Elastic Net with CV
  cat("Training Elastic Net...\n")
  cv_glmnet <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5, nfolds = 5)
  
  # Model 3: Random Forest with hyperparameter tuning
  cat("Training Random Forest with hyperparameter tuning...\n")
  
  # Determine mtry range based on number of features
  n_features <- ncol(train_data) - 1
  max_mtry <- min(n_features, 10)  # Cap at 10 for efficiency
  
  # Create tuning grid
  if(n_features >= 3) {
    mtry_values <- unique(c(1, floor(sqrt(n_features)), floor(n_features/3), floor(n_features/2), max_mtry))
    mtry_values <- mtry_values[mtry_values <= n_features & mtry_values >= 1]
  } else {
    mtry_values <- 1:n_features
  }
  
  # Set up cross-validation
  train_control <- trainControl(
    method = "cv",
    number = 3,  # Reduced for speed with small data
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  
  # Prepare data for caret (needs factor with proper levels)
  train_data_caret <- train_data
  train_data_caret[[target_col]] <- factor(train_data_caret[[target_col]], 
                                           levels = c(0, 1), 
                                           labels = c("Class0", "Class1"))
  
  # Create tuning grid
  tune_grid <- expand.grid(
    mtry = mtry_values,
    splitrule = "gini",
    min.node.size = c(1, 3, 5)
  )
  
  cat("  Tuning grid size:", nrow(tune_grid), "combinations\n")
  cat("  mtry values:", paste(mtry_values, collapse = ", "), "\n")
  
  # Train with tuning using ranger (faster than randomForest)
  rf_tuned <- train(
    as.formula(paste(target_col, "~ .")),
    data = train_data_caret,
    method = "ranger",
    trControl = train_control,
    tuneGrid = tune_grid,
    num.trees = 200,
    importance = 'impurity',
    metric = "ROC",
    num.threads = 1  # Avoid nested parallelization
  )
  
  cat("  Best parameters: mtry =", rf_tuned$bestTune$mtry, 
      ", min.node.size =", rf_tuned$bestTune$min.node.size, "\n")
  
  # Train final RF model with best parameters using randomForest for consistency
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
  
  cat("Models trained successfully\n\n")
  
  return(list(
    gbm = gbm_model,
    glmnet = cv_glmnet,
    rf = rf_model,
    rf_tuned = rf_tuned,  # Keep tuning results
    optimal_trees = optimal_trees,
    x_train = x_train,
    y_train = y_train,
    best_rf_params = list(mtry = best_mtry, nodesize = best_nodesize)
  ))
}

###########################################################################
# MODEL EVALUATION FUNCTION - POSITIVE CLASS = 1 (UNCHANGED)
###########################################################################

evaluate_on_test_set <- function(models, test_data, target_col = "gt") {
  
  cat("=== EVALUATING ON TEST SET ===\n")
  
  # Prepare test data
  x_test <- as.matrix(test_data[, !names(test_data) %in% target_col])
  y_test <- test_data[[target_col]]
  
  # Get predictions
  gbm_pred <- predict(models$gbm, test_data, n.trees = models$optimal_trees, type = "response")
  glmnet_pred <- predict(models$glmnet, x_test, s = "lambda.min", type = "response")[,1]
  rf_pred <- predict(models$rf, test_data, type = "prob")[,2]  # Class "1" probabilities
  
  # Ensemble prediction
  ensemble_pred <- (gbm_pred + glmnet_pred + rf_pred) / 3
  
  # Test multiple thresholds
  thresholds <- seq(0.3, 0.7, 0.1)
  threshold_results <- list()
  
  for(thresh in thresholds) {
    pred_class <- ifelse(ensemble_pred > thresh, 1, 0)
    pred_class <- factor(pred_class, levels = c(0, 1))
    actual <- factor(y_test, levels = c(0, 1))
    
    # Set positive class = "1" for confusionMatrix
    cm <- confusionMatrix(pred_class, actual, positive = "1")
    
    threshold_results[[as.character(thresh)]] <- list(
      threshold = thresh,
      accuracy = cm$overall['Accuracy'],
      kappa = cm$overall['Kappa'],
      sensitivity = cm$byClass['Sensitivity'],  # Now for class 1
      specificity = cm$byClass['Specificity'],  # Now for class 1
      precision = cm$byClass['Pos Pred Value'], # Precision for class 1
      confusion_matrix = cm
    )
  }
  
  # Find best threshold
  kappa_scores <- sapply(threshold_results, function(x) x$kappa)
  best_thresh <- names(threshold_results)[which.max(kappa_scores)]
  
  cat("Best threshold:", best_thresh, "\n")
  print(threshold_results[[best_thresh]]$confusion_matrix)
  
  # Individual model performance
  cat("\n=== INDIVIDUAL MODEL COMPARISON ===\n")
  thresh_val <- as.numeric(best_thresh)
  
  individual_results <- list()
  model_preds <- list(GBM = gbm_pred, ElasticNet = glmnet_pred, RandomForest = rf_pred, Ensemble = ensemble_pred)
  
  for(model_name in names(model_preds)) {
    pred_class <- ifelse(model_preds[[model_name]] > thresh_val, 1, 0)
    pred_class <- factor(pred_class, levels = c(0, 1))
    actual <- factor(y_test, levels = c(0, 1))
    
    # Set positive class = "1"
    cm <- confusionMatrix(pred_class, actual, positive = "1")
    individual_results[[model_name]] <- list(
      accuracy = cm$overall['Accuracy'],
      kappa = cm$overall['Kappa'],
      sensitivity = cm$byClass['Sensitivity'],  # Sensitivity for class 1
      specificity = cm$byClass['Specificity'],  # Specificity for class 1
      precision = cm$byClass['Pos Pred Value']  # Precision for class 1
    )
    
    cat(sprintf("%s - Accuracy: %.3f, Kappa: %.3f, Sensitivity: %.3f, Precision: %.3f\n", 
                model_name, cm$overall['Accuracy'], cm$overall['Kappa'],
                cm$byClass['Sensitivity'], cm$byClass['Pos Pred Value']))
  }
  
  # PC importance (instead of individual metabolites)
  cat("\n=== PRINCIPAL COMPONENT IMPORTANCE ===\n")
  gbm_imp <- summary(models$gbm, n.trees = models$optimal_trees, plotit = FALSE)
  print(head(gbm_imp, 10))
  
  # RF tuning results
  cat("\n=== RF TUNING RESULTS ===\n")
  cat("Best RF parameters: mtry =", models$best_rf_params$mtry, 
      ", nodesize =", models$best_rf_params$nodesize, "\n")
  
  return(list(
    threshold_results = threshold_results,
    best_threshold = thresh_val,
    individual_results = individual_results,
    predictions = model_preds,
    feature_importance = gbm_imp,
    rf_params = models$best_rf_params
  ))
}

###########################################################################
# MAIN EVALUATION FUNCTION WITH PCA (UPDATED)
###########################################################################

evaluate_metabolomics_pca <- function(data, target_col = "gt", test_prop = 0.3, seed = 123,
                                      variance_threshold = 0.8, max_components = 10,
                                      min_variance_filter = 0.01) {
  
  set.seed(seed)
  
  # Single train/test split
  train_idx <- createDataPartition(data[[target_col]], p = 1 - test_prop, list = FALSE)
  test_idx <- setdiff(1:nrow(data), train_idx)
  
  cat("Train samples:", length(train_idx), "Test samples:", length(test_idx), "\n")
  
  # Preprocess data with PCA
  preprocessed <- preprocess_metabolomics_pca(
    data, target_col, train_idx,
    variance_threshold = variance_threshold,
    max_components = max_components,
    min_variance_filter = min_variance_filter
  )
  
  processed_data <- preprocessed$data
  
  # Split processed data
  train_data <- processed_data[train_idx, ]
  test_data <- processed_data[test_idx, ]
  
  cat("Final training data dimensions:", dim(train_data), "\n")
  cat("Final test data dimensions:", dim(test_data), "\n\n")
  
  # Train models
  models <- train_multiple_models(train_data, target_col)
  
  # Evaluate on test set
  results <- evaluate_on_test_set(models, test_data, target_col)
  
  return(list(
    models = models,
    results = results,
    train_data = train_data,
    test_data = test_data,
    train_idx = train_idx,
    test_idx = test_idx,
    preprocessing_params = preprocessed$params,
    loadings = preprocessed$loadings,
    selected_features = preprocessed$selected_features
  ))
}

###########################################################################
# USAGE AND CROSS-VALIDATION
###########################################################################

# Convert gt to numeric (1 = GT_present, 0 = GT_absent) - POSITIVE CLASS = 1
boost_data <- boost_data %>%
  mutate(gt = ifelse(gt == "GT_present", 1,
              ifelse(gt == "GT_absent", 0, gt))) %>%
  mutate(gt = as.numeric(gt))

# Verify class distribution
cat("=== CLASS DISTRIBUTION ===\n")
table(boost_data$gt)
cat("0 = GT_absent (negative class)\n")
cat("1 = GT_present (positive class)\n\n")

# Run single example
cat("=== RUNNING SINGLE PCA-BASED ANALYSIS (EXAMPLE) ===\n")
pca_results <- evaluate_metabolomics_pca(
  boost_data, 
  target_col = "gt", 
  test_prop = 0.3,
  variance_threshold = 0.8,
  max_components = 10
)

# ====================================
# MANUAL 10-FOLD CV FOR PCA MODEL WITH TUNING
# ====================================
cat("\n=== MANUAL 10-FOLD CV FOR TUNED PCA MODEL ===\n")

# Function to extract balanced accuracy from PCA result
extract_pca_metrics <- function(pca_result) {
  # The confusion matrix is in the threshold results
  cm <- pca_result$results$threshold_results[["0.5"]]$confusion_matrix
  
  metrics <- list(
    accuracy = cm$overall['Accuracy'],
    balanced_accuracy = cm$byClass['Balanced Accuracy'],
    sensitivity = cm$byClass['Sensitivity'],
    specificity = cm$byClass['Specificity'],
    kappa = cm$overall['Kappa']
  )
  
  return(metrics)
}

# Run PCA 10 times with different seeds for CV
n_runs <- 10
pca_cv_results <- list()
important_pcs <- list()
rf_params_list <- list()

for(i in 1:n_runs) {
  cat("PCA Run", i, "...")
  
  # Use different seed for each run
  pca_result <- evaluate_metabolomics_pca(
    boost_data, 
    target_col = "gt", 
    test_prop = 0.3,
    seed = 123 + i * 10,  # Different seed each time
    variance_threshold = 0.8,
    max_components = 10
  )
  
  # Extract metrics
  metrics <- extract_pca_metrics(pca_result)
  pca_cv_results[[i]] <- metrics
  
  # Extract which PC was most important in GBM
  gbm_importance <- pca_result$results$feature_importance
  most_important_pc <- gbm_importance$var[1]
  important_pcs[[i]] <- most_important_pc
  
  # Extract RF tuning parameters
  rf_params_list[[i]] <- pca_result$results$rf_params
  
  cat(" Balanced Acc:", round(metrics$balanced_accuracy, 3), 
      " Most Important PC:", most_important_pc,
      " RF mtry:", rf_params_list[[i]]$mtry, "\n")
}

# Calculate summary statistics
accuracies <- sapply(pca_cv_results, function(x) x$accuracy)
balanced_accuracies <- sapply(pca_cv_results, function(x) x$balanced_accuracy)
sensitivities <- sapply(pca_cv_results, function(x) x$sensitivity)
specificities <- sapply(pca_cv_results, function(x) x$specificity)
kappas <- sapply(pca_cv_results, function(x) x$kappa)

# Extract RF parameters
mtry_values <- sapply(rf_params_list, function(x) x$mtry)
nodesize_values <- sapply(rf_params_list, function(x) x$nodesize)

# Print individual results
cat("\n=== INDIVIDUAL RUN RESULTS ===\n")
for(i in 1:n_runs) {
  cat(sprintf("Run %2d: Bal.Acc = %.3f, Accuracy = %.3f, Kappa = %.3f, Most Important PC = %s, RF mtry = %d, nodesize = %d\n",
              i, balanced_accuracies[i], accuracies[i], kappas[i], important_pcs[[i]], 
              mtry_values[i], nodesize_values[i]))
}

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat(sprintf("Mean Balanced Accuracy: %.3f ± %.3f\n", 
            mean(balanced_accuracies, na.rm = TRUE), 
            sd(balanced_accuracies, na.rm = TRUE)))
cat(sprintf("Mean Accuracy: %.3f ± %.3f\n", 
            mean(accuracies, na.rm = TRUE), 
            sd(accuracies, na.rm = TRUE)))
cat(sprintf("Mean Kappa: %.3f ± %.3f\n", 
            mean(kappas, na.rm = TRUE), 
            sd(kappas, na.rm = TRUE)))

# RF parameter variability
cat("\n=== RF PARAMETER VARIABILITY ===\n")
cat("mtry values used:", paste(sort(unique(mtry_values)), collapse = ", "), "\n")
cat("mtry distribution:\n")
print(table(mtry_values))
cat("nodesize values used:", paste(sort(unique(nodesize_values)), collapse = ", "), "\n")
cat("nodesize distribution:\n")
print(table(nodesize_values))

# PC importance variability
cat("\n=== PC IMPORTANCE VARIABILITY ===\n")
pc_importance_table <- table(unlist(important_pcs))
cat("Most important PC across runs:\n")
print(pc_importance_table)

# Compare to No Information Rate
nir <- max(table(boost_data$gt)) / length(boost_data$gt)
cat(sprintf("\nNo Information Rate: %.3f\n", nir))
cat(sprintf("Mean improvement over NIR: %.3f (%.1f percentage points)\n", 
            mean(balanced_accuracies, na.rm = TRUE) - nir,
            (mean(balanced_accuracies, na.rm = TRUE) - nir) * 100))

# Statistical test: is performance significantly better than NIR?
t_test_result <- t.test(balanced_accuracies, mu = nir, alternative = "greater")
cat(sprintf("T-test vs NIR: p-value = %.4f\n", t_test_result$p.value))
if(t_test_result$p.value < 0.05) {
  cat("Performance is significantly better than NIR (p < 0.05)\n")
} else {
  cat("Performance is NOT significantly better than NIR (p >= 0.05)\n")
}

# Create summary dataframe
pca_manual_cv_summary <- data.frame(
  Run = 1:n_runs,
  Accuracy = round(accuracies, 3),
  Balanced_Accuracy = round(balanced_accuracies, 3),
  Sensitivity = round(sensitivities, 3),
  Specificity = round(specificities, 3),
  Kappa = round(kappas, 3),
  Most_Important_PC = unlist(important_pcs),
  RF_mtry = mtry_values,
  RF_nodesize = nodesize_values,
  stringsAsFactors = FALSE
)

# Save results
write.csv(pca_manual_cv_summary, "../results/metab_rf/pca_tuned_manual_cv_10runs.csv", row.names = FALSE)
cat("\nResults saved to 'pca_tuned_manual_cv_10runs.csv'\n")

# Show the dataframe
print(pca_manual_cv_summary)

# Stop parallel processing
stopImplicitCluster()








# # Check the t-test manually
# balanced_accuracies <- c(0.452, 0.583, 0.432, 0.708, 0.542, 0.614, 0.571, 0.632, 0.583, 0.597)
# nir <- 0.521

# # Manual t-test
# mean_diff <- mean(balanced_accuracies) - nir
# se_diff <- sd(balanced_accuracies) / sqrt(length(balanced_accuracies))
# t_stat <- mean_diff / se_diff

# cat("Mean balanced accuracy:", mean(balanced_accuracies), "\n")
# cat("NIR:", nir, "\n") 
# cat("Mean difference:", mean_diff, "\n")
# cat("Standard error:", se_diff, "\n")
# cat("T-statistic:", t_stat, "\n")
# cat("Degrees of freedom:", length(balanced_accuracies) - 1, "\n")

# # Critical value for one-tailed test at α = 0.05, df = 9
# critical_value <- qt(0.95, df = 9)
# cat("Critical value (one-tailed, α=0.05):", critical_value, "\n")

# # P-value
# p_value <- 1 - pt(t_stat, df = 9)
# cat("P-value:", p_value, "\n")