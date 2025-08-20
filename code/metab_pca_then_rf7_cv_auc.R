# this script runs machine learning models on preprocessed metabolomic data
# POSITIVE CLASS = 1 (GT_present)
# UPDATED: Includes AUC calculations with standard error

# Load required libraries
library(gbm)
library(caret)
library(glmnet)
library(randomForest)
library(dplyr)
library(doParallel)
library(pROC)  # Added for AUC calculations

# Set up parallel processing for caret
registerDoParallel(cores = detectCores() - 1)

# Load preprocessed data (user will provide this)
boost_data <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")

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
# UPDATED MODEL EVALUATION FUNCTION - INCLUDES AUC CALCULATIONS
###########################################################################

evaluate_on_test_set <- function(models, test_data, target_col = "gt") {
  
  cat("=== EVALUATING ON TEST SET ===\n")
  
  # Prepare test data
  x_test <- as.matrix(test_data[, !names(test_data) %in% target_col])
  y_test <- test_data[[target_col]]
  
  # Get predictions (probabilities)
  gbm_pred <- predict(models$gbm, test_data, n.trees = models$optimal_trees, type = "response")
  glmnet_pred <- predict(models$glmnet, x_test, s = "lambda.min", type = "response")[,1]
  rf_pred <- predict(models$rf, test_data, type = "prob")[,2]  # Class "1" probabilities
  
  # Ensemble prediction
  ensemble_pred <- (gbm_pred + glmnet_pred + rf_pred) / 3
  
  # Calculate AUC for each model
  gbm_auc <- auc(y_test, gbm_pred, quiet = TRUE)
  glmnet_auc <- auc(y_test, glmnet_pred, quiet = TRUE)
  rf_auc <- auc(y_test, rf_pred, quiet = TRUE)
  ensemble_auc <- auc(y_test, ensemble_pred, quiet = TRUE)
  
  # Test multiple thresholds for ensemble (for other metrics)
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
      sensitivity = cm$byClass['Sensitivity'],
      specificity = cm$byClass['Specificity'],
      precision = cm$byClass['Pos Pred Value'],
      balanced_accuracy = cm$byClass['Balanced Accuracy'],
      confusion_matrix = cm
    )
  }
  
  # Find best threshold
  kappa_scores <- sapply(threshold_results, function(x) x$kappa)
  best_thresh <- names(threshold_results)[which.max(kappa_scores)]
  thresh_val <- as.numeric(best_thresh)
  
  cat("Best threshold:", best_thresh, "\n")
  
  # Individual model performance with best threshold
  cat("\n=== INDIVIDUAL MODEL COMPARISON (with AUC) ===\n")
  
  individual_results <- list()
  model_preds <- list(GBM = gbm_pred, ElasticNet = glmnet_pred, RandomForest = rf_pred, Ensemble = ensemble_pred)
  model_aucs <- list(GBM = gbm_auc, ElasticNet = glmnet_auc, RandomForest = rf_auc, Ensemble = ensemble_auc)
  
  for(model_name in names(model_preds)) {
    pred_class <- ifelse(model_preds[[model_name]] > thresh_val, 1, 0)
    pred_class <- factor(pred_class, levels = c(0, 1))
    actual <- factor(y_test, levels = c(0, 1))
    
    # Set positive class = "1"
    cm <- confusionMatrix(pred_class, actual, positive = "1")
    individual_results[[model_name]] <- list(
      accuracy = cm$overall['Accuracy'],
      balanced_accuracy = cm$byClass['Balanced Accuracy'],
      kappa = cm$overall['Kappa'],
      sensitivity = cm$byClass['Sensitivity'],
      specificity = cm$byClass['Specificity'],
      precision = cm$byClass['Pos Pred Value'],
      auc = as.numeric(model_aucs[[model_name]])
    )
    
    cat(sprintf("%s - AUC: %.3f, Accuracy: %.3f, Balanced Acc: %.3f\n", 
                model_name, as.numeric(model_aucs[[model_name]]), 
                cm$overall['Accuracy'], cm$byClass['Balanced Accuracy']))
  }
  
  return(list(
    threshold_results = threshold_results,
    best_threshold = thresh_val,
    individual_results = individual_results,
    predictions = model_preds,
    rf_params = models$best_rf_params
  ))
}

###########################################################################
# MAIN EVALUATION FUNCTION - SIMPLIFIED
###########################################################################

evaluate_models <- function(data, target_col = "gt", test_prop = 0.3, seed = 123) {
  
  set.seed(seed)
  
  # Single train/test split
  train_idx <- createDataPartition(data[[target_col]], p = 1 - test_prop, list = FALSE)
  test_idx <- setdiff(1:nrow(data), train_idx)
  
  # Split data
  train_data <- data[train_idx, ]
  test_data <- data[test_idx, ]
  
  cat("Train samples:", nrow(train_data), "Test samples:", nrow(test_data), "\n")
  cat("Features:", ncol(data) - 1, "\n\n")
  
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
    test_idx = test_idx
  ))
}

###########################################################################
# UPDATED CROSS-VALIDATION FUNCTION WITH AUC AND SE
###########################################################################

run_cv_all_models <- function(data, target_col = "gt", n_runs = 10, test_prop = 0.3) {
  
  cat("=== RUNNING", n_runs, "RUNS FOR ALL MODELS WITH AUC ===\n")
  
  # Initialize storage for all models
  all_results <- list(
    GBM = list(accuracy = c(), balanced_accuracy = c(), auc = c()),
    ElasticNet = list(accuracy = c(), balanced_accuracy = c(), auc = c()),
    RandomForest = list(accuracy = c(), balanced_accuracy = c(), auc = c()),
    Ensemble = list(accuracy = c(), balanced_accuracy = c(), auc = c())
  )
  
  rf_params_list <- list()
  
  for(i in 1:n_runs) {
    cat("Run", i, "...")
    
    # Use different seed for each run
    result <- evaluate_models(
      data, 
      target_col = target_col, 
      test_prop = test_prop,
      seed = 123 + i * 10
    )
    
    # Extract metrics for each model
    for(model_name in names(all_results)) {
      metrics <- result$results$individual_results[[model_name]]
      all_results[[model_name]]$accuracy[i] <- metrics$accuracy
      all_results[[model_name]]$balanced_accuracy[i] <- metrics$balanced_accuracy
      all_results[[model_name]]$auc[i] <- metrics$auc
    }
    
    # Store RF parameters
    rf_params_list[[i]] <- result$results$rf_params
    
    cat(" Ensemble AUC:", round(all_results$Ensemble$auc[i], 3),
        " Acc:", round(all_results$Ensemble$accuracy[i], 3), "\n")
  }
  
  # Calculate summary statistics with AUC and SE
  summary_stats <- data.frame(
    Model = character(),
    Mean_Accuracy = numeric(),
    SD_Accuracy = numeric(),
    SE_Accuracy = numeric(),
    Mean_Balanced_Accuracy = numeric(),
    SD_Balanced_Accuracy = numeric(),
    SE_Balanced_Accuracy = numeric(),
    Mean_AUC = numeric(),
    SD_AUC = numeric(),
    SE_AUC = numeric(),
    AUC_CI_Lower = numeric(),
    AUC_CI_Upper = numeric(),
    N_Features = integer(),
    stringsAsFactors = FALSE
  )
  
  for(model_name in names(all_results)) {
    # Calculate means, SDs, and SEs
    acc_vals <- all_results[[model_name]]$accuracy
    bal_acc_vals <- all_results[[model_name]]$balanced_accuracy
    auc_vals <- all_results[[model_name]]$auc
    
    # Remove NA values
    acc_vals <- acc_vals[!is.na(acc_vals)]
    bal_acc_vals <- bal_acc_vals[!is.na(bal_acc_vals)]
    auc_vals <- auc_vals[!is.na(auc_vals)]
    
    n_valid <- length(auc_vals)
    
    # AUC statistics
    mean_auc <- mean(auc_vals)
    sd_auc <- sd(auc_vals)
    se_auc <- sd_auc / sqrt(n_valid)
    ci_lower <- mean_auc - 1.96 * se_auc
    ci_upper <- mean_auc + 1.96 * se_auc
    
    summary_stats <- rbind(summary_stats, data.frame(
      Model = model_name,
      Mean_Accuracy = round(mean(acc_vals), 3),
      SD_Accuracy = round(sd(acc_vals), 3),
      SE_Accuracy = round(sd(acc_vals) / sqrt(length(acc_vals)), 3),
      Mean_Balanced_Accuracy = round(mean(bal_acc_vals), 3),
      SD_Balanced_Accuracy = round(sd(bal_acc_vals), 3),
      SE_Balanced_Accuracy = round(sd(bal_acc_vals) / sqrt(length(bal_acc_vals)), 3),
      Mean_AUC = round(mean_auc, 3),
      SD_AUC = round(sd_auc, 3),
      SE_AUC = round(se_auc, 3),
      AUC_CI_Lower = round(ci_lower, 3),
      AUC_CI_Upper = round(ci_upper, 3),
      N_Features = ncol(data) - 1,
      stringsAsFactors = FALSE
    ))
  }
  
  # Print summary with AUC
  cat("\n=== SUMMARY STATISTICS WITH AUC ± SE ===\n")
  print(summary_stats[, c("Model", "N_Features", "Mean_AUC", "SE_AUC", "AUC_CI_Lower", "AUC_CI_Upper",
                          "Mean_Balanced_Accuracy", "SE_Balanced_Accuracy")])
  
  # Create comparison table
  cat("\n=== COMPARISON TABLE (AUC ± SE) ===\n")
  comparison_table <- data.frame(
    Model = summary_stats$Model,
    N_Features = summary_stats$N_Features,
    AUC = paste0(summary_stats$Mean_AUC, " ± ", summary_stats$SE_AUC),
    AUC_CI = paste0("(", summary_stats$AUC_CI_Lower, "-", summary_stats$AUC_CI_Upper, ")"),
    Balanced_Accuracy = paste0(summary_stats$Mean_Balanced_Accuracy, " ± ", summary_stats$SE_Balanced_Accuracy),
    stringsAsFactors = FALSE
  )
  
  print(comparison_table)
  
  # Compare to chance (AUC = 0.5)
  cat("\n=== STATISTICAL SIGNIFICANCE (vs AUC = 0.5) ===\n")
  for(model_name in names(all_results)) {
    auc_vals <- all_results[[model_name]]$auc
    auc_vals <- auc_vals[!is.na(auc_vals)]
    
    if(length(auc_vals) > 1) {
      t_test <- t.test(auc_vals, mu = 0.5, alternative = "greater")
      sig_text <- ifelse(t_test$p.value < 0.05, "significant", "not significant")
      cat(sprintf("%s: AUC = %.3f ± %.3f (vs 0.5, p = %.4f, %s)\n", 
                  model_name, mean(auc_vals), sd(auc_vals), 
                  t_test$p.value, sig_text))
    }
  }
  
  # RF parameter analysis
  cat("\n=== RF PARAMETER VARIABILITY ===\n")
  mtry_values <- sapply(rf_params_list, function(x) x$mtry)
  nodesize_values <- sapply(rf_params_list, function(x) x$nodesize)
  
  cat("mtry distribution:\n")
  print(table(mtry_values))
  cat("nodesize distribution:\n")
  print(table(nodesize_values))
  
  # Create detailed results dataframe
  detailed_results <- data.frame(
    Run = rep(1:n_runs, 4),
    Model = rep(names(all_results), each = n_runs),
    Accuracy = c(all_results$GBM$accuracy, all_results$ElasticNet$accuracy, 
                all_results$RandomForest$accuracy, all_results$Ensemble$accuracy),
    Balanced_Accuracy = c(all_results$GBM$balanced_accuracy, all_results$ElasticNet$balanced_accuracy,
                         all_results$RandomForest$balanced_accuracy, all_results$Ensemble$balanced_accuracy),
    AUC = c(all_results$GBM$auc, all_results$ElasticNet$auc,
           all_results$RandomForest$auc, all_results$Ensemble$auc),
    RF_mtry = rep(mtry_values, 4),
    RF_nodesize = rep(nodesize_values, 4),
    stringsAsFactors = FALSE
  )
  
  return(list(
    summary_stats = summary_stats,
    detailed_results = detailed_results,
    all_results = all_results,
    rf_params = rf_params_list,
    comparison_table = comparison_table
  ))
}

###########################################################################
# USAGE EXAMPLE
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

# Run comprehensive evaluation for all models with AUC
cv_results <- run_cv_all_models(boost_data, target_col = "gt", n_runs = 10, test_prop = 0.3)

# Save results
write.csv(cv_results$summary_stats, "../results/metab_rf/all_models_summary_cv_with_auc.csv", row.names = FALSE)
write.csv(cv_results$detailed_results, "../results/metab_rf/all_models_detailed_cv_with_auc.csv", row.names = FALSE)
write.csv(cv_results$comparison_table, "../results/metab_rf/comparison_table_auc.csv", row.names = FALSE)

cat("\nResults with AUC ± SE saved to CSV files\n")

# Stop parallel processing
stopImplicitCluster()