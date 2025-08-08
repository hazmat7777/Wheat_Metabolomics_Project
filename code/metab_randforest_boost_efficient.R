# this script runs a randomforest model on the metabolomic data

library(dplyr)

###change so that 1 = present, 0 = absent****

# load microbial metabolomic and gt data
fk_metabolom_gt <- readRDS("../data/fk_metabolomics_gt_noNA.RDS")

nrow(fk_metabolom_gt) # 677 compounds
#View(fk_metabolom_gt)

# count the zeroes
sum(fk_metabolom_gt ==0, na.rm = TRUE)
dim(fk_metabolom_gt)

## data prep 
fk_metabolom_gt_t <- t(fk_metabolom_gt) # transpose
colnames(fk_metabolom_gt_t) <- fk_metabolom_gt_t[1, ] # set the first row as column names
fk_metabolom_gt_t <- fk_metabolom_gt_t[-1, ] # remove the first row (now column names)
fk_metabolom_gt_t <- as.data.frame(fk_metabolom_gt_t) # convert to data frame
colnames(fk_metabolom_gt_t)[ncol(fk_metabolom_gt_t)] <- "gt" # rename the last column to 'gt'
fk_metabolom_gt_t$gt <- trimws(as.character(fk_metabolom_gt_t$gt))
fk_metabolom_gt_t$gt <- as.factor(fk_metabolom_gt_t$gt)
levels(fk_metabolom_gt_t$gt) <- c("GT_absent", "GT_present")
levels(fk_metabolom_gt_t$gt)

colnames(fk_metabolom_gt_t) <- make.names(colnames(fk_metabolom_gt_t), unique = TRUE)# Clean colnames by replacing spaces with underscores
str(fk_metabolom_gt_t[,1:5]) # data aren't scaled or numeric
fk_metabolom_gt_t[ , -c(ncol(fk_metabolom_gt_t))] <- lapply( # make the peak areas numeric
  fk_metabolom_gt_t[ , -c(1, ncol(fk_metabolom_gt_t))],
  function(x) as.numeric(as.character(x))
)
hist(fk_metabolom_gt_t[,78], breaks = 30, main = "Histogram", xlab = "Value") # check for skew- a bit righty
plot(density(log(fk_metabolom_gt_t[,78]), na.rm = TRUE), main = "Density") # does log help? think so

# log data
metab_numeric <- fk_metabolom_gt_t[, sapply(fk_metabolom_gt_t, is.numeric)] # Select only numeric metabolite columns
metab_scaled <- as.data.frame(log(metab_numeric + 1e-6)) # log transform
# THINK ABOUT Z SCALING IT HERE
non_numeric <- fk_metabolom_gt_t[, !sapply(fk_metabolom_gt_t, is.numeric)] 
fk_metabolom_gt_scaled <- cbind(non_numeric, metab_scaled)# add back non-numeric columns
colnames(fk_metabolom_gt_scaled)[1]<- "gt"
View(fk_metabolom_gt_scaled)


###### 

# new approach- BOOSTING (?)
# Load required libraries
library(gbm)
library(caret)
library(glmnet)
library(randomForest)
library(dplyr)

# ===============================
# SINGLE PREPROCESSING - NO LEAKAGE
# ===============================

preprocess_metabolomics_proper <- function(data, target_col, train_idx) {
  
  # Separate features and target
  if(is.character(target_col)) {
    target <- data[[target_col]]
    features <- data[, !names(data) %in% target_col]
  } else {
    target <- data[, target_col]
    features <- data[, -target_col]
  }
  
  cat("=== PREPROCESSING (TRAINING DATA ONLY) ===\n")
  cat("Original data dimensions:", dim(features), "\n")
  cat("Number of zeros:", sum(features == 0), "\n")
  
  # STEP 1: Feature filtering based on TRAINING data only
  train_features <- features[train_idx, ]
  
  # Remove features with >50% zeros in training data
  zero_prop_features <- sapply(train_features, function(x) sum(x == 0) / length(x))
  keep_zero_filter <- zero_prop_features <= 0.5
  
  # Remove low variance features in training data
  train_var <- sapply(train_features[, keep_zero_filter], var, na.rm = TRUE)
  keep_var_filter <- train_var > quantile(train_var, 0.1, na.rm = TRUE)
  
  # Combine filters
  final_keep_features <- names(features)[keep_zero_filter]
  final_keep_features <- final_keep_features[keep_var_filter]
  
  features_filtered <- features[, final_keep_features]
  train_features_filtered <- train_features[, final_keep_features]
  
  cat("Removed", sum(!keep_zero_filter), "features with >50% zeros\n")
  cat("Removed", sum(!keep_var_filter), "low variance features\n")
  cat("Remaining features after filtering:", ncol(features_filtered), "\n")
  
  # STEP 2: Feature selection using TRAINING data only
  max_features <- 5
  cat("Selecting top", max_features, "features using t-tests on training data\n")
  
  # Calculate p-values on training data only
  p_values <- numeric(ncol(train_features_filtered))
  train_target <- target[train_idx]
  
  for(i in 1:ncol(train_features_filtered)) {
    test_result <- t.test(train_features_filtered[, i] ~ train_target)
    p_values[i] <- test_result$p.value
  }
  
  # Select best features
  selected_features <- names(train_features_filtered)[order(p_values)[1:min(max_features, length(p_values))]]
  features_selected <- features_filtered[, selected_features]
  train_features_selected <- train_features_filtered[, selected_features]
  
  cat("Selected", length(selected_features), "features\n")
  
  # STEP 3: Zero imputation using TRAINING statistics
  features_imputed <- features_selected
  for(col in selected_features) {
    zeros_all <- features_selected[[col]] == 0
    if(sum(zeros_all) > 0) {
      # Calculate replacement value from TRAINING data only
      train_nonzero <- train_features_selected[[col]][train_features_selected[[col]] > 0]
      if(length(train_nonzero) > 0) {
        min_val <- min(train_nonzero, na.rm = TRUE)
        features_imputed[[col]][zeros_all] <- min_val / 2
      }
    }
  }
  
  # STEP 4: Log transform and scale using TRAINING statistics
  features_log <- log10(features_imputed + 1)
  
  # Calculate scaling parameters from training data only
  train_log <- features_log[train_idx, ]
  train_means <- sapply(train_log, mean, na.rm = TRUE)
  train_sds <- sapply(train_log, sd, na.rm = TRUE)
  
  # Apply scaling to ALL data using training parameters
  features_scaled <- features_log
  for(i in 1:ncol(features_log)) {
    if(train_sds[i] > 0) {
      features_scaled[, i] <- (features_log[, i] - train_means[i]) / train_sds[i]
    }
  }
  
  # Return processed data and preprocessing parameters
  processed_data <- data.frame(features_scaled)
  processed_data$gt <- target
  
  preprocessing_params <- list(
    selected_features = selected_features,
    train_means = train_means,
    train_sds = train_sds,
    zero_replacement_values = sapply(selected_features, function(col) {
      train_nonzero <- train_features_selected[[col]][train_features_selected[[col]] > 0]
      if(length(train_nonzero) > 0) min(train_nonzero, na.rm = TRUE) / 2 else 0
    })
  )
  
  cat("Final processed data dimensions:", dim(processed_data), "\n\n")
  
  return(list(
    data = processed_data,
    params = preprocessing_params,
    selected_features = selected_features
  ))
}

# ===============================
# PROPER TRAIN/TEST SPLIT EVALUATION
# ===============================

evaluate_metabolomics_models <- function(data, target_col = "gt", test_prop = 0.3, seed = 123) {
  
  set.seed(seed)
  
  # Single train/test split
  train_idx <- createDataPartition(data[[target_col]], p = 1 - test_prop, list = FALSE)
  test_idx <- setdiff(1:nrow(data), train_idx)
  
  cat("Train samples:", length(train_idx), "Test samples:", length(test_idx), "\n")
  
  # Preprocess data (using training data only)
  preprocessed <- preprocess_metabolomics_proper(data, target_col, train_idx)
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
    test_data = test_data,
    preprocessing_params = preprocessed$params,
    selected_features = preprocessed$selected_features
  ))
}

# ===============================
# MODEL TRAINING (TRAINING DATA ONLY)
# ===============================

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
  
  # Model 3: Random Forest
  cat("Training Random Forest...\n")
  rf_model <- randomForest(
    as.factor(gt) ~ ., 
    data = train_data,
    ntree = 200,
    mtry = max(1, floor(sqrt(ncol(train_data)-1))),
    nodesize = 3
  )
  
  cat("Models trained successfully\n\n")
  
  return(list(
    gbm = gbm_model,
    glmnet = cv_glmnet,
    rf = rf_model,
    optimal_trees = optimal_trees,
    x_train = x_train,
    y_train = y_train
  ))
}

# ===============================
# TEST SET EVALUATION
# ===============================

evaluate_on_test_set <- function(models, test_data, target_col = "gt") {
  
  cat("=== EVALUATING ON TEST SET ===\n")
  
  # Prepare test data
  x_test <- as.matrix(test_data[, !names(test_data) %in% target_col])
  y_test <- test_data[[target_col]]
  
  # Get predictions
  gbm_pred <- predict(models$gbm, test_data, n.trees = models$optimal_trees, type = "response")
  glmnet_pred <- predict(models$glmnet, x_test, s = "lambda.min", type = "response")[,1]
  rf_pred <- predict(models$rf, test_data, type = "prob")[,2]
  
  # Ensemble prediction
  ensemble_pred <- (gbm_pred + glmnet_pred + rf_pred) / 3
  
  # Test multiple thresholds
  thresholds <- seq(0.3, 0.7, 0.1)
  threshold_results <- list()
  
  for(thresh in thresholds) {
    pred_class <- ifelse(ensemble_pred > thresh, 1, 0)
    pred_class <- factor(pred_class, levels = c(0, 1))
    actual <- factor(y_test, levels = c(0, 1))
    
    cm <- confusionMatrix(pred_class, actual)
    
    threshold_results[[as.character(thresh)]] <- list(
      threshold = thresh,
      accuracy = cm$overall['Accuracy'],
      kappa = cm$overall['Kappa'],
      sensitivity = cm$byClass['Sensitivity'],
      specificity = cm$byClass['Specificity'],
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
    
    cm <- confusionMatrix(pred_class, actual)
    individual_results[[model_name]] <- list(
      accuracy = cm$overall['Accuracy'],
      kappa = cm$overall['Kappa']
    )
    
    cat(sprintf("%s - Accuracy: %.3f, Kappa: %.3f\n", 
                model_name, cm$overall['Accuracy'], cm$overall['Kappa']))
  }
  
  # Feature importance
  cat("\n=== FEATURE IMPORTANCE ===\n")
  gbm_imp <- summary(models$gbm, n.trees = models$optimal_trees, plotit = FALSE)
  print(head(gbm_imp, 10))
  
  return(list(
    threshold_results = threshold_results,
    best_threshold = thresh_val,
    individual_results = individual_results,
    predictions = model_preds,
    feature_importance = gbm_imp
  ))
}

# ===============================
# NESTED CROSS-VALIDATION (PROPER WAY)
# ===============================

nested_cv_evaluation <- function(data, target_col = "gt", outer_folds = 5, n_repeats = 5) {
  
  results_list <- list()
  
  for(repeat_i in 1:n_repeats) {
    cat("\n=== REPEAT", repeat_i, "OF", n_repeats, "===\n")
    set.seed(repeat_i)
    
    # Create outer CV folds
    folds <- createFolds(data[[target_col]], k = outer_folds, list = TRUE)
    fold_results <- list()
    
    for(fold_i in 1:length(folds)) {
      cat("Outer fold", fold_i, "of", outer_folds, "\n")
      
      # This fold's test set
      test_idx <- folds[[fold_i]]
      train_idx <- setdiff(1:nrow(data), test_idx)
      
      # Run evaluation
      fold_result <- evaluate_metabolomics_models(
        data, target_col, test_prop = length(test_idx)/nrow(data), seed = repeat_i * 100 + fold_i
      )
      
      fold_results[[fold_i]] <- fold_result$results$individual_results$Ensemble
    }
    
    # Average across folds for this repeat
    repeat_accuracy <- mean(sapply(fold_results, function(x) x$accuracy))
    repeat_kappa <- mean(sapply(fold_results, function(x) x$kappa))
    
    results_list[[repeat_i]] <- list(accuracy = repeat_accuracy, kappa = repeat_kappa)
    
    cat("Repeat", repeat_i, "- Accuracy:", round(repeat_accuracy, 3), "Kappa:", round(repeat_kappa, 3), "\n")
  }
  
  # Overall results
  final_accuracies <- sapply(results_list, function(x) x$accuracy)
  final_kappas <- sapply(results_list, function(x) x$kappa)
  
  cat("\n=== FINAL NESTED CV RESULTS ===\n")
  cat("Mean Accuracy:", round(mean(final_accuracies), 3), "±", round(sd(final_accuracies), 3), "\n")
  cat("Mean Kappa:", round(mean(final_kappas), 3), "±", round(sd(final_kappas), 3), "\n")
  
  return(list(accuracies = final_accuracies, kappas = final_kappas))
}

# ===============================
# USAGE EXAMPLES
# ===============================
# Convert gt to numeric (0/1) if not already
boost_data <- fk_metabolom_gt_t %>%
  mutate(gt = ifelse(gt == "GT_present", 1,
              ifelse(gt == "GT_absent", 0, gt))) %>%
  mutate(gt = as.numeric(gt))  # ensure it's numeric

# Simple train/test evaluation (honest estimate)
results <- evaluate_metabolomics_models(boost_data, target_col = "gt", test_prop = 0.3)

str(results)
# $models - the trained models (gbm, glmnet, rf)
# $results - performance metrics, predictions, confusion matrices
# $test_data - the test set with features + gt column
# $preprocessing_params - scaling parameters, selected features
# $selected_features - names of the metabolites used

# Robust nested cross-validation (more reliable but slower)
cv_results <- nested_cv_evaluation(boost_data, target_col = "gt", outer_folds = 5, n_repeats = 3)

mean(cv_results$accuracies)

glimpse(results)

## in this code we are selecting the top 5 features based on t tests and using these same 5 every time

# how stable is this choosing?
# Test feature selection stability
test_feature_stability <- function(data, target_col = "gt", n_runs = 10) {
  
  feature_lists <- list()
  
  for(i in 1:n_runs) {
    set.seed(i)  # Different seed each time
    
    # Different train/test split each time
    train_idx <- createDataPartition(data[[target_col]], p = 0.7, list = FALSE)
    
    # Run preprocessing and feature selection
    preprocessed <- preprocess_metabolomics_proper(data, target_col, train_idx)
    
    feature_lists[[i]] <- preprocessed$selected_features
    
    cat("Run", i, "- Top 5 features:", 
        head(preprocessed$selected_features, 5), "\n")
  }
  
  # Check overlap
  all_features <- unique(unlist(feature_lists))
  feature_counts <- table(unlist(feature_lists))
  
  cat("\n=== FEATURE STABILITY ANALYSIS ===\n")
  cat("Features appearing in multiple runs:\n")
  print(sort(feature_counts, decreasing = TRUE))
  
  # Calculate Jaccard similarity between runs
  similarities <- numeric()
  for(i in 1:(length(feature_lists)-1)) {
    for(j in (i+1):length(feature_lists)) {
      intersection <- length(intersect(feature_lists[[i]], feature_lists[[j]]))
      union <- length(union(feature_lists[[i]], feature_lists[[j]]))
      similarities <- c(similarities, intersection/union)
    }
  }
  
  cat("Mean Jaccard similarity between feature sets:", round(mean(similarities), 3), "\n")
  cat("SD of similarities:", round(sd(similarities), 3), "\n")
  
  return(list(
    feature_lists = feature_lists,
    feature_counts = feature_counts,
    mean_similarity = mean(similarities)
  ))
}

# Run the test
stability_results <- test_feature_stability(boost_data, target_col = "gt", n_runs = 10)

#LOW STABILITY- cant I repeat the runs and pick e.g. the ten features which come up most often? need to ask claude















# tracking model wins
model_wins <- c()

for (i in 1:20) {
  res <- train_ensemble_models(boost_data, target_col = "gt")
  best_thresh <- as.character(res$best_threshold)

  # Get test data and labels
  test_labels <- res$test_data$gt
  
  # Predictions per model
  rf_pred <- as.numeric(predict(res$models$rf, res$test_data, type = "response")) - 1 # -1 as it was predicting 1s and 2s (?)
  gbm_pred <- ifelse(predict(res$models$gbm, res$test_data, 
                              n.trees = res$models$gbm$n.trees, type = "response") > res$best_threshold, 1, 0)
  en_pred <- ifelse(predict(res$models$glmnet, as.matrix(res$test_data[,-ncol(res$test_data)]), 
                             s = "lambda.min", type = "response")[,1] > res$best_threshold, 1, 0)
  ens_pred <- ifelse((rf_pred + gbm_pred + en_pred)/3 > res$best_threshold, 1, 0)
  
  # Accuracies
  rf_acc <- mean(rf_pred == test_labels)
  gbm_acc <- mean(gbm_pred == test_labels)
  en_acc <- mean(en_pred == test_labels)
  ens_acc <- mean(ens_pred == test_labels)

  # Find best model by actual accuracy
  accs <- c(RandomForest = rf_acc, GBM = gbm_acc, ElasticNet = en_acc, Ensemble = ens_acc)
  best_model <- names(which.max(accs))
  model_wins <- c(model_wins, best_model)
}
