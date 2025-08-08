# this script runs a randomforest model on the metabolomic data

library(dplyr)

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
library(glmnet)  # for regularized regression
library(randomForest)
library(dplyr)
library(VIM)

# ===============================
# IMPROVED PREPROCESSING FOR SMALL DATASETS
# ===============================

preprocess_small_dataset <- function(data, target_col, train_idx) {
  
  # Separate features and target
  if(is.character(target_col)) {
    target <- data[[target_col]]
    features <- data[, !names(data) %in% target_col]
  } else {
    target <- data[, target_col]
    features <- data[, -target_col]
  }
  
  cat("Original data dimensions:", dim(features), "\n")
  cat("Number of zeros:", sum(features == 0), "\n")
  
  # CRITICAL: Only use training data for preprocessing decisions
  train_features <- features[train_idx, ]
  
  # Step 1: More aggressive feature filtering for small datasets
  # Remove features with >50% zeros in TRAINING data only
  zero_prop_features <- sapply(train_features, function(x) sum(x == 0) / length(x))
  keep_features <- zero_prop_features <= 0.5
  
  features_filtered <- features[, keep_features]
  train_features_filtered <- train_features[, keep_features]
  
  cat("Removed", sum(!keep_features), "features with >50% zeros\n")
  cat("Remaining features:", ncol(features_filtered), "\n")
  
  # Step 2: Remove low variance features in training data
  train_var <- sapply(train_features_filtered, var, na.rm = TRUE)
  keep_var <- train_var > quantile(train_var, 0.1, na.rm = TRUE)  # Keep top 90%
  
  features_filtered <- features_filtered[, keep_var]
  train_features_filtered <- train_features_filtered[, keep_var]
  
  cat("Removed", sum(!keep_var), "low variance features\n")
  cat("Remaining features:", ncol(features_filtered), "\n")
  
  # Step 3: Handle zeros using training data statistics only
  features_imputed <- features_filtered
  for(col in names(features_filtered)) {
    zeros <- features_filtered[[col]] == 0
    if(sum(zeros) > 0) {
      # Use minimum from TRAINING data only
      min_val <- min(train_features_filtered[[col]][train_features_filtered[[col]] > 0], na.rm = TRUE)
      if(is.finite(min_val)) {
        features_imputed[[col]][zeros] <- min_val / 2
      }
    }
  }
  
  # Step 4: Log transform and scale using training statistics
  features_log <- log10(features_imputed + 1)
  
  # Calculate scaling parameters from training data only
  train_means <- sapply(features_log[train_idx, ], mean, na.rm = TRUE)
  train_sds <- sapply(features_log[train_idx, ], sd, na.rm = TRUE)
  
  # Apply scaling to all data using training statistics
  features_scaled <- features_log
  for(i in 1:ncol(features_log)) {
    if(train_sds[i] > 0) {
      features_scaled[, i] <- (features_log[, i] - train_means[i]) / train_sds[i]
    }
  }
  
  # Combine back
  processed_data <- data.frame(features_scaled)
  processed_data$gt <- target
  
  cat("Final data dimensions:", dim(processed_data), "\n")
  
  return(processed_data)
}

# ===============================
# FEATURE SELECTION FOR SMALL DATASETS
# ===============================

select_important_features <- function(data, target_col = "gt", train_idx, max_features = NULL) {
  
  train_data <- data[train_idx, ]
  
  # Determine max features based on sample size (rule of thumb: n_samples/10)
  if(is.null(max_features)) {
    max_features <- max(5, min(50, floor(nrow(train_data) / 10)))
  }
  
  cat("Selecting top", max_features, "features\n")
  
  # Method 1: Univariate feature selection (t-test or correlation)
  feature_cols <- !names(train_data) %in% target_col
  p_values <- numeric(sum(feature_cols))
  
  for(i in which(feature_cols)) {
    test_result <- t.test(train_data[, i] ~ train_data[[target_col]])
    p_values[i] <- test_result$p.value
  }
  
  # Select features with lowest p-values
  feature_names <- names(train_data)[feature_cols]
  selected_features <- feature_names[order(p_values[feature_cols])[1:min(max_features, length(feature_names))]]
  
  # Create reduced dataset
  reduced_data <- data[, c(selected_features, target_col)]
  
  cat("Selected features:", length(selected_features), "\n")
  return(reduced_data)
}

# ===============================
# ENSEMBLE APPROACH FOR SMALL DATASETS
# ===============================

train_ensemble_models <- function(data, target_col = "gt", train_prop = 0.7) {
  
  train_idx <- createDataPartition(data[[target_col]], p = train_prop, list = FALSE) # select training data
  
  # Preprocess data properly
  processed_data <- preprocess_small_dataset(data, target_col, train_idx)
  
  # Feature selection
  selected_data <- select_important_features(processed_data, target_col, train_idx, max_features = 15)
  
  train_data <- selected_data[train_idx, ]
  test_data <- selected_data[-train_idx, ]
  
  cat("\nFinal training data dimensions:", dim(train_data), "\n")
  cat("Final test data dimensions:", dim(test_data), "\n")
  
  # Model 1: Simplified GBM
  cat("\nTraining GBM...\n")
  gbm_model <- gbm(
    gt ~ ., 
    data = train_data,
    distribution = "bernoulli",
    n.trees = 100,  # Much fewer trees for small dataset
    interaction.depth = 2,  # Simpler interactions
    shrinkage = 0.1,
    n.minobsinnode = 3,  # Lower minimum
    bag.fraction = 0.8,
    cv.folds = 5
  )
  
  optimal_trees <- gbm.perf(gbm_model, method = "cv", plot.it = FALSE)
  gbm_pred <- predict(gbm_model, test_data, n.trees = optimal_trees, type = "response")
  
  # Model 2: Elastic Net (good for small datasets)
  cat("Training Elastic Net...\n")
  x_train <- as.matrix(train_data[, !names(train_data) %in% target_col])
  y_train <- train_data[[target_col]]
  x_test <- as.matrix(test_data[, !names(test_data) %in% target_col])
  
  cv_glmnet <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5)
  glmnet_pred <- predict(cv_glmnet, x_test, s = "lambda.min", type = "response")[,1]
  
  # Model 3: Random Forest
  cat("Training Random Forest...\n")
  rf_model <- randomForest(
    as.factor(gt) ~ ., 
    data = train_data,
    ntree = 100,
    mtry = max(1, floor(sqrt(ncol(train_data)-1))),
    nodesize = 3
  )
  
  #cat("Test data dimensions before prediction:", dim(test_data), "\n")  # Should be 21 rows
  
  rf_pred <- predict(rf_model, test_data, type = "prob")[,2]
  
  #cat("Length of Random Forest Predictions:", length(rf_pred), "\n")  # Should be 21

  # Ensemble prediction (average)
  ensemble_pred <- (gbm_pred + glmnet_pred + rf_pred) / 3
  
  # Test different thresholds
  results <- list()
  thresholds <- seq(0.3, 0.7, 0.1)
  
  for(thresh in thresholds) {
    pred_class <- ifelse(ensemble_pred > thresh, 1, 0)
    pred_class <- factor(pred_class, levels = c(0, 1))
    actual <- factor(test_data[[target_col]], levels = c(0, 1))
    
    cm <- confusionMatrix(pred_class, actual)
    
    results[[as.character(thresh)]] <- list(
      threshold = thresh,
      accuracy = cm$overall['Accuracy'],
      kappa = cm$overall['Kappa'],
      sensitivity = cm$byClass['Sensitivity'],
      specificity = cm$byClass['Specificity'],
      confusion_matrix = cm
    )
  }
  
  # Find best threshold
  kappa_scores <- sapply(results, function(x) x$kappa)
  best_thresh <- names(results)[which.max(kappa_scores)]
  
  cat("\n=== ENSEMBLE RESULTS ===\n")
  cat("Best threshold:", best_thresh, "\n")
  print(results[[best_thresh]]$confusion_matrix)
  
  # Individual model performance at best threshold
  cat("\n=== INDIVIDUAL MODEL COMPARISON ===\n")
  thresh_val <- as.numeric(best_thresh)
  
  models_pred <- list(
    GBM = gbm_pred,
    ElasticNet = glmnet_pred,
    RandomForest = rf_pred,
    Ensemble = ensemble_pred
  )
  
  for(model_name in names(models_pred)) {
    pred_class <- ifelse(models_pred[[model_name]] > thresh_val, 1, 0)
    pred_class <- factor(pred_class, levels = c(0, 1))
    actual <- factor(test_data[[target_col]], levels = c(0, 1))
    
    cm <- confusionMatrix(pred_class, actual)
    cat(sprintf("%s - Accuracy: %.3f, Kappa: %.3f\n", 
                model_name, cm$overall['Accuracy'], cm$overall['Kappa']))
  }
  
  # Feature importance from different models
  cat("\n=== FEATURE IMPORTANCE ===\n")
  
  # GBM importance
  gbm_imp <- summary(gbm_model, n.trees = optimal_trees, plotit = FALSE)
  cat("Top GBM features:\n")
  print(head(gbm_imp[gbm_imp$rel.inf > 0, ], 10))
  
  return(list(
    models = list(gbm = gbm_model, glmnet = cv_glmnet, rf = rf_model),
    results = results,
    best_threshold = thresh_val,
    test_data = test_data,
    feature_importance = gbm_imp
  ))
}

# ===============================
# CROSS-VALIDATION FOR SMALL DATASETS
# ===============================

repeated_cv_evaluation <- function(data, target_col = "gt", n_repeats = 10) {
  
  results <- list()
  
  for(i in 1:n_repeats) {
    cat("Repeat", i, "of", n_repeats, "\n")
    set.seed(i)
    
    # Different random split each time
    result <- train_ensemble_models(data, target_col)
    results[[i]] <- result$results[[as.character(result$best_threshold)]]
  }
  
  # Aggregate results
  accuracies <- sapply(results, function(x) x$accuracy)
  kappas <- sapply(results, function(x) x$kappa)
  
  cat("\n=== REPEATED CV RESULTS ===\n")
  cat("Mean Accuracy:", round(mean(accuracies), 3), "±", round(sd(accuracies), 3), "\n")
  cat("Mean Kappa:", round(mean(kappas), 3), "±", round(sd(kappas), 3), "\n")
  
  return(list(accuracies = accuracies, kappas = kappas))
}

# ===============================
# USAGE EXAMPLES
# ===============================
set.seed(123)

# Convert gt to numeric (0/1) if not already
boost_data <- fk_metabolom_gt_t %>%
  mutate(gt = ifelse(gt == "GT_present", 1,
              ifelse(gt == "GT_absent", 0, gt))) %>%
  mutate(gt = as.numeric(gt))  # ensure it's numeric


# Single run with ensemble
results <- train_ensemble_models(boost_data, target_col = "gt")

# More robust evaluation with repeated CV
cv_results <- repeated_cv_evaluation(boost_data, target_col = "gt", n_repeats = 20)

results
cv_results

# Save winner per run
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

# storing more results in a df

library(future)
library(future.apply)

# Use only a portion of available cores
available_cores <- future::availableCores()
cores_to_use <- max(1, available_cores - 2)  # Leave 2 cores free

plan(multisession, workers = cores_to_use)

# Run 100 simulations in parallel
results_list <- future_lapply(1:100, function(i) {
  res <- train_ensemble_models(boost_data, target_col = "gt")
  best_thresh <- as.numeric(res$best_threshold)
  test_labels <- factor(res$test_data$gt, levels = c(0, 1))

  # Random Forest
  rf_raw <- predict(res$models$rf, res$test_data, type = "response")
  rf_prob <- as.numeric(rf_raw) - 1
  rf_class <- factor(ifelse(rf_prob > best_thresh, 1, 0), levels = c(0, 1))
  cm_rf <- confusionMatrix(rf_class, test_labels)
  auc_rf <- tryCatch(auc(test_labels, rf_prob), error = function(e) NA)

  # GBM
  gbm_prob <- predict(res$models$gbm, res$test_data, n.trees = res$models$gbm$n.trees, type = "response")
  gbm_class <- factor(ifelse(gbm_prob > best_thresh, 1, 0), levels = c(0, 1))
  cm_gbm <- confusionMatrix(gbm_class, test_labels)
  auc_gbm <- tryCatch(auc(test_labels, gbm_prob), error = function(e) NA)

  # Elastic Net
  en_prob <- predict(res$models$glmnet,
                     as.matrix(res$test_data[, -ncol(res$test_data)]),
                     s = "lambda.min", type = "response")[, 1]
  en_class <- factor(ifelse(en_prob > best_thresh, 1, 0), levels = c(0, 1))
  cm_en <- confusionMatrix(en_class, test_labels)
  auc_en <- tryCatch(auc(test_labels, en_prob), error = function(e) NA)

  # Ensemble
  ensemble_prob <- (rf_prob + gbm_prob + en_prob) / 3
  ensemble_class <- factor(ifelse(ensemble_prob > best_thresh, 1, 0), levels = c(0, 1))
  cm_ens <- confusionMatrix(ensemble_class, test_labels)
  auc_ens <- tryCatch(auc(test_labels, ensemble_prob), error = function(e) NA)

  # Collect results
  data.frame(
    run = i,

    # Random Forest
    rf_acc = cm_rf$overall["Accuracy"],
    rf_kappa = cm_rf$overall["Kappa"],
    rf_sens = cm_rf$byClass["Sensitivity"],
    rf_spec = cm_rf$byClass["Specificity"],
    rf_auc = auc_rf,

    # GBM
    gbm_acc = cm_gbm$overall["Accuracy"],
    gbm_kappa = cm_gbm$overall["Kappa"],
    gbm_sens = cm_gbm$byClass["Sensitivity"],
    gbm_spec = cm_gbm$byClass["Specificity"],
    gbm_auc = auc_gbm,

    # Elastic Net
    en_acc = cm_en$overall["Accuracy"],
    en_kappa = cm_en$overall["Kappa"],
    en_sens = cm_en$byClass["Sensitivity"],
    en_spec = cm_en$byClass["Specificity"],
    en_auc = auc_en,

    # Ensemble
    ens_acc = cm_ens$overall["Accuracy"],
    ens_kappa = cm_ens$overall["Kappa"],
    ens_sens = cm_ens$byClass["Sensitivity"],
    ens_spec = cm_ens$byClass["Specificity"],
    ens_auc = auc_ens
  )
})

results_df <- do.call(rbind, results_list)

View(results_df)




#############################################################

res$models$rf





# Split the data
training <- sample(nrow(boost_data), round(nrow(boost_data)/2))
test <- setdiff(1:nrow(boost_data), training)


# Boosted tree
#install.packages("gbm")
library(gbm)

boostedmodel <- gbm(gt ~ ., 
             data = boost_data[training,], 
             distribution = "bernoulli", 
             n.trees = 1000, 
             bag.fraction = 0.9,   # how many samples used per tree
             n.minobsinnode = 5,  # Increase this value to ensure enough samples per node
             shrinkage = 0.1, # increase this to make it faster?
             interaction.depth = 3, 
             cv.folds = 5)

head(summary(boostedmodel), n = 20) # A plot of variable importance 

# Get predicted probabilities
pred_prob <- predict(boostedmodel, newdata = boost_data[test, ], type = "response")

# Convert to class labels
predicted_classes <- ifelse(pred_prob > 0.5, 1, 0)

# Convert both to factors
predicted_classes <- factor(predicted_classes, levels = c(0, 1))
actual_labels <- factor(boost_data$gt[test], levels = c(0, 1))

# Compute confusion matrix
library(caret)
conf_matrix <- confusionMatrix(predicted_classes, actual_labels)
print(conf_matrix)


#############################################################

## Attempt 2

#install.packages("caret") 
library(caret)

library(dplyr)

# ===============================
# 1. DATA PREPROCESSING & EXPLORATION
# ===============================

# Check data balance = class distribution
table(fk_metabolom_gt_scaled$gt) # somewhat balanced

# ===============================
# 2. FEATURE SCALING/TRANSFORMATION
# ===============================

fk_scaled <- fk_metabolom_gt_scaled # or could z scale it. dont think its necessary
#fk_scaled <- scale(fk_metabolom_gt_scaled)

# ===============================
# 3. TRAIN/TEST SPLIT
# ===============================

set.seed(123)
# Use stratified sampling to maintain class balance
train_idx <- createDataPartition(fk_scaled$gt, p = 0.7, list = FALSE) # train with 70% of the data
train_data <- fk_scaled[train_idx, ]
test_data <- fk_scaled[-train_idx, ]

table(train_data$gt) # training class distrib
table(test_data$gt) # test class distrib

# ===============================
# 4. HYPERPARAMETER TUNING
# ===============================

# Grid search for optimal parameters
set.seed(123)

# Define specific mtry values to test
mtry_values <- as.integer(c(
  sqrt(ncol(train_data) - 1),          # Rule of thumb
  ncol(train_data) / 3,                # Conservative
  ncol(train_data) / 2,                # Original
  ncol(train_data) * 2 / 3,            # More features
  281,                                 # Between 225 and 338
  394,                                  # Between 338 and 451
  150                                   # an extra one
))

# or let it test lots of diff mtry values

# Run tuneRF with a limited range
tuneRF(train_data[, -ncol(train_data)], train_data$gt,
       stepFactor = 1.5,
       improve = 0.01,
       ntreeTry = 100,
       plot = TRUE)  # Limit maximum mtry value

# Remove duplicates and sort
mtry_values <- sort(unique(mtry_values))
mtry_values
# Define the tuning grid
tuneGrid <- expand.grid(mtry = mtry_values)

# Cross-validation setup
ctrl <- trainControl(
  method = "cv",
  number = 10,                    
  classProbs = TRUE,              
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)


# Parallel training setup
library(doParallel)
cores <- parallel::detectCores()
cl <- makeCluster(cores - 1)
registerDoParallel(cl)

# Train Random Forest with tuning
rf_tuned <- train(
  gt ~ ., 
  data = train_data,
  method = "rf",
  tuneGrid = tuneGrid,
  trControl = ctrl,
  metric = "ROC",
  ntree = 1000,
  importance = TRUE
)

# Stop parallel backend
stopCluster(cl)
registerDoSEQ()

# Show results
print(rf_tuned)
plot(rf_tuned)

# ===============================
# 5. MODEL 2: ADVANCED MODEL WITH BEST PARAMETERS
# ===============================

# Build final model with best parameters
best_mtry <- rf_tuned$bestTune$mtry

set.seed(123)
rf_final <- randomForest(
  gt ~ ., 
  data = train_data,
  mtry = best_mtry,
  ntree = 1000,                   # More trees for stability
  importance = TRUE,
  nodesize = 5,                   # Prevent overfitting
  maxnodes = NULL,                # Let trees grow
  replace = TRUE,
  classwt = NULL                  # Could add class weights if imbalanced
)

rf_final
# save the final model
saveRDS(rf_final, file = "../results/rf_metab_log_tuned.rds")




# ===============================
# 6. MODEL3- FEATURE SELECTION
# ===============================

# Variable importance
importance_scores <- importance(rf_final)
varImpPlot(rf_final, n.var = 10)


# Select top 10* features (see below- "how many features" code)
top_features <- head(order(importance_scores[,1], decreasing = TRUE), 10)
train_selected <- train_data[, c(top_features, which(names(train_data) == "gt"))]
test_selected <- test_data[, c(top_features, which(names(test_data) == "gt"))]

# Retrain with selected features
set.seed(123)
rf_selected <- randomForest(
  gt ~ ., 
  data = train_selected,
  mtry = floor(sqrt(ncol(train_selected)-1)),
  ntree = 1000,
  importance = TRUE
)

# ===============================
# 7. MODEL EVALUATION
# ===============================

# Evaluate all models
models <- list(
  "Original" = model,
  "Tuned" = rf_final,
  "Feature_Selected" = rf_selected
)

results <- list()

for(model_name in names(models)) {
  current_model <- models[[model_name]]
  
  # Make predictions
  if(model_name == "Original") {
    pred <- predict(current_model, fk_metabolom_gt_scaled[-training, ])
    actual <- fk_metabolom_gt_scaled$gt[-training]
  } else if(model_name == "Feature_Selected") {
    pred <- predict(current_model, test_selected)
    actual <- test_selected$gt
  } else {
    pred <- predict(current_model, test_data)
    actual <- test_data$gt
  }
  
  # Calculate metrics
  accuracy <- mean(pred == actual)
  conf_matrix <- confusionMatrix(pred, actual)
  
  results[[model_name]] <- list(
    accuracy = accuracy,
    sensitivity = conf_matrix$byClass["Sensitivity"],
    specificity = conf_matrix$byClass["Specificity"],
    conf_matrix = conf_matrix$table
  )
  
  cat("\n", model_name, "Model Results:\n")
  print(conf_matrix)
}



# ===============================
# SYSTEMATIC FEATURE SELECTION (trying selecting more/fewer features)
# ===============================

# Get variable importance from your best model
importance_scores <- importance(rf_final)
ranked_features <- order(importance_scores[,1], decreasing = TRUE)
length(ranked_features)
# Rule of thumb for sample size
n_samples <- nrow(train_data)
n_positive <- sum(train_data$gt == "GT_present")

# dataset info
cat("Total training samples:", n_samples, "\n")
cat("Positive samples:", n_positive, "\n")
cat("Total features:", ncol(train_data) - 1, "\n\n")

# Conservative feature selection rules
max_features_conservative <- min(
  n_positive * 2,           # Max 2 features per positive sample
  sqrt(ncol(train_data)-1), # Square root rule
  20                        # Hard cap for small datasets
)

cat("Recommended maximum features:", max_features_conservative, "\n\n")

# ===============================
# TEST DIFFERENT FEATURE COUNTS
# ===============================

# Test different numbers of top features
feature_counts <- c(5, 10, 15, 20, 30, 50)
# Remove counts that exceed our data
feature_counts <- feature_counts[feature_counts <= max_features_conservative]

results_df <- data.frame()

for(n_features in feature_counts) {
  cat("Testing with", n_features, "features...\n")
  
  # Select top features
  top_features <- ranked_features[1:n_features]
  train_subset <- train_data[, c(top_features, which(names(train_data) == "gt"))]
  test_subset <- test_data[, c(top_features, which(names(test_data) == "gt"))]
  
  # Train model
  set.seed(123)
  rf_subset <- randomForest(
    gt ~ ., 
    data = train_subset,
    mtry = floor(sqrt(n_features)),
    ntree = 1000,
    nodesize = 3  # Smaller for small datasets
  )
  
  # Predict and evaluate
  pred <- predict(rf_subset, test_subset)
  actual <- test_subset$gt
  
  # Calculate metrics
  conf_mat <- confusionMatrix(pred, actual)
  accuracy <- conf_mat$overall["Accuracy"]
  sensitivity <- conf_mat$byClass["Sensitivity"]
  specificity <- conf_mat$byClass["Specificity"]
  balanced_acc <- (sensitivity + specificity) / 2
  
  # Calculate F1 score
  precision <- conf_mat$byClass["Precision"]
  if(is.na(precision)) precision <- 0
  f1 <- ifelse(is.na(precision) || is.na(sensitivity), 0,
               2 * (precision * sensitivity) / (precision + sensitivity))
  
  # Store results
  results_df <- rbind(results_df, data.frame(
    n_features = n_features,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = balanced_acc,
    f1_score = f1,
    oob_error = rf_subset$err.rate[nrow(rf_subset$err.rate), "OOB"]
  ))
}

# Display results
cat("\n=== FEATURE SELECTION RESULTS ===\n")
print(round(results_df, 3))

# Find optimal number
best_f1 <- which.max(results_df$f1_score)
best_balanced <- which.max(results_df$balanced_accuracy)
best_oob <- which.min(results_df$oob_error)

cat("\nOptimal feature counts:\n")
cat("Best F1 score:", results_df$n_features[best_f1], "features (F1 =", round(results_df$f1_score[best_f1], 3), ")\n")
cat("Best balanced accuracy:", results_df$n_features[best_balanced], "features (Bal.Acc =", round(results_df$balanced_accuracy[best_balanced], 3), ")\n")
cat("Best OOB error:", results_df$n_features[best_oob], "features (OOB =", round(results_df$oob_error[best_oob], 3), ")\n")

# Plot results
par(mfrow = c(2, 2))
plot(results_df$n_features, results_df$accuracy, type = "b", main = "Accuracy vs Features", xlab = "Features", ylab = "Accuracy")
plot(results_df$n_features, results_df$balanced_accuracy, type = "b", main = "Balanced Accuracy vs Features", xlab = "Features", ylab = "Balanced Accuracy")
plot(results_df$n_features, results_df$f1_score, type = "b", main = "F1 Score vs Features", xlab = "Features", ylab = "F1 Score")
plot(results_df$n_features, results_df$oob_error, type = "b", main = "OOB Error vs Features", xlab = "Features", ylab = "OOB Error")
par(mfrow = c(1, 1))

# ===============================
# SHOW TOP METABOLITES
# ===============================

optimal_n <- results_df$n_features[best_f1]
cat("\n=== TOP", optimal_n, "METABOLITES ===\n")

top_metabolites <- ranked_features[1:optimal_n]
feature_names <- names(train_data)[top_metabolites]
importance_values <- importance_scores[top_metabolites, 1]

metabolite_ranking <- data.frame(
  Rank = 1:optimal_n,
  Metabolite = feature_names,
  Importance = round(importance_values, 4)
)

print(metabolite_ranking)








# ===============================
# FINAL MODEL WITH OPTIMAL 10 FEATURES
# ===============================

# Set n_features to 10 (instead of looping)
n_features <- 10

# Select top features
top_features <- ranked_features[1:n_features]
train_subset <- train_data[, c(top_features, which(names(train_data) == "gt"))]
test_subset <- test_data[, c(top_features, which(names(test_data) == "gt"))]

# Train model
set.seed(123)
rf_subset <- randomForest(
  gt ~ .,
  data = train_subset,
  mtry = floor(sqrt(n_features)),
  ntree = 1000,
  nodesize = 3 # Smaller for small datasets
)

# Predict and evaluate
pred <- predict(rf_subset, test_subset)
actual <- test_subset$gt

# Calculate metrics
conf_mat <- confusionMatrix(pred, actual)
accuracy <- conf_mat$overall["Accuracy"]
sensitivity <- conf_mat$byClass["Sensitivity"]
specificity <- conf_mat$byClass["Specificity"]
balanced_acc <- (sensitivity + specificity) / 2

# Calculate F1 score
precision <- conf_mat$byClass["Precision"]

f1 <- ifelse(is.na(precision) || is.na(sensitivity), 0,
             2 * (precision * sensitivity) / (precision + sensitivity))

# Print results
cat("Accuracy:", round(accuracy, 3), "\n")
cat("F1 Score:", round(f1, 3), "\n")
cat("Balanced Accuracy:", round(balanced_acc, 3), "\n")

# Save this as your final model
rf_final_10 <- rf_subset

saveRDS(rf_final_10, "../results/metab_topten_log_tuned_rf")

conf_mat

























# Variable importance
importance_scores <- importance(rf_final)
varImpPlot(rf_final, n.var = 10)

# Select top 10 features based on importance
ranked_features <- names(sort(importance_scores[,1], decreasing = TRUE))
train_selected <- train_data[, c(ranked_features[1:10], "gt")]
test_selected <- test_data[, c(ranked_features[1:10], "gt")]

# Retrain with selected features
set.seed(123)
rf_selected <- randomForest(
  gt ~ ., 
  data = train_selected,
  mtry = floor(sqrt(ncol(train_selected)-1)),
  ntree = 1000,
  importance = TRUE
)

# ===============================
# 7. MODEL EVALUATION
# ===============================

# Evaluate all models
models <- list(
  "Original" = model,
  "Tuned" = rf_final,
  "Feature_Selected" = rf_selected
)

results <- list()

for(model_name in names(models)) {
  current_model <- models[[model_name]]
  
  # Make predictions
  if(model_name == "Original") {
    pred <- predict(current_model, fk_metabolom_gt_scaled[-training, ])
    actual <- fk_metabolom_gt_scaled$gt[-training]
  } else if(model_name == "Feature_Selected") {
    pred <- predict(current_model, test_selected)
    actual <- test_selected$gt
  } else {
    pred <- predict(current_model, test_data)
    actual <- test_data$gt
  }
  
  # Calculate metrics
  accuracy <- mean(pred == actual)
  conf_matrix <- confusionMatrix(pred, actual)
  
  results[[model_name]] <- list(
    accuracy = accuracy,
    sensitivity = conf_matrix$byClass["Sensitivity"],
    specificity = conf_matrix$byClass["Specificity"],
    conf_matrix = conf_matrix$table
  )
  
  cat("\n", model_name, "Model Results:\n")
  print(conf_matrix)
}

# ===============================
# SYSTEMATIC FEATURE SELECTION (trying selecting more/fewer features)
# ===============================

# Get variable importance from your best model
importance_scores <- importance(rf_final)
ranked_features <- names(sort(importance_scores[,1], decreasing = TRUE))

# Rule of thumb for sample size
n_samples <- nrow(train_data)
n_positive <- sum(train_data$gt == "GT_present")

# dataset info
cat("Total training samples:", n_samples, "\n")
cat("Positive samples:", n_positive, "\n")
cat("Total features:", ncol(train_data) - 1, "\n\n")

# Conservative feature selection rules
max_features_conservative <- min(
  n_positive * 2,           # Max 2 features per positive sample
  sqrt(ncol(train_data)-1), # Square root rule
  20                        # Hard cap for small datasets
)

cat("Recommended maximum features:", max_features_conservative, "\n\n")

# ===============================
# TEST DIFFERENT FEATURE COUNTS
# ===============================

# Test different numbers of top features
feature_counts <- c(5, 10, 15, 20, 30, 50)

# Ensure we don't exceed available features
feature_counts <- feature_counts[feature_counts <= length(ranked_features)]

results_df <- data.frame()

for(n_features in feature_counts) {
  cat("Testing with", n_features, "features...\n")
  
  # Select top features
  top_features <- ranked_features[1:n_features]
  train_subset <- train_data[, c(top_features, "gt")]
  test_subset <- test_data[, c(top_features, "gt")]
  
  # Train model
  set.seed(123)
  rf_subset <- randomForest(
    gt ~ ., 
    data = train_subset,
    mtry = floor(sqrt(n_features)),
    ntree = 1000,
    nodesize = 3  # Smaller for small datasets
  )
  
  # Predict and evaluate
  pred <- predict(rf_subset, test_subset)
  actual <- test_subset$gt
  
  # Calculate metrics
  conf_mat <- confusionMatrix(pred, actual)
  accuracy <- conf_mat$overall["Accuracy"]
  sensitivity <- conf_mat$byClass["Sensitivity"]
  specificity <- conf_mat$byClass["Specificity"]
  balanced_acc <- (sensitivity + specificity) / 2
  
  # Calculate F1 score
  precision <- conf_mat$byClass["Precision"]
  if(is.na(precision)) precision <- 0
  f1 <- ifelse(is.na(precision) || is.na(sensitivity), 0,
               2 * (precision * sensitivity) / (precision + sensitivity))
  
  # Store results
  results_df <- rbind(results_df, data.frame(
    n_features = n_features,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = balanced_acc,
    f1_score = f1,
    oob_error = rf_subset$err.rate[nrow(rf_subset$err.rate), "OOB"]
  ))
}

# Display results
cat("\n=== FEATURE SELECTION RESULTS ===\n")
print(round(results_df, 3))

# Find optimal number of features based on F1 score, balanced accuracy, and OOB error
best_f1 <- which.max(results_df$f1_score)
best_balanced <- which.max(results_df$balanced_accuracy)
best_oob <- which.min(results_df$oob_error)

cat("\nOptimal feature counts:\n")
cat("Best F1 score:", results_df$n_features[best_f1], "features (F1 =", round(results_df$f1_score[best_f1], 3), ")\n")
cat("Best balanced accuracy:", results_df$n_features[best_balanced], "features (Bal.Acc =", round(results_df$balanced_accuracy[best_balanced], 3), ")\n")
cat("Best OOB error:", results_df$n_features[best_oob], "features (OOB =", round(results_df$oob_error[best_oob], 3), ")\n")

# Plot results
par(mfrow = c(2, 2))
plot(results_df$n_features, results_df$accuracy, type = "b", main = "Accuracy vs Features", xlab = "Features", ylab = "Accuracy")
plot(results_df$n_features, results_df$balanced_accuracy, type = "b", main = "Balanced Accuracy vs Features", xlab = "Features", ylab = "Balanced Accuracy")
plot(results_df$n_features, results_df$f1_score, type = "b", main = "F1 Score vs Features", xlab = "Features", ylab = "F1 Score")
plot(results_df$n_features, results_df$oob_error, type = "b", main = "OOB Error vs Features", xlab = "Features", ylab = "OOB Error")
par(mfrow = c(1, 1))
