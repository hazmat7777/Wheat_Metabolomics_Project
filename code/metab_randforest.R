# this script runs a randomforest model on the metabolomic data

library(dplyr)

# load microbial metabolomic and gt data
fk_metabolom_gt <- readRDS("../data/fk_metabolomics_gt.RDS")

nrow(fk_metabolom_gt) # 677 compounds

#data prep 
fk_metabolom_gt_t <- t(fk_metabolom_gt) # transpose
colnames(fk_metabolom_gt_t) <- fk_metabolom_gt_t[1, ] # set the first row as column names
fk_metabolom_gt_t <- fk_metabolom_gt_t[-1, ] # remove the first row (now column names)

fk_metabolom_gt_t <- as.data.frame(fk_metabolom_gt_t) # convert to data frame
colnames(fk_metabolom_gt_t)[ncol(fk_metabolom_gt_t)] <- "gt" # rename the last column to 'gt'
fk_metabolom_gt_t$gt <- trimws(as.character(fk_metabolom_gt_t$gt))
fk_metabolom_gt_t$gt <- as.factor(fk_metabolom_gt_t$gt)
levels(fk_metabolom_gt_t$gt) <- c("GT_absent", "GT_present")
levels(fk_metabolom_gt_t$gt)


unique(fk_metabolom_gt_t$gt) # check unique values in gt column

colnames(fk_metabolom_gt_t) <- make.names(colnames(fk_metabolom_gt_t), unique = TRUE)# Clean colnames by replacing spaces with underscores

#View(fk_metabolom_gt_t)

# quick random forest
library(randomForest)
set.seed(123) # for reproducibility

training <- sample(nrow(fk_metabolom_gt_t), round(nrow(fk_metabolom_gt_t)/2)) # Pick some training data

model <- randomForest(gt ~ ., data = fk_metabolom_gt_t[training,], mtry = ncol(fk_metabolom_gt_t)/2)


# Examine the model

# accuracy
predicted <- predict(model, fk_metabolom_gt_t[-training, ]) # Predict classes on the test set
actual <- fk_metabolom_gt_t$gt[-training] # Actual classes
table(Predicted = predicted, Actual = actual) # Confusion matrix
mean(predicted == actual) # 0.606- basically as good as random guessing, maybe worse

# model summary
print(model)     # Shows OOB error, confusion matrix, etc.
plot(model)      # Shows OOB error vs. number of trees

#############################################################

## Attempt 2

library(randomForest)
#install.packages("caret") 
library(caret)

library(dplyr)

# ===============================
# 1. DATA PREPROCESSING & EXPLORATION
# ===============================

# Check data balance
cat("Class distribution:\n")
table(fk_metabolom_gt_t$gt)

# Check for missing values
cat("\nMissing values:\n")
sum(is.na(fk_metabolom_gt_t))

# Check for near-zero variance predictors
nzv <- nearZeroVar(fk_metabolom_gt_t[, -which(names(fk_metabolom_gt_t) == "gt")])

if(length(nzv) > 0) {
  cat("\nRemoving", length(nzv), "near-zero variance predictors\n")
  fk_clean <- fk_metabolom_gt_t[, -nzv]
} else {
  fk_clean <- fk_metabolom_gt_t
}

# ===============================
# 2. FEATURE SCALING/TRANSFORMATION
# ===============================

# Scale numeric predictors (often helps with metabolomics data)
predictors <- fk_clean[, !names(fk_clean) %in% "gt"]
predictors[] <- lapply(predictors, function(x) as.numeric(as.character(x)))
predictors_scaled <- scale(predictors)
fk_scaled <- data.frame(predictors_scaled, gt = fk_clean$gt)
View(fk_scaled)

# ===============================
# 3. BETTER TRAIN/TEST SPLIT
# ===============================

set.seed(123)
# Use stratified sampling to maintain class balance
train_idx <- createDataPartition(fk_scaled$gt, p = 0.7, list = FALSE)
train_data <- fk_scaled[train_idx, ]
test_data <- fk_scaled[-train_idx, ]

cat("Training set class distribution:\n")
table(train_data$gt)
cat("Test set class distribution:\n")
table(test_data$gt)

# ===============================
# 4. HYPERPARAMETER TUNING
# ===============================

# Grid search for optimal parameters
set.seed(123)
tuneGrid <- expand.grid(
  mtry = c(sqrt(ncol(train_data)-1),        # Rule of thumb for classification
           ncol(train_data)/3,               # Conservative
           ncol(train_data)/2,               # Your original
           ncol(train_data) * 2/3)           # More features
)

# Cross-validation setup
ctrl <- trainControl(
  method = "cv",
  number = 10,                    # 10-fold CV
  classProbs = TRUE,              # For ROC curves
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# rename the target variable to a valid name (needs a letter at the start)
levels(train_data$gt) <- c("GT_absent", "GT_present")
levels(train_data$gt)


# Train with cross-validation
rf_tuned <- train(
  gt ~ ., 
  data = train_data,
  method = "rf",
  tuneGrid = tuneGrid,
  trControl = ctrl,
  metric = "ROC",                 # Optimize for AUC
  ntree = 1000,                   # More trees
  importance = TRUE
)

print(rf_tuned)
plot(rf_tuned)

# ===============================
# 5. ADVANCED MODEL WITH BEST PARAMETERS
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

# ===============================
# 6. FEATURE SELECTION
# ===============================

# Variable importance
importance_scores <- importance(rf_final)
varImpPlot(rf_final, n.var = 20)


# Select top features (try with top 20, 50, etc.)
top_features <- head(order(importance_scores[,1], decreasing = TRUE), 50)
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
    pred <- predict(current_model, fk_metabolom_gt_t[-training, ])
    actual <- fk_metabolom_gt_t$gt[-training]
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
  cat("Accuracy:", round(accuracy, 3), "\n")
  cat("Sensitivity:", round(conf_matrix$byClass["Sensitivity"], 3), "\n")
  cat("Specificity:", round(conf_matrix$byClass["Specificity"], 3), "\n")
  print(conf_matrix$table)
}


# ===============================
# SYSTEMATIC FEATURE SELECTION (trying selecting more/fewer features)
# ===============================

# Get variable importance from your best model
importance_scores <- importance(rf_final)
ranked_features <- order(importance_scores[,1], decreasing = TRUE)

# Rule of thumb for sample size
n_samples <- nrow(train_data)
n_positive <- sum(train_data$gt == "GT_present")

cat("Dataset info:\n")
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
# FINAL MODEL WITH 5 FEATURES
# ===============================

# Variable importance
importance_scores <- importance(rf_final)
varImpPlot(rf_final, n.var = 20)


# Select top features (try with top 20, 50, etc.)
top_features <- head(order(importance_scores[,1], decreasing = TRUE), 5)
train_selected <- train_data[, c(top_features, which(names(train_data) == "gt"))]
test_selected <- test_data[, c(top_features, which(names(test_data) == "gt"))]

# Retrain with selected features
rf_selected <- randomForest(
  gt ~ ., 
  data = train_selected,
  mtry = floor(sqrt(ncol(train_selected)-1)),
  ntree = 1000,
  importance = TRUE
)

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
    pred <- predict(current_model, fk_metabolom_gt_t[-training, ])
    actual <- fk_metabolom_gt_t$gt[-training]
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
  cat("Accuracy:", round(accuracy, 3), "\n")
  cat("Sensitivity:", round(conf_matrix$byClass["Sensitivity"], 3), "\n")
  cat("Specificity:", round(conf_matrix$byClass["Specificity"], 3), "\n")
  print(conf_matrix$table)
}