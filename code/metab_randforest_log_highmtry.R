# this script runs a randomforest model on the metabolomic data

library(dplyr)

# load microbial metabolomic and gt data
fk_metabolom_gt <- readRDS("../data/fk_metabolomics_gt_noNA.RDS")

nrow(fk_metabolom_gt) # 677 compounds
#View(fk_metabolom_gt)

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
hist(fk_metabolom_gt_t[,3], breaks = 30, main = "Histogram", xlab = "Value") # check for skew- a bit righty
plot(density(log(fk_metabolom_gt_t[,3]), na.rm = TRUE), main = "Density") # does log help? think so

# log data
metab_numeric <- fk_metabolom_gt_t[, sapply(fk_metabolom_gt_t, is.numeric)] # Select only numeric metabolite columns
metab_scaled <- as.data.frame(log(metab_numeric + 1e-6)) # log transform
# THINK ABOUT Z SCALING IT HERE
non_numeric <- fk_metabolom_gt_t[, !sapply(fk_metabolom_gt_t, is.numeric)] # Optionally add back non-numeric columns
fk_metabolom_gt_scaled <- cbind(non_numeric, metab_scaled)
colnames(fk_metabolom_gt_scaled)[1]<- "gt"
View(fk_metabolom_gt_scaled)



# quick random forest
library(randomForest)
set.seed(123) # for reproducibility

training <- sample(nrow(fk_metabolom_gt_scaled), round(nrow(fk_metabolom_gt_scaled)/2)) # Pick some training data

model <- randomForest(gt ~ ., data = fk_metabolom_gt_scaled[training,], mtry = ncol(fk_metabolom_gt_scaled)/2)


# Examine the model

# accuracy
predicted <- predict(model, fk_metabolom_gt_scaled[-training, ]) # Predict classes on the test set
actual <- fk_metabolom_gt_scaled$gt[-training] # Actual classes
table(Predicted = predicted, Actual = actual) # Confusion matrix
mean(predicted == actual) # 0.54-  as random guessing worse

# model summary
print(model)     # Shows OOB error, confusion matrix, etc.
plot(model)      # Shows OOB error vs. number of trees

#############################################################

## Attempt 2

#install.packages("caret") 
library(caret)

library(dplyr)

# ===============================
# 1. DATA PREPROCESSING & EXPLORATION
# ===============================

# Check data balance = class distribution
table(fk_metabolom_gt_scaled$gt) # not balanced

fk_clean <- fk_metabolom_gt_scaled

# ===============================
# 2. FEATURE SCALING/TRANSFORMATION
# ===============================

fk_scaled <- fk_metabolom_gt_scaled
# ===============================
# 3. TRAIN/TEST SPLIT
# ===============================

set.seed(123)
# Use stratified sampling to maintain class balance
train_idx <- createDataPartition(fk_scaled$gt, p = 0.7, list = FALSE)
train_data <- fk_scaled[train_idx, ]
test_data <- fk_scaled[-train_idx, ]

table(train_data$gt) # training class distrib
table(test_data$gt) # test class distrib

# ===============================
# 4. HYPERPARAMETER TUNING
# ===============================

# Grid search for optimal parameters.
set.seed(123)
tuneGrid <- expand.grid( # mtry is the number of predictors you use to classify. Will try all these
  mtry = c(sqrt(ncol(train_data)-1),        # Rule of thumb for classification
           ncol(train_data)/3,               # Conservative
           ncol(train_data)/2,               # Your original
           ncol(train_data) * 2/3)           # More features
)

# Cross-validation setup
ctrl <- trainControl(
  method = "cv",
  number = 10,                    # 10-fold CV- model trains on 9 parts and tests on 1 part
  classProbs = TRUE,              # For ROC curves
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

## Train in parallel with cross-validation
# install.packages("doParallel")
library(doParallel)

# Detect number of cores and set up parallel backend
cores <- parallel::detectCores()
cl <- makeCluster(cores - 1)  # leave one core free
registerDoParallel(cl)

# Now run your training code
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

# Stop parallel cluster after training
stopCluster(cl)
registerDoSEQ()

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
# Test different numbers of top features WITH HIGHER MTRY
feature_counts <- c(5, 10, 15, 20)
results_df_highmtry <- data.frame()

for(n_features in feature_counts) {
  cat("Testing with", n_features, "features...\n")
  
  # Select top features
  top_features <- ranked_features[1:n_features]
  train_subset <- train_data[, c(top_features, which(names(train_data) == "gt"))]
  test_subset <- test_data[, c(top_features, which(names(test_data) == "gt"))]
  
  # Try HIGHER mtry values - test both conservative and aggressive
  mtry_values <- c(floor(sqrt(n_features)),           # Conservative (your original)
                   floor(n_features/2),               # Moderate  
                   min(n_features, floor(n_features * 0.7)))  # Aggressive
  
  for(mtry_val in mtry_values) {
    # Train model
    set.seed(123)
    rf_subset <- randomForest(
      gt ~ ., 
      data = train_subset,
      mtry = mtry_val,          # Higher mtry
      ntree = 1000,
      nodesize = 3
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
    results_df_highmtry <- rbind(results_df_highmtry, data.frame(
      n_features = n_features,
      mtry = mtry_val,
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      balanced_accuracy = balanced_acc,
      f1_score = f1,
      oob_error = rf_subset$err.rate[nrow(rf_subset$err.rate), "OOB"]
    ))
  }
}

# Display results
cat("\n=== HIGHER MTRY RESULTS ===\n")
print(round(results_df_highmtry, 3))

# Find best combination
best_combo <- which.max(results_df_highmtry$f1_score)
cat("\nBest combination:\n")
cat("Features:", results_df_highmtry$n_features[best_combo], "\n")
cat("mtry:", results_df_highmtry$mtry[best_combo], "\n") 
cat("F1 Score:", round(results_df_highmtry$f1_score[best_combo], 3), "\n")
cat("Accuracy:", round(results_df_highmtry$accuracy[best_combo], 3), "\n")


### trying again but more aggressively:
# Test much higher mtry values
feature_counts <- c(5, 10, 15, 20)
results_aggressive <- data.frame()

for(n_features in feature_counts) {
  top_features <- ranked_features[1:n_features]
  train_subset <- train_data[, c(top_features, which(names(train_data) == "gt"))]
  test_subset <- test_data[, c(top_features, which(names(test_data) == "gt"))]
  
  # Try MUCH higher mtry values
  mtry_values <- c(2,                    # Your current best (terrible)
                   floor(n_features/3),  # More reasonable
                   floor(n_features/2),  # Aggressive
                   n_features-1)         # Nearly all features
  
  for(mtry_val in mtry_values) {
    set.seed(123)
    rf_subset <- randomForest(
      gt ~ ., 
      data = train_subset,
      mtry = mtry_val,
      ntree = 1000,
      nodesize = 3
    )
    
    pred <- predict(rf_subset, test_subset)
    conf_mat <- confusionMatrix(pred, test_subset$gt)
    
    results_aggressive <- rbind(results_aggressive, data.frame(
      n_features = n_features,
      mtry = mtry_val,
      accuracy = conf_mat$overall["Accuracy"],
      f1_score = conf_mat$byClass["F1"]
    ))
    
    cat("Features:", n_features, "mtry:", mtry_val, 
        "Accuracy:", round(conf_mat$overall["Accuracy"], 3), "\n")
  }
}

print(results_aggressive)


## next code is from the old one
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














