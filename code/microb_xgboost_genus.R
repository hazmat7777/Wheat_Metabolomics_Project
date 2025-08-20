# Save this as a new clean file
library(tidyverse)
library(phyloseq)
#install.packages("xgboost")
library(xgboost)
library(caret)

ps_highdiv_relative <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")
ps_genus <- tax_glom(ps_highdiv_relative, "Genus", NArm = TRUE)
taxa_names(ps_genus) <- make.unique(as.character(tax_table(ps_genus)[, "Genus"]))

cat("ESVs after tax_glom:", ntaxa(ps_genus), "\n")

# Check the taxonomy table after aggregation
tax_after <- data.frame(tax_table(ps_genus))
genus_names_after <- as.character(tax_after$Genus)
cat("Unique genera after tax_glom:", length(unique(genus_names_after)), "\n")
cat("Duplicated genera:", sum(duplicated(genus_names_after)), "\n")

# Show some examples of what's duplicated
dups <- genus_names_after[duplicated(genus_names_after)]
cat("Example duplicated genera:", head(unique(dups), 5), "\n") # only incertae sedis is duplicated.

# Prepare data
sample_df <- data.frame(sample_data(ps_genus))
otu_df <- data.frame(t(as.matrix(otu_table(ps_genus))))
xgb_data <- cbind(sample_df, otu_df)
xgb_data <- xgb_data[!is.na(xgb_data$gt), ]

# Clean columns
xgb_data$sample_name <- NULL
xgb_data$og_sample_names <- NULL
if("tillage" %in% names(xgb_data)) xgb_data$tillage <- as.factor(xgb_data$tillage)

# Convert target variable to numeric for XGBoost (0-indexed)
xgb_data$gt_numeric <- as.numeric(as.factor(xgb_data$gt)) - 1
original_levels <- levels(as.factor(xgb_data$gt))
cat("Target variable levels:", paste(original_levels, collapse = ", "), "\n")
cat("Numeric encoding:", paste(0:(length(original_levels)-1), collapse = ", "), "\n")

# Cross-validation
set.seed(12550)
folds <- createFolds(xgb_data$gt, k = 5, returnTrain = FALSE)
accuracies <- c()

for(i in 1:5) {
  train_data <- xgb_data[-folds[[i]], ]
  test_data <- xgb_data[folds[[i]], ]
  
  # Prepare data for XGBoost
  # Remove the original factor gt column, keep only gt_numeric
  train_features <- train_data %>% select(-gt, -gt_numeric)
  test_features <- test_data %>% select(-gt, -gt_numeric)
  
  # Convert factor variables to numeric if any
  factor_cols <- sapply(train_features, is.factor)
  if(any(factor_cols)) {
    train_features[factor_cols] <- lapply(train_features[factor_cols], as.numeric)
    test_features[factor_cols] <- lapply(test_features[factor_cols], as.numeric)
  }
  
  # Create DMatrix objects
  dtrain <- xgb.DMatrix(data = as.matrix(train_features), 
                        label = train_data$gt_numeric)
  dtest <- xgb.DMatrix(data = as.matrix(test_features), 
                       label = test_data$gt_numeric)
  
  # Set parameters for multiclass classification
  num_classes <- length(original_levels)
  
  if(num_classes > 2) {
    # Multiclass classification
    params <- list(
      objective = "multi:softprob",
      eval_metric = "mlogloss",
      num_class = num_classes,
      eta = 0.1,              # learning rate
      max_depth = 6,          # maximum depth of trees
      subsample = 0.8,        # subsample ratio
      colsample_bytree = 0.8, # subsample ratio of columns
      min_child_weight = 1
    )
  } else {
    # Binary classification
    params <- list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      eta = 0.1,
      max_depth = 6,
      subsample = 0.8,
      colsample_bytree = 0.8,
      min_child_weight = 1
    )
  }
  
  # Train XGBoost model
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 100,          # number of boosting rounds (equivalent to num.trees in ranger)
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0,            # suppress output
    early_stopping_rounds = 10
  )
  
  # Make predictions
  if(num_classes > 2) {
    # For multiclass, predict returns probabilities for each class
    pred_probs <- predict(xgb_model, dtest, reshape = TRUE)
    predictions <- apply(pred_probs, 1, which.max) - 1  # Convert back to 0-indexed
  } else {
    # For binary, threshold at 0.5
    pred_probs <- predict(xgb_model, dtest)
    predictions <- ifelse(pred_probs > 0.5, 1, 0)
  }
  
  # Calculate accuracy
  accuracy <- mean(predictions == test_data$gt_numeric)
  accuracies[i] <- accuracy
  cat("Fold", i, "accuracy:", round(accuracy, 3), "\n")
}

cat("Mean accuracy:", round(mean(accuracies), 3), "\n")

# Optional: Feature importance analysis
cat("\nTraining final model for feature importance...\n")
final_features <- xgb_data %>% select(-gt, -gt_numeric)
factor_cols <- sapply(final_features, is.factor)
if(any(factor_cols)) {
  final_features[factor_cols] <- lapply(final_features[factor_cols], as.numeric)
}

final_dtrain <- xgb.DMatrix(data = as.matrix(final_features), 
                           label = xgb_data$gt_numeric)

final_model <- xgb.train(
  params = params,
  data = final_dtrain,
  nrounds = 100,
  verbose = 0
)

# Get feature importance
importance_matrix <- xgb.importance(colnames(final_features), model = final_model)
cat("Top 10 most important features:\n")
print(head(importance_matrix, 10))