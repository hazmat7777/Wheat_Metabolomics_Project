# this script runs a randomforest model on the microobial community data (vs gt baiting data)
# the real problem is the small sample size- 6 samples in test set.



library(dplyr)

# load microbial metabarcoding data
ps_16S_highdiv_absolute <- readRDS("../data/ps_16S_highdiv_absolute.rds")

sample_names(ps_16S_highdiv_absolute)

View(sample_data(ps_16S_highdiv_absolute))


ps_16S_rf <- prune_samples(!is.na(sample_data(ps_16S_highdiv_absolute)$gt), ps_16S_highdiv_absolute) # remove samples without GT data
ps_16S_rf <- prune_taxa(taxa_sums(ps_16S_rf) > 0, ps_16S_rf) #removes missing taxa

ps_16S_rf

#===================================================
# THE PLAN
#===================================================

# Prelim- do one by ESVs

# then do dimensionality reduction:
    # First 5 rfs- do a rf aggregated by each taxonomic level (phylum, class, order, genus)
    # 6th one- pcoa on jenson-shannon distance between samples (???)


# Load required libraries
library(phyloseq)
library(randomForest)
library(caret)
library(dplyr)

# Extract response variable from sample_data
sample_df <- as.data.frame(sample_data(ps_16S_rf))
gt <- sample_df$gt  # Extract the gt column
names(gt) <- rownames(sample_df)  # Assign sample names

# Convert to factor for classification
gt <- as.factor(gt)

# Check dimensions and levels
cat("Sample names in phyloseq:", length(sample_names(ps_16S_rf)), "\n")
cat("Response variable length:", length(gt), "\n")
cat("Response variable levels:", levels(gt), "\n")
cat("Response variable distribution:", table(gt), "\n\n")

# ===== TAXONOMIC AGGREGATION =====

# Function to aggregate OTUs at different taxonomic ranks
aggregate_taxa <- function(ps_obj, tax_rank) {
  # Aggregate taxa at specified rank
  ps_agg <- tax_glom(ps_obj, taxrank = tax_rank)
  
  # Get OTU table (samples as rows, taxa as columns)
  otu_table_agg <- as.data.frame(t(otu_table(ps_agg)))
  
  # Optional: Add taxonomic names as column names
  if (!is.null(tax_table(ps_agg))) {
    tax_names <- as.data.frame(tax_table(ps_agg))[, tax_rank]
    # Handle NAs and create unique names
    tax_names[is.na(tax_names)] <- paste0("Unknown_", tax_rank, "_", seq_len(sum(is.na(tax_names))))
    tax_names <- make.unique(as.character(tax_names))
    
    # Clean up names for R compatibility
    tax_names <- make.names(tax_names, unique = TRUE)
    
    colnames(otu_table_agg) <- tax_names
  } else {
    # If no taxonomy table, create generic names
    colnames(otu_table_agg) <- make.names(paste0("Taxa_", seq_len(ncol(otu_table_agg))), unique = TRUE)
  }
  
  return(otu_table_agg)
}

# Aggregate at different taxonomic ranks
taxonomic_ranks <- c("Genus", "Family", "Order", "Class", "Phylum")
aggregated_data <- list()

for (rank in taxonomic_ranks) {
  cat("Aggregating at", rank, "level...\n")
  aggregated_data[[rank]] <- aggregate_taxa(ps_16S_rf, rank)
  cat("Number of", rank, "features:", ncol(aggregated_data[[rank]]), "\n\n")
}

# ===== IMPROVED RANDOM FOREST WITH TUNING AND STRATIFIED SAMPLING =====

# Function to run improved random forest with tuning and stratified sampling
run_improved_rf <- function(feature_data, response_var, approach_name, 
                           train_prop = 0.7, cv_folds = 5, tune_length = 5) {
  
  cat("Running Improved Random Forest for", approach_name, "...\n")
  
  # Ensure sample names match
  feature_samples <- rownames(feature_data)
  response_samples <- names(response_var)
  
  if (is.null(response_samples)) {
    response_samples <- sample_names(ps_16S_rf)
    names(response_var) <- response_samples
  }
  
  # Find common samples
  common_samples <- intersect(feature_samples, response_samples)
  cat("Common samples:", length(common_samples), "\n")
  
  if (length(common_samples) == 0) {
    stop("No matching samples between feature data and response variable")
  }
  
  # Subset both to common samples
  feature_data_matched <- feature_data[common_samples, , drop = FALSE]
  response_var_matched <- response_var[common_samples]
  
  # Clean column names
  colnames(feature_data_matched) <- make.names(colnames(feature_data_matched), unique = TRUE)
  
  # Remove zero variance features
  zero_var_cols <- which(apply(feature_data_matched, 2, var, na.rm = TRUE) == 0)
  if (length(zero_var_cols) > 0) {
    feature_data_matched <- feature_data_matched[, -zero_var_cols]
    cat("Removed", length(zero_var_cols), "zero variance features\n")
  }
  
  # Create final dataset
  rf_data <- data.frame(feature_data_matched, response = as.factor(response_var_matched))
  
  # Remove rows with missing response
  complete_rows <- complete.cases(rf_data$response)
  rf_data <- rf_data[complete_rows, ]
  
  cat("Final dataset dimensions:", dim(rf_data), "\n")
  cat("Response distribution:", table(rf_data$response), "\n")
  
  # ===== STRATIFIED TRAIN/TEST SPLIT =====
  set.seed(123)  # For reproducibility
  
  # Create stratified split (automatically stratified by response variable)
  train_indices <- createDataPartition(rf_data$response, 
                                      p = train_prop, 
                                      list = FALSE)
  
  train_data <- rf_data[train_indices, ]
  test_data <- rf_data[-train_indices, ]
  
  cat("Training set size:", nrow(train_data), "\n")
  cat("Training set distribution:", table(train_data$response), "\n")
  cat("Test set size:", nrow(test_data), "\n")
  cat("Test set distribution:", table(test_data$response), "\n")
  
  # ===== HYPERPARAMETER TUNING =====
  
  # Define tuning grid for mtry
  max_features <- ncol(train_data) - 1
  mtry_values <- floor(seq(sqrt(max_features), max_features/2, length.out = tune_length))
  mtry_values <- unique(pmax(1, mtry_values))  # Ensure at least 1
  
  cat("Tuning mtry values:", mtry_values, "\n")
  
  # Set up cross-validation
  ctrl <- trainControl(
    method = "cv",
    number = cv_folds,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    sampling = "up"  # Handle class imbalance if present
  )
  
  # Tune random forest
  cat("Performing hyperparameter tuning...\n")
  rf_tune <- train(
    response ~ ., 
    data = train_data,
    method = "rf",
    trControl = ctrl,
    tuneGrid = data.frame(mtry = mtry_values),
    ntree = 500,
    importance = TRUE,
    metric = "ROC"
  )
  
  # Best parameters
  best_mtry <- rf_tune$bestTune$mtry
  cat("Best mtry:", best_mtry, "\n")
  
  # ===== FINAL MODEL TRAINING =====
  
  # Train final model with best parameters
  final_rf <- randomForest(
    response ~ ., 
    data = train_data,
    mtry = best_mtry,
    ntree = 500,
    importance = TRUE
  )
  
  # ===== MODEL EVALUATION =====
  
  # Predictions on test set
  test_pred <- predict(final_rf, test_data)
  test_prob <- predict(final_rf, test_data, type = "prob")
  
  # Confusion matrix
  cm <- confusionMatrix(test_pred, test_data$response)
  
  # Results summary
  cat("\n=== RESULTS SUMMARY ===\n")
  cat("Best mtry:", best_mtry, "\n")
  cat("Training OOB Error:", round(final_rf$err.rate[nrow(final_rf$err.rate), "OOB"] * 100, 2), "%\n")
  cat("Test Accuracy:", round(cm$overall["Accuracy"] * 100, 2), "%\n")
  cat("Test Sensitivity:", round(cm$byClass["Sensitivity"] * 100, 2), "%\n")
  cat("Test Specificity:", round(cm$byClass["Specificity"] * 100, 2), "%\n")
  
  cat("\nTest Set Confusion Matrix:\n")
  print(cm$table)
  
  # Variable importance
  cat("\nTop 10 Important Features:\n")
  importance_scores <- importance(final_rf)
  top_features <- head(
    importance_scores[order(importance_scores[, "MeanDecreaseAccuracy"], decreasing = TRUE), ],
    10
  )
  print(round(top_features, 4))
  
  cat("\n" , rep("=", 50), "\n\n")
  
  # Return results
  results <- list(
    model = final_rf,
    tuning_results = rf_tune,
    best_mtry = best_mtry,
    confusion_matrix = cm,
    train_data = train_data,
    test_data = test_data,
    predictions = test_pred,
    probabilities = test_prob,
    importance = importance_scores
  )
  
  return(results)
}

# ===== RUN IMPROVED RANDOM FORESTS =====

rf_results <- list()

# Run for each taxonomic rank
for (rank in taxonomic_ranks) {
  rf_results[[rank]] <- run_improved_rf(
    aggregated_data[[rank]], 
    gt,  
    paste("Taxonomic aggregation -", rank),
    train_prop = 0.7,  # 70% for training
    cv_folds = 5,      # 5-fold CV for tuning
    tune_length = 5    # Test 5 different mtry values
  )
}

# ===== COMPARE MODELS =====

# Extract performance metrics
performance_summary <- data.frame(
  Approach = names(rf_results),
  Best_mtry = sapply(rf_results, function(x) x$best_mtry),
  Train_OOB_Error = sapply(rf_results, function(x) 
    round(x$model$err.rate[nrow(x$model$err.rate), "OOB"] * 100, 2)),
  Test_Accuracy = sapply(rf_results, function(x) 
    round(x$confusion_matrix$overall["Accuracy"] * 100, 2)),
  Test_Sensitivity = sapply(rf_results, function(x) 
    round(x$confusion_matrix$byClass["Sensitivity"] * 100, 2)),
  Test_Specificity = sapply(rf_results, function(x) 
    round(x$confusion_matrix$byClass["Specificity"] * 100, 2)),
  stringsAsFactors = FALSE
)

cat("=== MODEL COMPARISON ===\n")
print(performance_summary)

# Plot comparison
par(mfrow = c(2, 2))

# Test Accuracy
barplot(performance_summary$Test_Accuracy, 
        names.arg = performance_summary$Approach,
        main = "Test Accuracy by Taxonomic Level",
        ylab = "Accuracy (%)",
        las = 2, ylim = c(0, 100))

# Sensitivity (True Positive Rate)
barplot(performance_summary$Test_Sensitivity, 
        names.arg = performance_summary$Approach,
        main = "Test Sensitivity by Taxonomic Level",
        ylab = "Sensitivity (%)",
        las = 2, ylim = c(0, 100))

# Specificity (True Negative Rate)
barplot(performance_summary$Test_Specificity, 
        names.arg = performance_summary$Approach,
        main = "Test Specificity by Taxonomic Level",
        ylab = "Specificity (%)",
        las = 2, ylim = c(0, 100))

# Best mtry values
barplot(performance_summary$Best_mtry, 
        names.arg = performance_summary$Approach,
        main = "Optimal mtry by Taxonomic Level",
        ylab = "mtry",
        las = 2)

par(mfrow = c(1, 1))