# Load required libraries
library(ranger)
library(caret)
library(doParallel)

# Set up parallel processing
n_cores <- detectCores() - 1
registerDoParallel(cores = n_cores)
cat("Using", n_cores, "cores for parallel processing\n")

# ====================================
# PREPARE DATA (your existing code)
# ====================================
genera_implicated <- c("Pseudomonas", "Bacillus") 

tax_df <- as.data.frame(tax_table(ps_highdiv_relative))
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]
ps_fs <- prune_taxa(taxa_to_keep, ps_highdiv_relative)

# View(tax_table(ps_fs))
# tax_df <- as.data.frame(tax_table(ps_fs))
# sort(unique(tax_df$Species))
# ps_fs

# Prepare data
sample_df <- data.frame(sample_data(ps_fs))
otu_df <- data.frame(t(as.matrix(otu_table(ps_fs))))
rf_data <- cbind(sample_df, otu_df)
rf_data <- rf_data[!is.na(rf_data$gt), ]
rf_data$gt <- as.factor(rf_data$gt)

# Clean columns
rf_data$sample_name <- NULL
rf_data$og_sample_names <- NULL
rf_data$tillage <- NULL

# Check dimensions
n_features <- ncol(rf_data) - 1  # Should be ~3700
n_samples <- nrow(rf_data)        # 65
cat("\nDataset:", n_samples, "samples x", n_features, "features\n")
cat("Target distribution:\n")
print(table(rf_data$gt))

# ====================================
# USING CARET (sophisticated)
# ====================================
cat("\n\n=== METHOD 2: CARET TUNING ===\n")

# Create training control
train_control <- trainControl(
  method = "cv",
  number = 5,
  search = "grid",
  verboseIter = TRUE,
  allowParallel = TRUE
)

# Define tuning grid for caret
# Note: caret uses different parameter names
tune_grid <- expand.grid(
  mtry = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
  splitrule = c("gini", "extratrees"),
  min.node.size = c(1, 3, 5)
)

# Train model with tuning
cat("Starting caret grid search...\n")
model_tuned <- train(
  gt ~ .,
  data = rf_data,
  method = "ranger",
  trControl = train_control,
  tuneGrid = tune_grid,
  num.trees = 2000,  # Fixed for caret
  importance = 'impurity',
  num.threads = n_cores
)

# Display results
print(model_tuned)
cat("\nBest parameters from caret:\n")
print(model_tuned$bestTune)

# Plot results if you want
plot(model_tuned)



# ====================================
# FINAL EVALUATION WITH BEST PARAMETERS
# ====================================
cat("\n\n=== FINAL EVALUATION WITH BEST PARAMETERS ===\n")

# Use best parameters from any method
best_mtry <- best_params$mtry
best_trees <- best_params$num.trees
best_nodesize <- best_params$min.node.size
best_fraction <- best_params$sample.fraction

cat("Using: trees =", best_trees, ", mtry =", best_mtry, 
    ", min.node.size =", best_nodesize, ", sample.fraction =", best_fraction, "\n\n")

# Final 5-fold CV with best parameters
set.seed(1250)
folds <- createFolds(rf_data$gt, k = 5, returnTrain = FALSE)
final_accuracies <- c()
importance_matrix <- matrix(0, nrow = n_features, ncol = 5)

for(i in 1:5) {
  train_data <- rf_data[-folds[[i]], ]
  test_data <- rf_data[folds[[i]], ]
  
  rf_model <- ranger(
    gt ~ ., 
    data = train_data, 
    num.trees = best_trees, 
    mtry = best_mtry,
    min.node.size = best_nodesize,
    sample.fraction = best_fraction,
    importance = 'impurity'
  )
  
  predictions <- predict(rf_model, test_data)
  accuracy <- mean(predictions$predictions == test_data$gt)
  final_accuracies[i] <- accuracy
  importance_matrix[,i] <- rf_model$variable.importance
  
  cat("Fold", i, "accuracy:", round(accuracy, 3), "\n")
}

cat("\nFinal Mean Accuracy:", round(mean(final_accuracies), 3), "±", round(sd(final_accuracies), 3), "\n")
cat("Range: [", round(min(final_accuracies), 3), ",", round(max(final_accuracies), 3), "]\n")

# Average importance across folds
mean_importance <- rowMeans(importance_matrix)
names(mean_importance) <- names(rf_data)[-which(names(rf_data) == "gt")]

# Top important features
top_features <- sort(mean_importance, decreasing = TRUE)[1:20]
cat("\nTop 20 Important Features:\n")
print(round(top_features, 4))

# Save results
write.csv(results, "rf_tuning_results_2genera.csv")
cat("\nResults saved to rf_tuning_results_2genera.csv\n")

# Stop parallel backend
stopCluster(cl = NULL)