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
genera_implicated <- c("Pseudomonas", "Nocardioides", "Bacillus") # think I'll try again with just p and b... or even just p fluorescens? have a look at tax table.
tax_df <- as.data.frame(tax_table(ps_highdiv_relative))
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]
ps_fs <- prune_taxa(taxa_to_keep, ps_highdiv_relative)

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
# METHOD 1: QUICK PARAMETER SCAN
# ====================================
cat("\n=== QUICK PARAMETER SCAN ===\n")

# For 3700 features, mtry should be tested around sqrt(3700) ≈ 60
param_grid <- expand.grid(
  num.trees = c(500, 1000, 2000),
  mtry = c(20, 40, 60, 100, 200, 500),  # sqrt(3700) ≈ 60, but test wide range
  min.node.size = c(1, 3, 5),
  sample.fraction = c(0.632, 0.8),
  stringsAsFactors = FALSE
)

cat("Testing", nrow(param_grid), "parameter combinations\n\n")

# Store results
results <- data.frame()

# Progress counter
pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)

for(i in 1:nrow(param_grid)) {
  set.seed(1250)
  params <- param_grid[i,]
  
  # 5-fold CV for this parameter set
  folds <- createFolds(rf_data$gt, k = 5, returnTrain = FALSE)
  fold_accuracies <- numeric(5)
  
  for(j in 1:5) {
    train_data <- rf_data[-folds[[j]], ]
    test_data <- rf_data[folds[[j]], ]
    
    rf_model <- ranger(
      gt ~ .,
      data = train_data,
      num.trees = params$num.trees,
      mtry = params$mtry,
      min.node.size = params$min.node.size,
      sample.fraction = params$sample.fraction,
      num.threads = 1  # Use 1 thread per model since we're already parallel
    )
    
    predictions <- predict(rf_model, test_data)
    fold_accuracies[j] <- mean(predictions$predictions == test_data$gt)
  }
  
  # Store results
  results <- rbind(results, data.frame(
    num.trees = params$num.trees,
    mtry = params$mtry,
    min.node.size = params$min.node.size,
    sample.fraction = params$sample.fraction,
    mean_accuracy = mean(fold_accuracies),
    sd_accuracy = sd(fold_accuracies),
    min_accuracy = min(fold_accuracies),
    max_accuracy = max(fold_accuracies)
  ))
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Sort by mean accuracy
results <- results[order(-results$mean_accuracy),]

# Display top 10
cat("\n\nTOP 10 PARAMETER COMBINATIONS:\n")
print(head(results, 10))

# Get best parameters
best_params <- results[1,]
cat("\n=== BEST PARAMETERS ===\n")
print(best_params)

# ====================================
# METHOD 2: USING CARET (More sophisticated)
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
  mtry = c(20, 40, 60, 100, 200),
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
  num.trees = 1000,  # Fixed for caret
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
# METHOD 3: RANDOM SEARCH (Faster for many parameters)
# ====================================
cat("\n\n=== METHOD 3: RANDOM SEARCH ===\n")

# Random search - test 50 random combinations
n_random <- 50
random_params <- data.frame(
  num.trees = sample(c(500, 1000, 1500, 2000, 3000), n_random, replace = TRUE),
  mtry = sample(10:500, n_random, replace = TRUE),
  min.node.size = sample(1:10, n_random, replace = TRUE),
  sample.fraction = runif(n_random, 0.5, 1.0),
  max.depth = sample(c(0, 5, 10, 15, 20, 30), n_random, replace = TRUE)
)

# Test each random combination
random_results <- foreach(i = 1:n_random, .combine = rbind) %dopar% {
  set.seed(1250 + i)
  params <- random_params[i,]
  
  # 5-fold CV
  folds <- createFolds(rf_data$gt, k = 5, returnTrain = FALSE)
  fold_accuracies <- numeric(5)
  
  for(j in 1:5) {
    train_data <- rf_data[-folds[[j]], ]
    test_data <- rf_data[folds[[j]], ]
    
    rf_model <- ranger(
      gt ~ .,
      data = train_data,
      num.trees = params$num.trees,
      mtry = params$mtry,
      min.node.size = params$min.node.size,
      sample.fraction = params$sample.fraction,
      max.depth = ifelse(params$max.depth == 0, NULL, params$max.depth)
    )
    
    predictions <- predict(rf_model, test_data)
    fold_accuracies[j] <- mean(predictions$predictions == test_data$gt)
  }
  
  c(params, mean_acc = mean(fold_accuracies), sd_acc = sd(fold_accuracies))
}

# Convert to data frame and sort
random_results <- as.data.frame(random_results)
random_results <- random_results[order(-random_results$mean_acc),]

cat("\nTop 10 from random search:\n")
print(head(random_results, 10))

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
write.csv(results, "rf_tuning_results_3genera.csv")
cat("\nResults saved to rf_tuning_results_3genera.csv\n")

# Stop parallel backend
stopCluster(cl = NULL)