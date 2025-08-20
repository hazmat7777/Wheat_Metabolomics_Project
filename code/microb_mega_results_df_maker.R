# Load required libraries
library(ranger)
library(caret)
library(doParallel)
library(phyloseq)

# Set up parallel processing
n_cores <- detectCores() - 1
registerDoParallel(cores = n_cores)
cat("Using", n_cores, "cores for parallel processing\n")

# ====================================
# PREPARE FULL DATASET
# ====================================
cat("\n=== PREPARING FULL DATASET ===\n")

# Prepare full data
sample_df_full <- data.frame(sample_data(ps_highdiv_relative))
otu_df_full <- data.frame(t(as.matrix(otu_table(ps_highdiv_relative))))
rf_data_full <- cbind(sample_df_full, otu_df_full)
rf_data_full <- rf_data_full[!is.na(rf_data_full$gt), ]
rf_data_full$gt <- as.factor(rf_data_full$gt)

# Clean columns
rf_data_full$sample_name <- NULL
rf_data_full$og_sample_names <- NULL
rf_data_full$tillage <- NULL

cat("Full dataset:", nrow(rf_data_full), "samples x", ncol(rf_data_full)-1, "features\n")

# ====================================
# PREPARE FEATURE-SELECTED DATASET
# ====================================
cat("\n=== PREPARING FEATURE-SELECTED DATASET (2 GENERA) ===\n")

genera_implicated <- c("Pseudomonas", "Bacillus") 

tax_df <- as.data.frame(tax_table(ps_highdiv_relative))
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]
ps_fs <- prune_taxa(taxa_to_keep, ps_highdiv_relative)

# Prepare feature-selected data
sample_df_fs <- data.frame(sample_data(ps_fs))
otu_df_fs <- data.frame(t(as.matrix(otu_table(ps_fs))))
rf_data_fs <- cbind(sample_df_fs, otu_df_fs)
rf_data_fs <- rf_data_fs[!is.na(rf_data_fs$gt), ]
rf_data_fs$gt <- as.factor(rf_data_fs$gt)

# Clean columns
rf_data_fs$sample_name <- NULL
rf_data_fs$og_sample_names <- NULL
rf_data_fs$tillage <- NULL

cat("Feature-selected dataset:", nrow(rf_data_fs), "samples x", ncol(rf_data_fs)-1, "features\n")
cat("Target distribution:\n")
print(table(rf_data_fs$gt))

# ====================================
# INITIALIZE RESULTS DATAFRAME
# ====================================
results_df <- data.frame(
  Model = character(),
  Dataset = character(),
  N_Features = integer(),
  Method = character(),
  num.trees = integer(),
  mtry = integer(),
  min.node.size = integer(),
  sample.fraction = numeric(),
  splitrule = character(),
  Mean_Accuracy = numeric(),
  SD_Accuracy = numeric(),
  Min_Accuracy = numeric(),
  Max_Accuracy = numeric(),
  stringsAsFactors = FALSE
)




















#################### this is the bit that made the code break

# ====================================
# FUNCTION TO RUN TUNING AND EVALUATION
# ====================================
run_rf_tuning <- function(rf_data, dataset_name, model_name) {
  
  cat("\n\n=== RUNNING", model_name, "-", dataset_name, "===\n")
  cat("Dataset dimensions:", nrow(rf_data), "x", ncol(rf_data)-1, "\n\n")
  
  n_features <- ncol(rf_data) - 1
  
  # Create training control
  train_control <- trainControl(
    method = "cv",
    number = 5,
    search = "grid",
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  
  # Define tuning grid for caret
  tune_grid <- expand.grid(
    mtry = c(5, 6, 7, 8, 9, 12),
    splitrule = c("gini", "extratrees"),
    min.node.size = c(1, 3, 5)
  )
  
  # Train model with tuning
  cat("Starting caret grid search for", dataset_name, "...\n")
  model_tuned <- train(
    gt ~ .,
    data = rf_data,
    method = "ranger",
    trControl = train_control,
    tuneGrid = tune_grid,
    num.trees = 2000,
    importance = 'impurity',
    num.threads = n_cores
  )
  
  # Display results
  cat("\nBest parameters for", dataset_name, ":\n")
  print(model_tuned$bestTune)
  cat("Best CV Accuracy:", round(max(model_tuned$results$Accuracy), 3), "\n")
  
  # Extract best parameters
  best_params <- model_tuned$bestTune
  best_mtry <- best_params$mtry
  best_nodesize <- best_params$min.node.size
  best_splitrule <- as.character(best_params$splitrule)
  
  # Final evaluation with best parameters
  cat("\nFinal evaluation with best parameters...\n")
  set.seed(1250)
  folds <- createFolds(rf_data$gt, k = 5, returnTrain = FALSE)
  final_accuracies <- c()
  
  for(i in 1:5) {
    train_data <- rf_data[-folds[[i]], ]
    test_data <- rf_data[folds[[i]], ]
    
    rf_model <- ranger(
      gt ~ ., 
      data = train_data, 
      num.trees = 2000,
      mtry = best_mtry,
      min.node.size = best_nodesize,
      splitrule = best_splitrule,
      importance = 'impurity'
    )
    
    predictions <- predict(rf_model, test_data)
    accuracy <- mean(predictions$predictions == test_data$gt)
    final_accuracies[i] <- accuracy
    
    cat("  Fold", i, "accuracy:", round(accuracy, 3), "\n")
  }
  
  # Store results
  result_row <- data.frame(
    Model = model_name,
    Dataset = dataset_name,
    N_Features = n_features,
    Method = "Caret_GridSearch",
    num.trees = 2000,
    mtry = best_mtry,
    min.node.size = best_nodesize,
    sample.fraction = 0.632,  # ranger default
    splitrule = best_splitrule,
    Mean_Accuracy = round(mean(final_accuracies), 3),
    SD_Accuracy = round(sd(final_accuracies), 3),
    Min_Accuracy = round(min(final_accuracies), 3),
    Max_Accuracy = round(max(final_accuracies), 3),
    stringsAsFactors = FALSE
  )
  
  cat("\nMean Accuracy:", result_row$Mean_Accuracy, "±", result_row$SD_Accuracy, "\n")
  
  return(list(
    result_row = result_row,
    model = model_tuned,
    accuracies = final_accuracies
  ))
}










########## code got to here before breaking

# ====================================
# RUN EXPERIMENTS
# ====================================

# 1. Run on FULL dataset
cat("\n", paste(rep("=", 60), collapse=""), "\n", sep="")
cat("EXPERIMENT 1: FULL DATASET (ALL FEATURES)\n")
cat(paste(rep("=", 60), collapse=""), "\n", sep="")

full_results <- run_rf_tuning(
  rf_data = rf_data_full,
  dataset_name = "Full_Dataset",
  model_name = "RandomForest_AllFeatures"
)
results_df <- rbind(results_df, full_results$result_row)

# 2. Run on FEATURE-SELECTED dataset
cat("\n", paste(rep("=", 60), collapse=""), "\n", sep="")
cat("EXPERIMENT 2: FEATURE-SELECTED (2 GENERA)\n")
cat(paste(rep("=", 60), collapse=""), "\n", sep="")

fs_results <- run_rf_tuning(
  rf_data = rf_data_fs,
  dataset_name = "Feature_Selected_2Genera",
  model_name = "RandomForest_Pseudomonas_Bacillus"
)
results_df <- rbind(results_df, fs_results$result_row)

# ====================================
# OPTIONAL: Add more feature selection methods
# ====================================
# You can add more feature selection approaches here
# For example, top variance features, correlation-based selection, etc.

# Example: Top 100 variance features
if(FALSE) {  # Set to TRUE to include this
  cat("\n", paste(rep("=", 60), collapse=""), "\n", sep="")
  cat("EXPERIMENT 3: TOP 100 VARIANCE FEATURES\n")
  cat(paste(rep("=", 60), collapse=""), "\n", sep="")
  
  # Calculate variance for each feature
  feature_vars <- apply(rf_data_full[, -which(names(rf_data_full) == "gt")], 2, var)
  top_100_features <- names(sort(feature_vars, decreasing = TRUE)[1:100])
  rf_data_topvar <- rf_data_full[, c("gt", top_100_features)]
  
  topvar_results <- run_rf_tuning(
    rf_data = rf_data_topvar,
    dataset_name = "Top100_Variance",
    model_name = "RandomForest_TopVariance"
  )
  results_df <- rbind(results_df, topvar_results$result_row)
}

# ====================================
# DISPLAY FINAL COMPARISON
# ====================================
cat("\n\n", paste(rep("=", 60), collapse=""), "\n", sep="")
cat("FINAL COMPARISON OF ALL MODELS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n", sep="")

# Sort by mean accuracy
results_df <- results_df[order(results_df$Mean_Accuracy, decreasing = TRUE), ]
print(results_df)

# Create summary comparison
cat("\n=== PERFORMANCE SUMMARY ===\n")
for(i in 1:nrow(results_df)) {
  cat(sprintf("%-30s: %d features, Accuracy = %.3f ± %.3f\n",
              results_df$Dataset[i],
              results_df$N_Features[i],
              results_df$Mean_Accuracy[i],
              results_df$SD_Accuracy[i]))
}

# Identify best model
best_idx <- which.max(results_df$Mean_Accuracy)
cat("\n=== BEST PERFORMING MODEL ===\n")
cat("Dataset:", results_df$Dataset[best_idx], "\n")
cat("Features:", results_df$N_Features[best_idx], "\n")
cat("Accuracy:", results_df$Mean_Accuracy[best_idx], "±", results_df$SD_Accuracy[best_idx], "\n")
cat("Parameters: mtry =", results_df$mtry[best_idx], 
    ", min.node.size =", results_df$min.node.size[best_idx],
    ", splitrule =", results_df$splitrule[best_idx], "\n")

# Calculate performance difference
if(nrow(results_df) >= 2) {
  diff_accuracy <- results_df$Mean_Accuracy[1] - results_df$Mean_Accuracy[2]
  cat("\n=== FEATURE SELECTION IMPACT ===\n")
  cat(sprintf("Difference in accuracy: %.3f\n", abs(diff_accuracy)))
  
  if(which(results_df$Dataset == "Feature_Selected_2Genera") < 
     which(results_df$Dataset == "Full_Dataset")) {
    cat("Feature selection IMPROVED accuracy!\n")
    cat(sprintf("Reduced features from %d to %d (%.1f%% reduction)\n",
                results_df$N_Features[results_df$Dataset == "Full_Dataset"],
                results_df$N_Features[results_df$Dataset == "Feature_Selected_2Genera"],
                (1 - results_df$N_Features[results_df$Dataset == "Feature_Selected_2Genera"] / 
                   results_df$N_Features[results_df$Dataset == "Full_Dataset"]) * 100))
  } else {
    cat("Full dataset performed better than feature selection\n")
  }
}

# Save results
write.csv(results_df, "../results/microb/rf_feature_selection_comparison.csv", row.names = FALSE)
cat("\nResults saved to rf_feature_selection_comparison.csv\n")

# Create visualization if ggplot2 is available
if(require(ggplot2, quietly = TRUE)) {
  
  # Plot comparing accuracies
  p1 <- ggplot(results_df, aes(x = reorder(Dataset, Mean_Accuracy), 
                                y = Mean_Accuracy, fill = Dataset)) +
    geom_col(alpha = 0.7) +
    geom_errorbar(aes(ymin = Mean_Accuracy - SD_Accuracy,
                      ymax = Mean_Accuracy + SD_Accuracy),
                  width = 0.2) +
    coord_flip() +
    labs(title = "Random Forest Performance: Full vs Feature-Selected",
         x = "Dataset",
         y = "Mean Accuracy (5-fold CV)",
         subtitle = paste("Error bars show ±1 SD across folds")) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Set2")
  
  print(p1)
  ggsave("../results/microb/rf_feature_selection_comparison.png", p1, width = 10, height = 6)
  
  # Plot features vs accuracy
  p2 <- ggplot(results_df, aes(x = N_Features, y = Mean_Accuracy)) +
    geom_point(size = 4, aes(color = Dataset)) +
    geom_text(aes(label = Dataset), vjust = -1, size = 3) +
    labs(title = "Number of Features vs Model Performance",
         x = "Number of Features",
         y = "Mean Accuracy") +
    theme_minimal() +
    scale_x_log10() +
    theme(legend.position = "none")
  
  print(p2)
  ggsave("../results/microb/rf_features_vs_performance.png", p2, width = 8, height = 6)
}

# Stop parallel backend
stopImplicitCluster()