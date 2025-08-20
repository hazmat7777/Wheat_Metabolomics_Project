# NOTE TO USER- this script takes a long time to run (c. 15 min with 12 cores)

# Load required libraries
library(ranger)
library(caret)
library(doParallel)
library(phyloseq)
library(dplyr)


# ====================================
# SETUP
# ====================================

# Full dataset
ps_highdiv_relative <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")

# Set up parallel processing
n_cores <- detectCores() - 1
registerDoParallel(cores = n_cores)
cat("Using", n_cores, "cores for parallel processing\n")

# ====================================
# INITIALIZE RESULTS DATAFRAME
# ====================================
results_df <- data.frame(
  Model = character(),
  Dataset = character(),
  Taxonomic_Level = character(),
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
  Mean_Balanced_Accuracy = numeric(),
  SD_Balanced_Accuracy = numeric(),
  stringsAsFactors = FALSE
)

# ====================================
# FUNCTION TO RUN RF WITH TUNING
# ====================================
run_rf_tuning <- function(rf_data, dataset_name, model_name, tax_level = NA) {
  
  cat("\n\n=== RUNNING", model_name, "-", dataset_name, "===\n")
  cat("Dataset dimensions:", nrow(rf_data), "x", ncol(rf_data)-1, "\n\n")
  
  n_features <- ncol(rf_data) - 1
  
  # Adjust mtry range based on number of features
  max_mtry <- min(12, floor(sqrt(n_features)))
  mtry_range <- unique(c(3, 5, 7, 10, max_mtry))
  mtry_range <- mtry_range[mtry_range <= n_features]
  
  # Create training control
  train_control <- trainControl(
    method = "cv",
    number = 5,
    search = "grid",
    verboseIter = FALSE,
    allowParallel = TRUE,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  # Define tuning grid
  tune_grid <- expand.grid(
    mtry = mtry_range,
    splitrule = c("gini", "extratrees"),
    min.node.size = c(1, 3, 5)
  )
  
  # Train model with tuning
  cat("Starting grid search (", nrow(tune_grid), "combinations)...\n")
  model_tuned <- train(
    gt ~ .,
    data = rf_data,
    method = "ranger",
    trControl = train_control,
    tuneGrid = tune_grid,
    num.trees = 1000,  # Reduced for faster computation
    importance = 'impurity',
    num.threads = n_cores,
    metric = "ROC"
  )
  
  # Display results
  cat("\nBest parameters:\n")
  print(model_tuned$bestTune)
  
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
  final_balanced_accuracies <- c()
  
  for(i in 1:5) {
    train_data <- rf_data[-folds[[i]], ]
    test_data <- rf_data[folds[[i]], ]
    
    rf_model <- ranger(
      gt ~ ., 
      data = train_data, 
      num.trees = 1000,
      mtry = best_mtry,
      min.node.size = best_nodesize,
      splitrule = best_splitrule,
      importance = 'impurity'
    )
    
    predictions <- predict(rf_model, test_data)
    cm <- confusionMatrix(predictions$predictions, test_data$gt)
    
    accuracy <- cm$overall['Accuracy']
    balanced_accuracy <- cm$byClass['Balanced Accuracy']
    
    final_accuracies[i] <- accuracy
    final_balanced_accuracies[i] <- balanced_accuracy
    
    cat("  Fold", i, "- Acc:", round(accuracy, 3), 
        ", Balanced Acc:", round(balanced_accuracy, 3), "\n")
  }
  
  # Store results
  result_row <- data.frame(
    Model = model_name,
    Dataset = dataset_name,
    Taxonomic_Level = ifelse(is.na(tax_level), "ASV", tax_level),
    N_Features = n_features,
    Method = "Caret_GridSearch",
    num.trees = 1000,
    mtry = best_mtry,
    min.node.size = best_nodesize,
    sample.fraction = 0.632,
    splitrule = best_splitrule,
    Mean_Accuracy = round(mean(final_accuracies), 3),
    SD_Accuracy = round(sd(final_accuracies), 3),
    Min_Accuracy = round(min(final_accuracies), 3),
    Max_Accuracy = round(max(final_accuracies), 3),
    Mean_Balanced_Accuracy = round(mean(final_balanced_accuracies, na.rm = TRUE), 3),
    SD_Balanced_Accuracy = round(sd(final_balanced_accuracies, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
  
  cat("\nMean Accuracy:", result_row$Mean_Accuracy, "±", result_row$SD_Accuracy, "\n")
  cat("Mean Balanced Accuracy:", result_row$Mean_Balanced_Accuracy, 
      "±", result_row$SD_Balanced_Accuracy, "\n")
  
  return(list(
    result_row = result_row,
    model = model_tuned,
    accuracies = final_accuracies
  ))
}

# ====================================
# 1. TAXONOMIC LEVEL MODELS
# ====================================
cat("\n============================================================\n")
cat("EXPERIMENT SET 1: TAXONOMIC AGGLOMERATION MODELS\n")
cat("============================================================\n")

taxonomic_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

for(tax_level in taxonomic_levels) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Processing", tax_level, "level\n")
  cat("------------------------------------------------------------\n")
  
  # Tax agglomeration
  ps_agg <- tax_glom(ps_highdiv_relative, tax_level, NArm = TRUE)
  
  # Set taxa names to the taxonomic level
  taxa_names(ps_agg) <- make.unique(as.character(tax_table(ps_agg)[, tax_level]))
  
  cat("Features after aggregation to", tax_level, ":", ntaxa(ps_agg), "\n")
  
  # Prepare data
  sample_df <- data.frame(sample_data(ps_agg))
  otu_df <- data.frame(t(as.matrix(otu_table(ps_agg))))
  rf_data <- cbind(sample_df, otu_df)
  rf_data <- rf_data[!is.na(rf_data$gt), ]
  rf_data$gt <- as.factor(rf_data$gt)
  
  # Clean columns
  rf_data$sample_name <- NULL
  rf_data$og_sample_names <- NULL
  if("tillage" %in% names(rf_data)) rf_data$tillage <- NULL
  
  # Run RF with tuning
  tax_results <- run_rf_tuning(
    rf_data = rf_data,
    dataset_name = paste0("TaxAgg_", tax_level),
    model_name = paste0("RandomForest_", tax_level),
    tax_level = tax_level
  )
  
  results_df <- rbind(results_df, tax_results$result_row)
}

# ====================================
# 2. FEATURE-SELECTED MODEL (2 GENERA)
# ====================================
cat("\n============================================================\n")
cat("EXPERIMENT 2: FEATURE-SELECTED (2 GENERA)\n")
cat("============================================================\n")

genera_implicated <- c("Pseudomonas", "Bacillus")

tax_df <- as.data.frame(tax_table(ps_highdiv_relative))
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]
ps_fs <- prune_taxa(taxa_to_keep, ps_highdiv_relative)

cat("Selected genera:", paste(genera_implicated, collapse = ", "), "\n")
cat("Number of ASVs from selected genera:", ntaxa(ps_fs), "\n")

# Prepare feature-selected data
sample_df_fs <- data.frame(sample_data(ps_fs))
otu_df_fs <- data.frame(t(as.matrix(otu_table(ps_fs))))
rf_data_fs <- cbind(sample_df_fs, otu_df_fs)
rf_data_fs <- rf_data_fs[!is.na(rf_data_fs$gt), ]
rf_data_fs$gt <- as.factor(rf_data_fs$gt)

# Clean columns
rf_data_fs$sample_name <- NULL
rf_data_fs$og_sample_names <- NULL
if("tillage" %in% names(rf_data_fs)) rf_data_fs$tillage <- NULL

fs_results <- run_rf_tuning(
  rf_data = rf_data_fs,
  dataset_name = "Feature_Selected_2Genera",
  model_name = "RandomForest_Pseudomonas_Bacillus",
  tax_level = "Selected_ASVs"
)

results_df <- rbind(results_df, fs_results$result_row)

# ====================================
# 3. TOP VARIANCE ASVs MODEL (OPTIONAL)
# ====================================
if(FALSE) {  # Set to FALSE to skip
  cat("\n============================================================\n")
  cat("EXPERIMENT 3: TOP VARIANCE ASVs\n")
  cat("============================================================\n")
  
  # Get all ASV data
  sample_df_all <- data.frame(sample_data(ps_highdiv_relative))
  otu_df_all <- data.frame(t(as.matrix(otu_table(ps_highdiv_relative))))
  rf_data_all <- cbind(sample_df_all, otu_df_all)
  rf_data_all <- rf_data_all[!is.na(rf_data_all$gt), ]
  rf_data_all$gt <- as.factor(rf_data_all$gt)
  
  # Clean columns
  rf_data_all$sample_name <- NULL
  rf_data_all$og_sample_names <- NULL
  if("tillage" %in% names(rf_data_all)) rf_data_all$tillage <- NULL
  
  # Calculate variance for each ASV
  feature_cols <- setdiff(names(rf_data_all), "gt")
  feature_vars <- apply(rf_data_all[, feature_cols], 2, var)
  
  # Select top 1000 variance ASVs
  n_top <- min(1000, length(feature_vars))
  top_features <- names(sort(feature_vars, decreasing = TRUE)[1:n_top])
  rf_data_topvar <- rf_data_all[, c("gt", top_features)]
  
  cat("Selected top", n_top, "variance ASVs\n")
  
  topvar_results <- run_rf_tuning(
    rf_data = rf_data_topvar,
    dataset_name = paste0("Top", n_top, "_Variance_ASVs"),
    model_name = paste0("RandomForest_Top", n_top, "Var"),
    tax_level = "Selected_ASVs"
  )
  
  results_df <- rbind(results_df, topvar_results$result_row)
}

# ====================================
# DISPLAY FINAL COMPARISON
# ====================================
cat("\n\n============================================================\n")
cat("FINAL COMPARISON OF ALL MODELS\n")
cat("============================================================\n\n")

# Sort by balanced accuracy (more robust for imbalanced data)
results_df <- results_df[order(results_df$Mean_Balanced_Accuracy, decreasing = TRUE), ]

# Display results
print(results_df[, c("Model", "Taxonomic_Level", "N_Features", 
                     "Mean_Accuracy", "SD_Accuracy",
                     "Mean_Balanced_Accuracy", "SD_Balanced_Accuracy")])

# # Create detailed summary
# cat("\n=== PERFORMANCE SUMMARY ===\n")
# cat("Sorted by Balanced Accuracy:\n\n")
# for(i in 1:nrow(results_df)) {
#   cat(sprintf("%2d. %-30s: %5d features, Acc = %.3f±%.3f, Balanced = %.3f±%.3f\n",
#               i,
#               results_df$Dataset[i],
#               results_df$N_Features[i],
#               results_df$Mean_Accuracy[i],
#               results_df$SD_Accuracy[i],
#               results_df$Mean_Balanced_Accuracy[i],
#               results_df$SD_Balanced_Accuracy[i]))
# }

# # Identify best model
# best_idx <- which.max(results_df$Mean_Balanced_Accuracy)
# cat("\n=== BEST PERFORMING MODEL (by Balanced Accuracy) ===\n")
# cat("Dataset:", results_df$Dataset[best_idx], "\n")
# cat("Taxonomic Level:", results_df$Taxonomic_Level[best_idx], "\n")
# cat("Features:", results_df$N_Features[best_idx], "\n")
# cat("Accuracy:", results_df$Mean_Accuracy[best_idx], "±", results_df$SD_Accuracy[best_idx], "\n")
# cat("Balanced Accuracy:", results_df$Mean_Balanced_Accuracy[best_idx], 
#     "±", results_df$SD_Balanced_Accuracy[best_idx], "\n")
# cat("Best Parameters: mtry =", results_df$mtry[best_idx], 
#     ", min.node.size =", results_df$min.node.size[best_idx],
#     ", splitrule =", results_df$splitrule[best_idx], "\n")

# # Save results
write.csv(results_df, "../results/microb/rf_taxonomic_comparison_results.csv", row.names = FALSE)
# cat("\nResults saved to 'rf_taxonomic_comparison_results.csv'\n")

# # Create visualization if ggplot2 is available
# if(require(ggplot2, quietly = TRUE)) {
#   library(ggplot2)
  
#   # Plot comparing balanced accuracies across taxonomic levels
#   p1 <- ggplot(results_df, aes(x = reorder(Dataset, Mean_Balanced_Accuracy), 
#                                 y = Mean_Balanced_Accuracy, 
#                                 fill = Taxonomic_Level)) +
#     geom_col(alpha = 0.7) +
#     geom_errorbar(aes(ymin = Mean_Balanced_Accuracy - SD_Balanced_Accuracy,
#                       ymax = Mean_Balanced_Accuracy + SD_Balanced_Accuracy),
#                   width = 0.2) +
#     coord_flip() +
#     labs(title = "Random Forest Performance Across Taxonomic Levels",
#          x = "Dataset",
#          y = "Mean Balanced Accuracy (5-fold CV)",
#          subtitle = "Error bars show ±1 SD across folds") +
#     theme_minimal() +
#     scale_fill_brewer(palette = "Set2")
  
#   ggsave("rf_taxonomic_comparison.png", p1, width = 12, height = 8)
  
#   # Plot features vs accuracy
#   p2 <- ggplot(results_df, aes(x = N_Features, y = Mean_Balanced_Accuracy)) +
#     geom_point(size = 4, aes(color = Taxonomic_Level)) +
#     geom_text(aes(label = Taxonomic_Level), vjust = -1, size = 3) +
#     geom_errorbar(aes(ymin = Mean_Balanced_Accuracy - SD_Balanced_Accuracy,
#                       ymax = Mean_Balanced_Accuracy + SD_Balanced_Accuracy),
#                   width = 0.01, alpha = 0.3) +
#     labs(title = "Number of Features vs Model Performance",
#          x = "Number of Features (log scale)",
#          y = "Mean Balanced Accuracy") +
#     theme_minimal() +
#     scale_x_log10() +
#     scale_color_brewer(palette = "Set2")
  
#   ggsave("rf_features_vs_performance.png", p2, width = 10, height = 6)
  
#   cat("Visualizations saved.\n")
# }

# # Stop parallel backend
# stopImplicitCluster()

# cat("\n=== ANALYSIS COMPLETE ===\n")