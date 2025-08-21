# NOTE TO USER- this script takes a long time to run (c. 15 min with 12 cores)

# Load required libraries
library(ranger)
library(caret)
library(doParallel)
library(phyloseq)
library(dplyr)
library(pROC)  # Added for AUC calculations

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
# INITIALIZE RESULTS DATAFRAME WITH AUC AND SE
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
  SE_Accuracy = numeric(),
  Mean_Balanced_Accuracy = numeric(),
  SD_Balanced_Accuracy = numeric(),
  SE_Balanced_Accuracy = numeric(),
  Mean_AUC = numeric(),
  SD_AUC = numeric(),
  SE_AUC = numeric(),
  AUC_CI_Lower = numeric(),
  AUC_CI_Upper = numeric(),
  stringsAsFactors = FALSE
)

# ====================================
# UPDATED FUNCTION TO RUN RF WITH AUC CALCULATIONS
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
  
  # Final evaluation with best parameters INCLUDING AUC
  cat("\nFinal evaluation with best parameters (including AUC)...\n")
  set.seed(1250)
  folds <- createFolds(rf_data$gt, k = 5, returnTrain = FALSE)
  
  # Storage for all metrics
  final_accuracies <- c()
  final_balanced_accuracies <- c()
  final_aucs <- c()
  
  for(i in 1:5) {
    train_data <- rf_data[-folds[[i]], ]
    test_data <- rf_data[folds[[i]], ]
    
    # Ensure both datasets have the same factor levels
    standard_levels <- c("GT_absent", "GT_present")
    train_data$gt <- factor(train_data$gt, levels = standard_levels)
    test_data$gt <- factor(test_data$gt, levels = standard_levels)
    
    # Check if both classes are present in test data
    cat("  Fold", i, "- Test classes:", paste(table(test_data$gt), collapse = "/"), "\n")
    
    # Skip fold if only one class present
    if(length(unique(test_data$gt)) < 2) {
      cat("    Skipping fold - only one class present\n")
      final_accuracies[i] <- NA
      final_balanced_accuracies[i] <- NA
      final_aucs[i] <- NA
      next
    }
    
    # Train ranger model with probability predictions
    rf_model <- ranger(
      gt ~ ., 
      data = train_data, 
      num.trees = 1000,
      mtry = best_mtry,
      min.node.size = best_nodesize,
      splitrule = best_splitrule,
      importance = 'impurity',
      probability = TRUE  # Enable probability predictions for AUC
    )
    
    # Get class predictions for accuracy (predict without probabilities for class prediction)
    rf_model_class <- ranger(
      gt ~ ., 
      data = train_data, 
      num.trees = 1000,
      mtry = best_mtry,
      min.node.size = best_nodesize,
      splitrule = best_splitrule,
      importance = 'impurity',
      probability = FALSE  # Get class predictions
    )
    
    predictions_class <- predict(rf_model_class, test_data)
    
    # Get probability predictions for AUC
    predictions_prob <- predict(rf_model, test_data)
    
    # Debug: Check dimensions
    cat("    Debug - Test data rows:", nrow(test_data), "\n")
    cat("    Debug - Class predictions length:", length(predictions_class$predictions), "\n")
    cat("    Debug - Prob predictions dims:", dim(predictions_prob$predictions), "\n")
    
    # Ensure predictions have same factor levels
    predicted_classes <- factor(predictions_class$predictions, levels = standard_levels)
    actual_classes <- factor(test_data$gt, levels = standard_levels)
    
    cat("    Final - Predicted length:", length(predicted_classes), 
        ", Actual length:", length(actual_classes), "\n")
    
    # Calculate confusion matrix
    accuracy <- NA
    balanced_accuracy <- NA
    
    tryCatch({
      cm <- confusionMatrix(predicted_classes, actual_classes, positive = "GT_present")
      accuracy <- as.numeric(cm$overall['Accuracy'])
      balanced_accuracy <- as.numeric(cm$byClass['Balanced Accuracy'])
      
      # Handle NA balanced accuracy (can happen with extreme imbalance)
      if(is.na(balanced_accuracy)) {
        balanced_accuracy <- accuracy
      }
    }, error = function(e) {
      cat("    Error in confusion matrix:", e$message, "\n")
      cat("    Predicted classes length:", length(predicted_classes), "\n")
      cat("    Actual classes length:", length(actual_classes), "\n")
      cat("    Predicted levels:", paste(levels(predicted_classes), collapse = ", "), "\n")
      cat("    Actual levels:", paste(levels(actual_classes), collapse = ", "), "\n")
      
      # Manual accuracy calculation as fallback
      if(length(predicted_classes) == length(actual_classes)) {
        accuracy <- sum(predicted_classes == actual_classes, na.rm = TRUE) / length(actual_classes)
        balanced_accuracy <- accuracy
        cat("    Manual accuracy calculation:", round(accuracy, 3), "\n")
      } else {
        accuracy <- 0.5
        balanced_accuracy <- 0.5
        cat("    Using default accuracy (0.5) due to length mismatch\n")
      }
    })
    
    # Extract probabilities for GT_present 
    if(ncol(predictions_prob$predictions) == 2) {
      # Check if columns are named
      if("GT_present" %in% colnames(predictions_prob$predictions)) {
        prob_positive <- predictions_prob$predictions[, "GT_present"]
      } else {
        # Assume second column is positive class
        prob_positive <- predictions_prob$predictions[, 2]
      }
    } else {
      cat("    Warning: Expected 2 probability columns, got", ncol(predictions_prob$predictions), "\n")
      prob_positive <- predictions_prob$predictions[, ncol(predictions_prob$predictions)]
    }
    
    cat("    Debug - Probability vector length:", length(prob_positive), "\n")
    cat("    Debug - Test gt length:", length(test_data$gt), "\n")
    
    # Calculate AUC
    auc_value <- 0.5  # Default value
    
    tryCatch({
      # Convert factor to numeric for pROC (GT_absent=1, GT_present=2)
      numeric_gt <- as.numeric(test_data$gt) - 1  # Convert to 0/1
      
      # Check that we have valid probabilities
      if(length(prob_positive) == length(numeric_gt) && !any(is.na(prob_positive))) {
        roc_obj <- roc(numeric_gt, prob_positive, quiet = TRUE)
        auc_value <- as.numeric(auc(roc_obj))
      } else {
        cat("    Warning: Invalid probabilities for AUC calculation\n")
        cat("    Probability length:", length(prob_positive), ", GT length:", length(numeric_gt), "\n")
        auc_value <- 0.5
      }
    }, error = function(e) {
      cat("    Error calculating AUC:", e$message, "\n")
      auc_value <- 0.5  # Default to random performance
    })
    
    # Ensure all metrics are valid numbers
    if(is.na(accuracy) || !is.finite(accuracy)) accuracy <- 0.5
    if(is.na(balanced_accuracy) || !is.finite(balanced_accuracy)) balanced_accuracy <- 0.5
    if(is.na(auc_value) || !is.finite(auc_value)) auc_value <- 0.5
    
    # Store metrics
    final_accuracies[i] <- accuracy
    final_balanced_accuracies[i] <- balanced_accuracy
    final_aucs[i] <- auc_value
    
    cat("    AUC:", round(auc_value, 3),
        ", Acc:", round(accuracy, 3), 
        ", Balanced Acc:", round(balanced_accuracy, 3), "\n")
  }
  
  # Calculate summary statistics with SE
  n_folds <- length(final_aucs)
  
  # Remove any NA values
  final_accuracies <- final_accuracies[!is.na(final_accuracies)]
  final_balanced_accuracies <- final_balanced_accuracies[!is.na(final_balanced_accuracies)]
  final_aucs <- final_aucs[!is.na(final_aucs)]
  
  # AUC statistics
  mean_auc <- mean(final_aucs)
  sd_auc <- sd(final_aucs)
  se_auc <- sd_auc / sqrt(length(final_aucs))
  ci_lower <- mean_auc - 1.96 * se_auc
  ci_upper <- mean_auc + 1.96 * se_auc
  
  # Accuracy statistics
  mean_acc <- mean(final_accuracies)
  sd_acc <- sd(final_accuracies)
  se_acc <- sd_acc / sqrt(length(final_accuracies))
  
  # Balanced Accuracy statistics
  mean_bal_acc <- mean(final_balanced_accuracies)
  sd_bal_acc <- sd(final_balanced_accuracies)
  se_bal_acc <- sd_bal_acc / sqrt(length(final_balanced_accuracies))
  
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
    Mean_Accuracy = round(mean_acc, 3),
    SD_Accuracy = round(sd_acc, 3),
    SE_Accuracy = round(se_acc, 3),
    Mean_Balanced_Accuracy = round(mean_bal_acc, 3),
    SD_Balanced_Accuracy = round(sd_bal_acc, 3),
    SE_Balanced_Accuracy = round(se_bal_acc, 3),
    Mean_AUC = round(mean_auc, 3),
    SD_AUC = round(sd_auc, 3),
    SE_AUC = round(se_auc, 3),
    AUC_CI_Lower = round(ci_lower, 3),
    AUC_CI_Upper = round(ci_upper, 3),
    stringsAsFactors = FALSE
  )
  
  cat("\nSUMMARY STATISTICS:\n")
  cat("AUC:", result_row$Mean_AUC, "±", result_row$SE_AUC, 
      "(95% CI:", result_row$AUC_CI_Lower, "-", result_row$AUC_CI_Upper, ")\n")
  cat("Accuracy:", result_row$Mean_Accuracy, "±", result_row$SE_Accuracy, "\n")
  cat("Balanced Accuracy:", result_row$Mean_Balanced_Accuracy, "±", result_row$SE_Balanced_Accuracy, "\n")
  
  return(list(
    result_row = result_row,
    model = model_tuned,
    aucs = final_aucs,
    accuracies = final_accuracies,
    balanced_accuracies = final_balanced_accuracies
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
  
  # STANDARDIZE GT FACTOR LEVELS
  rf_data$gt <- as.factor(rf_data$gt)
  cat("Original gt levels:", paste(levels(rf_data$gt), collapse = ", "), "\n")
  cat("gt distribution:", paste(table(rf_data$gt), collapse = "/"), "\n")
  
  # Ensure consistent factor levels (GT_absent, GT_present)
  if(all(levels(rf_data$gt) %in% c("0", "1"))) {
    # Convert numeric to proper labels
    rf_data$gt <- factor(rf_data$gt, levels = c("0", "1"), labels = c("GT_absent", "GT_present"))
  } else if(all(levels(rf_data$gt) %in% c("GT_absent", "GT_present"))) {
    # Already correct, but ensure proper ordering
    rf_data$gt <- factor(rf_data$gt, levels = c("GT_absent", "GT_present"))
  } else {
    cat("Warning: Unexpected gt levels. Converting to standard format.\n")
    # Assume first level is negative, second is positive
    current_levels <- levels(rf_data$gt)
    rf_data$gt <- factor(rf_data$gt, levels = current_levels, labels = c("GT_absent", "GT_present"))
  }
  
  cat("Standardized gt levels:", paste(levels(rf_data$gt), collapse = ", "), "\n")
  cat("Final gt distribution:", paste(table(rf_data$gt), collapse = "/"), "\n")
  
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
# DISPLAY FINAL COMPARISON WITH AUC
# ====================================
cat("\n\n============================================================\n")
cat("FINAL COMPARISON OF ALL MODELS (with AUC ± SE)\n")
cat("============================================================\n\n")

# Sort by AUC
results_df <- results_df[order(results_df$Mean_AUC, decreasing = TRUE), ]

# Display results with AUC
print(results_df[, c("Model", "Taxonomic_Level", "N_Features", 
                     "Mean_AUC", "SE_AUC", "AUC_CI_Lower", "AUC_CI_Upper",
                     "Mean_Balanced_Accuracy", "SE_Balanced_Accuracy")])

# Create comparison table
cat("\n=== COMPARISON TABLE (AUC ± SE) ===\n")
comparison_table <- data.frame(
  Model = results_df$Model,
  Taxonomic_Level = results_df$Taxonomic_Level,
  N_Features = results_df$N_Features,
  AUC = paste0(results_df$Mean_AUC, " ± ", results_df$SE_AUC),
  AUC_CI = paste0("(", results_df$AUC_CI_Lower, "-", results_df$AUC_CI_Upper, ")"),
  Balanced_Accuracy = paste0(results_df$Mean_Balanced_Accuracy, " ± ", results_df$SE_Balanced_Accuracy),
  stringsAsFactors = FALSE
)

print(comparison_table)

# Statistical significance testing (vs AUC = 0.5)
cat("\n=== STATISTICAL SIGNIFICANCE (vs AUC = 0.5) ===\n")
for(i in 1:nrow(results_df)) {
  model_name <- results_df$Model[i]
  mean_auc <- results_df$Mean_AUC[i]
  se_auc <- results_df$SE_AUC[i]
  
  # Simple t-test approximation
  t_stat <- (mean_auc - 0.5) / se_auc
  df <- 4  # 5-fold CV = 4 degrees of freedom
  p_value <- pt(t_stat, df, lower.tail = FALSE)
  
  sig_text <- ifelse(p_value < 0.05, "significant", "not significant")
  cat(sprintf("%-40s: AUC = %.3f ± %.3f (vs 0.5, p ≈ %.4f, %s)\n", 
              model_name, mean_auc, se_auc, p_value, sig_text))
}

# Save results with all metrics
write.csv(results_df, "../results/microb/rf_taxonomic_comparison_results_with_auc.csv", row.names = FALSE)
write.csv(comparison_table, "../results/microb/rf_comparison_table_auc.csv", row.names = FALSE)

cat("\nResults with AUC ± SE saved to CSV files\n")

# Stop parallel backend
stopImplicitCluster()

cat("\n=== ANALYSIS COMPLETE ===\n")