# Feature Selection Stability Analysis and Final Model

library(Boruta)
library(car)
library(pROC)
library(ggplot2)

# Load data and create train/test split
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")

# Create reproducible train/test split
set.seed(123)
train_indices <- sample(1:nrow(fk_metabolom_gt_scaled), size = 0.7 * nrow(fk_metabolom_gt_scaled))
train_data <- fk_metabolom_gt_scaled[train_indices, ]
test_data <- fk_metabolom_gt_scaled[-train_indices, ]

cat("Dataset loaded successfully\n")
cat("Total samples:", nrow(fk_metabolom_gt_scaled), "\n")
cat("Training samples:", nrow(train_data), "\n")
cat("Test samples:", nrow(test_data), "\n")
cat("Total features:", ncol(fk_metabolom_gt_scaled) - 1, "\n\n")

# 1. Feature Selection Stability Test
stability_test <- function(n_runs = 20) {
  all_selections <- list()
  
  cat("Running", n_runs, "Boruta iterations for stability analysis...\n")
  
  for(i in 1:n_runs) {
    set.seed(100 + i)
    boruta_result <- Boruta(gt ~ ., data = train_data, maxRuns = 100, doTrace = 0)
    selected <- getSelectedAttributes(boruta_result, withTentative = TRUE)
    all_selections[[i]] <- selected
    
    if(i %% 5 == 0) cat("Completed", i, "runs\n")
  }
  
  # Count how often each feature is selected
  all_features <- unique(unlist(all_selections))
  selection_frequency <- sapply(all_features, function(feat) {
    sum(sapply(all_selections, function(x) feat %in% x))
  })
  
  stability_df <- data.frame(
    feature = names(selection_frequency),
    selected_times = selection_frequency,
    selection_rate = selection_frequency / n_runs,
    stringsAsFactors = FALSE
  )
  
  return(stability_df[order(stability_df$selection_rate, decreasing = TRUE), ])
}

# Run stability analysis
stability_results <- stability_test(20)

# Display results
cat("\n=== FEATURE SELECTION STABILITY RESULTS ===\n")
print(stability_results)

# 2. Define stable features (selected in ≥50% of runs)
stable_features <- stability_results$feature[stability_results$selection_rate >= 0.5]
cat("\nStable features (selected ≥50% of time):", length(stable_features), "\n")
cat("Features:", paste(stable_features, collapse = ", "), "\n")

# 3. Build final model with stable features
if(length(stable_features) > 0) {
  
  # Remove 'predicted_prob' if it appears (data leakage)
  stable_features_clean <- stable_features[stable_features != "predicted_prob"]
  
  if(length(stable_features_clean) > 0) {
    cat("\n=== FINAL MODEL WITH STABLE FEATURES ===\n")
    
    # Build model
    formula_stable <- as.formula(paste("gt ~", paste(stable_features_clean, collapse = " + ")))
    stable_model <- glm(formula_stable, data = train_data, family = binomial)
    
    print(summary(stable_model))
    
    # Check VIF
    if(length(stable_features_clean) > 1) {
      cat("\nVariance Inflation Factors:\n")
      print(vif(stable_model))
    }
    
    # Test set performance
    test_predictions_stable <- predict(stable_model, newdata = test_data, type = "response")
    test_auc_stable <- auc(test_data$gt, test_predictions_stable)
    cat("\nTest AUC for stable feature model:", round(test_auc_stable, 3), "\n")
    
    # Create a simple model with top 3 most stable features (if you want)
    if(length(stable_features_clean) >= 3) {
      top3_features <- stable_features_clean[1:3]
      cat("\n=== SIMPLE MODEL WITH TOP 3 STABLE FEATURES ===\n")
      cat("Features:", paste(top3_features, collapse = ", "), "\n")
      
      formula_top3 <- as.formula(paste("gt ~", paste(top3_features, collapse = " + ")))
      top3_model <- glm(formula_top3, data = train_data, family = binomial)
      
      print(summary(top3_model))
      
      # VIF for top 3
      cat("\nVIF for top 3 model:\n")
      print(vif(top3_model))
      
      # Test performance
      test_pred_top3 <- predict(top3_model, newdata = test_data, type = "response")
      test_auc_top3 <- auc(test_data$gt, test_pred_top3)
      cat("Test AUC for top 3 stable features:", round(test_auc_top3, 3), "\n")
    }
    
  } else {
    cat("No stable features found after removing data leakage variables\n")
  }
  
} else {
  cat("No features selected in ≥50% of runs. Consider lowering threshold.\n")
  # Show features selected in ≥25% of runs
  somewhat_stable <- stability_results$feature[stability_results$selection_rate >= 0.25]
  cat("Features selected ≥25% of time:", paste(somewhat_stable, collapse = ", "), "\n")
}

# 4. Create visualization with stable features (like your original figure)
if(length(stable_features_clean) >= 3) {
  cat("\n=== CREATING COMBINED VISUALIZATION ===\n")
  
  library(cowplot)
  library(ggsignif)
  library(scales)
  
  # Create predictions for plots
  train_data$predicted_prob <- predict(stable_model, type = "response")
  
  # Get the top 3 stable features for individual plots
  top3_stable <- stable_features_clean[1:3]
  
  # Plot A: Protective Score Boxplot using stable features
  # Calculate protective score using ONLY training data scaling for stable features
  train_stable_metabolites <- scale(train_data[, top3_stable])
  train_stable_protective_score <- rowSums(train_stable_metabolites)
  
  # Add to training data
  train_data$stable_protective_score <- train_stable_protective_score
  
  pA_stable <- ggplot(train_data, aes(x = factor(gt), y = stable_protective_score, fill = factor(gt))) +
    geom_boxplot(alpha = 0.7) +
    geom_signif(comparisons = list(c("GT_absent", "GT_present")), 
                test = "t.test",
                map_signif_level = TRUE,  # This will show *, **, *** based on p-value
                y_position = max(train_data$stable_protective_score) + 0.5) +
    labs(
      title = "D",
      x = "Disease Status",
      y = "Stable Protective Score (Standardized)",
      fill = "Outcome"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  # Plot B: First stable feature prediction
  p1_stable <- ggplot(train_data, aes(x = get(top3_stable[1]), y = predicted_prob)) +
    geom_point(aes(color = factor(gt)), alpha = 0.9) +
    geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    labs(title = "A",
         x = paste(gsub("\\.", " ", gsub("^X", "", top3_stable[1])), "(log peak area)"),
         y = "Predicted Disease Probability",
         color = "Disease Status") +
    theme_bw() +
    theme(legend.position = "none")
  
  # Plot C: Second stable feature prediction  
  p2_stable <- ggplot(train_data, aes(x = get(top3_stable[2]), y = predicted_prob)) +
    geom_point(aes(color = factor(gt)), alpha = 0.9) +
    geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    labs(title = "B",
         x = paste(gsub("\\.", " ", gsub("^X", "", top3_stable[2])), "(log peak area)"),
         y = "Predicted Probability",
         color = "Disease Status") +
    theme_bw() +
    theme(legend.position = "none")
  
  # Plot D: Third stable feature prediction
  p3_stable <- ggplot(train_data, aes(x = get(top3_stable[3]), y = predicted_prob)) +
    geom_point(aes(color = factor(gt)), alpha = 0.9) +
    geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    labs(title = "C",
         x = paste(gsub("\\.", " ", gsub("^X", "", top3_stable[3])), "(log peak area)"),
         y = "Predicted Probability",
         color = "Disease Status") +
    scale_y_continuous(labels = percent_format()) +
    theme_bw() +
    theme(legend.position = "none")
  
  # Create legend from one plot
  p2_legend <- p2_stable + theme(legend.position = "bottom")
  shared_legend <- get_legend(p2_legend)
  
  # Combine plots
  combined_plots_stable <- plot_grid(
    p1_stable, p2_stable,
    p3_stable, pA_stable,
    labels = NULL,
    ncol = 2,
    align = "hv"
  )
  
  # Add shared legend
  final_plot_stable <- plot_grid(
    combined_plots_stable,
    shared_legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )
  
  print(final_plot_stable)
}

# save the plot
ggsave("../results/plots/stable_metabolite_visualization_final.png", 
       final_plot_stable, width = 12, height = 10, dpi = 300)

# 5. Create clean feature names and export stability table
cat("\n=== CREATING CSV EXPORT ===\n")

# Function to clean metabolite names
clean_metabolite_names <- function(names) {
  cleaned <- names
  
  # Remove leading X and numbers followed by dots
  cleaned <- gsub("^X\\d+\\.", "", cleaned)
  cleaned <- gsub("^X", "", cleaned)
  
  # Replace dots with spaces
  cleaned <- gsub("\\.", " ", cleaned)
  
  # Handle specific patterns
  cleaned <- gsub("Spectral Match to (.+) from NIST14", "\\1", cleaned)
  cleaned <- gsub("Massbank\\.[A-Z0-9]+\\.", "", cleaned)
  
  # Clean up common chemical notation
  cleaned <- gsub("sn glycero", "sn-glycero", cleaned)
  cleaned <- gsub("  +", " ", cleaned)  # Remove multiple spaces
  cleaned <- trimws(cleaned)  # Remove leading/trailing spaces
  
  # Capitalize first letter
  cleaned <- paste0(toupper(substr(cleaned, 1, 1)), substr(cleaned, 2, nchar(cleaned)))
  
  return(cleaned)
}

# Create clean stability table
stability_export <- stability_results
stability_export$clean_name <- clean_metabolite_names(stability_results$feature)

# Reorder columns
stability_export <- stability_export[, c("clean_name", "feature", "selected_times", "selection_rate")]
colnames(stability_export) <- c("Metabolite_Name", "Original_Feature_Name", "Times_Selected", "Selection_Rate")

# Export to CSV
write.csv(stability_export, "metabolite_stability_results.csv", row.names = FALSE)
cat("Stability table exported to: metabolite_stability_results.csv\n")

# Show preview of clean names
cat("\nPreview of cleaned names:\n")
print(head(stability_export[, c("Metabolite_Name", "Times_Selected", "Selection_Rate")], 10))

# 6. Cross-validation with stable features for robust AUC ± SE WITH ACCURACY METRICS
cat("\n=== CROSS-VALIDATION WITH STABLE FEATURES ===\n")

library(caret)

# Function to calculate balanced accuracy
calculate_balanced_accuracy <- function(actual, predicted_class) {
  cm <- table(Actual = actual, Predicted = predicted_class)
  # Handle case where one class might be missing in predictions
  if(nrow(cm) == 1 || ncol(cm) == 1) {
    return(0)
  }
  sensitivity <- cm[2,2] / (cm[2,2] + cm[2,1])
  specificity <- cm[1,1] / (cm[1,1] + cm[1,2])
  balanced_acc <- (sensitivity + specificity) / 2
  return(balanced_acc)
}

# Use the top 3 stable features consistently across all CV folds
if(length(stable_features_clean) >= 3) {
  
  cv_features <- stable_features_clean[1:3]
  cat("Using stable features for CV:", paste(cv_features, collapse = ", "), "\n")
  
  # 10-fold cross-validation repeated 5 times for robust estimates
  set.seed(123)
  
  # Store all metrics
  cv_results_detailed <- replicate(5, {
    # Create 10 folds
    folds <- createFolds(fk_metabolom_gt_scaled$gt, k = 10, returnTrain = FALSE)
    
    fold_metrics <- t(sapply(1:10, function(i) {
      # Split data
      test_fold_indices <- folds[[i]]
      train_fold <- fk_metabolom_gt_scaled[-test_fold_indices, ]
      test_fold <- fk_metabolom_gt_scaled[test_fold_indices, ]
      
      # Build model with stable features only
      formula_cv <- as.formula(paste("gt ~", paste(cv_features, collapse = " + ")))
      model_fold <- glm(formula_cv, data = train_fold, family = binomial)
      
      # Predict on test fold
      pred_prob <- predict(model_fold, test_fold, type = "response")
      
      # Calculate AUC
      auc_fold <- as.numeric(auc(test_fold$gt, pred_prob))
      
      # Calculate accuracy metrics using 0.5 threshold
      pred_class <- ifelse(pred_prob > 0.5, "GT_present", "GT_absent")
      actual_class <- test_fold$gt
      
      # Accuracy
      accuracy <- mean(pred_class == actual_class)
      
      # Balanced Accuracy
      balanced_acc <- calculate_balanced_accuracy(actual_class, pred_class)
      
      return(c(auc = auc_fold, accuracy = accuracy, balanced_accuracy = balanced_acc))
    }))
    
    return(fold_metrics)
  }, simplify = FALSE)
  
  # Flatten results - combine all 5 repetitions
  all_metrics <- do.call(rbind, cv_results_detailed)
  
  # Calculate statistics for each metric
  metrics_summary <- data.frame(
    Metric = c("AUC", "Accuracy", "Balanced_Accuracy"),
    Mean = c(mean(all_metrics[,"auc"]), 
             mean(all_metrics[,"accuracy"]), 
             mean(all_metrics[,"balanced_accuracy"])),
    SE = c(sd(all_metrics[,"auc"]) / sqrt(nrow(all_metrics)),
           sd(all_metrics[,"accuracy"]) / sqrt(nrow(all_metrics)),
           sd(all_metrics[,"balanced_accuracy"]) / sqrt(nrow(all_metrics))),
    SD = c(sd(all_metrics[,"auc"]),
           sd(all_metrics[,"accuracy"]),
           sd(all_metrics[,"balanced_accuracy"])),
    CI_Lower = c(mean(all_metrics[,"auc"]) - 1.96 * sd(all_metrics[,"auc"]) / sqrt(nrow(all_metrics)),
                 mean(all_metrics[,"accuracy"]) - 1.96 * sd(all_metrics[,"accuracy"]) / sqrt(nrow(all_metrics)),
                 mean(all_metrics[,"balanced_accuracy"]) - 1.96 * sd(all_metrics[,"balanced_accuracy"]) / sqrt(nrow(all_metrics))),
    CI_Upper = c(mean(all_metrics[,"auc"]) + 1.96 * sd(all_metrics[,"auc"]) / sqrt(nrow(all_metrics)),
                 mean(all_metrics[,"accuracy"]) + 1.96 * sd(all_metrics[,"accuracy"]) / sqrt(nrow(all_metrics)),
                 mean(all_metrics[,"balanced_accuracy"]) + 1.96 * sd(all_metrics[,"balanced_accuracy"]) / sqrt(nrow(all_metrics)))
  )
  
  # Results
  cat("\n=== CROSS-VALIDATION RESULTS (INDIVIDUAL METABOLITES) ===\n")
  cat("Number of CV iterations:", nrow(all_metrics), "\n")
  print(metrics_summary)
  cat("\nFormatted for reporting:\n")
  cat("AUC =", round(metrics_summary$Mean[1], 3), "±", round(metrics_summary$SE[1], 3), "\n")
  cat("Accuracy =", round(metrics_summary$Mean[2], 3), "±", round(metrics_summary$SE[2], 3), "\n")
  cat("Balanced Accuracy =", round(metrics_summary$Mean[3], 3), "±", round(metrics_summary$SE[3], 3), "\n")
  
  # Store for later comparison
  individual_metrics_summary <- metrics_summary
  
  # Also calculate with protective score approach
  cat("\n=== CV WITH PROTECTIVE SCORE ===\n")
  
  cv_protective_detailed <- replicate(5, {
    folds <- createFolds(fk_metabolom_gt_scaled$gt, k = 10, returnTrain = FALSE)
    
    fold_metrics_protective <- t(sapply(1:10, function(i) {
      test_fold_indices <- folds[[i]]
      train_fold <- fk_metabolom_gt_scaled[-test_fold_indices, ]
      test_fold <- fk_metabolom_gt_scaled[test_fold_indices, ]
      
      # Calculate protective score using training data scaling
      train_metabolites_cv <- scale(train_fold[, cv_features])
      train_protective_score_cv <- rowSums(train_metabolites_cv)
      
      # Apply same scaling to test fold
      test_means_cv <- attr(train_metabolites_cv, "scaled:center")
      test_sds_cv <- attr(train_metabolites_cv, "scaled:scale")
      test_metabolites_cv <- scale(test_fold[, cv_features], 
                                  center = test_means_cv, scale = test_sds_cv)
      test_protective_score_cv <- rowSums(test_metabolites_cv)
      
      # Build protective score model
      train_fold$protective_score <- train_protective_score_cv
      test_fold$protective_score <- test_protective_score_cv
      
      protective_model_cv <- glm(gt ~ protective_score, data = train_fold, family = binomial)
      pred_prob_protective <- predict(protective_model_cv, test_fold, type = "response")
      
      # AUC
      auc_protective <- as.numeric(auc(test_fold$gt, pred_prob_protective))
      
      # Accuracy metrics
      pred_class_protective <- ifelse(pred_prob_protective > 0.5, "GT_present", "GT_absent")
      actual_class <- test_fold$gt
      
      accuracy_protective <- mean(pred_class_protective == actual_class)
      balanced_acc_protective <- calculate_balanced_accuracy(actual_class, pred_class_protective)
      
      return(c(auc = auc_protective, 
               accuracy = accuracy_protective, 
               balanced_accuracy = balanced_acc_protective))
    }))
    
    return(fold_metrics_protective)
  }, simplify = FALSE)
  
  # Flatten protective score results
  all_metrics_protective <- do.call(rbind, cv_protective_detailed)
  
  # Calculate statistics for protective score
  protective_metrics_summary <- data.frame(
    Metric = c("AUC", "Accuracy", "Balanced_Accuracy"),
    Mean = c(mean(all_metrics_protective[,"auc"]), 
             mean(all_metrics_protective[,"accuracy"]), 
             mean(all_metrics_protective[,"balanced_accuracy"])),
    SE = c(sd(all_metrics_protective[,"auc"]) / sqrt(nrow(all_metrics_protective)),
           sd(all_metrics_protective[,"accuracy"]) / sqrt(nrow(all_metrics_protective)),
           sd(all_metrics_protective[,"balanced_accuracy"]) / sqrt(nrow(all_metrics_protective))),
    SD = c(sd(all_metrics_protective[,"auc"]),
           sd(all_metrics_protective[,"accuracy"]),
           sd(all_metrics_protective[,"balanced_accuracy"])),
    CI_Lower = c(mean(all_metrics_protective[,"auc"]) - 1.96 * sd(all_metrics_protective[,"auc"]) / sqrt(nrow(all_metrics_protective)),
                 mean(all_metrics_protective[,"accuracy"]) - 1.96 * sd(all_metrics_protective[,"accuracy"]) / sqrt(nrow(all_metrics_protective)),
                 mean(all_metrics_protective[,"balanced_accuracy"]) - 1.96 * sd(all_metrics_protective[,"balanced_accuracy"]) / sqrt(nrow(all_metrics_protective))),
    CI_Upper = c(mean(all_metrics_protective[,"auc"]) + 1.96 * sd(all_metrics_protective[,"auc"]) / sqrt(nrow(all_metrics_protective)),
                 mean(all_metrics_protective[,"accuracy"]) + 1.96 * sd(all_metrics_protective[,"accuracy"]) / sqrt(nrow(all_metrics_protective)),
                 mean(all_metrics_protective[,"balanced_accuracy"]) + 1.96 * sd(all_metrics_protective[,"balanced_accuracy"]) / sqrt(nrow(all_metrics_protective)))
  )
  
  cat("\n=== PROTECTIVE SCORE RESULTS ===\n")
  print(protective_metrics_summary)
  cat("\nFormatted for reporting:\n")
  cat("AUC =", round(protective_metrics_summary$Mean[1], 3), "±", round(protective_metrics_summary$SE[1], 3), "\n")
  cat("Accuracy =", round(protective_metrics_summary$Mean[2], 3), "±", round(protective_metrics_summary$SE[2], 3), "\n")
  cat("Balanced Accuracy =", round(protective_metrics_summary$Mean[3], 3), "±", round(protective_metrics_summary$SE[3], 3), "\n")
  
  # Create comparison table for your paper
  cat("\n=== FINAL COMPARISON TABLE ===\n")
  cat("Model Comparison:\n")
  cat("----------------------------------------\n")
  cat("Individual metabolites:\n")
  cat("  AUC:", round(individual_metrics_summary$Mean[1], 3), "±", round(individual_metrics_summary$SE[1], 3), "\n")
  cat("  Accuracy:", round(individual_metrics_summary$Mean[2], 3), "±", round(individual_metrics_summary$SE[2], 3), "\n")
  cat("  Balanced Accuracy:", round(individual_metrics_summary$Mean[3], 3), "±", round(individual_metrics_summary$SE[3], 3), "\n")
  cat("\nProtective score:\n")
  cat("  AUC:", round(protective_metrics_summary$Mean[1], 3), "±", round(protective_metrics_summary$SE[1], 3), "\n")
  cat("  Accuracy:", round(protective_metrics_summary$Mean[2], 3), "±", round(protective_metrics_summary$SE[2], 3), "\n")
  cat("  Balanced Accuracy:", round(protective_metrics_summary$Mean[3], 3), "±", round(protective_metrics_summary$SE[3], 3), "\n")
  cat("\nNumber of features: 3\n")
  cat("Method: Stability-based feature selection (≥50% selection rate)\n")
  
} else {
  cat("Not enough stable features for cross-validation\n")
}

# 7. Summary table for paper
cat("\n=== SUMMARY FOR PAPER ===\n")
cat("Feature selection performed across", 20, "independent Boruta runs\n")
cat("Most stable feature:", stability_export$Metabolite_Name[1], 
    "(selected", stability_results$selected_times[1], "/20 times)\n")
cat("Final model based on features selected ≥50% of runs\n")

# Create comprehensive results dataframe with all metrics
performance_results <- data.frame(
  Model = c("Stability_Logistic", "Protective_Score"),
  Method = c("Individual metabolites", "Combined score"),
  N_Features = c(3, 3),
  AUC = c(round(individual_metrics_summary$Mean[1], 3), 
          round(protective_metrics_summary$Mean[1], 3)),
  AUC_SE = c(round(individual_metrics_summary$SE[1], 3), 
             round(protective_metrics_summary$SE[1], 3)),
  AUC_SD = c(round(individual_metrics_summary$SD[1], 3), 
             round(protective_metrics_summary$SD[1], 3)),
  AUC_CI_Lower = c(round(individual_metrics_summary$CI_Lower[1], 3), 
                   round(protective_metrics_summary$CI_Lower[1], 3)),
  AUC_CI_Upper = c(round(individual_metrics_summary$CI_Upper[1], 3), 
                   round(protective_metrics_summary$CI_Upper[1], 3)),
  Accuracy = c(round(individual_metrics_summary$Mean[2], 3), 
               round(protective_metrics_summary$Mean[2], 3)),
  Accuracy_SE = c(round(individual_metrics_summary$SE[2], 3), 
                  round(protective_metrics_summary$SE[2], 3)),
  Accuracy_SD = c(round(individual_metrics_summary$SD[2], 3), 
                  round(protective_metrics_summary$SD[2], 3)),
  Balanced_Accuracy = c(round(individual_metrics_summary$Mean[3], 3), 
                        round(protective_metrics_summary$Mean[3], 3)),
  Balanced_Accuracy_SE = c(round(individual_metrics_summary$SE[3], 3), 
                           round(protective_metrics_summary$SE[3], 3)),
  Balanced_Accuracy_SD = c(round(individual_metrics_summary$SD[3], 3), 
                           round(protective_metrics_summary$SD[3], 3)),
  Features = c(paste(cv_features, collapse = "; "), 
               paste(cv_features, collapse = "; ")),
  Feature_Selection_Method = c("Stability-based (≥50% selection)", 
                               "Stability-based (≥50% selection)"),
  CV_Iterations = c(nrow(all_metrics), nrow(all_metrics_protective))
)

# Print results
print(performance_results)

# Save to CSV
write.csv(performance_results, "../results/boruta_stable_metabolite_CV_performance_complete.csv", row.names = FALSE)
cat("\nPerformance results saved to: results/boruta_stable_metabolite_CV_performance_complete.csv\n")

# Also save detailed CV results for supplementary material
detailed_cv_results <- data.frame(
  Iteration = rep(1:nrow(all_metrics), 2),
  Model = c(rep("Stability_Logistic", nrow(all_metrics)), 
            rep("Protective_Score", nrow(all_metrics_protective))),
  AUC = c(all_metrics[,"auc"], all_metrics_protective[,"auc"]),
  Accuracy = c(all_metrics[,"accuracy"], all_metrics_protective[,"accuracy"]),
  Balanced_Accuracy = c(all_metrics[,"balanced_accuracy"], 
                       all_metrics_protective[,"balanced_accuracy"])
)

write.csv(detailed_cv_results, "../results/boruta_detailed_cv_results_complete.csv", row.names = FALSE)
cat("Detailed CV results saved to: results/boruta_detailed_cv_results_complete.csv\n")

# Create a summary metrics table for easy comparison
summary_metrics_table <- rbind(
  cbind(Model = "Individual_Metabolites", individual_metrics_summary),
  cbind(Model = "Protective_Score", protective_metrics_summary)
)

write.csv(summary_metrics_table, "../results/boruta_metab_summary_allmetrics_comparison.csv", row.names = FALSE)
cat("Summary metrics comparison saved to: results/summary_metrics_comparison.csv\n")