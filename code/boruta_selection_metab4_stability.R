# Feature Selection Stability Analysis and Final Model

library(Boruta)
library(car)
library(pROC)

# Load data and create train/test split
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")

# Create reproducible train/test split
set.seed(123)
train_indices <- sample(1:nrow(fk_metabolom_gt_scaled), size = 0.7 * nrow(fk_metabolom_gt_scaled))
train_data <- fk_metabolom_gt_scaled[train_indices, ]
test_data <- fk_metabolom_gt_scaled[-train_indices, ]

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
                y_position = max(train_data$stable_protective_score) + 0.5,
                map_signif_level =TRUE) +
    labs(
      title = "A",
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
    labs(title = "B",
         x = paste(gsub("\\.", " ", gsub("^X", "", top3_stable[1])), "(log peak area)"),
         y = "Predicted Disease Probability",
         color = "Disease Status") +
    theme_bw() +
    theme(legend.position = "none")
  
  # Plot C: Second stable feature prediction  
  p2_stable <- ggplot(train_data, aes(x = get(top3_stable[2]), y = predicted_prob)) +
    geom_point(aes(color = factor(gt)), alpha = 0.9) +
    geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    labs(title = "C",
         x = paste(gsub("\\.", " ", gsub("^X", "", top3_stable[2])), "(log peak area)"),
         y = "Predicted Probability",
         color = "Disease Status") +
    theme_bw() +
    theme(legend.position = "none")
  
  # Plot D: Third stable feature prediction
  p3_stable <- ggplot(train_data, aes(x = get(top3_stable[3]), y = predicted_prob)) +
    geom_point(aes(color = factor(gt)), alpha = 0.9) +
    geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    labs(title = "D",
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
    pA_stable, p1_stable,
    p2_stable, p3_stable,
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
  


ggsave("../results/plots/stable_metab_visualization.png", 
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
write.csv(stability_export, "results/metabolite_stability_results.csv", row.names = FALSE)
cat("Stability table exported to: metabolite_stability_results.csv\n")

# Show preview of clean names
cat("\nPreview of cleaned names:\n")
print(head(stability_export[, c("Metabolite_Name", "Times_Selected", "Selection_Rate")], 10))

# 6. Summary table for paper
cat("\n=== SUMMARY FOR PAPER ===\n")
cat("Feature selection performed across", 20, "independent Boruta runs\n")
cat("Most stable feature:", stability_export$Metabolite_Name[1], 
    "(selected", stability_results$selected_times[1], "/20 times)\n")
cat("Final model based on features selected ≥50% of runs\n")