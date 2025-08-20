# Feature Selection Stability Analysis and Final Model

library(Boruta)
library(car)
library(pROC)

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

# 4. Summary table for paper
cat("\n=== SUMMARY FOR PAPER ===\n")
cat("Feature selection performed across", 20, "independent Boruta runs\n")
cat("Most stable feature:", stability_results$feature[1], 
    "(selected", stability_results$selected_times[1], "/20 times)\n")
cat("Final model based on features selected ≥50% of runs\n")