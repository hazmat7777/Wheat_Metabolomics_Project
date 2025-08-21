# Microbiome Stability-Based Feature Selection Analysis
# Following the same approach as metabolomics: Boruta stability → GLM with CV

library(phyloseq)
library(Boruta)
library(car)
library(pROC)
library(caret)
library(ggplot2)
library(cowplot)
library(ggsignif)

# ====================================
# DATA PREPARATION
# ====================================

# Load and prepare data
ps <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")
ps_genus <- readRDS("../data/metabarcoding/ps_16S_highdiv_genus_relative.rds")
# taxa_names(ps_genus) <- make.unique(as.character(tax_table(ps_genus)[, "Genus"])) 
# sample_names(ps_genus) <- as.data.frame(sample_data(ps_genus))$sample_name

cat("ESVs after tax_glom:", ntaxa(ps_genus), "\n")

# Calculate prevalence
otu_tab <- otu_table(ps_genus)
prevalence_df <- data.frame(
  Prevalence = taxa_sums(otu_tab > 0),
  TotalAbundance = taxa_sums(otu_tab),
  Genus = taxa_names(ps_genus)
)
prevalence_df$PrevalenceProp <- prevalence_df$Prevalence / nsamples(ps_genus)

cat("Total samples:", nsamples(ps_genus), "\n")
cat("Total genera:", ntaxa(ps_genus), "\n")
cat("Prevalence summary:\n")
summary(prevalence_df$PrevalenceProp)

# Filter to genera present in at least 10% of samples
ps_genus_filtered <- prune_taxa(prevalence_df$PrevalenceProp >= 0.10, ps_genus)
cat("After 10% prevalence filtering:", ntaxa(ps_genus_filtered), "genera remain\n")

# Convert to data frame
genus_data <- as.data.frame(t(otu_table(ps_genus_filtered)))
outcome_var <- as.factor(sample_data(ps_genus_filtered)$gt)
full_data <- cbind(gt = outcome_var, genus_data)

cat("Final dataset dimensions:", nrow(full_data), "x", ncol(full_data)-1, "\n")
cat("Class distribution:\n")
print(table(full_data$gt))

# ====================================
# STABILITY-BASED FEATURE SELECTION
# ====================================

stability_test_microbiome <- function(data, n_runs = 20) {
  all_selections <- list()
  
  cat("Running", n_runs, "Boruta iterations for microbiome stability analysis...\n")
  
  for(i in 1:n_runs) {
    set.seed(100 + i)
    boruta_result <- Boruta(gt ~ ., data = data, maxRuns = 100, doTrace = 0)
    selected <- getSelectedAttributes(boruta_result, withTentative = TRUE)
    all_selections[[i]] <- selected
    
    if(i %% 5 == 0) cat("Completed", i, "runs\n")
  }
  
  # Count selection frequency
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

# Run stability analysis on full dataset
microbiome_stability_results <- stability_test_microbiome(full_data, 20)

cat("\n=== MICROBIOME FEATURE SELECTION STABILITY RESULTS ===\n")
print(microbiome_stability_results)

# Define stable features (≥50% selection rate)
stable_genera <- microbiome_stability_results$feature[microbiome_stability_results$selection_rate >= 0.5]
cat("\nStable genera (selected ≥50% of time):", length(stable_genera), "\n")

if(length(stable_genera) > 0) {
  cat("Stable genera:", paste(stable_genera, collapse = ", "), "\n")
  
  # Show taxonomic information for stable genera
  if(length(stable_genera) <= 10) {  # Only show if reasonable number
    cat("\n=== TAXONOMIC INFORMATION FOR STABLE GENERA ===\n")
    stable_tax <- tax_table(ps_genus_filtered)[stable_genera, ]
    print(stable_tax)
  }
} else {
  cat("No genera selected ≥50% of time. Trying 25% threshold...\n")
  stable_genera <- microbiome_stability_results$feature[microbiome_stability_results$selection_rate >= 0.25]
  cat("Genera selected ≥25% of time:", length(stable_genera), "\n")
  if(length(stable_genera) > 0) {
    cat("Features:", paste(stable_genera, collapse = ", "), "\n")
  }
}

# ====================================
# CROSS-VALIDATION WITH STABLE FEATURES AND ALL METRICS
# ====================================

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

if(length(stable_genera) > 0) {
  cat("\n=== CROSS-VALIDATION WITH STABLE MICROBIOME FEATURES ===\n")
  
  # Determine how many features to use (max 5 for stability)
  n_features_to_use <- min(length(stable_genera), 5)
  cv_genera <- stable_genera[1:n_features_to_use]
  
  cat("Using top", n_features_to_use, "stable genera for CV:", paste(cv_genera, collapse = ", "), "\n")
  
  # 10-fold cross-validation repeated 5 times with all metrics
  set.seed(123)
  
  cv_results_microbiome_detailed <- replicate(5, {
    folds <- createFolds(full_data$gt, k = 10, returnTrain = FALSE)
    
    fold_metrics <- t(sapply(1:10, function(i) {
      test_fold_indices <- folds[[i]]
      train_fold <- full_data[-test_fold_indices, ]
      test_fold <- full_data[test_fold_indices, ]
      
      # Build model with stable genera only
      if(n_features_to_use == 1) {
        formula_cv <- as.formula(paste("gt ~", paste0("`", cv_genera, "`")))
      } else {
        genus_names_clean <- paste0("`", cv_genera, "`")
        formula_cv <- as.formula(paste("gt ~", paste(genus_names_clean, collapse = " + ")))
      }
      
      model_fold <- glm(formula_cv, data = train_fold, family = binomial)
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
  all_metrics_microbiome <- do.call(rbind, cv_results_microbiome_detailed)
  
  # Calculate statistics for each metric
  microbiome_metrics_summary <- data.frame(
    Metric = c("AUC", "Accuracy", "Balanced_Accuracy"),
    Mean = c(mean(all_metrics_microbiome[,"auc"]), 
             mean(all_metrics_microbiome[,"accuracy"]), 
             mean(all_metrics_microbiome[,"balanced_accuracy"])),
    SE = c(sd(all_metrics_microbiome[,"auc"]) / sqrt(nrow(all_metrics_microbiome)),
           sd(all_metrics_microbiome[,"accuracy"]) / sqrt(nrow(all_metrics_microbiome)),
           sd(all_metrics_microbiome[,"balanced_accuracy"]) / sqrt(nrow(all_metrics_microbiome))),
    SD = c(sd(all_metrics_microbiome[,"auc"]),
           sd(all_metrics_microbiome[,"accuracy"]),
           sd(all_metrics_microbiome[,"balanced_accuracy"])),
    CI_Lower = c(mean(all_metrics_microbiome[,"auc"]) - 1.96 * sd(all_metrics_microbiome[,"auc"]) / sqrt(nrow(all_metrics_microbiome)),
                 mean(all_metrics_microbiome[,"accuracy"]) - 1.96 * sd(all_metrics_microbiome[,"accuracy"]) / sqrt(nrow(all_metrics_microbiome)),
                 mean(all_metrics_microbiome[,"balanced_accuracy"]) - 1.96 * sd(all_metrics_microbiome[,"balanced_accuracy"]) / sqrt(nrow(all_metrics_microbiome))),
    CI_Upper = c(mean(all_metrics_microbiome[,"auc"]) + 1.96 * sd(all_metrics_microbiome[,"auc"]) / sqrt(nrow(all_metrics_microbiome)),
                 mean(all_metrics_microbiome[,"accuracy"]) + 1.96 * sd(all_metrics_microbiome[,"accuracy"]) / sqrt(nrow(all_metrics_microbiome)),
                 mean(all_metrics_microbiome[,"balanced_accuracy"]) + 1.96 * sd(all_metrics_microbiome[,"balanced_accuracy"]) / sqrt(nrow(all_metrics_microbiome)))
  )
  
  cat("\n=== MICROBIOME CROSS-VALIDATION RESULTS ===\n")
  cat("Number of CV iterations:", nrow(all_metrics_microbiome), "\n")
  print(microbiome_metrics_summary)
  cat("\nFormatted for reporting:\n")
  cat("AUC =", round(microbiome_metrics_summary$Mean[1], 3), "±", round(microbiome_metrics_summary$SE[1], 3), "\n")
  cat("Accuracy =", round(microbiome_metrics_summary$Mean[2], 3), "±", round(microbiome_metrics_summary$SE[2], 3), "\n")
  cat("Balanced Accuracy =", round(microbiome_metrics_summary$Mean[3], 3), "±", round(microbiome_metrics_summary$SE[3], 3), "\n")
  
  # Store values for later use
  mean_auc_microbiome <- microbiome_metrics_summary$Mean[1]
  se_auc_microbiome <- microbiome_metrics_summary$SE[1]
  sd_auc_microbiome <- microbiome_metrics_summary$SD[1]
  ci_lower_microbiome <- microbiome_metrics_summary$CI_Lower[1]
  ci_upper_microbiome <- microbiome_metrics_summary$CI_Upper[1]
  
} else {
  cat("No stable genera found for cross-validation\n")
}

# ====================================
# SAVE RESULTS AND CREATE COMPARISON
# ====================================

# Create comprehensive results summary with all metrics
if(length(stable_genera) > 0 && exists("microbiome_metrics_summary")) {
  microbiome_performance_results <- data.frame(
    Model = "Microbiome_Stability_GLM",
    Method = "Stability-based genus selection",
    N_Features = n_features_to_use,
    AUC = round(microbiome_metrics_summary$Mean[1], 3),
    AUC_SE = round(microbiome_metrics_summary$SE[1], 3),
    AUC_SD = round(microbiome_metrics_summary$SD[1], 3),
    AUC_CI_Lower = round(microbiome_metrics_summary$CI_Lower[1], 3),
    AUC_CI_Upper = round(microbiome_metrics_summary$CI_Upper[1], 3),
    Accuracy = round(microbiome_metrics_summary$Mean[2], 3),
    Accuracy_SE = round(microbiome_metrics_summary$SE[2], 3),
    Accuracy_SD = round(microbiome_metrics_summary$SD[2], 3),
    Balanced_Accuracy = round(microbiome_metrics_summary$Mean[3], 3),
    Balanced_Accuracy_SE = round(microbiome_metrics_summary$SE[3], 3),
    Balanced_Accuracy_SD = round(microbiome_metrics_summary$SD[3], 3),
    Features = paste(cv_genera, collapse = "; "),
    Feature_Selection_Method = "Stability-based (≥50% selection)",
    CV_Iterations = nrow(all_metrics_microbiome)
  )
  
  # Save comprehensive results
  write.csv(microbiome_performance_results, "../results/microbiome_stability_performance_complete_balacc.csv", row.names = FALSE)
  cat("\nMicrobiome performance results saved to: results/microbiome_stability_performance_complete.csv\n")
  
  # Save detailed CV results for supplementary material
  detailed_cv_results_microbiome <- data.frame(
    Iteration = 1:nrow(all_metrics_microbiome),
    Model = "Microbiome_Stability_GLM",
    AUC = all_metrics_microbiome[,"auc"],
    Accuracy = all_metrics_microbiome[,"accuracy"],
    Balanced_Accuracy = all_metrics_microbiome[,"balanced_accuracy"]
  )
  
  write.csv(detailed_cv_results_microbiome, "../results/microbiome_detailed_cv_results_complete.csv", row.names = FALSE)
  cat("Detailed CV results saved to: results/microbiome_detailed_cv_results_complete.csv\n")
  
  # Also save stability table with clean names
  microbiome_stability_export <- microbiome_stability_results
  microbiome_stability_export$clean_name <- gsub("\\.", " ", microbiome_stability_results$feature)
  microbiome_stability_export <- microbiome_stability_export[, c("clean_name", "feature", "selected_times", "selection_rate")]
  colnames(microbiome_stability_export) <- c("Genus_Name", "Original_Feature_Name", "Times_Selected", "Selection_Rate")
  
  write.csv(microbiome_stability_export, "../results/microbiome_stability_results.csv", row.names = FALSE)
  cat("Stability table exported to: results/microbiome_stability_results.csv\n")
  
  # Save summary metrics table for easy viewing
  write.csv(microbiome_metrics_summary, "../results/microbiome_metrics_summary.csv", row.names = FALSE)
  cat("Metrics summary saved to: results/microbiome_metrics_summary.csv\n")
  
  # Print comparison format
  cat("\n=== FOR COMPARISON WITH METABOLOMICS ===\n")
  cat("Microbiome (Stability GLM):\n")
  cat("  Features:", paste(cv_genera, collapse = ", "), "\n")
  cat("  AUC: ", round(microbiome_metrics_summary$Mean[1], 3), " ± ", round(microbiome_metrics_summary$SE[1], 3), 
      " (95% CI: ", round(microbiome_metrics_summary$CI_Lower[1], 3), "-", round(microbiome_metrics_summary$CI_Upper[1], 3), ")\n", sep = "")
  cat("  Accuracy: ", round(microbiome_metrics_summary$Mean[2], 3), " ± ", round(microbiome_metrics_summary$SE[2], 3),
      " (95% CI: ", round(microbiome_metrics_summary$CI_Lower[2], 3), "-", round(microbiome_metrics_summary$CI_Upper[2], 3), ")\n", sep = "")
  cat("  Balanced Accuracy: ", round(microbiome_metrics_summary$Mean[3], 3), " ± ", round(microbiome_metrics_summary$SE[3], 3),
      " (95% CI: ", round(microbiome_metrics_summary$CI_Lower[3], 3), "-", round(microbiome_metrics_summary$CI_Upper[3], 3), ")\n", sep = "")
  cat("  Method: Stability-based feature selection (≥50% selection rate)\n")
  cat("  CV: 10-fold repeated 5 times (n =", nrow(all_metrics_microbiome), "iterations)\n")
}

# ====================================
# STATISTICAL SIGNIFICANCE TEST
# ====================================

# if(exists("all_metrics_microbiome")) {
#   cat("\n=== STATISTICAL SIGNIFICANCE (vs baseline) ===\n")
  
#   # Test AUC vs 0.5
#   t_test_auc <- t.test(all_metrics_microbiome[,"auc"], mu = 0.5, alternative = "greater")
#   sig_text_auc <- ifelse(t_test_auc$p.value < 0.05, "significant", "not significant")
#   cat(sprintf("AUC: %.3f