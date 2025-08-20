# ====================================
# METABOLOMICS 5-FOLD CROSS-VALIDATION WITH AUC ± SE
# ====================================

library(randomForest)
library(caret)
library(dplyr)
library(doParallel)
library(pROC)
library(Boruta)
library(glmnet)

# Set up parallel processing
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# ====================================
# LOAD DATA
# ====================================
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
cat("Loaded data:", nrow(fk_metabolom_gt_scaled), "samples x", 
    ncol(fk_metabolom_gt_scaled)-1, "features\n")

# Initialize results with AUC SE
results_df <- data.frame(
  Model = character(),
  Method = character(),
  N_Features = character(),
  Mean_Accuracy = numeric(),
  SD_Accuracy = numeric(),
  SE_Accuracy = numeric(),
  Mean_Balanced_Accuracy = numeric(),
  SD_Balanced_Accuracy = numeric(),
  SE_Balanced_Accuracy = numeric(),
  Mean_Kappa = numeric(),
  SD_Kappa = numeric(),
  SE_Kappa = numeric(),
  Mean_AUC = numeric(),
  SD_AUC = numeric(),
  SE_AUC = numeric(),
  AUC_CI_Lower = numeric(),
  AUC_CI_Upper = numeric(),
  stringsAsFactors = FALSE
)

# ====================================
# CREATE 5-FOLD SPLITS
# ====================================
set.seed(123)  # Same seed as original for reproducibility
folds <- createFolds(fk_metabolom_gt_scaled$gt, k = 5, returnTrain = FALSE)

cat("Created 5-fold CV splits:\n")
for(i in 1:5) {
  fold_data <- fk_metabolom_gt_scaled[folds[[i]], ]
  cat("Fold", i, ":", nrow(fold_data), "samples,", 
      table(fold_data$gt)[1], "GT_absent,", table(fold_data$gt)[2], "GT_present\n")
}

# ====================================
# HELPER FUNCTION FOR SINGLE FOLD EVALUATION
# ====================================
evaluate_single_fold <- function(model, test_data) {
  pred <- predict(model, test_data)
  cm <- confusionMatrix(pred, test_data$gt)
  
  pred_prob <- predict(model, test_data, type = "prob")
  roc_obj <- roc(test_data$gt, pred_prob[,2], quiet = TRUE)
  auc_val <- auc(roc_obj)
  
  # Handle NA values
  sensitivity <- ifelse(is.na(cm$byClass['Sensitivity']), 0, cm$byClass['Sensitivity'])
  specificity <- ifelse(is.na(cm$byClass['Specificity']), 0, cm$byClass['Specificity'])
  precision <- ifelse(is.na(cm$byClass['Pos Pred Value']), 0, cm$byClass['Pos Pred Value'])
  f1 <- ifelse(is.na(cm$byClass['F1']), 0, cm$byClass['F1'])
  balanced_acc <- ifelse(is.na(cm$byClass['Balanced Accuracy']), 0, cm$byClass['Balanced Accuracy'])
  
  return(list(
    Accuracy = as.numeric(cm$overall['Accuracy']),
    Balanced_Accuracy = as.numeric(balanced_acc),
    Sensitivity = as.numeric(sensitivity),
    Specificity = as.numeric(specificity),
    Precision = as.numeric(precision),
    F1_Score = as.numeric(f1),
    Kappa = as.numeric(cm$overall['Kappa']),
    AUC = as.numeric(auc_val)
  ))
}

# ====================================
# UPDATED FUNCTION TO AGGREGATE RESULTS WITH SE
# ====================================
aggregate_cv_results <- function(fold_results, model_name, method_name, n_features_info) {
  metrics <- c("Accuracy", "Balanced_Accuracy", "Kappa", "AUC")
  
  aggregated <- list()
  for(metric in metrics) {
    values <- sapply(fold_results, function(x) x[[metric]])
    n_folds <- length(values)
    
    # Calculate mean, SD, and SE
    mean_val <- mean(values, na.rm = TRUE)
    sd_val <- sd(values, na.rm = TRUE)
    se_val <- sd_val / sqrt(n_folds)
    
    aggregated[[paste0("Mean_", metric)]] <- round(mean_val, 3)
    aggregated[[paste0("SD_", metric)]] <- round(sd_val, 3)
    aggregated[[paste0("SE_", metric)]] <- round(se_val, 3)
    
    # 95% Confidence intervals for AUC
    if(metric == "AUC") {
      ci_lower <- mean_val - 1.96 * se_val
      ci_upper <- mean_val + 1.96 * se_val
      aggregated[["AUC_CI_Lower"]] <- round(ci_lower, 3)
      aggregated[["AUC_CI_Upper"]] <- round(ci_upper, 3)
    }
  }
  
  result_row <- data.frame(
    Model = model_name,
    Method = method_name,
    N_Features = n_features_info,
    Mean_Accuracy = aggregated$Mean_Accuracy,
    SD_Accuracy = aggregated$SD_Accuracy,
    SE_Accuracy = aggregated$SE_Accuracy,
    Mean_Balanced_Accuracy = aggregated$Mean_Balanced_Accuracy,
    SD_Balanced_Accuracy = aggregated$SD_Balanced_Accuracy,
    SE_Balanced_Accuracy = aggregated$SE_Balanced_Accuracy,
    Mean_Kappa = aggregated$Mean_Kappa,
    SD_Kappa = aggregated$SD_Kappa,
    SE_Kappa = aggregated$SE_Kappa,
    Mean_AUC = aggregated$Mean_AUC,
    SD_AUC = aggregated$SD_AUC,
    SE_AUC = aggregated$SE_AUC,
    AUC_CI_Lower = aggregated$AUC_CI_Lower,
    AUC_CI_Upper = aggregated$AUC_CI_Upper,
    stringsAsFactors = FALSE
  )
  
  return(result_row)
}

# ====================================
# MODEL 1: BORUTA FEATURE SELECTION + RF
# ====================================
cat("\n=== MODEL 1: BORUTA FEATURE SELECTION + RF ===\n")

boruta_fold_results <- list()
boruta_features_per_fold <- list()

for(fold in 1:5) {
  cat("\nFold", fold, ":\n")
  
  # Create train/test split
  test_indices <- folds[[fold]]
  train_data <- fk_metabolom_gt_scaled[-test_indices, ]
  test_data <- fk_metabolom_gt_scaled[test_indices, ]
  
  tryCatch({
    # Run Boruta on training data only
    set.seed(123 + fold)
    boruta_result <- Boruta(gt ~ ., data = train_data, 
                           maxRuns = 100, 
                           doTrace = 0)
    
    important_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
    boruta_features_per_fold[[fold]] <- important_features
    
    cat("  Boruta selected", length(important_features), "features\n")
    if(length(important_features) > 0) {
      cat("  Features:", paste(important_features[1:min(3, length(important_features))], collapse = ", "))
      if(length(important_features) > 3) cat(" ...")
      cat("\n")
    }
    
    if(length(important_features) > 0) {
      # Train RF with selected features
      train_boruta <- train_data[, c(important_features, "gt")]
      test_boruta <- test_data[, c(important_features, "gt")]
      
      set.seed(123 + fold)
      rf_boruta <- randomForest(
        gt ~ ., 
        data = train_boruta,
        mtry = max(1, floor(sqrt(length(important_features)))),
        ntree = 2000,
        importance = TRUE,
        nodesize = 3,
        replace = TRUE,
        classwt = c(1, 1)
      )
      
      boruta_fold_results[[fold]] <- evaluate_single_fold(rf_boruta, test_boruta)
      cat("  AUC:", round(boruta_fold_results[[fold]]$AUC, 3), 
          "| Balanced Accuracy:", round(boruta_fold_results[[fold]]$Balanced_Accuracy, 3), "\n")
      
    } else {
      # No features selected
      boruta_fold_results[[fold]] <- list(
        Accuracy = 0.5, Balanced_Accuracy = 0.5, Sensitivity = 0.5,
        Specificity = 0.5, Precision = 0.5, F1_Score = 0.5,
        Kappa = 0, AUC = 0.5
      )
      cat("  No features selected - using dummy results\n")
    }
    
  }, error = function(e) {
    cat("  Error in fold", fold, ":", e$message, "\n")
    boruta_fold_results[[fold]] <- list(
      Accuracy = 0.5, Balanced_Accuracy = 0.5, Sensitivity = 0.5,
      Specificity = 0.5, Precision = 0.5, F1_Score = 0.5,
      Kappa = 0, AUC = 0.5
    )
    boruta_features_per_fold[[fold]] <- character(0)
  })
}

# Aggregate Boruta results
feature_counts <- sapply(boruta_features_per_fold, length)
n_features_boruta <- paste0(round(mean(feature_counts), 1), " ± ", round(sd(feature_counts), 1), 
                           " (", min(feature_counts), "-", max(feature_counts), ")")

boruta_summary <- aggregate_cv_results(boruta_fold_results, "Boruta_RF", "5fold_CV", n_features_boruta)
results_df <- rbind(results_df, boruta_summary)

cat("\nBoruta RF Summary: AUC =", boruta_summary$Mean_AUC, "±", boruta_summary$SE_AUC, 
    "| Balanced Accuracy =", boruta_summary$Mean_Balanced_Accuracy, "±", boruta_summary$SE_Balanced_Accuracy, "\n")

# ====================================
# MODEL 2: PROPERLY TUNED RF
# ====================================
cat("\n=== MODEL 2: PROPERLY TUNED RF ===\n")

tuned_fold_results <- list()

for(fold in 1:5) {
  cat("\nFold", fold, ":\n")
  
  # Create train/test split
  test_indices <- folds[[fold]]
  train_data <- fk_metabolom_gt_scaled[-test_indices, ]
  test_data <- fk_metabolom_gt_scaled[test_indices, ]
  
  tryCatch({
    # Conservative tuning grid
    tuneGrid <- expand.grid(
      mtry = c(floor(sqrt(ncol(train_data)-1)),
               floor(log2(ncol(train_data)-1)),
               floor((ncol(train_data)-1)/10),
               floor((ncol(train_data)-1)/20))
    )
    tuneGrid <- unique(tuneGrid[tuneGrid$mtry > 0, , drop = FALSE])
    
    # Cross-validation within training set only
    ctrl <- trainControl(
      method = "cv",
      number = 3,  # Reduced due to smaller training set in each fold
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      savePredictions = "final",
      sampling = "up"
    )
    
    set.seed(123 + fold)
    rf_tuned <- train(
      gt ~ ., 
      data = train_data,
      method = "rf",
      tuneGrid = tuneGrid,
      trControl = ctrl,
      metric = "ROC",
      ntree = 1000,
      importance = TRUE,
      nodesize = 3,
      replace = TRUE
    )
    
    cat("  Best mtry:", rf_tuned$bestTune$mtry, "\n")
    
    # Evaluate on test set
    tuned_fold_results[[fold]] <- evaluate_single_fold(rf_tuned, test_data)
    cat("  AUC:", round(tuned_fold_results[[fold]]$AUC, 3), 
        "| Balanced Accuracy:", round(tuned_fold_results[[fold]]$Balanced_Accuracy, 3), "\n")
    
  }, error = function(e) {
    cat("  Error in fold", fold, ":", e$message, "\n")
    tuned_fold_results[[fold]] <- list(
      Accuracy = 0.5, Balanced_Accuracy = 0.5, Sensitivity = 0.5,
      Specificity = 0.5, Precision = 0.5, F1_Score = 0.5,
      Kappa = 0, AUC = 0.5
    )
  })
}

# Aggregate tuned RF results
tuned_summary <- aggregate_cv_results(tuned_fold_results, "Properly_Tuned_RF", "5fold_CV", "676")
results_df <- rbind(results_df, tuned_summary)

cat("\nTuned RF Summary: AUC =", tuned_summary$Mean_AUC, "±", tuned_summary$SE_AUC, 
    "| Balanced Accuracy =", tuned_summary$Mean_Balanced_Accuracy, "±", tuned_summary$SE_Balanced_Accuracy, "\n")

# ====================================
# MODEL 3: IMPROVED ANTIFUNGAL RF
# ====================================
cat("\n=== MODEL 3: IMPROVED ANTIFUNGAL RF ===\n")

# Define antifungal compounds
interesting_compounds <- c("Salicylate", "ISOLEUCINE", "Ornithine", "THREONINE", 
"ASPARTATE", "NICOTINAMIDE", "TRYPTOPHAN", "VALINE", "succinic acid", 
"LACTIC ACID", "Linalool", "Geranyl acetate", "trans-Nerolidol", "Nerol", 
"FERULATE", "CAFFEATE", "caffeine", "Cinnamic acid - 40.0 eV", 
"Vanillic acid", "BETAINE", "CHOLINE", "SARCOSINE", "N,N-Dimethylglycine", 
"3-HYDROXYBENZOATE", "4-HYDROXYBENZOATE", "Sinapic acid", 
"Hydroquinone", "TYROSINE", "trans-3-Hydroxycinnamic acid", 
"3,4-DIHYDROXYBENZOATE", "COUMARIN", "L-Arginine", "L-Proline", 
"PHENYLALANINE")

interesting_compounds <- make.names(interesting_compounds, unique = TRUE)
available_antifungal <- intersect(interesting_compounds, colnames(fk_metabolom_gt_scaled))
cat("Using", length(available_antifungal), "antifungal metabolites\n")

antifungal_fold_results <- list()

for(fold in 1:5) {
  cat("\nFold", fold, ":\n")
  
  # Create train/test split with antifungal features
  test_indices <- folds[[fold]]
  
  antifungal_data <- fk_metabolom_gt_scaled[, c(available_antifungal, "gt")]
  train_antifungal <- antifungal_data[-test_indices, ]
  test_antifungal <- antifungal_data[test_indices, ]
  
  # Convert to numeric
  predictor_cols <- names(train_antifungal)[names(train_antifungal) != "gt"]
  train_antifungal[predictor_cols] <- lapply(train_antifungal[predictor_cols], as.numeric)
  test_antifungal[predictor_cols] <- lapply(test_antifungal[predictor_cols], as.numeric)
  
  tryCatch({
    set.seed(123 + fold)
    rf_antifungal <- randomForest(
      gt ~ ., 
      data = train_antifungal,
      mtry = max(1, floor(sqrt(length(available_antifungal)))),
      ntree = 2000,
      importance = TRUE,
      nodesize = 2,
      replace = TRUE,
      classwt = c(1, 1)
    )
    
    antifungal_fold_results[[fold]] <- evaluate_single_fold(rf_antifungal, test_antifungal)
    cat("  AUC:", round(antifungal_fold_results[[fold]]$AUC, 3), 
        "| Balanced Accuracy:", round(antifungal_fold_results[[fold]]$Balanced_Accuracy, 3), "\n")
    
  }, error = function(e) {
    cat("  Error in fold", fold, ":", e$message, "\n")
    antifungal_fold_results[[fold]] <- list(
      Accuracy = 0.5, Balanced_Accuracy = 0.5, Sensitivity = 0.5,
      Specificity = 0.5, Precision = 0.5, F1_Score = 0.5,
      Kappa = 0, AUC = 0.5
    )
  })
}

# Aggregate antifungal RF results
antifungal_summary <- aggregate_cv_results(antifungal_fold_results, "Improved_Antifungal_RF", 
                                          "5fold_CV", as.character(length(available_antifungal)))
results_df <- rbind(results_df, antifungal_summary)

cat("\nAntifungal RF Summary: AUC =", antifungal_summary$Mean_AUC, "±", antifungal_summary$SE_AUC, 
    "| Balanced Accuracy =", antifungal_summary$Mean_Balanced_Accuracy, "±", antifungal_summary$SE_Balanced_Accuracy, "\n")

# ====================================
# DISPLAY FINAL RESULTS WITH AUC ± SE
# ====================================
cat("\n============================================================\n")
cat("SUMMARY OF 5-FOLD CROSS-VALIDATION RESULTS WITH AUC ± SE\n")
cat("============================================================\n\n")

# Sort by AUC
results_df <- results_df[order(results_df$Mean_AUC, decreasing = TRUE), ]

# Display results with AUC
print(results_df[, c("Model", "N_Features", "Mean_AUC", "SE_AUC", "AUC_CI_Lower", "AUC_CI_Upper",
                     "Mean_Balanced_Accuracy", "SE_Balanced_Accuracy")])

# Create comparison table with AUC ± SE
cat("\n=== COMPARISON TABLE WITH AUC ± SE ===\n")
comparison_table <- data.frame(
  Model = results_df$Model,
  N_Features = results_df$N_Features,
  AUC = paste0(results_df$Mean_AUC, " ± ", results_df$SE_AUC),
  AUC_CI = paste0("(", results_df$AUC_CI_Lower, "-", results_df$AUC_CI_Upper, ")"),
  Balanced_Accuracy = paste0(results_df$Mean_Balanced_Accuracy, " ± ", results_df$SE_Balanced_Accuracy),
  stringsAsFactors = FALSE
)

print(comparison_table)

# Save results
write.csv(results_df, "../results/metabolomics_5fold_cv_results_with_auc_se.csv", row.names = FALSE)

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

cat("\nResults saved to 'metabolomics_5fold_cv_results_with_auc_se.csv'\n")
cat("\n=== ANALYSIS COMPLETE ===\n")

#install.packages("patchwork")
library(patchwork)

# Create heatmap with all diagonal labels
p1 <- ggplot(cor_melted, aes(x = Metabolite, y = Genus, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 3)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "A",
       x = "Metabolites", y = "Genera")

# Create boxplots of combined protection scores
p2 <- ggplot(merged_df_combinedprotection, aes(x = gt.x, y = combined_protection, fill = gt.x)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_fill_manual(values = c("GT_absent" = "lightblue", "GT_present" = "lightcoral")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "B",
       x = "Take-all Status", 
       y = "Combined Protection Score")

# Combine plots
combined_plot <- p1 + p2 + plot_layout(widths = c(1, 1))

# Save plot
ggsave("../results/plots/heatmap_and_combined_prediction.png", 
       plot = combined_plot, 
       width = 12, height = 6, 
       dpi = 300)

summary(glm_combined)

# Calculate pseudo-R² values
#install.packages("DescTools")
library(DescTools)

# McFadden's pseudo-R²
mcfadden_r2 <- PseudoR2(glm_combined, which = "McFadden")

# Nagelkerke's pseudo-R²  
nagelkerke_r2 <- PseudoR2(glm_combined, which = "Nagelkerke")

# Or calculate McFadden's manually:
null_deviance <- glm_combined$null.deviance
residual_deviance <- glm_combined$deviance
mcfadden_manual <- 1 - (residual_deviance / null_deviance)

cat("McFadden's pseudo-R²:", round(mcfadden_manual, 3), "\n")


# Create summary dataframe of model comparisons
model_summary <- data.frame(
 Model = c("Combined Protection", "Metabolite Only", "Microbe Only", "Both Separate"),
 AIC = c(AIC(glm_combined), AIC(glm_metab), AIC(glm_microb), AIC(glm_both)),
 Delta_AIC = c(0, NA, NA, NA)
)

# Calculate delta AIC (difference from best model)
best_aic <- min(model_summary$AIC)
model_summary$Delta_AIC <- model_summary$AIC - best_aic

# Add pseudo-R² (McFadden's)
model_summary$Pseudo_R2 <- c(
 1 - (glm_combined$deviance / glm_combined$null.deviance),
 1 - (glm_metab$deviance / glm_metab$null.deviance),
 1 - (glm_microb$deviance / glm_microb$null.deviance),
 1 - (glm_both$deviance / glm_both$null.deviance)
)

# Round values
model_summary$AIC <- round(model_summary$AIC, 2)
model_summary$Delta_AIC <- round(model_summary$Delta_AIC, 2)
model_summary$Pseudo_R2 <- round(model_summary$Pseudo_R2, 3)

# Display table
print(model_summary)

write.csv(model_summary, "../supplementary/borutamodel_AIC_comparison")