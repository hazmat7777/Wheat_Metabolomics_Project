# ====================================
# METABOLOMICS 5-FOLD CROSS-VALIDATION
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

# Initialize results
results_df <- data.frame(
  Model = character(),
  Method = character(),
  N_Features = character(),  # Changed to character to handle variable numbers
  Mean_Accuracy = numeric(),
  SD_Accuracy = numeric(),
  Min_Accuracy = numeric(),
  Max_Accuracy = numeric(),
  Mean_Balanced_Accuracy = numeric(),
  SD_Balanced_Accuracy = numeric(),
  Min_Balanced_Accuracy = numeric(),
  Max_Balanced_Accuracy = numeric(),
  Mean_Kappa = numeric(),
  SD_Kappa = numeric(),
  Mean_AUC = numeric(),
  SD_AUC = numeric(),
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
# FUNCTION TO AGGREGATE RESULTS
# ====================================
aggregate_cv_results <- function(fold_results, model_name, method_name, n_features_info) {
  metrics <- c("Accuracy", "Balanced_Accuracy", "Kappa", "AUC")
  
  aggregated <- list()
  for(metric in metrics) {
    values <- sapply(fold_results, function(x) x[[metric]])
    aggregated[[paste0("Mean_", metric)]] <- round(mean(values, na.rm = TRUE), 3)
    aggregated[[paste0("SD_", metric)]] <- round(sd(values, na.rm = TRUE), 3)
    aggregated[[paste0("Min_", metric)]] <- round(min(values, na.rm = TRUE), 3)
    aggregated[[paste0("Max_", metric)]] <- round(max(values, na.rm = TRUE), 3)
  }
  
  result_row <- data.frame(
    Model = model_name,
    Method = method_name,
    N_Features = n_features_info,
    Mean_Accuracy = aggregated$Mean_Accuracy,
    SD_Accuracy = aggregated$SD_Accuracy,
    Min_Accuracy = aggregated$Min_Accuracy,
    Max_Accuracy = aggregated$Max_Accuracy,
    Mean_Balanced_Accuracy = aggregated$Mean_Balanced_Accuracy,
    SD_Balanced_Accuracy = aggregated$SD_Balanced_Accuracy,
    Min_Balanced_Accuracy = aggregated$Min_Balanced_Accuracy,
    Max_Balanced_Accuracy = aggregated$Max_Balanced_Accuracy,
    Mean_Kappa = aggregated$Mean_Kappa,
    SD_Kappa = aggregated$SD_Kappa,
    Mean_AUC = aggregated$Mean_AUC,
    SD_AUC = aggregated$SD_AUC,
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
      cat("  Balanced Accuracy:", round(boruta_fold_results[[fold]]$Balanced_Accuracy, 3), "\n")
      
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

cat("\nBoruta RF Summary:", boruta_summary$Mean_Balanced_Accuracy, "±", boruta_summary$SD_Balanced_Accuracy, "\n")

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
    cat("  Balanced Accuracy:", round(tuned_fold_results[[fold]]$Balanced_Accuracy, 3), "\n")
    
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

cat("\nTuned RF Summary:", tuned_summary$Mean_Balanced_Accuracy, "±", tuned_summary$SD_Balanced_Accuracy, "\n")

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
    cat("  Balanced Accuracy:", round(antifungal_fold_results[[fold]]$Balanced_Accuracy, 3), "\n")
    
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

cat("\nAntifungal RF Summary:", antifungal_summary$Mean_Balanced_Accuracy, "±", antifungal_summary$SD_Balanced_Accuracy, "\n")

# ====================================
# ANALYZE BORUTA FEATURE CONSISTENCY
# ====================================
cat("\n============================================================\n")
cat("BORUTA FEATURE SELECTION ANALYSIS\n")
cat("============================================================\n\n")

# Collect all confirmed features across folds
all_confirmed_features <- unlist(boruta_features_per_fold)

if(length(all_confirmed_features) > 0) {
  cat("=== CONFIRMED FEATURES FREQUENCY ===\n")
  confirmed_freq <- table(all_confirmed_features)
  confirmed_freq <- sort(confirmed_freq, decreasing = TRUE)
  print(confirmed_freq)
  
  cat("\n=== TOP CONFIRMED FEATURES ===\n")
  for(i in 1:min(5, length(confirmed_freq))) {
    feature_name <- names(confirmed_freq)[i]
    frequency <- confirmed_freq[i]
    cat(sprintf("%-40s: %d/5 folds (%.1f%%)\n", 
                feature_name, frequency, frequency/5*100))
  }
  
  # Check specifically for X3.4.Dihydroxybenzaldehyde
  target_feature <- "X3.4.Dihydroxybenzaldehyde"
  if(target_feature %in% names(confirmed_freq)) {
    cat(sprintf("\n*** %s was selected in %d/5 folds (%.1f%%) ***\n", 
                target_feature, confirmed_freq[target_feature], 
                confirmed_freq[target_feature]/5*100))
  } else {
    cat(sprintf("\n*** %s was NOT selected in any fold ***\n", target_feature))
  }
} else {
  cat("No confirmed features found in any fold.\n")
}

# Show feature selection summary for each fold
cat("\n=== FEATURE SELECTION PER FOLD ===\n")
for(i in 1:5) {
  features <- boruta_features_per_fold[[i]]
  cat(sprintf("Fold %d: %d features", i, length(features)))
  
  if(length(features) > 0) {
    cat(" | Features: ", paste(features[1:min(2, length(features))], collapse = ", "))
    if(length(features) > 2) cat(" ...")
  }
  cat("\n")
}

# ====================================
# DISPLAY FINAL RESULTS
# ====================================
cat("\n============================================================\n")
cat("SUMMARY OF 5-FOLD CROSS-VALIDATION RESULTS\n")
cat("============================================================\n\n")

# Sort by balanced accuracy
results_df <- results_df[order(results_df$Mean_Balanced_Accuracy, decreasing = TRUE), ]

# Display results
print(results_df[, c("Model", "N_Features", "Mean_Balanced_Accuracy", "SD_Balanced_Accuracy", 
                     "Mean_Accuracy", "SD_Accuracy", "Mean_Kappa", "SD_Kappa")])

# Create comparison table
cat("\n=== COMPARISON TABLE (Mean ± SD) ===\n")
comparison_table <- data.frame(
  Model = results_df$Model,
  N_Features = results_df$N_Features,
  Balanced_Accuracy = paste0(results_df$Mean_Balanced_Accuracy, " ± ", results_df$SD_Balanced_Accuracy),
  Accuracy = paste0(results_df$Mean_Accuracy, " ± ", results_df$SD_Accuracy),
  Kappa = paste0(results_df$Mean_Kappa, " ± ", results_df$SD_Kappa),
  stringsAsFactors = FALSE
)

print(comparison_table)

# Save results
write.csv(results_df, "metabolomics_5fold_cv_results.csv", row.names = FALSE)

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

cat("\nResults saved to 'metabolomics_5fold_cv_results.csv'\n")
cat("\n=== ANALYSIS COMPLETE ===\n")



#### this is its own script
# run boruta algo
boruta_result <- Boruta(gt ~ ., data = train_data, 
                        maxRuns = 100, 
                        doTrace = 0)
# Use Boruta features in logistic regression
selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
selected_features

# More liberal statistical threshold
all_selected_features <- c()

for(run in 1:10) {
  set.seed(100 + run)  # Different seed each time
  boruta_result <- Boruta(gt ~ ., data = train_data, 
                         maxRuns = 100, 
                         doTrace = 0)
  
  selected_features <- getSelectedAttributes(boruta_result, withTentative = TRUE)
  all_selected_features <- c(all_selected_features, selected_features)
  
  cat("Run", run, ":", length(selected_features), "features\n")
}



formula <- as.formula(paste("gt ~", paste(all_selected_features, collapse = " + ")))
logistic_model <- glm(formula, data = train_data, family = binomial)

simple_model <- glm(gt ~ Haematommic.acid + Spectral.Match.to.Pinolenic.acid.from.NIST14, 
                   data = train_data, family = binomial)
summary(simple_model)

# boxplot for this



summary(logistic_model)
# This would tell you:
# - Which features have positive/negative associations
# - Relative importance (coefficients)
# - Statistical significance of each