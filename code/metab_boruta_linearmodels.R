# ====================================
# TEST: BORUTA FEATURES + LINEAR MODELS
# ====================================

library(randomForest)
library(caret)
library(Boruta)
library(pROC)
library(MASS)

# Load data
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")

# Create same 5-fold splits as before
set.seed(123)
folds <- createFolds(fk_metabolom_gt_scaled$gt, k = 5, returnTrain = FALSE)

# Initialize results
linear_results <- list()

# ====================================
# TEST LINEAR MODELS WITH BORUTA FEATURES
# ====================================

for(fold in 1:5) {
  cat("\n=== FOLD", fold, "===\n")
  
  # Create train/test split
  test_indices <- folds[[fold]]
  train_data <- fk_metabolom_gt_scaled[-test_indices, ]
  test_data <- fk_metabolom_gt_scaled[test_indices, ]
  
  # Run Boruta on training data
  set.seed(123 + fold)
  boruta_result <- Boruta(gt ~ ., data = train_data, 
                         maxRuns = 100, 
                         doTrace = 0)
  
  important_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
  cat("Boruta selected", length(important_features), "features:", 
      paste(important_features[1:min(3, length(important_features))], collapse = ", "))
  if(length(important_features) > 3) cat(" ...")
  cat("\n")
  
  if(length(important_features) > 0) {
    # Prepare data with selected features
    train_selected <- train_data[, c(important_features, "gt")]
    test_selected <- test_data[, c(important_features, "gt")]
    
    # Create formula
    formula_str <- paste("gt ~", paste(important_features, collapse = " + "))
    formula_obj <- as.formula(formula_str)
    
    fold_results <- list()
    
    # ====================================
    # MODEL 1: LOGISTIC REGRESSION
    # ====================================
    tryCatch({
      logistic_model <- glm(formula_obj, data = train_selected, family = binomial)
      
      # Predict on test set
      logistic_pred_prob <- predict(logistic_model, test_selected, type = "response")
      logistic_pred <- factor(ifelse(logistic_pred_prob > 0.5, "GT_present", "GT_absent"),
                             levels = levels(test_selected$gt))
      
      # Evaluate
      cm_logistic <- confusionMatrix(logistic_pred, test_selected$gt)
      auc_logistic <- auc(roc(test_selected$gt, logistic_pred_prob, quiet = TRUE))
      
      fold_results$logistic <- list(
        accuracy = as.numeric(cm_logistic$overall['Accuracy']),
        balanced_accuracy = as.numeric(cm_logistic$byClass['Balanced Accuracy']),
        kappa = as.numeric(cm_logistic$overall['Kappa']),
        auc = as.numeric(auc_logistic),
        coefficients = coef(logistic_model)
      )
      
      cat("Logistic Regression - Balanced Acc:", round(fold_results$logistic$balanced_accuracy, 3), "\n")
      
    }, error = function(e) {
      cat("Logistic regression failed:", e$message, "\n")
      fold_results$logistic <- NULL
    })
    
    # ====================================
    # MODEL 2: LINEAR DISCRIMINANT ANALYSIS
    # ====================================
    tryCatch({
      lda_model <- lda(formula_obj, data = train_selected)
      
      # Predict on test set
      lda_pred <- predict(lda_model, test_selected)
      
      # Evaluate
      cm_lda <- confusionMatrix(lda_pred$class, test_selected$gt)
      auc_lda <- auc(roc(test_selected$gt, lda_pred$posterior[,2], quiet = TRUE))
      
      fold_results$lda <- list(
        accuracy = as.numeric(cm_lda$overall['Accuracy']),
        balanced_accuracy = as.numeric(cm_lda$byClass['Balanced Accuracy']),
        kappa = as.numeric(cm_lda$overall['Kappa']),
        auc = as.numeric(auc_lda)
      )
      
      cat("LDA - Balanced Acc:", round(fold_results$lda$balanced_accuracy, 3), "\n")
      
    }, error = function(e) {
      cat("LDA failed:", e$message, "\n")
      fold_results$lda <- NULL
    })
    
    # ====================================
    # MODEL 3: SIMPLE METABOLITE SCORE
    # ====================================
    tryCatch({
      # Calculate simple mean of standardized features
      train_means <- sapply(train_selected[important_features], mean)
      train_sds <- sapply(train_selected[important_features], sd)
      
      # Standardize features
      train_std <- train_selected
      test_std <- test_selected
      
      for(feature in important_features) {
        train_std[[feature]] <- (train_selected[[feature]] - train_means[feature]) / train_sds[feature]
        test_std[[feature]] <- (test_selected[[feature]] - train_means[feature]) / train_sds[feature]
      }
      
      # Create simple score (mean of standardized features)
      train_score <- rowMeans(train_std[important_features])
      test_score <- rowMeans(test_std[important_features])
      
      # Find optimal threshold on training set
      train_roc <- roc(train_selected$gt, train_score, quiet = TRUE)
      optimal_threshold <- coords(train_roc, "best", ret = "threshold")$threshold
      
      # Apply threshold to test set
      score_pred <- factor(ifelse(test_score > optimal_threshold, "GT_present", "GT_absent"),
                          levels = levels(test_selected$gt))
      
      # Evaluate
      cm_score <- confusionMatrix(score_pred, test_selected$gt)
      auc_score <- auc(roc(test_selected$gt, test_score, quiet = TRUE))
      
      fold_results$score <- list(
        accuracy = as.numeric(cm_score$overall['Accuracy']),
        balanced_accuracy = as.numeric(cm_score$byClass['Balanced Accuracy']),
        kappa = as.numeric(cm_score$overall['Kappa']),
        auc = as.numeric(auc_score),
        threshold = optimal_threshold
      )
      
      cat("Metabolite Score - Balanced Acc:", round(fold_results$score$balanced_accuracy, 3), "\n")
      
    }, error = function(e) {
      cat("Metabolite score failed:", e$message, "\n")
      fold_results$score <- NULL
    })
    
    # ====================================
    # MODEL 4: RANDOM FOREST (for comparison)
    # ====================================
    tryCatch({
      set.seed(123 + fold)
      rf_model <- randomForest(
        formula_obj, 
        data = train_selected,
        mtry = max(1, floor(sqrt(length(important_features)))),
        ntree = 2000,
        importance = TRUE,
        nodesize = 3,
        replace = TRUE,
        classwt = c(1, 1)
      )
      
      # Predict on test set
      rf_pred <- predict(rf_model, test_selected)
      rf_pred_prob <- predict(rf_model, test_selected, type = "prob")
      
      # Evaluate
      cm_rf <- confusionMatrix(rf_pred, test_selected$gt)
      auc_rf <- auc(roc(test_selected$gt, rf_pred_prob[,2], quiet = TRUE))
      
      fold_results$rf <- list(
        accuracy = as.numeric(cm_rf$overall['Accuracy']),
        balanced_accuracy = as.numeric(cm_rf$byClass['Balanced Accuracy']),
        kappa = as.numeric(cm_rf$overall['Kappa']),
        auc = as.numeric(auc_rf)
      )
      
      cat("Random Forest - Balanced Acc:", round(fold_results$rf$balanced_accuracy, 3), "\n")
      
    }, error = function(e) {
      cat("Random Forest failed:", e$message, "\n")
      fold_results$rf <- NULL
    })
    
  } else {
    cat("No features selected - skipping this fold\n")
    fold_results <- list(logistic = NULL, lda = NULL, score = NULL, rf = NULL)
  }
  
  linear_results[[fold]] <- fold_results
}

# ====================================
# AGGREGATE RESULTS
# ====================================
cat("\n============================================================\n")
cat("COMPARISON: BORUTA FEATURES + DIFFERENT MODELS\n")
cat("============================================================\n\n")

models <- c("logistic", "lda", "score", "rf")
model_names <- c("Logistic Regression", "Linear Discriminant Analysis", 
                "Metabolite Score", "Random Forest")

comparison_results <- data.frame(
  Model = model_names,
  Mean_Balanced_Accuracy = numeric(4),
  SD_Balanced_Accuracy = numeric(4),
  Mean_AUC = numeric(4),
  SD_AUC = numeric(4),
  N_Successful_Folds = integer(4),
  stringsAsFactors = FALSE
)

for(i in 1:4) {
  model <- models[i]
  
  # Extract results across folds
  balanced_accs <- c()
  aucs <- c()
  
  for(fold in 1:5) {
    if(!is.null(linear_results[[fold]][[model]])) {
      balanced_accs <- c(balanced_accs, linear_results[[fold]][[model]]$balanced_accuracy)
      aucs <- c(aucs, linear_results[[fold]][[model]]$auc)
    }
  }
  
  # Calculate summary statistics
  if(length(balanced_accs) > 0) {
    comparison_results$Mean_Balanced_Accuracy[i] <- round(mean(balanced_accs, na.rm = TRUE), 3)
    comparison_results$SD_Balanced_Accuracy[i] <- round(sd(balanced_accs, na.rm = TRUE), 3)
    comparison_results$Mean_AUC[i] <- round(mean(aucs, na.rm = TRUE), 3)
    comparison_results$SD_AUC[i] <- round(sd(aucs, na.rm = TRUE), 3)
    comparison_results$N_Successful_Folds[i] <- length(balanced_accs)
  } else {
    comparison_results$Mean_Balanced_Accuracy[i] <- NA
    comparison_results$SD_Balanced_Accuracy[i] <- NA
    comparison_results$Mean_AUC[i] <- NA
    comparison_results$SD_AUC[i] <- NA
    comparison_results$N_Successful_Folds[i] <- 0
  }
}

# Display results
print(comparison_results)

cat("\n=== SUMMARY ===\n")
for(i in 1:4) {
  if(!is.na(comparison_results$Mean_Balanced_Accuracy[i])) {
    cat(sprintf("%-25s: %.3f ± %.3f (AUC: %.3f ± %.3f) [%d/5 folds]\n",
                comparison_results$Model[i],
                comparison_results$Mean_Balanced_Accuracy[i],
                comparison_results$SD_Balanced_Accuracy[i],
                comparison_results$Mean_AUC[i],
                comparison_results$SD_AUC[i],
                comparison_results$N_Successful_Folds[i]))
  }
}

# Save results
write.csv(comparison_results, "boruta_linear_models_comparison.csv", row.names = FALSE)
cat("\nResults saved to 'boruta_linear_models_comparison.csv'\n")
cat("\n=== ANALYSIS COMPLETE ===\n")












# Simple t-tests to confirm directions
wilcox.test(Haematommic.acid ~ gt, data = fk_metabolom_gt_scaled)
wilcox.test(Spectral.Match.to.Pinolenic.acid.from.NIST14 ~ gt, data = fk_metabolom_gt_scaled)
wilcox.test(X2..4.Hydroxyphenyl.propionic.Acid ~ gt, data = fk_metabolom_gt_scaled)




# Plot the distributions
boxplot(Haematommic.acid ~ gt, data = fk_metabolom_gt_scaled, 
        main = "Haematommic acid levels by GT status")
boxplot(Spectral.Match.to.Pinolenic.acid.from.NIST14 ~ gt, data = fk_metabolom_gt_scaled, 
        main = "Spectral Match to Pinolenic acid by GT status")
boxplot(X2..4.Hydroxyphenyl.propionic.Acid ~ gt, data = fk_metabolom_gt_scaled, 
        main = "4-Hydroxyphenyl propionic Acid by GT status")

install.packages("effsize")
library(effsize)
cohen.d(fk_metabolom_gt_scaled$Haematommic.acid, fk_metabolom_gt_scaled$gt)


# First create the protective score
protective_metabolites <- c("Haematommic.acid", 
                           "Spectral.Match.to.Pinolenic.acid.from.NIST14", 
                           "X2..4.Hydroxyphenyl.propionic.Acid")

protective_score <- rowSums(scale(fk_metabolom_gt_scaled[, protective_metabolites]))

# Add it to your dataframe
fk_metabolom_gt_scaled$protective_score <- protective_score

# For binary classification, use logistic regression (glm)
logistic_model <- glm(gt ~ protective_score, 
                     data = fk_metabolom_gt_scaled, 
                     family = binomial)

summary(logistic_model)

# make a boxplot of the protective score by GT status
ggplot(fk_metabolom_gt_scaled, aes(x = gt, y = protective_score, fill = gt)) +
  geom_boxplot() +
  labs(title = "Protective Metabolite Score by GT Status",
       x = "GT Status",
       y = "Protective Metabolite Score") +
  theme_minimal() +
  scale_fill_manual(values = c("GT_absent" = "#E69F00", "GT_present" = "#56B4E9")) +
  theme(legend.position = "none")