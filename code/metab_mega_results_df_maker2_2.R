# ====================================
# IMPROVED METABOLOMICS RANDOM FOREST MODELS
# ====================================

library(randomForest)
library(caret)
library(dplyr)
library(doParallel)
library(pROC)
#install.packages("Boruta", dependencies = TRUE)
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
  N_Features = integer(),
  Accuracy = numeric(),
  Balanced_Accuracy = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  Precision = numeric(),
  F1_Score = numeric(),
  Kappa = numeric(),
  AUC = numeric(),
  stringsAsFactors = FALSE
)

# ====================================
# BETTER TRAIN/TEST SPLIT WITH STRATIFICATION
# ====================================
set.seed(123)

# Use stratified sampling to ensure balanced classes in both sets
train_idx <- createDataPartition(fk_metabolom_gt_scaled$gt, p = 0.75, list = FALSE)  # Increased to 75%
train_data <- fk_metabolom_gt_scaled[train_idx, ]
test_data <- fk_metabolom_gt_scaled[-train_idx, ]

cat("Train/Test split: ", nrow(train_data), "/", nrow(test_data), "\n")
cat("Train class dist:", table(train_data$gt), "\n")
cat("Test class dist:", table(test_data$gt), "\n")

# ====================================
# HELPER FUNCTION FOR EVALUATION
# ====================================
evaluate_model <- function(model, test_data, model_name, method_name, n_features) {
  pred <- predict(model, test_data)
  cm <- confusionMatrix(pred, test_data$gt)
  
  pred_prob <- predict(model, test_data, type = "prob")
  roc_obj <- roc(test_data$gt, pred_prob[,2])
  auc_val <- auc(roc_obj)
  
  # Handle NA values in confusion matrix
  sensitivity <- ifelse(is.na(cm$byClass['Sensitivity']), 0, cm$byClass['Sensitivity'])
  specificity <- ifelse(is.na(cm$byClass['Specificity']), 0, cm$byClass['Specificity'])
  precision <- ifelse(is.na(cm$byClass['Pos Pred Value']), 0, cm$byClass['Pos Pred Value'])
  f1 <- ifelse(is.na(cm$byClass['F1']), 0, cm$byClass['F1'])
  balanced_acc <- ifelse(is.na(cm$byClass['Balanced Accuracy']), 0, cm$byClass['Balanced Accuracy'])
  
  data.frame(
    Model = model_name,
    Method = method_name,
    N_Features = n_features,
    Accuracy = round(cm$overall['Accuracy'], 3),
    Balanced_Accuracy = round(balanced_acc, 3),
    Sensitivity = round(sensitivity, 3),
    Specificity = round(specificity, 3),
    Precision = round(precision, 3),
    F1_Score = round(f1, 3),
    Kappa = round(cm$overall['Kappa'], 3),
    AUC = round(auc_val, 3),
    stringsAsFactors = FALSE
  )
}

# ====================================
# MODEL 1: FEATURE SELECTION WITH BORUTA
# ====================================
cat("\n=== MODEL 1: BORUTA FEATURE SELECTION + RF ===\n")

# Run Boruta for feature selection (on training data only)
set.seed(123)
boruta_result <- Boruta(gt ~ ., data = train_data, 
                       maxRuns = 100, 
                       doTrace = 2)

# Get confirmed important features
important_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
cat("Boruta selected", length(important_features), "important features\n")

if(length(important_features) > 0) {
  # Create dataset with selected features
  train_boruta <- train_data[, c(important_features, "gt")]
  test_boruta <- test_data[, c(important_features, "gt")]
  
  # Train RF with selected features
  set.seed(123)
  rf_boruta <- randomForest(
    gt ~ ., 
    data = train_boruta,
    mtry = max(1, floor(sqrt(length(important_features)))),  # Much more conservative
    ntree = 2000,  # More trees for stability
    importance = TRUE,
    nodesize = 3,  # Smaller node size for small dataset
    replace = TRUE,
    classwt = c(1, 1)  # Balanced class weights
  )
  
  # Evaluate
  result <- evaluate_model(rf_boruta, test_boruta, "Boruta_RF", 
                          paste0("Boruta_", length(important_features), "_features"), 
                          length(important_features))
  results_df <- rbind(results_df, result)
  cat("Boruta RF - Balanced Accuracy:", result$Balanced_Accuracy, "\n")
}

# ====================================
# MODEL 2: LASSO FEATURE SELECTION + RF
# ====================================
cat("\n=== MODEL 2: LASSO FEATURE SELECTION + RF ===\n")

# Prepare data for LASSO
X_train <- as.matrix(train_data[, !names(train_data) %in% "gt"])
y_train <- as.factor(train_data$gt)

# LASSO with cross-validation
set.seed(123)
lasso_cv <- cv.glmnet(X_train, y_train, 
                     family = "binomial", 
                     alpha = 1, 
                     nfolds = 10,
                     type.measure = "class")

# Get features selected by LASSO
lasso_coef <- coef(lasso_cv, s = "lambda.1se")
selected_features <- rownames(lasso_coef)[which(lasso_coef[,1] != 0)]
selected_features <- selected_features[selected_features != "(Intercept)"]

cat("LASSO selected", length(selected_features), "features\n")

if(length(selected_features) > 0) {
  # Create dataset with LASSO-selected features
  train_lasso <- train_data[, c(selected_features, "gt")]
  test_lasso <- test_data[, c(selected_features, "gt")]
  
  # Train RF with LASSO-selected features
  set.seed(123)
  rf_lasso <- randomForest(
    gt ~ ., 
    data = train_lasso,
    mtry = max(1, floor(sqrt(length(selected_features)))),
    ntree = 2000,
    importance = TRUE,
    nodesize = 3,
    replace = TRUE,
    classwt = c(1, 1)
  )
  
  # Evaluate
  result <- evaluate_model(rf_lasso, test_lasso, "LASSO_RF", 
                          paste0("LASSO_", length(selected_features), "_features"), 
                          length(selected_features))
  results_df <- rbind(results_df, result)
  cat("LASSO RF - Balanced Accuracy:", result$Balanced_Accuracy, "\n")
}

# ====================================
# MODEL 3: PROPERLY TUNED RF
# ====================================
cat("\n=== MODEL 3: PROPERLY TUNED RF ===\n")

# Much more conservative tuning grid
tuneGrid <- expand.grid(
  mtry = c(floor(sqrt(ncol(train_data)-1)),          # ~26 features
           floor(log2(ncol(train_data)-1)),           # ~9 features  
           floor((ncol(train_data)-1)/10),            # ~67 features
           floor((ncol(train_data)-1)/20))            # ~34 features
)

tuneGrid <- unique(tuneGrid[tuneGrid$mtry > 0, , drop = FALSE])

# Cross-validation with more folds for small dataset
ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,           # Fewer folds due to small sample size
  repeats = 3,          # Multiple repeats for stability
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  sampling = "up"       # Upsample minority class
)

cat("Running conservative hyperparameter tuning...\n")
set.seed(123)
rf_tuned_proper <- train(
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

best_mtry <- rf_tuned_proper$bestTune$mtry
cat("Best mtry:", best_mtry, "(vs", ncol(train_data)-1, "total features)\n")

# Evaluate
pred_tuned <- predict(rf_tuned_proper, test_data)
result <- evaluate_model(rf_tuned_proper, test_data, "Properly_Tuned_RF", 
                        paste0("CV_Tuned_mtry", best_mtry), 
                        ncol(train_data)-1)
results_df <- rbind(results_df, result)
cat("Properly Tuned RF - Balanced Accuracy:", result$Balanced_Accuracy, "\n")

# ====================================
# MODEL 4: IMPROVED ANTIFUNGAL RF
# ====================================
cat("\n=== MODEL 4: IMPROVED ANTIFUNGAL RF ===\n")

# Your antifungal compounds
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

# Create improved antifungal model
antifungal_data <- fk_metabolom_gt_scaled[, c(available_antifungal, "gt")]
train_antifungal <- antifungal_data[train_idx, ]
test_antifungal <- antifungal_data[-train_idx, ]

# Convert to numeric
predictor_cols <- names(train_antifungal)[names(train_antifungal) != "gt"]
train_antifungal[predictor_cols] <- lapply(train_antifungal[predictor_cols], as.numeric)
test_antifungal[predictor_cols] <- lapply(test_antifungal[predictor_cols], as.numeric)

set.seed(123)
rf_antifungal_improved <- randomForest(
  gt ~ ., 
  data = train_antifungal,
  mtry = max(1, floor(sqrt(length(available_antifungal)))),  # Conservative mtry
  ntree = 2000,
  importance = TRUE,
  nodesize = 2,  # Even smaller for the smaller feature set
  replace = TRUE,
  classwt = c(1, 1)
)

result <- evaluate_model(rf_antifungal_improved, test_antifungal, "Improved_Antifungal_RF", 
                        "Selected_Compounds_Tuned", 
                        length(available_antifungal))
results_df <- rbind(results_df, result)
cat("Improved Antifungal RF - Balanced Accuracy:", result$Balanced_Accuracy, "\n")

# ====================================
# MODEL 5: ENSEMBLE OF TOP FEATURES
# ====================================
cat("\n=== MODEL 5: ENSEMBLE OF TOP IMPORTANCE FEATURES ===\n")

# Get feature importance from the properly tuned model
importance_scores <- importance(rf_tuned_proper$finalModel)
top_features <- rownames(importance_scores)[order(importance_scores[, "MeanDecreaseGini"], decreasing = TRUE)][1:50]

# Create ensemble dataset
train_ensemble <- train_data[, c(top_features, "gt")]
test_ensemble <- test_data[, c(top_features, "gt")]

set.seed(123)
rf_ensemble <- randomForest(
  gt ~ ., 
  data = train_ensemble,
  mtry = max(1, floor(sqrt(50))),  # ~7 features
  ntree = 2000,
  importance = TRUE,
  nodesize = 2,
  replace = TRUE,
  classwt = c(1, 1)
)

result <- evaluate_model(rf_ensemble, test_ensemble, "Top_Features_RF", 
                        "Top50_by_Importance", 50)
results_df <- rbind(results_df, result)
cat("Top Features RF - Balanced Accuracy:", result$Balanced_Accuracy, "\n")

# ====================================
# DISPLAY RESULTS
# ====================================
cat("\n============================================================\n")
cat("IMPROVED METABOLOMICS RF MODELS COMPARISON\n")
cat("============================================================\n\n")

results_df <- results_df[order(results_df$Balanced_Accuracy, decreasing = TRUE), ]
print(results_df)

# Summary
cat("\n=== PERFORMANCE SUMMARY ===\n")
for(i in 1:nrow(results_df)) {
  cat(sprintf("%d. %-25s: %3d features, Bal.Acc = %.3f, AUC = %.3f\n",
              i,
              results_df$Model[i],
              results_df$N_Features[i],
              results_df$Balanced_Accuracy[i],
              results_df$AUC[i]))
}

# Save results
write.csv(results_df, "improved_metabolomics_rf_results.csv", row.names = FALSE)

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

cat("\n=== ANALYSIS COMPLETE ===\n")



# After running Boruta, check what feature was selected
important_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
print(important_features)

# Or if you want to see the full Boruta results with scores
print(boruta_result)

# To see a summary of feature importance decisions
plot(boruta_result, las = 2, cex.axis = 0.7)