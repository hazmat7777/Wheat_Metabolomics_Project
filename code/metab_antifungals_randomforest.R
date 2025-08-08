
# ONLY LOOKING AT THE ANTIFUNGAL METABOLITES

library(randomForest)
library(caret)
library(dplyr)

# ===============================
# ANTIFUNGAL METABOLITES MODEL
# ===============================

# Define the antifungal metabolites identified from your list
# Based on the analysis, these are the compounds with documented antifungal properties
antifungal_metabolites <- c(
  "Salicylate",                    # [1] - keratolytic antifungal
  "berberine",                     # [630] - broad spectrum antifungal
  "VANILLIN",                     # [339] - demonstrated antifungal activity
  "caffeine",                     # [485] - direct antifungal vs Candida
  "COUMARIN",                     # [674] - antifungal vs dermatophytes
  "TAURINE",                      # [619] - antifungal properties
  "NICOTINAMIDE",                 # [81] - B-vitamin with antifungal activity
  "THYMIDINE",                    # [214] - nucleoside with antifungal properties
  "QUINOLINE",                    # [498] - scaffold in antifungal medications
  "Usnic.acid",                   # [573] - lichen metabolite, antifungal
  "Chrysin",                      # [477] - flavonoid with antifungal activity
  "CAFFEATE",                     # [426] - phenolic with antifungal properties
  "trans.3.Hydroxycinnamic.acid"  # [267] - related to caffeic acid
)

# Check which antifungal metabolites are present in your dataset
available_antifungal <- intersect(antifungal_metabolites, names(fk_metabolom_gt_t))

# ===============================
# DATA PREPARATION WITH ANTIFUNGAL SUBSET
# ===============================

# Create subset with only antifungal metabolites + target variable
antifungal_data <- fk_metabolom_gt_t[, c(available_antifungal, "gt")]

# Scale the predictors
predictor_cols <- names(antifungal_data)[names(antifungal_data) != "gt"]
antifungal_data[predictor_cols] <- lapply(antifungal_data[predictor_cols], function(x) {
as.numeric(as.character(x))
})

# Scale the data
antifungal_data[predictor_cols] <- scale(antifungal_data[predictor_cols])

# ===============================
# TRAIN/TEST SPLIT
# ===============================

set.seed(123)
train_idx_af <- createDataPartition(antifungal_data$gt, p = 0.7, list = FALSE)
train_antifungal <- antifungal_data[train_idx_af, ]
test_antifungal <- antifungal_data[-train_idx_af, ]

cat("\nAntifungal model - Training set class distribution:\n")
print(table(train_antifungal$gt))
cat("Antifungal model - Test set class distribution:\n")
print(table(test_antifungal$gt))

# ===============================
# HYPERPARAMETER TUNING FOR ANTIFUNGAL MODEL
# ===============================

# Rename target variable levels
levels(train_antifungal$gt) <- c("GT_absent", "GT_present")
levels(test_antifungal$gt) <- c("GT_absent", "GT_present")

# Grid search for optimal mtry (fewer features now)
n_features <- ncol(train_antifungal) - 1
tuneGrid_af <- expand.grid(
mtry = c(
    max(1, floor(sqrt(n_features))),      # Rule of thumb
    max(1, floor(n_features/3)),          # Conservative
    max(1, floor(n_features/2)),          # Middle ground
    max(1, n_features)                    # All features (if small number)
)
)

# Remove duplicates
tuneGrid_af <- tuneGrid_af[!duplicated(tuneGrid_af$mtry), , drop = FALSE]

cat("\nTuning grid for antifungal model:\n")
print(tuneGrid_af)

# Cross-validation setup
ctrl_af <- trainControl(
method = "cv",
number = 10,
classProbs = TRUE,
summaryFunction = twoClassSummary,
savePredictions = "final"
)

# Train with cross-validation
set.seed(123)
rf_antifungal_tuned <- train(
    gt ~ ., 
    data = train_antifungal,
    method = "rf",
    tuneGrid = tuneGrid_af,
    trControl = ctrl_af,
    metric = "ROC",
    ntree = 1000,
    importance = TRUE
)

cat("\nAntifungal Random Forest Tuning Results:\n")
print(rf_antifungal_tuned)

# ===============================
# FINAL ANTIFUNGAL MODEL
# ===============================

best_mtry_af <- rf_antifungal_tuned$bestTune$mtry

set.seed(123)
rf_antifungal_final <- randomForest(
    gt ~ ., 
    data = train_antifungal,
    mtry = best_mtry_af,
    ntree = 1000,
    importance = TRUE,
    nodesize = 5
)

# ===============================
# MODEL EVALUATION
# ===============================

# Make predictions
pred_antifungal <- predict(rf_antifungal_final, test_antifungal)
actual_antifungal <- test_antifungal$gt

# Calculate metrics
accuracy_af <- mean(pred_antifungal == actual_antifungal)
conf_matrix_af <- confusionMatrix(pred_antifungal, actual_antifungal)

cat("\n", paste(rep("=", 50), collapse=""))
cat("\nANTIFUNGAL METABOLITES MODEL RESULTS:\n")
cat(paste(rep("=", 50), collapse=""), "\n")
cat("Number of features used:", length(available_antifungal), "\n")
cat("Features:", paste(available_antifungal, collapse = ", "), "\n")
cat("Accuracy:", round(accuracy_af, 3), "\n")
cat("Sensitivity:", round(conf_matrix_af$byClass["Sensitivity"], 3), "\n")
cat("Specificity:", round(conf_matrix_af$byClass["Specificity"], 3), "\n")
cat("AUC from CV:", round(max(rf_antifungal_tuned$results$ROC), 3), "\n")

print(conf_matrix_af$table)

# Variable importance for antifungal model
cat("\nVariable Importance (Antifungal Metabolites):\n")
importance_af <- importance(rf_antifungal_final)
importance_df_af <- data.frame(
Metabolite = rownames(importance_af),
MeanDecreaseAccuracy = importance_af[,1],
MeanDecreaseGini = importance_af[,2]
)
importance_df_af <- importance_df_af[order(importance_df_af$MeanDecreaseAccuracy, decreasing = TRUE), ]
print(importance_df_af)

# Plot variable importance
varImpPlot(rf_antifungal_final, main = "Variable Importance - Antifungal Metabolites Model")

# RESULTS SUMMARY
cat("Antifungal subset model - Accuracy:", round(accuracy_af, 3), "\n")
cat("Number of features - Full model: 676, Antifungal model:", length(available_antifungal), "\n")

# Save the antifungal model results
antifungal_results <- list(
model = rf_antifungal_final,
accuracy = accuracy_af,
sensitivity = conf_matrix_af$byClass["Sensitivity"],
specificity = conf_matrix_af$byClass["Specificity"],
conf_matrix = conf_matrix_af$table,
features_used = available_antifungal,
importance = importance_df_af,
cv_results = rf_antifungal_tuned
)
                                                                                                                                                                                                                    