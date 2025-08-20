# ====================================
# METABOLOMICS RANDOM FOREST MODELS COMPARISON
# ====================================
# This script compares 4 different RF approaches on metabolomics data

library(randomForest)
library(caret)
library(dplyr)
library(doParallel)

# Set up parallel processing
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel processing\n")

# ====================================
# LOAD DATA
# ====================================
cat("\n=== LOADING DATA ===\n")

# Load the logged metabolomics data (used by most models)
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
cat("Loaded logged metabolomics data:", nrow(fk_metabolom_gt_scaled), "samples x", 
    ncol(fk_metabolom_gt_scaled)-1, "features\n")

# Load the boost data for PCA model
boost_data <- readRDS("../data/fk_metabolomics_boost_data.RDS")
cat("Loaded boost data for PCA:", nrow(boost_data), "samples x", 
    ncol(boost_data)-1, "features\n")

# Check class distribution
cat("\nClass distribution in logged data:\n")
print(table(fk_metabolom_gt_scaled$gt))

# ====================================
# INITIALIZE RESULTS DATAFRAME
# ====================================
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
# SET SEED AND SPLIT DATA
# ====================================
set.seed(123)

# Create consistent train/test split for fair comparison
train_idx <- createDataPartition(fk_metabolom_gt_scaled$gt, p = 0.7, list = FALSE)
train_data <- fk_metabolom_gt_scaled[train_idx, ]
test_data <- fk_metabolom_gt_scaled[-train_idx, ]

cat("\nTrain/Test split: ", nrow(train_data), "/", nrow(test_data), "\n")

# ====================================
# MODEL 1: BASIC RANDOM FOREST
# ====================================
cat("\n=== MODEL 1: BASIC RANDOM FOREST ===\n")

set.seed(123)
rf_basic <- randomForest(
  gt ~ ., 
  data = train_data,
  mtry = floor(ncol(train_data)/2),
  ntree = 1000
)

# Predict and evaluate
pred_basic <- predict(rf_basic, test_data)
cm_basic <- confusionMatrix(pred_basic, test_data$gt)

# Calculate AUC
pred_basic_prob <- predict(rf_basic, test_data, type = "prob")
library(pROC)
roc_basic <- roc(test_data$gt, pred_basic_prob[,2])
auc_basic <- auc(roc_basic)

# Store results
results_df <- rbind(results_df, data.frame(
  Model = "Basic_RF",
  Method = "Default_Parameters",
  N_Features = ncol(train_data) - 1,
  Accuracy = round(cm_basic$overall['Accuracy'], 3),
  Balanced_Accuracy = round(cm_basic$byClass['Balanced Accuracy'], 3),
  Sensitivity = round(cm_basic$byClass['Sensitivity'], 3),
  Specificity = round(cm_basic$byClass['Specificity'], 3),
  Precision = round(cm_basic$byClass['Pos Pred Value'], 3),
  F1_Score = round(cm_basic$byClass['F1'], 3),
  Kappa = round(cm_basic$overall['Kappa'], 3),
  AUC = round(auc_basic, 3),
  stringsAsFactors = FALSE
))

cat("Basic RF - Accuracy:", round(cm_basic$overall['Accuracy'], 3), "\n")

# ====================================
# MODEL 2: TUNED RANDOM FOREST
# ====================================
cat("\n=== MODEL 2: TUNED RANDOM FOREST ===\n")

# Grid search for optimal parameters
tuneGrid <- expand.grid(
  mtry = c(floor(sqrt(ncol(train_data)-1)),
           floor(ncol(train_data)/3),
           floor(ncol(train_data)/2),
           floor(ncol(train_data) * 2/3))
)

# Remove duplicates
tuneGrid <- unique(tuneGrid)

# Cross-validation setup
ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Train with tuning
cat("Running hyperparameter tuning...\n")
rf_tuned_cv <- train(
  gt ~ ., 
  data = train_data,
  method = "rf",
  tuneGrid = tuneGrid,
  trControl = ctrl,
  metric = "ROC",
  ntree = 1000,
  importance = TRUE
)

# Get best parameters
best_mtry <- rf_tuned_cv$bestTune$mtry
cat("Best mtry:", best_mtry, "\n")

# Train final model with best parameters
set.seed(123)
rf_tuned <- randomForest(
  gt ~ ., 
  data = train_data,
  mtry = best_mtry,
  ntree = 1000,
  importance = TRUE,
  nodesize = 5
)

# Predict and evaluate
pred_tuned <- predict(rf_tuned, test_data)
cm_tuned <- confusionMatrix(pred_tuned, test_data$gt)

# Calculate AUC
pred_tuned_prob <- predict(rf_tuned, test_data, type = "prob")
roc_tuned <- roc(test_data$gt, pred_tuned_prob[,2])
auc_tuned <- auc(roc_tuned)

# Store results
results_df <- rbind(results_df, data.frame(
  Model = "Tuned_RF",
  Method = paste0("CV_Tuned_mtry", best_mtry),
  N_Features = ncol(train_data) - 1,
  Accuracy = round(cm_tuned$overall['Accuracy'], 3),
  Balanced_Accuracy = round(cm_tuned$byClass['Balanced Accuracy'], 3),
  Sensitivity = round(cm_tuned$byClass['Sensitivity'], 3),
  Specificity = round(cm_tuned$byClass['Specificity'], 3),
  Precision = round(cm_tuned$byClass['Pos Pred Value'], 3),
  F1_Score = round(cm_tuned$byClass['F1'], 3),
  Kappa = round(cm_tuned$overall['Kappa'], 3),
  AUC = round(auc_tuned, 3),
  stringsAsFactors = FALSE
))

cat("Tuned RF - Accuracy:", round(cm_tuned$overall['Accuracy'], 3), "\n")

# ====================================
# MODEL 3: ANTIFUNGAL METABOLITES
# ====================================
cat("\n=== MODEL 3: ANTIFUNGAL METABOLITES RF ===\n")

# Define interesting/antifungal compounds from your script
interesting_compounds <- c("Salicylate", "ISOLEUCINE", "Ornithine", "THREONINE", 
"ASPARTATE", "NICOTINAMIDE", "TRYPTOPHAN", "VALINE", "succinic acid", 
"LACTIC ACID", "Linalool", "Geranyl acetate", "trans-Nerolidol", "Nerol", 
"FERULATE", "CAFFEATE", "caffeine", "Cinnamic acid - 40.0 eV", 
"Vanillic acid", "BETAINE", "CHOLINE", "SARCOSINE", "N,N-Dimethylglycine", 
"3-HYDROXYBENZOATE", "4-HYDROXYBENZOATE", "Sinapic acid", 
"Hydroquinone", "TYROSINE", "trans-3-Hydroxycinnamic acid", 
"3,4-DIHYDROXYBENZOATE", "COUMARIN", "L-Arginine", "L-Proline", 
"PHENYLALANINE")

# Convert to valid column names and check availability
interesting_compounds <- make.names(interesting_compounds, unique = TRUE)
available_antifungal <- intersect(interesting_compounds, colnames(fk_metabolom_gt_scaled))

cat("Using", length(available_antifungal), "antifungal metabolites\n")

# Create antifungal subset
antifungal_data <- fk_metabolom_gt_scaled[, c(available_antifungal, "gt")]

# Split using same indices for consistency
train_antifungal <- antifungal_data[train_idx, ]
test_antifungal <- antifungal_data[-train_idx, ]

# Convert predictors to numeric
predictor_cols <- names(train_antifungal)[names(train_antifungal) != "gt"]
train_antifungal[predictor_cols] <- lapply(train_antifungal[predictor_cols], function(x) {
  as.numeric(as.character(x))
})
test_antifungal[predictor_cols] <- lapply(test_antifungal[predictor_cols], function(x) {
  as.numeric(as.character(x))
})

# Train antifungal model
set.seed(123)
rf_antifungal <- randomForest(
  gt ~ ., 
  data = train_antifungal,
  mtry = max(1, floor(sqrt(length(available_antifungal)))),
  ntree = 1000,
  importance = TRUE,
  nodesize = 5
)

# Predict and evaluate
pred_antifungal <- predict(rf_antifungal, test_antifungal)
cm_antifungal <- confusionMatrix(pred_antifungal, test_antifungal$gt)

# Calculate AUC
pred_antifungal_prob <- predict(rf_antifungal, test_antifungal, type = "prob")
roc_antifungal <- roc(test_antifungal$gt, pred_antifungal_prob[,2])
auc_antifungal <- auc(roc_antifungal)

# Store results
results_df <- rbind(results_df, data.frame(
  Model = "Antifungal_RF",
  Method = "Selected_Compounds",
  N_Features = length(available_antifungal),
  Accuracy = round(cm_antifungal$overall['Accuracy'], 3),
  Balanced_Accuracy = round(cm_antifungal$byClass['Balanced Accuracy'], 3),
  Sensitivity = round(cm_antifungal$byClass['Sensitivity'], 3),
  Specificity = round(cm_antifungal$byClass['Specificity'], 3),
  Precision = round(cm_antifungal$byClass['Pos Pred Value'], 3),
  F1_Score = round(cm_antifungal$byClass['F1'], 3),
  Kappa = round(cm_antifungal$overall['Kappa'], 3),
  AUC = round(auc_antifungal, 3),
  stringsAsFactors = FALSE
))

cat("Antifungal RF - Accuracy:", round(cm_antifungal$overall['Accuracy'], 3), "\n")

# ====================================
# MODEL 4: PCA-BASED RF
# ====================================
cat("\n=== MODEL 4: PCA-BASED RF ===\n")

# Convert gt to numeric for PCA model (1 = GT_present, 0 = GT_absent)
boost_data_numeric <- boost_data %>%
  mutate(gt = ifelse(gt == "GT_present", 1,
              ifelse(gt == "GT_absent", 0, gt))) %>%
  mutate(gt = as.numeric(gt))

# Use same train/test indices for consistency
train_boost <- boost_data_numeric[train_idx, ]
test_boost <- boost_data_numeric[-train_idx, ]

# Perform PCA on training data
cat("Performing PCA on training data...\n")

# Separate features and target
train_features <- train_boost[, names(train_boost) != "gt"]
test_features <- test_boost[, names(test_boost) != "gt"]

# Log transform and scale
train_features_log <- log10(train_features + 1)
test_features_log <- log10(test_features + 1)

# Calculate scaling parameters from training data
train_means <- colMeans(train_features_log)
train_sds <- apply(train_features_log, 2, sd)

# Scale both train and test using training parameters
train_features_scaled <- scale(train_features_log, center = train_means, scale = train_sds)
test_features_scaled <- scale(test_features_log, center = train_means, scale = train_sds)

# Perform PCA on training data
pca_model <- prcomp(train_features_scaled, center = FALSE, scale. = FALSE)

# Determine number of components (80% variance or max 10)
explained_var <- pca_model$sdev^2 / sum(pca_model$sdev^2)
cumulative_var <- cumsum(explained_var)
n_components <- min(which(cumulative_var >= 0.8), 10)

cat("Using", n_components, "principal components (", 
    round(cumulative_var[n_components] * 100, 1), "% variance explained)\n")

# Transform both train and test data
train_pca <- predict(pca_model, train_features_scaled)[, 1:n_components]
test_pca <- predict(pca_model, test_features_scaled)[, 1:n_components]

# Create dataframes for modeling
train_pca_df <- as.data.frame(train_pca)
colnames(train_pca_df) <- paste0("PC", 1:n_components)
train_pca_df$gt <- as.factor(train_boost$gt)
levels(train_pca_df$gt) <- c("GT_absent", "GT_present")

test_pca_df <- as.data.frame(test_pca)
colnames(test_pca_df) <- paste0("PC", 1:n_components)
test_pca_df$gt <- as.factor(test_boost$gt)
levels(test_pca_df$gt) <- c("GT_absent", "GT_present")

# Train PCA-based RF
set.seed(123)
rf_pca <- randomForest(
  gt ~ ., 
  data = train_pca_df,
  mtry = max(1, floor(sqrt(n_components))),
  ntree = 1000,
  nodesize = 3
)

# Predict and evaluate
pred_pca <- predict(rf_pca, test_pca_df)
cm_pca <- confusionMatrix(pred_pca, test_pca_df$gt)

# Calculate AUC
pred_pca_prob <- predict(rf_pca, test_pca_df, type = "prob")
roc_pca <- roc(test_pca_df$gt, pred_pca_prob[,2])
auc_pca <- auc(roc_pca)

# Store results
results_df <- rbind(results_df, data.frame(
  Model = "PCA_RF",
  Method = paste0(n_components, "_PCs"),
  N_Features = n_components,
  Accuracy = round(cm_pca$overall['Accuracy'], 3),
  Balanced_Accuracy = round(cm_pca$byClass['Balanced Accuracy'], 3),
  Sensitivity = round(cm_pca$byClass['Sensitivity'], 3),
  Specificity = round(cm_pca$byClass['Specificity'], 3),
  Precision = round(cm_pca$byClass['Pos Pred Value'], 3),
  F1_Score = round(cm_pca$byClass['F1'], 3),
  Kappa = round(cm_pca$overall['Kappa'], 3),
  AUC = round(auc_pca, 3),
  stringsAsFactors = FALSE
))

cat("PCA RF - Accuracy:", round(cm_pca$overall['Accuracy'], 3), "\n")

# ====================================
# DISPLAY FINAL RESULTS
# ====================================
cat("\n============================================================\n")
cat("FINAL COMPARISON OF ALL METABOLOMICS RF MODELS\n")
cat("============================================================\n\n")

# Sort by balanced accuracy
results_df <- results_df[order(results_df$Balanced_Accuracy, decreasing = TRUE), ]

# Display results
print(results_df)

# Summary statistics
cat("\n=== PERFORMANCE SUMMARY ===\n")
for(i in 1:nrow(results_df)) {
  cat(sprintf("%d. %-15s: %3d features, Acc = %.3f, Balanced = %.3f, AUC = %.3f\n",
              i,
              results_df$Model[i],
              results_df$N_Features[i],
              results_df$Accuracy[i],
              results_df$Balanced_Accuracy[i],
              results_df$AUC[i]))
}

# Identify best model
best_idx <- which.max(results_df$Balanced_Accuracy)
cat("\n=== BEST MODEL (by Balanced Accuracy) ===\n")
cat("Model:", results_df$Model[best_idx], "\n")
cat("Method:", results_df$Method[best_idx], "\n")
cat("Features:", results_df$N_Features[best_idx], "\n")
cat("Balanced Accuracy:", results_df$Balanced_Accuracy[best_idx], "\n")
cat("AUC:", results_df$AUC[best_idx], "\n")

# Save results
write.csv(results_df, "metabolomics_rf_comparison_results.csv", row.names = FALSE)
cat("\nResults saved to 'metabolomics_rf_comparison_results.csv'\n")

# Create visualization if ggplot2 is available
if(require(ggplot2, quietly = TRUE)) {
  library(ggplot2)
  
  # Reshape data for plotting
  results_long <- results_df %>%
    select(Model, Accuracy, Balanced_Accuracy, Sensitivity, Specificity, AUC) %>%
    tidyr::pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value")
  
  # Create comparison plot
  p <- ggplot(results_long, aes(x = Model, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    coord_flip() +
    labs(title = "Metabolomics Random Forest Models Comparison",
         x = "Model",
         y = "Performance Score",
         fill = "Metric") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    ylim(0, 1)
  
  ggsave("metabolomics_rf_comparison.png", p, width = 12, height = 6)
  cat("Visualization saved to 'metabolomics_rf_comparison.png'\n")
}

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

cat("\n=== ANALYSIS COMPLETE ===\n")