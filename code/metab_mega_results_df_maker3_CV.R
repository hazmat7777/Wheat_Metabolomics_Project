# ====================================
# METABOLOMICS RANDOM FOREST MODELS COMPARISON WITH CROSS-VALIDATION
# ====================================
# This script compares 4 different RF approaches on metabolomics data using 5-fold CV

library(randomForest)
library(caret)
library(dplyr)
library(doParallel)
library(pROC)

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
  Mean_Accuracy = numeric(),
  SD_Accuracy = numeric(),
  Mean_Balanced_Accuracy = numeric(),
  SD_Balanced_Accuracy = numeric(),
  Mean_Sensitivity = numeric(),
  SD_Sensitivity = numeric(),
  Mean_Specificity = numeric(),
  SD_Specificity = numeric(),
  Mean_AUC = numeric(),
  SD_AUC = numeric(),
  stringsAsFactors = FALSE
)

# ====================================
# FUNCTION TO RUN 5-FOLD CV FOR A MODEL
# ====================================
run_cv_evaluation <- function(data, model_name, method_name, train_function, predict_function = NULL) {
  
  cat("\n--- Running", model_name, "---\n")
  
  # Create 5-fold CV
  set.seed(123)
  folds <- createFolds(data$gt, k = 5, returnTrain = FALSE)
  
  # Initialize vectors to store results
  accuracies <- numeric(5)
  balanced_accuracies <- numeric(5)
  sensitivities <- numeric(5)
  specificities <- numeric(5)
  aucs <- numeric(5)
  
  for(i in 1:5) {
    cat("  Fold", i, "...")
    
    # Split data
    train_data <- data[-folds[[i]], ]
    test_data <- data[folds[[i]], ]
    
    # Train model using provided function
    model <- train_function(train_data)
    
    # Make predictions
    if(is.null(predict_function)) {
      pred <- predict(model, test_data)
      pred_prob <- predict(model, test_data, type = "prob")
    } else {
      pred_results <- predict_function(model, test_data)
      pred <- pred_results$pred
      pred_prob <- pred_results$pred_prob
    }
    
    # Calculate metrics
    cm <- confusionMatrix(pred, test_data$gt)
    
    accuracies[i] <- cm$overall['Accuracy']
    balanced_accuracies[i] <- cm$byClass['Balanced Accuracy']
    sensitivities[i] <- cm$byClass['Sensitivity']
    specificities[i] <- cm$byClass['Specificity']
    
    # Calculate AUC
    if(ncol(pred_prob) >= 2) {
      roc_obj <- roc(test_data$gt, pred_prob[,2], quiet = TRUE)
      aucs[i] <- auc(roc_obj)
    } else {
      aucs[i] <- NA
    }
    
    cat(" Acc:", round(accuracies[i], 3), 
        " Bal:", round(balanced_accuracies[i], 3), "\n")
  }
  
  # Calculate summary statistics
  result_row <- data.frame(
    Model = model_name,
    Method = method_name,
    N_Features = ncol(data) - 1,
    Mean_Accuracy = round(mean(accuracies, na.rm = TRUE), 3),
    SD_Accuracy = round(sd(accuracies, na.rm = TRUE), 3),
    Mean_Balanced_Accuracy = round(mean(balanced_accuracies, na.rm = TRUE), 3),
    SD_Balanced_Accuracy = round(sd(balanced_accuracies, na.rm = TRUE), 3),
    Mean_Sensitivity = round(mean(sensitivities, na.rm = TRUE), 3),
    SD_Sensitivity = round(sd(sensitivities, na.rm = TRUE), 3),
    Mean_Specificity = round(mean(specificities, na.rm = TRUE), 3),
    SD_Specificity = round(sd(specificities, na.rm = TRUE), 3),
    Mean_AUC = round(mean(aucs, na.rm = TRUE), 3),
    SD_AUC = round(sd(aucs, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
  
  cat("  Summary - Mean Balanced Accuracy:", result_row$Mean_Balanced_Accuracy, 
      "±", result_row$SD_Balanced_Accuracy, "\n")
  
  return(result_row)
}

# ====================================
# MODEL 1: BASIC RANDOM FOREST
# ====================================
cat("\n=== MODEL 1: BASIC RANDOM FOREST ===\n")

# Training function
train_basic_rf <- function(train_data) {
  randomForest(
    gt ~ ., 
    data = train_data,
    mtry = floor(ncol(train_data)/2),
    ntree = 1000
  )
}

# Run CV evaluation
basic_results <- run_cv_evaluation(
  data = fk_metabolom_gt_scaled,
  model_name = "Basic_RF",
  method_name = "Default_Parameters",
  train_function = train_basic_rf
)

results_df <- rbind(results_df, basic_results)

# ====================================
# MODEL 2: TUNED RANDOM FOREST
# ====================================
cat("\n=== MODEL 2: TUNED RANDOM FOREST ===\n")

# First, find optimal mtry using grid search on full dataset
tuneGrid <- expand.grid(
  mtry = c(floor(sqrt(ncol(fk_metabolom_gt_scaled)-1)),
           floor(ncol(fk_metabolom_gt_scaled)/3),
           floor(ncol(fk_metabolom_gt_scaled)/2),
           floor(ncol(fk_metabolom_gt_scaled) * 2/3))
)
tuneGrid <- unique(tuneGrid)

# Quick tuning to find best mtry
ctrl <- trainControl(method = "cv", number = 3, classProbs = TRUE, summaryFunction = twoClassSummary)
cat("Finding optimal mtry...\n")
rf_tune <- train(
  gt ~ ., 
  data = fk_metabolom_gt_scaled,
  method = "rf",
  tuneGrid = tuneGrid,
  trControl = ctrl,
  metric = "ROC",
  ntree = 500  # Reduced for speed during tuning
)

best_mtry <- rf_tune$bestTune$mtry
cat("Best mtry found:", best_mtry, "\n")

# Training function with optimal mtry
train_tuned_rf <- function(train_data) {
  randomForest(
    gt ~ ., 
    data = train_data,
    mtry = best_mtry,
    ntree = 1000,
    nodesize = 5
  )
}

# Run CV evaluation
tuned_results <- run_cv_evaluation(
  data = fk_metabolom_gt_scaled,
  model_name = "Tuned_RF",
  method_name = paste0("CV_Tuned_mtry", best_mtry),
  train_function = train_tuned_rf
)

results_df <- rbind(results_df, tuned_results)

# ====================================
# MODEL 3: ANTIFUNGAL METABOLITES RF
# ====================================
cat("\n=== MODEL 3: ANTIFUNGAL METABOLITES RF ===\n")

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

# Convert to valid column names and check availability
interesting_compounds <- make.names(interesting_compounds, unique = TRUE)
available_antifungal <- intersect(interesting_compounds, colnames(fk_metabolom_gt_scaled))

cat("Using", length(available_antifungal), "antifungal metabolites\n")

# Create antifungal subset
antifungal_data <- fk_metabolom_gt_scaled[, c(available_antifungal, "gt")]

# Convert predictors to numeric (data preparation function)
prepare_antifungal_data <- function(data) {
  predictor_cols <- names(data)[names(data) != "gt"]
  data[predictor_cols] <- lapply(data[predictor_cols], function(x) {
    as.numeric(as.character(x))
  })
  return(data)
}

# Apply data preparation
antifungal_data <- prepare_antifungal_data(antifungal_data)

# Training function
train_antifungal_rf <- function(train_data) {
  train_data <- prepare_antifungal_data(train_data)
  randomForest(
    gt ~ ., 
    data = train_data,
    mtry = max(1, floor(sqrt(length(available_antifungal)))),
    ntree = 1000,
    nodesize = 5
  )
}

# Custom prediction function to handle data preparation
predict_antifungal_rf <- function(model, test_data) {
  test_data <- prepare_antifungal_data(test_data)
  pred <- predict(model, test_data)
  pred_prob <- predict(model, test_data, type = "prob")
  return(list(pred = pred, pred_prob = pred_prob))
}

# Run CV evaluation
antifungal_results <- run_cv_evaluation(
  data = antifungal_data,
  model_name = "Antifungal_RF",
  method_name = "Selected_Compounds",
  train_function = train_antifungal_rf,
  predict_function = predict_antifungal_rf
)

results_df <- rbind(results_df, antifungal_results)

# ====================================
# MODEL 4: PCA-BASED RF
# ====================================
cat("\n=== MODEL 4: PCA-BASED RF ===\n")

# Source the PCA functions
source("metab_pca_then_rf_just_the_fns.R")

# Prepare boost data
boost_data <- readRDS("../data/fk_metabolomics_boost_data.RDS")
boost_data$gt <- factor(boost_data$gt, levels = c(0, 1), labels = c("GT_absent", "GT_present"))

# Training function for PCA-based RF
train_pca_rf <- function(train_data) {
  
  # Convert gt to numeric for PCA function (1 = GT_present, 0 = GT_absent)
  train_data_numeric <- train_data %>%
    mutate(gt = ifelse(gt == "GT_present", 1, 0))
  
  # Run PCA analysis on training data only
  pca_result <- evaluate_metabolomics_pca(
    train_data_numeric, 
    target_col = "gt", 
    test_prop = 0.01,  # Use almost all data for training
    variance_threshold = 0.8,
    max_components = 10
  )
  
  return(pca_result)
}

# Custom prediction function for PCA-based RF
predict_pca_rf <- function(pca_result, test_data) {
  
  # Convert test data to numeric format
  test_data_numeric <- test_data %>%
    mutate(gt_numeric = ifelse(gt == "GT_present", 1, 0))
  
  # Apply the same PCA transformation to test data
  # Extract PCA object and RF model
  pca_obj <- pca_result$preprocessing_params$pca_object
  rf_model <- pca_result$models$rf
  
  # Transform test data using the PCA from training
  test_features <- test_data_numeric[, !names(test_data_numeric) %in% c("gt", "gt_numeric")]
  test_pca <- predict(pca_obj, test_features)
  test_pca_df <- data.frame(test_pca, gt = test_data_numeric$gt_numeric)
  
  # Make predictions
  pred_numeric <- predict(rf_model, test_pca_df)
  pred_prob <- predict(rf_model, test_pca_df, type = "prob")
  
  # Convert back to factor
  pred <- factor(ifelse(pred_numeric == 1, "GT_present", "GT_absent"), 
                levels = c("GT_absent", "GT_present"))
  
  return(list(pred = pred, pred_prob = pred_prob))
}

# Create a simple training function that returns necessary components
train_pca_simple <- function(train_data) {
  # Convert to numeric
  train_numeric <- train_data %>%
    mutate(gt = ifelse(gt == "GT_present", 1, 0))
  
  # Get features only
  features <- train_numeric[, !names(train_numeric) %in% "gt"]
  
  # Perform PCA
  pca_obj <- prcomp(features, scale. = TRUE, center = TRUE)
  
  # Determine number of components for 80% variance
  var_explained <- cumsum(pca_obj$sdev^2) / sum(pca_obj$sdev^2)
  n_comp <- min(which(var_explained >= 0.8), 10)
  
  # Create training data with PCs
  train_pca <- data.frame(pca_obj$x[, 1:n_comp], gt = train_numeric$gt)
  
  # Train RF on PCA data
  rf_model <- randomForest(
    factor(gt) ~ ., 
    data = train_pca,
    ntree = 1000
  )
  
  return(list(pca = pca_obj, rf = rf_model, n_comp = n_comp))
}

# Simple prediction function
predict_pca_simple <- function(model_list, test_data) {
  # Convert test data
  test_numeric <- test_data %>%
    mutate(gt_numeric = ifelse(gt == "GT_present", 1, 0))
  
  features <- test_numeric[, !names(test_numeric) %in% c("gt", "gt_numeric")]
  
  # Transform using PCA
  test_pca <- predict(model_list$pca, features)
  test_pca_df <- data.frame(test_pca[, 1:model_list$n_comp, drop = FALSE], 
                           gt = test_numeric$gt_numeric)
  
  # Predict
  pred_numeric <- predict(model_list$rf, test_pca_df)
  pred_prob <- predict(model_list$rf, test_pca_df, type = "prob")
  
  # Convert predictions back
  pred <- factor(ifelse(as.numeric(as.character(pred_numeric)) == 1, "GT_present", "GT_absent"), 
                levels = c("GT_absent", "GT_present"))
  
  return(list(pred = pred, pred_prob = pred_prob))
}

# Run CV evaluation for PCA model
pca_results <- run_cv_evaluation(
  data = boost_data,
  model_name = "PCA_RF",
  method_name = "PCA_Components",
  train_function = train_pca_simple,
  predict_function = predict_pca_simple
)

results_df <- rbind(results_df, pca_results)

# ====================================
# DISPLAY FINAL RESULTS
# ====================================
cat("\n============================================================\n")
cat("FINAL COMPARISON OF ALL METABOLOMICS RF MODELS\n")
cat("============================================================\n\n")

# Sort by mean balanced accuracy
results_df <- results_df[order(results_df$Mean_Balanced_Accuracy, decreasing = TRUE), ]

# Display key results
print(results_df[, c("Model", "Method", "N_Features", "Mean_Accuracy", "SD_Accuracy", 
                     "Mean_Balanced_Accuracy", "SD_Balanced_Accuracy", "Mean_AUC", "SD_AUC")])

# Summary statistics
cat("\n=== PERFORMANCE SUMMARY ===\n")
cat("Ranked by Mean Balanced Accuracy:\n\n")
for(i in 1:nrow(results_df)) {
  cat(sprintf("%d. %-15s: %3d features, Bal.Acc = %.3f±%.3f, AUC = %.3f±%.3f\n",
              i,
              results_df$Model[i],
              results_df$N_Features[i],
              results_df$Mean_Balanced_Accuracy[i],
              results_df$SD_Balanced_Accuracy[i],
              results_df$Mean_AUC[i],
              results_df$SD_AUC[i]))
}

# Identify best model
best_idx <- which.max(results_df$Mean_Balanced_Accuracy)
cat("\n=== BEST MODEL (by Mean Balanced Accuracy) ===\n")
cat("Model:", results_df$Model[best_idx], "\n")
cat("Method:", results_df$Method[best_idx], "\n")
cat("Features:", results_df$N_Features[best_idx], "\n")
cat("Mean Balanced Accuracy:", results_df$Mean_Balanced_Accuracy[best_idx], 
    "±", results_df$SD_Balanced_Accuracy[best_idx], "\n")
cat("Mean AUC:", results_df$Mean_AUC[best_idx], 
    "±", results_df$SD_AUC[best_idx], "\n")

# Save results
write.csv(results_df, "../results/metabolomics_rf_cv_comparison_results2.csv", row.names = FALSE)
cat("\nResults saved to 'metabolomics_rf_cv_comparison_results.csv'\n")

# Create visualization if ggplot2 is available
if(require(ggplot2, quietly = TRUE)) {
  library(ggplot2)
  
  # Plot with error bars
  p1 <- ggplot(results_df, aes(x = reorder(Model, Mean_Balanced_Accuracy), 
                               y = Mean_Balanced_Accuracy)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(ymin = Mean_Balanced_Accuracy - SD_Balanced_Accuracy,
                      ymax = Mean_Balanced_Accuracy + SD_Balanced_Accuracy),
                  width = 0.2) +
    coord_flip() +
    labs(title = "Metabolomics Random Forest Models Comparison",
         subtitle = "5-fold Cross-Validation Results",
         x = "Model",
         y = "Mean Balanced Accuracy ± SD") +
    theme_minimal()
  
  ggsave("../results/plots/metabolomics_rf_cv_comparison2.png", p1, width = 10, height = 6)
  cat("Visualization saved to 'metabolomics_rf_cv_comparison.png'\n")
}

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

cat("\n=== ANALYSIS COMPLETE ===\n")








# Basic structure
class(boost_data$gt)
typeof(boost_data$gt)
str(boost_data$gt)

# Values and distribution
unique(boost_data$gt)
table(boost_data$gt)
summary(boost_data$gt)

# Check for factors
is.factor(boost_data$gt)
if(is.factor(boost_data$gt)) levels(boost_data$gt)

# Check first few rows
head(boost_data$gt, 10)
tail(boost_data$gt, 10)

# Check for any NAs or special values
sum(is.na(boost_data$gt))
any(is.infinite(boost_data$gt))

# Compare with the working data
class(fk_metabolom_gt_scaled$gt)
table(fk_metabolom_gt_scaled$gt)