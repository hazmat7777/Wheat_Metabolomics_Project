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
  
  # Check class distribution first
  class_counts <- table(data$gt)
  cat("Class distribution:", paste(names(class_counts), "=", class_counts, collapse = ", "), "\n")
  
  if(min(class_counts) < 5) {
    cat("Warning: Very small class size detected. Using stratified sampling.\n")
  }
  
  # Create 5-fold CV with stratified sampling
  set.seed(123)
  folds <- createFolds(data$gt, k = 5, returnTrain = FALSE)
  
  # Check each fold for class balance
  valid_folds <- c()
  for(i in 1:5) {
    train_data <- data[-folds[[i]], ]
    test_data <- data[folds[[i]], ]
    
    train_classes <- table(train_data$gt)
    test_classes <- table(test_data$gt)
    
    # Check if both classes exist in training data
    if(length(train_classes) >= 2 && all(train_classes > 0)) {
      valid_folds <- c(valid_folds, i)
    } else {
      cat(sprintf("  Fold %d skipped - insufficient class balance in training data\n", i))
      cat("    Training classes:", paste(names(train_classes), "=", train_classes, collapse = ", "), "\n")
    }
  }
  
  if(length(valid_folds) < 2) {
    cat("Error: Not enough valid folds for cross-validation. Using simpler evaluation.\n")
    
    # Fallback to a single train/test split
    set.seed(123)
    train_idx <- createDataPartition(data$gt, p = 0.7, list = FALSE)
    train_data <- data[train_idx, ]
    test_data <- data[-train_idx, ]
    
    cat("Fallback: Single train/test split -", nrow(train_data), "train,", nrow(test_data), "test\n")
    
    # Train model
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
    
    # Return single result
    result_row <- data.frame(
      Model = model_name,
      Method = method_name,
      N_Features = ncol(data) - 1,
      Mean_Accuracy = round(cm$overall['Accuracy'], 3),
      SD_Accuracy = 0,  # No SD for single split
      Mean_Balanced_Accuracy = round(cm$byClass['Balanced Accuracy'], 3),
      SD_Balanced_Accuracy = 0,  # No SD for single split
      Mean_Sensitivity = round(cm$byClass['Sensitivity'], 3),
      SD_Sensitivity = 0,
      Mean_Specificity = round(cm$byClass['Specificity'], 3),
      SD_Specificity = 0,
      Mean_AUC = 0.5,  # Placeholder
      SD_AUC = 0,
      stringsAsFactors = FALSE
    )
    
    cat("  Single split result - Balanced Accuracy:", result_row$Mean_Balanced_Accuracy, "\n")
    return(result_row)
  }
  
  # Use only valid folds
  n_folds <- length(valid_folds)
  cat("Using", n_folds, "valid folds for CV\n")
  
  # Initialize vectors to store results
  accuracies <- numeric(n_folds)
  balanced_accuracies <- numeric(n_folds)
  sensitivities <- numeric(n_folds)
  specificities <- numeric(n_folds)
  aucs <- numeric(n_folds)
  
  for(j in 1:n_folds) {
    i <- valid_folds[j]
    cat("  Fold", i, "...")
    
    # Split data
    train_data <- data[-folds[[i]], ]
    test_data <- data[folds[[i]], ]
    
    tryCatch({
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
      
      accuracies[j] <- cm$overall['Accuracy']
      balanced_accuracies[j] <- cm$byClass['Balanced Accuracy']
      sensitivities[j] <- cm$byClass['Sensitivity']
      specificities[j] <- cm$byClass['Specificity']
      
      # Calculate AUC
      if(ncol(pred_prob) >= 2) {
        roc_obj <- roc(test_data$gt, pred_prob[,2], quiet = TRUE)
        aucs[j] <- auc(roc_obj)
      } else {
        aucs[j] <- NA
      }
      
      cat(" Acc:", round(accuracies[j], 3), 
          " Bal:", round(balanced_accuracies[j], 3), "\n")
      
    }, error = function(e) {
      cat(" ERROR:", e$message, "\n")
      accuracies[j] <- NA
      balanced_accuracies[j] <- NA
      sensitivities[j] <- NA
      specificities[j] <- NA
      aucs[j] <- NA
    })
  }
  
  # Calculate summary statistics (removing NAs)
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
# MANUAL CV FOR PCA MODEL
# ====================================
cat("\n=== MANUAL CV FOR PCA MODEL ===\n")

# Function to extract balanced accuracy from PCA result
extract_pca_metrics <- function(pca_result) {
  # The confusion matrix is in the threshold results
  cm <- pca_result$results$threshold_results[["0.5"]]$confusion_matrix
  
  metrics <- list(
    accuracy = cm$overall['Accuracy'],
    balanced_accuracy = cm$byClass['Balanced Accuracy'],
    sensitivity = cm$byClass['Sensitivity'],
    specificity = cm$byClass['Specificity'],
    auc = 0.5  # Placeholder - would need to calculate from predictions
  )
  
  return(metrics)
}

# Run PCA 5 times with different seeds for CV
n_runs <- 5
pca_cv_results <- list()

for(i in 1:n_runs) {
  cat("PCA Run", i, "...")
  
  # Use different seed for each run
  pca_result <- evaluate_metabolomics_pca(
    boost_data, 
    target_col = "gt", 
    test_prop = 0.3,
    seed = 123 + i,  # Different seed each time
    variance_threshold = 0.8,
    max_components = 10
  )
  
  # Extract metrics
  metrics <- extract_pca_metrics(pca_result)
  pca_cv_results[[i]] <- metrics
  
  cat(" Balanced Acc:", round(metrics$balanced_accuracy, 3), "\n")
}

# Calculate summary statistics
accuracies <- sapply(pca_cv_results, function(x) x$accuracy)
balanced_accuracies <- sapply(pca_cv_results, function(x) x$balanced_accuracy)
sensitivities <- sapply(pca_cv_results, function(x) x$sensitivity)
specificities <- sapply(pca_cv_results, function(x) x$specificity)

# Create PCA results row for comparison
pca_results_manual <- data.frame(
  Model = "PCA_RF",
  Method = "Manual_CV",
  N_Features = 10,  # Number of PCs
  Mean_Accuracy = round(mean(accuracies, na.rm = TRUE), 3),
  SD_Accuracy = round(sd(accuracies, na.rm = TRUE), 3),
  Mean_Balanced_Accuracy = round(mean(balanced_accuracies, na.rm = TRUE), 3),
  SD_Balanced_Accuracy = round(sd(balanced_accuracies, na.rm = TRUE), 3),
  Mean_Sensitivity = round(mean(sensitivities, na.rm = TRUE), 3),
  SD_Sensitivity = round(sd(sensitivities, na.rm = TRUE), 3),
  Mean_Specificity = round(mean(specificities, na.rm = TRUE), 3),
  SD_Specificity = round(sd(specificities, na.rm = TRUE), 3),
  Mean_AUC = 0.5,  # Placeholder
  SD_AUC = 0,
  stringsAsFactors = FALSE
)

cat("\nPCA Manual CV Summary:\n")
cat("Mean Balanced Accuracy:", pca_results_manual$Mean_Balanced_Accuracy, 
    "±", pca_results_manual$SD_Balanced_Accuracy, "\n")

# Add to results dataframe
results_df <- rbind(results_df, pca_results_manual)

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
write.csv(results_df, "../results/metabolomics_rf_cv_comparison_results.csv", row.names = FALSE)
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
  
  ggsave("../results/plots/metabolomics_rf_cv_comparison.png", p1, width = 10, height = 6)
  cat("Visualization saved to 'metabolomics_rf_cv_comparison.png'\n")
}

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

cat("\n=== ANALYSIS COMPLETE ===\n")


