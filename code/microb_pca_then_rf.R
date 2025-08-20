# this script runs machine learning models on microbial abundance data using PCA
# POSITIVE CLASS = 1 (your target condition)

# Load required libraries
library(gbm)
library(caret)
library(glmnet)
library(randomForest)
library(dplyr)
library(phyloseq)

# Load data
ps <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")

ps
###########################################################################
# PHYLOSEQ TO DATAFRAME CONVERSION FUNCTION
###########################################################################

phyloseq_to_dataframe <- function(ps, target_variable) {
  
  cat("=== CONVERTING PHYLOSEQ TO DATAFRAME ===\n")
  
  # Extract OTU/ASV abundance table (samples as rows, taxa as columns)
  otu_table <- as.data.frame(t(otu_table(ps)))
  
  # Extract sample metadata
  sample_data <- as.data.frame(sample_data(ps))
  
  # Check if target variable exists
  if(!target_variable %in% colnames(sample_data)) {
    stop(paste("Target variable", target_variable, "not found in sample data. Available variables:", 
               paste(colnames(sample_data), collapse = ", ")))
  }
  
  # Combine abundance data with target variable
  combined_data <- cbind(otu_table, gt = sample_data[[target_variable]])
  
  cat("Original dimensions - Samples:", nrow(combined_data), "Taxa:", ncol(otu_table), "\n")
  cat("Target variable distribution:\n")
  print(table(combined_data$gt))
  
  return(combined_data)
}

###########################################################################
# PCA-BASED PREPROCESSING FUNCTION FOR MICROBIAL DATA
###########################################################################

preprocess_microbial_pca <- function(data, target_col, train_idx, 
                                     variance_threshold = 0.8,
                                     max_components = 10,
                                     min_abundance_filter = 0.0001,
                                     prevalence_threshold = 0.1) {
  
  # Separate features and target
  if(is.character(target_col)) {
    target <- data[[target_col]]
    features <- data[, !names(data) %in% target_col]
  } else {
    target <- data[, target_col]
    features <- data[, -target_col]
  }
  
  cat("=== PCA PREPROCESSING FOR MICROBIAL DATA (TRAINING DATA ONLY) ===\n")
  cat("Original data dimensions:", dim(features), "\n")
  cat("Number of zeros:", sum(features == 0), "\n")
  
  # STEP 1: Filter low abundance and low prevalence taxa using TRAINING data
  cat("Filtering low abundance and low prevalence taxa...\n")
  
  train_features <- features[train_idx, ]
  
  # Filter by abundance (mean relative abundance > threshold)
  mean_abundance <- colMeans(train_features)
  abundance_keep <- mean_abundance > min_abundance_filter
  
  # Filter by prevalence (present in > threshold proportion of training samples)
  prevalence <- colSums(train_features > 0) / nrow(train_features)
  prevalence_keep <- prevalence > prevalence_threshold
  
  # Combine filters
  keep_taxa <- abundance_keep & prevalence_keep
  
  cat("Abundance filter: keeping", sum(abundance_keep), "of", length(abundance_keep), "taxa\n")
  cat("Prevalence filter: keeping", sum(prevalence_keep), "of", length(prevalence_keep), "taxa\n")
  cat("Combined filter: keeping", sum(keep_taxa), "of", length(keep_taxa), "taxa\n")
  
  if(sum(keep_taxa) == 0) {
    cat("Warning: No taxa pass filters. Using all taxa.\n")
    features_filtered <- features
  } else {
    features_filtered <- features[, keep_taxa]
  }
  
  # STEP 2: Handle zeros - for relative abundance data, add small pseudocount
  cat("Adding pseudocount for zero values...\n")
  
  # Add small pseudocount (1/10 of minimum non-zero value in training data)
  train_nonzero <- features_filtered[train_idx, ][features_filtered[train_idx, ] > 0]
  if(length(train_nonzero) > 0) {
    min_nonzero <- min(train_nonzero)
    pseudocount <- min_nonzero / 10
  } else {
    pseudocount <- 1e-6
  }
  
  features_pseudocount <- features_filtered + pseudocount
  
  # STEP 3: CLR transformation (recommended for compositional data)
  cat("Performing CLR (Centered Log-Ratio) transformation...\n")
  
  # CLR transformation: log(x_i / geometric_mean(x))
  features_clr <- features_pseudocount
  for(i in 1:nrow(features_pseudocount)) {
    row_data <- as.numeric(features_pseudocount[i, ])
    geom_mean <- exp(mean(log(row_data)))
    features_clr[i, ] <- log(row_data / geom_mean)
  }
  
  # STEP 4: Robust scaling using TRAINING statistics
  cat("Performing robust scaling...\n")
  
  # Calculate robust scaling parameters from training data only
  train_clr <- features_clr[train_idx, ]
  train_medians <- sapply(train_clr, median, na.rm = TRUE)
  train_mads <- sapply(train_clr, mad, na.rm = TRUE)
  
  # Apply scaling to ALL data using training parameters
  features_scaled <- features_clr
  for(i in 1:ncol(features_clr)) {
    if(train_mads[i] > 0) {
      features_scaled[, i] <- (features_clr[, i] - train_medians[i]) / train_mads[i]
    }
  }
  
  # STEP 5: Additional variance filtering after transformation
  cat("Performing post-transformation variance filtering...\n")
  
  train_features_scaled <- features_scaled[train_idx, ]
  train_var <- sapply(train_features_scaled, var, na.rm = TRUE)
  
  # Keep features with reasonable variance
  valid_features <- train_var > 0.01 & is.finite(train_var)
  
  if(sum(valid_features) == 0) {
    cat("Warning: No features with sufficient variance. Using all features.\n")
    features_for_pca <- features_scaled
  } else {
    features_for_pca <- features_scaled[, valid_features]
    cat("Removed", sum(!valid_features), "features with very low variance after CLR\n")
  }
  
  cat("Features going into PCA:", ncol(features_for_pca), "\n")
  
  # STEP 6: Apply PCA using TRAINING data only
  cat("Applying PCA to training data...\n")
  
  train_features_for_pca <- features_for_pca[train_idx, ]
  
  # Perform PCA on training data
  pca_result <- tryCatch({
    prcomp(train_features_for_pca, center = FALSE, scale. = FALSE)  # Already scaled
  }, error = function(e) {
    cat("PCA failed with error:", e$message, "\n")
    cat("Trying with centering...\n")
    prcomp(train_features_for_pca, center = TRUE, scale. = FALSE)
  })
  
  # Calculate explained variance
  explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_var <- cumsum(explained_var)
  
  # Determine number of components
  if(variance_threshold < 1) {
    # Select components explaining desired variance threshold
    n_components <- min(which(cumulative_var >= variance_threshold))
    n_components <- min(n_components, max_components, ncol(train_features_for_pca))
  } else {
    # Use fixed number of components
    n_components <- min(max_components, ncol(train_features_for_pca))
  }
  
  cat("Selected", n_components, "principal components\n")
  cat("Explained variance:", round(cumulative_var[n_components], 3), "\n")
  cat("Individual PC variances:", round(explained_var[1:n_components], 3), "\n")
  
  # STEP 7: Transform ALL data using the PCA fitted on training data
  cat("Transforming all data using training-fitted PCA...\n")
  
  # Transform all data (both training and test)
  features_pca <- predict(pca_result, features_for_pca)[, 1:n_components, drop = FALSE]
  
  # Create column names for PCs
  colnames(features_pca) <- paste0("PC", 1:n_components)
  
  # Create final dataset
  processed_data <- as.data.frame(features_pca)
  processed_data$gt <- target
  
  # Store preprocessing parameters
  preprocessing_params <- list(
    pca_model = pca_result,
    n_components = n_components,
    explained_variance = explained_var[1:n_components],
    cumulative_variance = cumulative_var[n_components],
    train_medians = train_medians,
    train_mads = train_mads,
    pseudocount = pseudocount,
    kept_taxa = colnames(features_for_pca),
    min_abundance_filter = min_abundance_filter,
    prevalence_threshold = prevalence_threshold
  )
  
  # Feature loadings for interpretation
  loadings <- pca_result$rotation[, 1:n_components, drop = FALSE]
  
  cat("Final processed data dimensions:", dim(processed_data), "\n")
  cat("Principal components:", colnames(features_pca), "\n")
  
  return(list(
    data = processed_data,
    params = preprocessing_params,
    loadings = loadings,
    selected_features = paste0("PC", 1:n_components)
  ))
}

###########################################################################
# MODEL TRAINING FUNCTION (SAME AS ORIGINAL)
###########################################################################

train_multiple_models <- function(train_data, target_col = "gt") {
  
  cat("=== TRAINING MODELS ===\n")
  
  # Prepare data
  x_train <- as.matrix(train_data[, !names(train_data) %in% target_col])
  y_train <- train_data[[target_col]]
  
  # Model 1: GBM with CV for optimal trees
  cat("Training GBM...\n")
  gbm_model <- gbm(
    gt ~ ., 
    data = train_data,
    distribution = "bernoulli",
    n.trees = 200,
    interaction.depth = 2,
    shrinkage = 0.1,
    n.minobsinnode = 3,
    bag.fraction = 0.8,
    cv.folds = 5,
    verbose = FALSE
  )
  optimal_trees <- gbm.perf(gbm_model, method = "cv", plot.it = FALSE)
  
  # Model 2: Elastic Net with CV
  cat("Training Elastic Net...\n")
  cv_glmnet <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5, nfolds = 5)
  
  # Model 3: Random Forest
  cat("Training Random Forest...\n")
  rf_model <- randomForest(
    as.factor(gt) ~ ., 
    data = train_data,
    ntree = 200,
    mtry = max(1, floor(sqrt(ncol(train_data)-1))),
    nodesize = 3
  )
  
  cat("Models trained successfully\n\n")
  
  return(list(
    gbm = gbm_model,
    glmnet = cv_glmnet,
    rf = rf_model,
    optimal_trees = optimal_trees,
    x_train = x_train,
    y_train = y_train
  ))
}

###########################################################################
# MODEL EVALUATION FUNCTION (SAME AS ORIGINAL)
###########################################################################

evaluate_on_test_set <- function(models, test_data, target_col = "gt") {
  
  cat("=== EVALUATING ON TEST SET ===\n")
  
  # Prepare test data
  x_test <- as.matrix(test_data[, !names(test_data) %in% target_col])
  y_test <- test_data[[target_col]]
  
  # Get predictions
  gbm_pred <- predict(models$gbm, test_data, n.trees = models$optimal_trees, type = "response")
  glmnet_pred <- predict(models$glmnet, x_test, s = "lambda.min", type = "response")[,1]
  rf_pred <- predict(models$rf, test_data, type = "prob")[,2]  # Class "1" probabilities
  
  # Ensemble prediction
  ensemble_pred <- (gbm_pred + glmnet_pred + rf_pred) / 3
  
  # Test multiple thresholds
  thresholds <- seq(0.3, 0.7, 0.1)
  threshold_results <- list()
  
  for(thresh in thresholds) {
    pred_class <- ifelse(ensemble_pred > thresh, 1, 0)
    pred_class <- factor(pred_class, levels = c(0, 1))
    actual <- factor(y_test, levels = c(0, 1))
    
    # Set positive class = "1" for confusionMatrix
    cm <- confusionMatrix(pred_class, actual, positive = "1")
    
    threshold_results[[as.character(thresh)]] <- list(
      threshold = thresh,
      accuracy = cm$overall['Accuracy'],
      kappa = cm$overall['Kappa'],
      sensitivity = cm$byClass['Sensitivity'],  # Now for class 1
      specificity = cm$byClass['Specificity'],  # Now for class 1
      precision = cm$byClass['Pos Pred Value'], # Precision for class 1
      confusion_matrix = cm
    )
  }
  
  # Find best threshold
  kappa_scores <- sapply(threshold_results, function(x) x$kappa)
  best_thresh <- names(threshold_results)[which.max(kappa_scores)]
  
  cat("Best threshold:", best_thresh, "\n")
  print(threshold_results[[best_thresh]]$confusion_matrix)
  
  # Individual model performance
  cat("\n=== INDIVIDUAL MODEL COMPARISON ===\n")
  thresh_val <- as.numeric(best_thresh)
  
  individual_results <- list()
  model_preds <- list(GBM = gbm_pred, ElasticNet = glmnet_pred, RandomForest = rf_pred, Ensemble = ensemble_pred)
  
  for(model_name in names(model_preds)) {
    pred_class <- ifelse(model_preds[[model_name]] > thresh_val, 1, 0)
    pred_class <- factor(pred_class, levels = c(0, 1))
    actual <- factor(y_test, levels = c(0, 1))
    
    # Set positive class = "1"
    cm <- confusionMatrix(pred_class, actual, positive = "1")
    individual_results[[model_name]] <- list(
      accuracy = cm$overall['Accuracy'],
      kappa = cm$overall['Kappa'],
      sensitivity = cm$byClass['Sensitivity'],  # Sensitivity for class 1
      specificity = cm$byClass['Specificity'],  # Specificity for class 1
      precision = cm$byClass['Pos Pred Value']  # Precision for class 1
    )
    
    cat(sprintf("%s - Accuracy: %.3f, Kappa: %.3f, Sensitivity: %.3f, Precision: %.3f\n", 
                model_name, cm$overall['Accuracy'], cm$overall['Kappa'],
                cm$byClass['Sensitivity'], cm$byClass['Pos Pred Value']))
  }
  
  # PC importance (instead of individual taxa)
  cat("\n=== PRINCIPAL COMPONENT IMPORTANCE ===\n")
  gbm_imp <- summary(models$gbm, n.trees = models$optimal_trees, plotit = FALSE)
  print(head(gbm_imp, 10))
  
  return(list(
    threshold_results = threshold_results,
    best_threshold = thresh_val,
    individual_results = individual_results,
    predictions = model_preds,
    feature_importance = gbm_imp
  ))
}

###########################################################################
# MAIN EVALUATION FUNCTION FOR MICROBIAL DATA
###########################################################################

evaluate_microbial_pca <- function(ps, target_variable, test_prop = 0.3, seed = 123,
                                  variance_threshold = 0.8, max_components = 10,
                                  min_abundance_filter = 0.0001, prevalence_threshold = 0.1) {
  
  set.seed(seed)
  
  # Convert phyloseq to dataframe
  data <- phyloseq_to_dataframe(ps, target_variable)
  
  # Convert target to binary numeric (adjust these categories to match your data)
  # MODIFY THESE LINES TO MATCH YOUR TARGET VARIABLE VALUES
  data <- data %>%
    mutate(gt = case_when(
      gt == "your_positive_condition" ~ 1,  # CHANGE THIS
      gt == "your_negative_condition" ~ 0,  # CHANGE THIS
      TRUE ~ as.numeric(as.factor(gt)) - 1  # Fallback: convert to 0/1
    ))
  
  # Verify class distribution
  cat("=== CLASS DISTRIBUTION ===\n")
  print(table(data$gt))
  cat("0 = negative class\n")
  cat("1 = positive class\n\n")
  
  # Single train/test split
  train_idx <- createDataPartition(data$gt, p = 1 - test_prop, list = FALSE)
  test_idx <- setdiff(1:nrow(data), train_idx)
  
  cat("Train samples:", length(train_idx), "Test samples:", length(test_idx), "\n")
  
  # Preprocess data with PCA
  preprocessed <- preprocess_microbial_pca(
    data, "gt", train_idx,
    variance_threshold = variance_threshold,
    max_components = max_components,
    min_abundance_filter = min_abundance_filter,
    prevalence_threshold = prevalence_threshold
  )
  
  processed_data <- preprocessed$data
  
  # Split processed data
  train_data <- processed_data[train_idx, ]
  test_data <- processed_data[test_idx, ]
  
  cat("Final training data dimensions:", dim(train_data), "\n")
  cat("Final test data dimensions:", dim(test_data), "\n\n")
  
  # Train models
  models <- train_multiple_models(train_data, "gt")
  
  # Evaluate on test set
  results <- evaluate_on_test_set(models, test_data, "gt")
  
  return(list(
    models = models,
    results = results,
    train_data = train_data,
    test_data = test_data,
    train_idx = train_idx,
    test_idx = test_idx,
    preprocessing_params = preprocessed$params,
    loadings = preprocessed$loadings,
    selected_features = preprocessed$selected_features,
    original_data = data
  ))
}

###########################################################################
# MICROBIAL-SPECIFIC INTERPRETATION FUNCTIONS
###########################################################################

interpret_microbial_loadings <- function(pca_results, ps, top_n = 10) {
  cat("=== TOP TAXA CONTRIBUTING TO EACH PC ===\n")
  
  loadings <- pca_results$loadings
  
  # Get taxonomic information if available
  if(!is.null(tax_table(ps))) {
    tax_info <- as.data.frame(tax_table(ps))
  } else {
    tax_info <- NULL
  }
  
  for(pc in 1:ncol(loadings)) {
    cat("\n--- PC", pc, "---\n")
    
    # Get absolute loadings for this PC
    pc_loadings <- abs(loadings[, pc])
    names(pc_loadings) <- rownames(loadings)
    
    # Sort and get top contributors
    top_taxa <- sort(pc_loadings, decreasing = TRUE)[1:min(top_n, length(pc_loadings))]
    
    for(i in 1:length(top_taxa)) {
      taxa_name <- names(top_taxa)[i]
      loading_val <- top_taxa[i]
      
      # Add taxonomic info if available
      if(!is.null(tax_info) && taxa_name %in% rownames(tax_info)) {
        # Find most specific taxonomic level
        tax_row <- tax_info[taxa_name, ]
        tax_string <- ""
        for(tax_level in rev(colnames(tax_info))) {
          if(!is.na(tax_row[[tax_level]]) && tax_row[[tax_level]] != "") {
            tax_string <- paste0(" (", tax_row[[tax_level]], ")")
            break
          }
        }
        cat(sprintf("%2d. %s%s: %.3f\n", i, taxa_name, tax_string, loading_val))
      } else {
        cat(sprintf("%2d. %s: %.3f\n", i, taxa_name, loading_val))
      }
    }
  }
}

plot_microbial_pca <- function(pca_results) {
  # Simple PCA plot
  pca_data <- pca_results$train_data
  
  if(ncol(pca_data) >= 3) {  # At least PC1 and PC2 plus target
    plot(pca_data$PC1, pca_data$PC2, 
         col = ifelse(pca_data$gt == 1, "red", "blue"),
         pch = 16, 
         xlab = "PC1", ylab = "PC2",
         main = "PCA of Microbial Abundance Data")
    legend("topright", legend = c("Negative", "Positive"), 
           col = c("blue", "red"), pch = 16)
  }
}

###########################################################################
# USAGE EXAMPLE
###########################################################################

# IMPORTANT: MODIFY THE TARGET VARIABLE TO MATCH YOUR DATA
# Check what variables are available in your phyloseq object:
cat("=== AVAILABLE SAMPLE VARIABLES ===\n")
if(!is.null(sample_data(ps))) {
  print(colnames(sample_data(ps)))
} else {
  cat("No sample data found in phyloseq object\n")
}

# MODIFY THIS LINE - replace "your_target_variable" with the actual column name
target_var <- "your_target_variable"  # CHANGE THIS TO YOUR ACTUAL TARGET VARIABLE

# Run the analysis
cat("\n=== RUNNING PCA-BASED MICROBIAL ANALYSIS ===\n")
microbial_results <- evaluate_microbial_pca(
  ps, 
  target_variable = target_var,
  test_prop = 0.3,
  variance_threshold = 0.8,    # Explain 80% of variance
  max_components = 10,         # Maximum PCs to use
  min_abundance_filter = 0.0001,   # Filter taxa with <0.01% mean abundance
  prevalence_threshold = 0.1       # Keep taxa present in >10% of samples
)

# Interpret the loadings with taxonomic information
interpret_microbial_loadings(microbial_results, ps, top_n = 10)

# Plot PCA
plot_microbial_pca(microbial_results)

###########################################################################
# RANDOM FOREST SPECIFIC ANALYSIS
###########################################################################

cat("\n=== RANDOM FOREST SPECIFIC ANALYSIS ===\n")

# Extract RF model and test data
rf_model <- microbial_results$models$rf
test_data <- microbial_results$test_data

# Get RF importance
rf_importance <- importance(rf_model)
cat("Random Forest Variable Importance:\n")
print(rf_importance)

# Make RF predictions (probabilities for class 1)
rf_predictions <- predict(rf_model, test_data, type = "prob")[,2]
cat("\nRF Predicted Probabilities for positive class (class 1):\n")
print(rf_predictions)

# Get RF confusion matrix with positive class = 1
rf_pred_class <- predict(rf_model, test_data, type = "class")
actual <- factor(test_data$gt, levels = c(0, 1))
rf_cm <- confusionMatrix(rf_pred_class, actual, positive = "1")

cat("\n=== RANDOM FOREST CONFUSION MATRIX (Positive Class = 1) ===\n")
print(rf_cm)

# Summary of RF performance
cat("\n=== RANDOM FOREST PERFORMANCE SUMMARY ===\n")
cat(sprintf("Accuracy: %.3f\n", rf_cm$overall['Accuracy']))
cat(sprintf("Kappa: %.3f\n", rf_cm$overall['Kappa']))
cat(sprintf("Sensitivity (True Positive Rate): %.3f\n", rf_cm$byClass['Sensitivity']))
cat(sprintf("Specificity (True Negative Rate): %.3f\n", rf_cm$byClass['Specificity']))
cat(sprintf("Precision (Positive Predictive Value): %.3f\n", rf_cm$byClass['Pos Pred Value']))
cat(sprintf("F1-Score: %.3f\n", rf_cm$byClass['F1']))

###########################################################################
# ADDITIONAL MICROBIAL-SPECIFIC ANALYSES
###########################################################################

# Get most important taxa from the PCA loadings
cat("\n=== MOST INFLUENTIAL TAXA ACROSS ALL PCs ===\n")
all_loadings <- abs(microbial_results$loadings)
taxa_importance <- rowSums(all_loadings)
top_taxa <- sort(taxa_importance, decreasing = TRUE)[1:20]

if(!is.null(tax_table(ps))) {
  tax_info <- as.data.frame(tax_table(ps))
  for(i in 1:length(top_taxa)) {
    taxa_name <- names(top_taxa)[i]
    importance_val <- top_taxa[i]
    
    if(taxa_name %in% rownames(tax_info)) {
      tax_row <- tax_info[taxa_name, ]
      # Get most specific taxonomy
      for(tax_level in rev(colnames(tax_info))) {
        if(!is.na(tax_row[[tax_level]]) && tax_row[[tax_level]] != "") {
          cat(sprintf("%2d. %s (%s): %.3f\n", i, taxa_name, tax_row[[tax_level]], importance_val))
          break
        }
      }
    } else {
      cat(sprintf("%2d. %s: %.3f\n", i, taxa_name, importance_val))
    }
  }
} else {
  print(top_taxa)
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved in 'microbial_results' object\n")