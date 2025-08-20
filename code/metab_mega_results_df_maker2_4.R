# ====================================
# REPEATED MODEL EVALUATION (10 RUNS)
# ====================================

library(randomForest)
library(caret)
library(dplyr)
library(doParallel)
library(pROC)
library(Boruta)

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

# ====================================
# HELPER FUNCTION FOR SINGLE MODEL EVALUATION
# ====================================
evaluate_single_model <- function(model, test_data) {
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
# FUNCTION TO RUN MODELS FOR ONE ITERATION
# ====================================
run_single_iteration <- function(seed_val) {
  
  cat("Running iteration", seed_val, "\n")
  
  # ====================================
  # TRAIN/TEST SPLIT
  # ====================================
  set.seed(seed_val)
  train_idx <- createDataPartition(fk_metabolom_gt_scaled$gt, p = 0.75, list = FALSE)
  train_data <- fk_metabolom_gt_scaled[train_idx, ]
  test_data <- fk_metabolom_gt_scaled[-train_idx, ]
  
  results_iteration <- list()
  
  # ====================================
  # MODEL 1: BORUTA FEATURE SELECTION + RF
  # ====================================
  tryCatch({
    set.seed(seed_val)
    boruta_result <- Boruta(gt ~ ., data = train_data, 
                           maxRuns = 100,  # Back to original
                           doTrace = 0)    # No output
    
    important_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
    
    if(length(important_features) > 0) {
      train_boruta <- train_data[, c(important_features, "gt")]
      test_boruta <- test_data[, c(important_features, "gt")]
      
      set.seed(seed_val)
      rf_boruta <- randomForest(
        gt ~ ., 
        data = train_boruta,
        mtry = max(1, floor(sqrt(length(important_features)))),
        ntree = 2000,  # Back to original
        importance = TRUE,
        nodesize = 3,
        replace = TRUE,
        classwt = c(1, 1)
      )
      
      results_iteration$Boruta_RF <- evaluate_single_model(rf_boruta, test_boruta)
      results_iteration$Boruta_RF$N_Features <- length(important_features)
    } else {
      # If no features selected, create dummy results
      results_iteration$Boruta_RF <- list(
        Accuracy = 0.5, Balanced_Accuracy = 0.5, Sensitivity = 0.5,
        Specificity = 0.5, Precision = 0.5, F1_Score = 0.5,
        Kappa = 0, AUC = 0.5, N_Features = 0
      )
    }
  }, error = function(e) {
    cat("Error in Boruta for iteration", seed_val, ":", e$message, "\n")
    results_iteration$Boruta_RF <- list(
      Accuracy = 0.5, Balanced_Accuracy = 0.5, Sensitivity = 0.5,
      Specificity = 0.5, Precision = 0.5, F1_Score = 0.5,
      Kappa = 0, AUC = 0.5, N_Features = 0
    )
  })
  
  # ====================================
  # MODEL 2: PROPERLY TUNED RF
  # ====================================
  tryCatch({
    tuneGrid <- expand.grid(
      mtry = c(floor(sqrt(ncol(train_data)-1)),
               floor(log2(ncol(train_data)-1)),
               floor((ncol(train_data)-1)/10),
               floor((ncol(train_data)-1)/20))
    )
    tuneGrid <- unique(tuneGrid[tuneGrid$mtry > 0, , drop = FALSE])
    
    ctrl <- trainControl(
      method = "repeatedcv",
      number = 5,  # Back to original  
      repeats = 3, # Back to original
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      savePredictions = "final",
      sampling = "up"
    )
    
    set.seed(seed_val)
    rf_tuned_proper <- train(
      gt ~ ., 
      data = train_data,
      method = "rf",
      tuneGrid = tuneGrid,
      trControl = ctrl,
      metric = "ROC",
      ntree = 1000,  # Back to original
      importance = TRUE,
      nodesize = 3,
      replace = TRUE
    )
    
    results_iteration$Properly_Tuned_RF <- evaluate_single_model(rf_tuned_proper, test_data)
    results_iteration$Properly_Tuned_RF$N_Features <- ncol(train_data) - 1
    results_iteration$Properly_Tuned_RF$Best_mtry <- rf_tuned_proper$bestTune$mtry
    
  }, error = function(e) {
    cat("Error in Tuned RF for iteration", seed_val, ":", e$message, "\n")
    results_iteration$Properly_Tuned_RF <- list(
      Accuracy = 0.5, Balanced_Accuracy = 0.5, Sensitivity = 0.5,
      Specificity = 0.5, Precision = 0.5, F1_Score = 0.5,
      Kappa = 0, AUC = 0.5, N_Features = ncol(train_data) - 1, Best_mtry = NA
    )
  })
  
  # ====================================
  # MODEL 3: IMPROVED ANTIFUNGAL RF
  # ====================================
  tryCatch({
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
    
    antifungal_data <- fk_metabolom_gt_scaled[, c(available_antifungal, "gt")]
    train_antifungal <- antifungal_data[train_idx, ]
    test_antifungal <- antifungal_data[-train_idx, ]
    
    predictor_cols <- names(train_antifungal)[names(train_antifungal) != "gt"]
    train_antifungal[predictor_cols] <- lapply(train_antifungal[predictor_cols], as.numeric)
    test_antifungal[predictor_cols] <- lapply(test_antifungal[predictor_cols], as.numeric)
    
    set.seed(seed_val)
    rf_antifungal_improved <- randomForest(
      gt ~ ., 
      data = train_antifungal,
      mtry = max(1, floor(sqrt(length(available_antifungal)))),
      ntree = 2000,  # Back to original
      importance = TRUE,
      nodesize = 2,
      replace = TRUE,
      classwt = c(1, 1)
    )
    
    results_iteration$Improved_Antifungal_RF <- evaluate_single_model(rf_antifungal_improved, test_antifungal)
    results_iteration$Improved_Antifungal_RF$N_Features <- length(available_antifungal)
    
  }, error = function(e) {
    cat("Error in Antifungal RF for iteration", seed_val, ":", e$message, "\n")
    results_iteration$Improved_Antifungal_RF <- list(
      Accuracy = 0.5, Balanced_Accuracy = 0.5, Sensitivity = 0.5,
      Specificity = 0.5, Precision = 0.5, F1_Score = 0.5,
      Kappa = 0, AUC = 0.5, N_Features = 33
    )
  })
  
  return(results_iteration)
}

# ====================================
# RUN 10 ITERATIONS
# ====================================
cat("\n=== RUNNING 10 ITERATIONS ===\n")

all_results <- list()
boruta_features_per_run <- list()  # Track Boruta features
seeds <- 100:109  # Different seeds for each iteration

for(i in 1:10) {
  cat("\n--- Iteration", i, "---\n")
  all_results[[i]] <- run_single_iteration(seeds[i])
}

# ====================================
# AGGREGATE RESULTS
# ====================================
cat("\n=== AGGREGATING RESULTS ===\n")

models <- c("Boruta_RF", "Properly_Tuned_RF", "Improved_Antifungal_RF")
metrics <- c("Accuracy", "Balanced_Accuracy", "Kappa", "AUC", "Sensitivity", "Specificity", "F1_Score")

# Initialize summary dataframe
summary_results <- data.frame(
  Model = character(),
  Metric = character(),
  Mean = numeric(),
  SD = numeric(),
  Min = numeric(),
  Max = numeric(),
  stringsAsFactors = FALSE
)

for(model in models) {
  cat("\nProcessing", model, "...\n")
  
  # Extract all values for each metric
  for(metric in metrics) {
    values <- sapply(1:10, function(i) {
      if(model %in% names(all_results[[i]]) && metric %in% names(all_results[[i]][[model]])) {
        return(all_results[[i]][[model]][[metric]])
      } else {
        return(NA)
      }
    })
    
    # Remove NA values
    values <- values[!is.na(values)]
    
    if(length(values) > 0) {
      summary_results <- rbind(summary_results, data.frame(
        Model = model,
        Metric = metric,
        Mean = round(mean(values), 3),
        SD = round(sd(values), 3),
        Min = round(min(values), 3),
        Max = round(max(values), 3),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Also get number of features (if available)
  n_features_values <- sapply(1:10, function(i) {
    if(model %in% names(all_results[[i]]) && "N_Features" %in% names(all_results[[i]][[model]])) {
      return(all_results[[i]][[model]][["N_Features"]])
    } else {
      return(NA)
    }
  })
  n_features_values <- n_features_values[!is.na(n_features_values)]
  
  if(length(n_features_values) > 0) {
    summary_results <- rbind(summary_results, data.frame(
      Model = model,
      Metric = "N_Features",
      Mean = round(mean(n_features_values), 1),
      SD = round(sd(n_features_values), 1),
      Min = min(n_features_values),
      Max = max(n_features_values),
      stringsAsFactors = FALSE
    ))
  }
}

# ====================================
# DISPLAY RESULTS
# ====================================
cat("\n============================================================\n")
cat("SUMMARY OF 10 ITERATIONS - MEAN ± SD\n")
cat("============================================================\n\n")

# Display in a nice format
for(model in models) {
  cat("=== ", model, " ===\n")
  model_results <- summary_results[summary_results$Model == model, ]
  
  for(i in 1:nrow(model_results)) {
    metric <- model_results$Metric[i]
    mean_val <- model_results$Mean[i]
    sd_val <- model_results$SD[i]
    min_val <- model_results$Min[i]
    max_val <- model_results$Max[i]
    
    cat(sprintf("%-18s: %6.3f ± %5.3f  (range: %6.3f - %6.3f)\n", 
                metric, mean_val, sd_val, min_val, max_val))
  }
  cat("\n")
}

# Create a comparison table for key metrics
cat("=== COMPARISON TABLE (Mean ± SD) ===\n")
key_metrics <- c("Accuracy", "Balanced_Accuracy", "Kappa")

comparison_table <- data.frame(Model = models)
for(metric in key_metrics) {
  metric_values <- character(length(models))
  for(j in 1:length(models)) {
    model_data <- summary_results[summary_results$Model == models[j] & summary_results$Metric == metric, ]
    if(nrow(model_data) > 0) {
      metric_values[j] <- paste0(sprintf("%.3f", model_data$Mean), " ± ", sprintf("%.3f", model_data$SD))
    } else {
      metric_values[j] <- "NA"
    }
  }
  comparison_table[[metric]] <- metric_values
}

print(comparison_table)

# Save detailed results
write.csv(summary_results, "repeated_evaluation_summary.csv", row.names = FALSE)

# Save individual iteration results for further analysis
detailed_results <- data.frame()
for(i in 1:10) {
  for(model in models) {
    if(model %in% names(all_results[[i]])) {
      model_result <- all_results[[i]][[model]]
      row_data <- data.frame(
        Iteration = i,
        Model = model,
        Accuracy = model_result$Accuracy,
        Balanced_Accuracy = model_result$Balanced_Accuracy,
        Kappa = model_result$Kappa,
        AUC = model_result$AUC,
        Sensitivity = model_result$Sensitivity,
        Specificity = model_result$Specificity,
        F1_Score = model_result$F1_Score,
        N_Features = model_result$N_Features,
        stringsAsFactors = FALSE
      )
      detailed_results <- rbind(detailed_results, row_data)
    }
  }
}

write.csv(detailed_results, "repeated_evaluation_detailed.csv", row.names = FALSE)

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

cat("\nResults saved to:\n")
cat("- repeated_evaluation_summary.csv (summary statistics)\n")
cat("- repeated_evaluation_detailed.csv (all individual results)\n")
cat("\n=== ANALYSIS COMPLETE ===\n")


# Read the detailed results file
detailed <- read.csv("repeated_evaluation_detailed.csv")

# Look at Boruta results
boruta_results <- detailed[detailed$Model == "Boruta_RF", ]
print("Boruta results per iteration:")
print(boruta_results[, c("Iteration", "N_Features", "Balanced_Accuracy", "AUC")])

# Check which runs had features vs performance
cat("\n=== ANALYSIS ===\n")
good_runs <- boruta_results[boruta_results$N_Features > 0, ]
bad_runs <- boruta_results[boruta_results$N_Features == 0, ]

cat("Runs with features found:\n")
if(nrow(good_runs) > 0) {
  print(good_runs[, c("Iteration", "N_Features", "Balanced_Accuracy")])
  cat("Average performance when features found:", round(mean(good_runs$Balanced_Accuracy), 3), "\n")
}

cat("\nRuns with NO features found:\n")
if(nrow(bad_runs) > 0) {
  print(bad_runs[, c("Iteration", "N_Features", "Balanced_Accuracy")])  
  cat("Average performance when no features found:", round(mean(bad_runs$Balanced_Accuracy), 3), "\n")
}