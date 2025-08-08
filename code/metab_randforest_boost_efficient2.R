# this script runs a randomforest model on the metabolomic data

library(dplyr)

###change so that 1 = present, 0 = absent****

# load microbial metabolomic and gt data
fk_metabolom_gt <- readRDS("../data/fk_metabolomics_gt_noNA.RDS")

nrow(fk_metabolom_gt) # 677 compounds
#View(fk_metabolom_gt)

# count the zeroes
sum(fk_metabolom_gt ==0, na.rm = TRUE)
dim(fk_metabolom_gt)

## data prep 
fk_metabolom_gt_t <- t(fk_metabolom_gt) # transpose
colnames(fk_metabolom_gt_t) <- fk_metabolom_gt_t[1, ] # set the first row as column names
fk_metabolom_gt_t <- fk_metabolom_gt_t[-1, ] # remove the first row (now column names)
fk_metabolom_gt_t <- as.data.frame(fk_metabolom_gt_t) # convert to data frame
colnames(fk_metabolom_gt_t)[ncol(fk_metabolom_gt_t)] <- "gt" # rename the last column to 'gt'
fk_metabolom_gt_t$gt <- trimws(as.character(fk_metabolom_gt_t$gt))
fk_metabolom_gt_t$gt <- as.factor(fk_metabolom_gt_t$gt)
levels(fk_metabolom_gt_t$gt) <- c("GT_absent", "GT_present")
levels(fk_metabolom_gt_t$gt)

colnames(fk_metabolom_gt_t) <- make.names(colnames(fk_metabolom_gt_t), unique = TRUE)# Clean colnames by replacing spaces with underscores
str(fk_metabolom_gt_t[,1:5]) # data aren't scaled or numeric
fk_metabolom_gt_t[ , -c(ncol(fk_metabolom_gt_t))] <- lapply( # make the peak areas numeric
  fk_metabolom_gt_t[ , -c(1, ncol(fk_metabolom_gt_t))],
  function(x) as.numeric(as.character(x))
)
hist(fk_metabolom_gt_t[,78], breaks = 30, main = "Histogram", xlab = "Value") # check for skew- a bit righty
plot(density(log(fk_metabolom_gt_t[,78]), na.rm = TRUE), main = "Density") # does log help? think so

# log data
metab_numeric <- fk_metabolom_gt_t[, sapply(fk_metabolom_gt_t, is.numeric)] # Select only numeric metabolite columns
metab_scaled <- as.data.frame(log(metab_numeric + 1e-6)) # log transform
# THINK ABOUT Z SCALING IT HERE
non_numeric <- fk_metabolom_gt_t[, !sapply(fk_metabolom_gt_t, is.numeric)] 
fk_metabolom_gt_scaled <- cbind(non_numeric, metab_scaled)# add back non-numeric columns
colnames(fk_metabolom_gt_scaled)[1]<- "gt"
View(fk_metabolom_gt_scaled)






















###### 

# new approach- BOOSTING (?)
# Load required libraries
library(gbm)
library(caret)
library(glmnet)
library(randomForest)
library(dplyr)

preprocess_metabolomics_stable <- function(data, target_col, train_idx, 
                                         stability_runs = 10, 
                                         stability_threshold = 0.7,
                                         max_features = 5,
                                         p_value_threshold = 0.01,
                                         effect_size_threshold = 0.5) {
  
  # Separate features and target
  if(is.character(target_col)) {
    target <- data[[target_col]]
    features <- data[, !names(data) %in% target_col]
  } else {
    target <- data[, target_col]
    features <- data[, -target_col]
  }
  
  cat("=== STABLE PREPROCESSING (TRAINING DATA ONLY) ===\n")
  cat("Original data dimensions:", dim(features), "\n")
  cat("Number of zeros:", sum(features == 0), "\n")
  
  # STEP 1: Adaptive feature filtering based on TRAINING data only
  train_features <- features[train_idx, ]
  n_features_original <- ncol(train_features)
  
  # Adaptive zero filtering - start strict, relax if needed
  zero_thresholds <- c(0.3, 0.5, 0.7, 0.8)  # Try increasingly relaxed thresholds
  
  for(zero_thresh in zero_thresholds) {
    zero_prop_features <- sapply(train_features, function(x) sum(x == 0) / length(x))
    keep_zero_filter <- zero_prop_features <= zero_thresh
    
    if(sum(keep_zero_filter) >= 50) {  # Keep at least 50 features for next step
      cat("Using zero threshold:", zero_thresh, "- kept", sum(keep_zero_filter), "features\n")
      break
    }
  }
  
  # Safety check
  if(sum(keep_zero_filter) == 0) {
    cat("Warning: All features have high zero proportion. Keeping top 100 features by zero proportion.\n")
    keep_zero_filter <- order(zero_prop_features)[1:min(100, length(zero_prop_features))]
    zero_filtered_names <- names(train_features)[keep_zero_filter]
  } else {
    zero_filtered_names <- names(train_features)[keep_zero_filter]
  }
  
  # Adaptive variance filtering
  train_features_zero_filtered <- train_features[, zero_filtered_names]
  train_var <- sapply(train_features_zero_filtered, var, na.rm = TRUE)
  
  # Remove features with zero or very low variance
  valid_var <- train_var > 0 & !is.na(train_var) & is.finite(train_var)
  
  if(sum(valid_var) == 0) {
    cat("Warning: No features with valid variance. Keeping all zero-filtered features.\n")
    var_filtered_names <- zero_filtered_names
  } else {
    train_var_valid <- train_var[valid_var]
    var_threshold <- quantile(train_var_valid, 0.1, na.rm = TRUE)  # Start with 10th percentile
    keep_var_filter <- train_var > var_threshold & valid_var
    
    # Ensure we keep enough features
    if(sum(keep_var_filter) < 20) {
      var_threshold <- quantile(train_var_valid, 0.05, na.rm = TRUE)  # Relax to 5th percentile
      keep_var_filter <- train_var > var_threshold & valid_var
    }
    
    var_filtered_names <- names(train_features_zero_filtered)[keep_var_filter]
  }
  
  # Adaptive SNR filtering - only if we have enough features
  if(length(var_filtered_names) >= 30) {
    train_features_var_filtered <- train_features[, var_filtered_names]
    train_means <- sapply(train_features_var_filtered, mean, na.rm = TRUE)
    train_sds <- sapply(train_features_var_filtered, sd, na.rm = TRUE)
    
    # Calculate SNR, handling division by zero
    snr <- ifelse(train_sds > 0, abs(train_means) / train_sds, 0)
    snr[!is.finite(snr)] <- 0
    
    if(sum(snr > 0) >= 10) {  # Only apply SNR filter if we have enough valid SNR values
      snr_threshold <- quantile(snr[snr > 0], 0.3, na.rm = TRUE)  # Keep top 70%
      keep_snr_filter <- snr > snr_threshold
      final_keep_features <- var_filtered_names[keep_snr_filter]
    } else {
      final_keep_features <- var_filtered_names
    }
  } else {
    final_keep_features <- var_filtered_names
  }
  
  # Final safety check
  if(length(final_keep_features) == 0) {
    cat("Warning: All filtering removed all features. Keeping top 50 features by variance.\n")
    all_vars <- sapply(train_features, var, na.rm = TRUE)
    all_vars[!is.finite(all_vars)] <- 0
    final_keep_features <- names(sort(all_vars, decreasing = TRUE))[1:min(50, length(all_vars))]
  }
  
  features_filtered <- features[, final_keep_features, drop = FALSE]
  train_features_filtered <- train_features[, final_keep_features, drop = FALSE]
  
  cat("Removed", n_features_original - length(zero_filtered_names), "features with high zero proportion\n")
  cat("Removed", length(zero_filtered_names) - length(var_filtered_names), "low variance features\n")
  cat("Removed", length(var_filtered_names) - length(final_keep_features), "low SNR features\n")
  cat("Remaining features after filtering:", ncol(features_filtered), "\n")
  
  # STEP 2: STABLE feature selection using bootstrap stability
  cat("Performing stable feature selection with", stability_runs, "bootstrap runs\n")
  
  train_target <- target[train_idx]
  n_train <- length(train_idx)
  
  # Storage for stability analysis
  feature_selection_matrix <- matrix(0, nrow = stability_runs, ncol = ncol(train_features_filtered))
  colnames(feature_selection_matrix) <- names(train_features_filtered)
  
  # Storage for effect sizes and p-values across runs
  all_p_values <- matrix(NA, nrow = stability_runs, ncol = ncol(train_features_filtered))
  all_effect_sizes <- matrix(NA, nrow = stability_runs, ncol = ncol(train_features_filtered))
  colnames(all_p_values) <- names(train_features_filtered)
  colnames(all_effect_sizes) <- names(train_features_filtered)
  
  # Bootstrap stability selection
  for(run in 1:stability_runs) {
    set.seed(run + 1000)  # Consistent but different seeds
    
    # Bootstrap sample from training data
    boot_idx <- sample(1:n_train, size = floor(n_train * 0.8), replace = TRUE)
    boot_features <- train_features_filtered[boot_idx, ]
    boot_target <- train_target[boot_idx]
    
    # Calculate statistics for this bootstrap
    run_p_values <- numeric(ncol(boot_features))
    run_effect_sizes <- numeric(ncol(boot_features))
    
    for(i in 1:ncol(boot_features)) {
      # T-test with error handling
      tryCatch({
        if(length(unique(boot_target)) > 1 && var(boot_features[, i], na.rm = TRUE) > 0) {
          test_result <- t.test(boot_features[, i] ~ boot_target)
          run_p_values[i] <- test_result$p.value
        } else {
          run_p_values[i] <- 1  # No difference if no variance or only one group
        }
      }, error = function(e) {
        run_p_values[i] <- 1
      })
      
      # Effect size (Cohen's d) with error handling
      tryCatch({
        group0 <- boot_features[boot_target == 0, i]
        group1 <- boot_features[boot_target == 1, i]
        
        if(length(group0) > 1 && length(group1) > 1) {
          var0 <- var(group0, na.rm = TRUE)
          var1 <- var(group1, na.rm = TRUE)
          
          if(var0 > 0 || var1 > 0) {
            pooled_sd <- sqrt(((length(group0)-1)*var0 + (length(group1)-1)*var1) / 
                             (length(group0) + length(group1) - 2))
            if(pooled_sd > 0) {
              cohens_d <- abs(mean(group1, na.rm = TRUE) - mean(group0, na.rm = TRUE)) / pooled_sd
              run_effect_sizes[i] <- cohens_d
            } else {
              run_effect_sizes[i] <- 0
            }
          } else {
            run_effect_sizes[i] <- 0
          }
        } else {
          run_effect_sizes[i] <- 0
        }
      }, error = function(e) {
        run_effect_sizes[i] <- 0
      })
    }
    
    # Store results
    all_p_values[run, ] <- run_p_values
    all_effect_sizes[run, ] <- run_effect_sizes
    
    # Select features for this run (both p-value and effect size criteria)
    significant_features <- (run_p_values < p_value_threshold) & 
                           (run_effect_sizes > effect_size_threshold)
    
    if(sum(significant_features) > max_features) {
      # If too many significant features, select top ones by combined score
      combined_score <- -log10(run_p_values) * run_effect_sizes
      top_features_idx <- order(combined_score, decreasing = TRUE)[1:max_features]
      significant_features <- rep(FALSE, length(significant_features))
      significant_features[top_features_idx] <- TRUE
    }
    
    feature_selection_matrix[run, significant_features] <- 1
  }
  
  # Calculate stability scores for each feature
  stability_scores <- colMeans(feature_selection_matrix)
  
  # Select features that appear in at least stability_threshold of runs
  stable_features <- names(stability_scores)[stability_scores >= stability_threshold]
  
  # If no features meet stability threshold, select most stable ones
  if(length(stable_features) == 0) {
    n_to_select <- min(max_features, length(stability_scores))
    stable_features <- names(sort(stability_scores, decreasing = TRUE))[1:n_to_select]
    cat("Warning: No features met stability threshold of", stability_threshold, 
        ". Selected", length(stable_features), "most stable features.\n")
  } else if(length(stable_features) > max_features) {
    # If too many stable features, select the most stable ones
    stable_features <- names(sort(stability_scores[stable_features], decreasing = TRUE))[1:max_features]
    cat("Selected top", length(stable_features), "most stable features from", 
        sum(stability_scores >= stability_threshold), "candidates.\n")
  }
  
  # Final validation: ensure selected features have reasonable average statistics
  if(length(stable_features) > 0) {
    avg_p_values <- colMeans(all_p_values, na.rm = TRUE)
    avg_effect_sizes <- colMeans(all_effect_sizes, na.rm = TRUE)
    
    # More lenient final filtering to avoid empty selection
    p_thresh_relaxed <- min(p_value_threshold * 3, 0.1)  # Relax p-value threshold
    effect_thresh_relaxed <- max(effect_size_threshold * 0.5, 0.2)  # Relax effect size threshold
    
    valid_features <- stable_features[
      !is.na(avg_p_values[stable_features]) & 
      !is.na(avg_effect_sizes[stable_features]) &
      is.finite(avg_p_values[stable_features]) &
      is.finite(avg_effect_sizes[stable_features])
    ]
    
    if(length(valid_features) > 0) {
      final_selected <- valid_features[
        (avg_p_values[valid_features] < p_thresh_relaxed) &
        (avg_effect_sizes[valid_features] > effect_thresh_relaxed)
      ]
      
      if(length(final_selected) == 0) {
        final_selected <- valid_features[1:min(max_features, length(valid_features))]
        cat("Warning: Statistical filtering too strict. Using most stable features.\n")
      }
    } else {
      final_selected <- stable_features[1:min(max_features, length(stable_features))]
      cat("Warning: Statistical validation failed. Using most stable features.\n")
    }
  } else {
    # Emergency fallback - select features with best individual statistics
    cat("Warning: No stable features found. Using emergency feature selection.\n")
    avg_p_values <- colMeans(all_p_values, na.rm = TRUE)
    avg_effect_sizes <- colMeans(all_effect_sizes, na.rm = TRUE)
    
    # Remove NAs and infinite values
    valid_idx <- !is.na(avg_p_values) & !is.na(avg_effect_sizes) & 
                 is.finite(avg_p_values) & is.finite(avg_effect_sizes)
    
    if(sum(valid_idx) > 0) {
      combined_score <- -log10(avg_p_values + 1e-10) * avg_effect_sizes
      combined_score[!valid_idx] <- -Inf
      final_selected <- names(sort(combined_score, decreasing = TRUE))[1:min(max_features, sum(valid_idx))]
    } else {
      # Ultimate fallback - select first few features
      final_selected <- names(train_features_filtered)[1:min(max_features, ncol(train_features_filtered))]
    }
  }
  
  cat("Selected", length(final_selected), "stable features\n")
  cat("Feature stability scores:\n")
  selected_stability <- stability_scores[final_selected]
  print(sort(selected_stability, decreasing = TRUE))
  
  # Apply selection to full dataset
  features_selected <- features_filtered[, final_selected]
  train_features_selected <- train_features_filtered[, final_selected]
  
  # STEP 3: Robust zero imputation using TRAINING statistics
  features_imputed <- features_selected
  zero_replacement_values <- numeric(length(final_selected))
  names(zero_replacement_values) <- final_selected
  
  for(col in final_selected) {
    zeros_all <- features_selected[[col]] == 0
    if(sum(zeros_all) > 0) {
      # Calculate replacement value from TRAINING data only
      train_nonzero <- train_features_selected[[col]][train_features_selected[[col]] > 0]
      if(length(train_nonzero) > 0) {
        # Use median of smallest 10% instead of minimum for robustness
        min_val <- quantile(train_nonzero, 0.1, na.rm = TRUE)
        replacement_val <- min_val / 2
        features_imputed[[col]][zeros_all] <- replacement_val
        zero_replacement_values[col] <- replacement_val
      }
    }
  }
  
  # STEP 4: Log transform and robust scaling using TRAINING statistics
  features_log <- log10(features_imputed + 1)
  
  # Calculate robust scaling parameters from training data only
  train_log <- features_log[train_idx, ]
  train_medians <- sapply(train_log, median, na.rm = TRUE)  # Use median instead of mean
  train_mads <- sapply(train_log, mad, na.rm = TRUE)       # Use MAD instead of SD
  
  # Apply scaling to ALL data using training parameters
  features_scaled <- features_log
  for(i in 1:ncol(features_log)) {
    if(train_mads[i] > 0) {
      features_scaled[, i] <- (features_log[, i] - train_medians[i]) / train_mads[i]
    }
  }
  
  # Return processed data and preprocessing parameters
  processed_data <- data.frame(features_scaled)
  processed_data$gt <- target
  
  preprocessing_params <- list(
    selected_features = final_selected,
    train_medians = train_medians,
    train_mads = train_mads,
    zero_replacement_values = zero_replacement_values,
    stability_scores = stability_scores[final_selected],
    avg_p_values = avg_p_values[final_selected],
    avg_effect_sizes = avg_effect_sizes[final_selected]
  )
  
  cat("Final processed data dimensions:", dim(processed_data), "\n")
  cat("Average p-values of selected features:", round(avg_p_values[final_selected], 4), "\n")
  cat("Average effect sizes of selected features:", round(avg_effect_sizes[final_selected], 3), "\n\n")
  
  return(list(
    data = processed_data,
    params = preprocessing_params,
    selected_features = final_selected
  ))
}

# Updated evaluation function using stable preprocessing
evaluate_metabolomics_models_stable <- function(data, target_col = "gt", test_prop = 0.3, seed = 123,
                                               stability_runs = 10, stability_threshold = 0.7,
                                               max_features = 5, p_value_threshold = 0.01,
                                               effect_size_threshold = 0.5) {
  
  set.seed(seed)
  
  # Single train/test split
  train_idx <- createDataPartition(data[[target_col]], p = 1 - test_prop, list = FALSE)
  test_idx <- setdiff(1:nrow(data), train_idx)
  
  cat("Train samples:", length(train_idx), "Test samples:", length(test_idx), "\n")
  
  # Preprocess data with stable feature selection
  preprocessed <- preprocess_metabolomics_stable(
    data, target_col, train_idx,
    stability_runs = stability_runs,
    stability_threshold = stability_threshold,
    max_features = max_features,
    p_value_threshold = p_value_threshold,
    effect_size_threshold = effect_size_threshold
  )
  
  processed_data <- preprocessed$data
  
  # Split processed data
  train_data <- processed_data[train_idx, ]
  test_data <- processed_data[test_idx, ]
  
  cat("Final training data dimensions:", dim(train_data), "\n")
  cat("Final test data dimensions:", dim(test_data), "\n\n")
  
  # Train models
  models <- train_multiple_models(train_data, target_col)
  
  # Evaluate on test set
  results <- evaluate_on_test_set(models, test_data, target_col)
  
  return(list(
    models = models,
    results = results,
    test_data = test_data,
    preprocessing_params = preprocessed$params,
    selected_features = preprocessed$selected_features
  ))
}

# test feature stability across multiple runs
test_feature_stability_improved <- function(data, target_col = "gt", n_runs = 10, 
                                           stability_threshold = 0.7, max_features = 5) {
  feature_lists <- list()
  stability_scores_all <- list()
  
  for(i in 1:n_runs) {
    set.seed(i * 42)  # More separated seeds
    
    # Different train/test split each time
    train_idx <- createDataPartition(data[[target_col]], p = 0.7, list = FALSE)
    
    # Run stable preprocessing and feature selection
    preprocessed <- preprocess_metabolomics_stable(
      data, target_col, train_idx,
      stability_runs = 10,  # Inner stability runs
      stability_threshold = stability_threshold,
      max_features = max_features
    )
    
    selected_features <- preprocessed$selected_features
    stability_scores <- preprocessed$params$stability_scores
    
    # Check if fewer features than expected were selected
    if (length(selected_features) < max_features) {
      cat("Warning: Run", i, "selected", length(selected_features), "features (less than expected).\n")
      cat("    Selected features:", selected_features, "\n")
    }
    
    # Check to ensure we don't access out of bounds
    if (length(selected_features) > 1) {
      feature_lists[[i]] <- selected_features
      stability_scores_all[[i]] <- stability_scores
    } else {
      cat("Warning: Run", i, "selected 0 or 1 feature. Skipping this run.\n")
    }
    
    # Debugging logs for selected features and stability scores
    cat("Run", i, "- Selected features:", selected_features, "\n")
    cat("    Stability scores:", round(stability_scores, 3), "\n")
  }
  
  # Check overlap between feature sets across runs
  all_features <- unique(unlist(feature_lists))
  feature_counts <- table(unlist(feature_lists))
  
  cat("\n=== IMPROVED FEATURE STABILITY ANALYSIS ===\n")
  cat("Features appearing in multiple runs:\n")
  print(sort(feature_counts, decreasing = TRUE))
  
  # Calculate Jaccard similarity between feature sets from different runs
  similarities <- numeric()
  for(i in 1:(length(feature_lists)-1)) {
    for(j in (i+1):length(feature_lists)) {
      intersection <- length(intersect(feature_lists[[i]], feature_lists[[j]]))
      union <- length(union(feature_lists[[i]], feature_lists[[j]]))
      similarities <- c(similarities, intersection/union)
    }
  }
  
  # Log the average similarity and standard deviation
  cat("Mean Jaccard similarity between feature sets:", round(mean(similarities), 3), "\n")
  cat("SD of similarities:", round(sd(similarities), 3), "\n")
  
  # Calculate consistency rate (how often features appear across runs)
  consistency_rate <- feature_counts / n_runs
  cat("\nMost consistent features (selection rate):\n")
  print(sort(consistency_rate, decreasing = TRUE))
  
  # Return a comprehensive list of results
  return(list(
    feature_lists = feature_lists,
    feature_counts = feature_counts,
    consistency_rates = consistency_rate,
    mean_similarity = mean(similarities),
    stability_scores_all = stability_scores_all
  ))
}



# Usage example with stable preprocessing
# Convert gt to numeric (0/1) if not already
boost_data <- fk_metabolom_gt_t %>%
  mutate(gt = ifelse(gt == "GT_present", 1,
              ifelse(gt == "GT_absent", 0, gt))) %>%
  mutate(gt = as.numeric(gt))

# Stable evaluation with more stringent criteria



# more lenient evaluation
stable_results <- evaluate_metabolomics_models_stable(
  boost_data, 
  target_col = "gt", 
  test_prop = 0.3,
  stability_runs = 20,        # More runs
  stability_threshold = 0.6,   # Lower threshold
  max_features = 3,           # Fewer features
  p_value_threshold = 0.05,   # Less stringent p-value
  effect_size_threshold = 0.4  # Smaller effect size
)

stable_results <- evaluate_metabolomics_models_stable(
  boost_data, 
  target_col = "gt", 
  test_prop = 0.3,
  stability_runs = 50,        # More runs for better stability
  stability_threshold = 0.4,   # Lower threshold but require consistency
  max_features = 4,           # Optimal for your sample size
  p_value_threshold = 0.05,   # Standard significance
  effect_size_threshold = 0.5  # Medium effect size
)

# Test the improved stability
stability_results_improved <- test_feature_stability_improved(
  boost_data, 
  target_col = "gt", 
  n_runs = 10,
  stability_threshold = 0.8,
  max_features = 5
)
