# Script to run random forest classification on phyloseq object to predict gt status

# Packages ----------------------------------
library(tidyverse)
library(phyloseq)
library(ranger)
library(mlr)
library(tuneRanger)
library(vegan)

# Functions ----------------------------------
preparePhyloseqForRF <- function(phyloseq_obj, target_variable) {
    # Extract sample data
    sample_df <- data.frame(sample_data(phyloseq_obj))
    
    # Extract OTU/ESV table and transpose (samples as rows)
    otu_df <- data.frame(t(otu_table(phyloseq_obj)))
    
    # Ensure sample names match
    if(!all(rownames(sample_df) == rownames(otu_df))) {
        common_samples <- intersect(rownames(sample_df), rownames(otu_df))
        sample_df <- sample_df[common_samples, ]
        otu_df <- otu_df[common_samples, ]
    }
    
    # Combine sample data with OTU abundances
    combined_df <- cbind(sample_df, otu_df)
    
    # Remove samples with missing target variable
    combined_df <- combined_df[!is.na(combined_df[[target_variable]]), ]
    
    # Convert target to factor if it isn't already
    combined_df[[target_variable]] <- as.factor(combined_df[[target_variable]])
    
    return(combined_df)
}

getBalancedSplit <- function(data, target_col, test_proportion = 0.2) {
    # Stratified sampling to maintain class balance
    library(caret)
    
    target_factor <- data[[target_col]]
    train_indices <- createDataPartition(target_factor, p = 1 - test_proportion, list = FALSE)
    
    return(list(
        "train" = data[train_indices, ],
        "test" = data[-train_indices, ]
    ))
}

# Custom metrics for classification
accuracy <- function(preds, actual) {
    mean(preds == actual)
}

# Prepare data ----------------------------------
cat("Preparing phyloseq data for Random Forest...\n")

# Load your phyloseq object (assuming it's already loaded as ps_highdiv_absolute)
ps_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")

# Check the phyloseq object
cat("Phyloseq object summary:\n")
print(ps_highdiv_absolute)

# Prepare data for Random Forest
rf_data <- preparePhyloseqForRF(ps_highdiv_absolute, "gt")

cat("Data dimensions:", dim(rf_data), "\n")
cat("GT status distribution:\n")
print(table(rf_data$gt))

# Optional: Filter low abundance ESVs (recommended for large datasets)
# Keep ESVs present in at least 5% of samples or with total abundance > 10
esv_cols <- grep("^ASV|^OTU|^ESV", colnames(rf_data), value = TRUE)
if(length(esv_cols) == 0) {
    # If no standard prefixes, assume all numeric columns after sample data are ESVs
    sample_data_cols <- ncol(data.frame(sample_data(ps_highdiv_absolute)))
    esv_cols <- colnames(rf_data)[(sample_data_cols + 1):ncol(rf_data)]
}

cat("Number of ESVs:", length(esv_cols), "\n")

# Filter ESVs (optional - remove if you want to keep all)
esv_prevalence <- colSums(rf_data[esv_cols] > 0) / nrow(rf_data)
esv_abundance <- colSums(rf_data[esv_cols])
keep_esvs <- names(esv_prevalence)[esv_prevalence >= 0.05 | esv_abundance >= 10]

cat("ESVs after filtering:", length(keep_esvs), "\n")

# Create final dataset with filtered ESVs
sample_cols <- setdiff(colnames(rf_data), esv_cols)
rf_data_filtered <- rf_data[, c(sample_cols, keep_esvs)]

# Remove any columns with zero variance or near-zero variance
nzv_cols <- nearZeroVar(rf_data_filtered[, -which(colnames(rf_data_filtered) == "gt")])
if(length(nzv_cols) > 0) {
    rf_data_filtered <- rf_data_filtered[, -nzv_cols]
    cat("Removed", length(nzv_cols), "near-zero variance features\n")
}

# Set parameters ----------------------------------
tuneRandomForest <- FALSE  # Set to TRUE if you want to tune hyperparameters
n_trees <- 1000
test_proportion <- 0.2
n_cv_folds <- 5

# Prepare cross-validation splits ----------------------------------
cat("Creating cross-validation splits...\n")

cv_splits <- list()
for(i in 1:n_cv_folds) {
    set.seed(123 + i)  # For reproducibility
    cv_splits[[i]] <- getBalancedSplit(rf_data_filtered, "gt", test_proportion)
}

# Hyperparameter tuning (optional) ----------------------------------
if (tuneRandomForest) {
    cat("Tuning Random Forest hyperparameters...\n")
    
    # Create task for tuning
    task <- makeClassifTask(data = cv_splits[[1]]$train, target = "gt")
    
    # Estimate tuning time
    estimateTimeTuneRanger(task)
    
    # Tune hyperparameters
    tuning_result <- tuneRanger(
        task, 
        measure = list(acc),  # Using accuracy for classification
        num.trees = n_trees,
        num.threads = 4, 
        iters = 50,  # Reduced iterations for faster tuning
        save.file.path = NULL
    )
    
    print("Tuning results:")
    print(tuning_result)
    
    # Extract optimal parameters
    optimal_mtry <- tuning_result$recommended.pars$mtry
    optimal_min_node_size <- tuning_result$recommended.pars$min.node.size
    optimal_sample_fraction <- tuning_result$recommended.pars$sample.fraction
    
} else {
    # Use default parameters (you can adjust these based on your data)
    optimal_mtry <- floor(sqrt(ncol(rf_data_filtered) - 1))  # Square root of features
    optimal_min_node_size <- 1  # Default for classification
    optimal_sample_fraction <- 0.632  # Default bootstrap sample size
}

cat("Using parameters: mtry =", optimal_mtry, 
    ", min.node.size =", optimal_min_node_size,
    ", sample.fraction =", optimal_sample_fraction, "\n")

# Fit Random Forest models with cross-validation ----------------------------------
cat("Fitting Random Forest models...\n")

predictions <- c()
observations <- c()
cv_accuracies <- c()
feature_importance_list <- list()

for (i in 1:length(cv_splits)) {
    cat("Processing CV fold", i, "of", length(cv_splits), "\n")
    
    # Fit Random Forest
    rf_model <- ranger(
        gt ~ .,
        data = cv_splits[[i]]$train, 
        importance = 'permutation',
        probability = TRUE,  # For classification probabilities
        mtry = optimal_mtry,
        min.node.size = optimal_min_node_size,
        sample.fraction = optimal_sample_fraction,
        num.trees = n_trees,
        seed = 123 + i
    )
    
    # Make predictions
    rf_predictions <- predict(rf_model, cv_splits[[i]]$test)
    pred_classes <- rf_predictions$predictions
    
    # Store predictions and observations
    predictions <- c(predictions, pred_classes)
    observations <- c(observations, cv_splits[[i]]$test$gt)
    
    # Calculate accuracy for this fold
    fold_accuracy <- accuracy(pred_classes, cv_splits[[i]]$test$gt)
    cv_accuracies <- c(cv_accuracies, fold_accuracy)
    
    # Store feature importance
    feature_importance_list[[i]] <- rf_model$variable.importance
    
    cat("Fold", i, "accuracy:", round(fold_accuracy, 3), "\n")
}

# Overall results ----------------------------------
cat("\n=== CROSS-VALIDATION RESULTS ===\n")
cat("Mean CV Accuracy:", round(mean(cv_accuracies), 3), "±", round(sd(cv_accuracies), 3), "\n")
cat("Individual fold accuracies:", round(cv_accuracies, 3), "\n")

# Confusion matrix
conf_matrix <- table(Predicted = predictions, Actual = observations)
print("Confusion Matrix:")
print(conf_matrix)

# Calculate additional metrics
precision <- diag(conf_matrix) / rowSums(conf_matrix)
recall <- diag(conf_matrix) / colSums(conf_matrix)
f1_score <- 2 * (precision * recall) / (precision + recall)

cat("\nPer-class metrics:\n")
for(class in names(precision)) {
    cat(class, "- Precision:", round(precision[class], 3), 
        ", Recall:", round(recall[class], 3), 
        ", F1:", round(f1_score[class], 3), "\n")
}

# Feature importance analysis ----------------------------------
cat("\n=== FEATURE IMPORTANCE ANALYSIS ===\n")

# Combine importance scores across CV folds
importance_df <- bind_rows(lapply(feature_importance_list, function(x) data.frame(t(x)))) %>%
    pivot_longer(cols = everything(), values_to = "Importance", names_to = "Feature")

# Calculate mean importance and create plot
importance_summary <- importance_df %>%
    group_by(Feature) %>%
    summarise(
        mean_importance = mean(Importance),
        sd_importance = sd(Importance),
        .groups = 'drop'
    ) %>%
    arrange(desc(mean_importance))

# Show top 20 most important features
cat("Top 20 most important features:\n")
print(head(importance_summary, 20))

# Create importance plot
library(ggplot2)

# Plot top 30 features
top_features <- head(importance_summary, 30)

importance_plot <- importance_summary %>%
    slice_max(mean_importance, n = 30) %>%
    mutate(Feature = fct_reorder(Feature, mean_importance)) %>%
    ggplot(aes(x = mean_importance, y = Feature)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(xmin = mean_importance - sd_importance, 
                      xmax = mean_importance + sd_importance),
                  width = 0.2, alpha = 0.7) +
    labs(
        title = "Top 30 Feature Importances for GT Status Prediction",
        x = "Mean Permutation Importance",
        y = "Feature"
    ) +
    theme_classic() +
    theme(
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 12, hjust = 0.5)
    )

print(importance_plot)

# Save results ----------------------------------
cat("\n=== SAVING RESULTS ===\n")

# Save importance results
write_csv(importance_summary, "gt_prediction_feature_importance.csv")

# Save prediction results
results_df <- data.frame(
    Sample = names(observations),
    Predicted = predictions,
    Actual = observations,
    Correct = predictions == observations
)
write_csv(results_df, "gt_prediction_results.csv")

# Save model summary
model_summary <- data.frame(
    Metric = c("Mean_CV_Accuracy", "SD_CV_Accuracy", "Number_of_Features", 
               "Number_of_Samples", "Number_of_Trees", "mtry", "min_node_size"),
    Value = c(mean(cv_accuracies), sd(cv_accuracies), ncol(rf_data_filtered) - 1,
              nrow(rf_data_filtered), n_trees, optimal_mtry, optimal_min_node_size)
)
write_csv(model_summary, "gt_prediction_model_summary.csv")

cat("Results saved to CSV files\n")
cat("Script completed successfully!\n")

# Optional: Identify taxonomic patterns in important ESVs ----------------------------------
if(exists("ps_highdiv_absolute") && any(grepl("^ASV|^OTU|^ESV", top_features$Feature))) {
    cat("\n=== TAXONOMIC ANALYSIS OF IMPORTANT ESVs ===\n")
    
    # Get taxonomy table
    tax_table_df <- data.frame(tax_table(ps_highdiv_absolute))
    
    # Find important ESVs that are in the taxonomy table
    important_esvs <- top_features$Feature[top_features$Feature %in% rownames(tax_table_df)]
    
    if(length(important_esvs) > 0) {
        important_taxa <- tax_table_df[important_esvs, ] %>%
            rownames_to_column("ESV") %>%
            left_join(top_features, by = c("ESV" = "Feature")) %>%
            arrange(desc(mean_importance))
        
        cat("Top taxonomically-identified ESVs:\n")
        print(head(important_taxa, 15))
        
        # Save taxonomic results
        write_csv(important_taxa, "gt_prediction_important_taxa.csv")
    }
}