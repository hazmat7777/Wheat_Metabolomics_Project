# this script runs a randomforest model on the microobial community data (vs ITS baiting data)
# there is no split of train/test, just OOB error- see explanation on rf_explainer

library(dplyr)

# load microbial metabarcoding data

ps_16S_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")

sample_names(ps_16S_highdiv_absolute)

ps_16S_rf <- prune_samples(!is.na(sample_data(ps_16S_highdiv_absolute)$gt), ps_16S_highdiv_absolute) # remove samples without GT data
ps_16S_rf <- prune_taxa(taxa_sums(ps_16S_rf) > 0, ps_16S_rf) #removes missing taxa

sample_names(ps_16S_rf)

#===================================================
# THE PLAN
#===================================================

# Prelim- do one by ESVs

# then do dimensionality reduction:
    # First 5 rfs- do a rf aggregated by each taxonomic level (phylum, class, order, genus)
    # 6th one- pcoa on jenson-shannon distance between samples (???)

# Load required libraries
library(phyloseq)
library(randomForest)
library(vegan)
library(ade4)
library(dplyr)

# ===== APPROACH 1-5: TAXONOMIC AGGREGATION =====

# Function to aggregate OTUs at different taxonomic ranks
# Function to aggregate OTUs at different taxonomic ranks
aggregate_taxa <- function(ps_obj, tax_rank) {
  # Aggregate taxa at specified rank
  ps_agg <- tax_glom(ps_obj, taxrank = tax_rank)
  
  # Get OTU table (samples as rows, taxa as columns)
  otu_table_agg <- as.data.frame(t(otu_table(ps_agg)))
  
  # Optional: Add taxonomic names as column names
  if (!is.null(tax_table(ps_agg))) {
    tax_names <- as.data.frame(tax_table(ps_agg))[, tax_rank]
    # Handle NAs and create unique names
    tax_names[is.na(tax_names)] <- paste0("Unknown_", tax_rank, "_", seq_len(sum(is.na(tax_names))))
    tax_names <- make.unique(as.character(tax_names))
    
    # Clean up names for R compatibility
    tax_names <- make.names(tax_names, unique = TRUE)
    
    colnames(otu_table_agg) <- tax_names
  } else {
    # If no taxonomy table, create generic names
    colnames(otu_table_agg) <- make.names(paste0("Taxa_", seq_len(ncol(otu_table_agg))), unique = TRUE)
  }
  
  return(otu_table_agg)
}

# Aggregate at different taxonomic ranks
taxonomic_ranks <- c("Genus", "Family", "Order", "Class", "Phylum")
aggregated_data <- list()

for (rank in taxonomic_ranks) {
  cat("Aggregating at", rank, "level...\n")
  aggregated_data[[rank]] <- aggregate_taxa(ps_16S_rf, rank)
  cat("Number of", rank, "features:", ncol(aggregated_data[[rank]]), "\n\n")
}
# ===== RANDOM FOREST MODELS =====

# Get sample data (make sure your response variable is in sample_data)
sample_df <- as.data.frame(sample_data(ps_16S_rf))

# Function to run random forest
run_random_forest <- function(feature_data, response_var, approach_name) {
  cat("Running Random Forest for", approach_name, "...\n")
  
  # Debug: Check dimensions
  cat("Feature data dimensions:", dim(feature_data), "\n")
  cat("Response variable length:", length(response_var), "\n")
  cat("Response variable type:", class(response_var), "\n")
  
  # Ensure sample names match between feature data and response variable
  feature_samples <- rownames(feature_data)
  response_samples <- names(response_var)
  
  # If response_var doesn't have names, try to get them from sample_data
  if (is.null(response_samples)) {
    response_samples <- sample_names(ps_16S_rf)
    names(response_var) <- response_samples
  }
  
  # Find common samples
  common_samples <- intersect(feature_samples, response_samples)
  cat("Common samples:", length(common_samples), "\n")
  
  if (length(common_samples) == 0) {
    stop("No matching samples between feature data and response variable")
  }
  
  # Subset both to common samples
  feature_data_matched <- feature_data[common_samples, , drop = FALSE]
  response_var_matched <- response_var[common_samples]
  
  # Clean column names to be R-compatible
  colnames(feature_data_matched) <- make.names(colnames(feature_data_matched), unique = TRUE)
  
  # Combine features with response variable
  rf_data <- cbind(feature_data_matched, response = response_var_matched)
  
  # Remove any columns with zero variance
  zero_var_cols <- which(apply(rf_data[, -ncol(rf_data)], 2, var, na.rm = TRUE) == 0)
  if (length(zero_var_cols) > 0) {
    rf_data <- rf_data[, -zero_var_cols]
    cat("Removed", length(zero_var_cols), "zero variance features\n")
  }
  
  # Remove rows with NA in response variable
  complete_rows <- complete.cases(rf_data$response)
  if (sum(complete_rows) < nrow(rf_data)) {
    cat("Removing", sum(!complete_rows), "samples with missing response values\n")
    rf_data <- rf_data[complete_rows, ]
  }
  
  # Ensure response is factor for classification
  rf_data$response <- as.factor(rf_data$response)
  cat("Response levels:", levels(rf_data$response), "\n")
  cat("Response distribution:", table(rf_data$response), "\n")
  
  # Run random forest classification
  rf_model <- randomForest(response ~ ., data = rf_data, 
                          importance = TRUE, ntree = 500)
  
  # Print results
  cat("Final dataset dimensions:", dim(rf_data), "\n")
  cat("OOB Error Rate:", round(rf_model$err.rate[nrow(rf_model$err.rate), "OOB"] * 100, 2), "%\n")
  cat("Number of features used:", ncol(rf_data) - 1, "\n")
  
  # Print confusion matrix
  cat("Confusion Matrix:\n")
  print(rf_model$confusion)
  cat("\n")
  
  return(rf_model)
}

# Extract response variable from sample_data
sample_df <- as.data.frame(sample_data(ps_16S_rf))
gt <- sample_df$gt  # Extract the gt column
names(gt) <- rownames(sample_df)  # Assign sample names
gt <- as.factor(gt) # Convert to factor for classification


# Approaches 1-5: Taxonomic aggregation
rf_models <- list()
for (rank in taxonomic_ranks) {
  rf_models[[rank]] <- run_random_forest(
    aggregated_data[[rank]], 
    gt,  
    paste("Taxonomic aggregation -", rank)
  )
}

saveRDS(rf_models, file = "../results/microb_rf_oob_error")
rf_models

# ===== COMPARE MODELS =====

# Extract OOB error rates
oob_errors <- sapply(rf_models, function(x) tail(x$err.rate[, "OOB"], 1))
model_comparison <- data.frame(
  Approach = names(oob_errors),
  OOB_Error = round(oob_errors * 100, 2),
  stringsAsFactors = FALSE
)

print("Model Comparison:")
print(model_comparison)

# Plot comparison
barplot(model_comparison$OOB_Error, 
        names.arg = model_comparison$Approach,
        main = "Random Forest OOB Error by Approach",
        ylab = "OOB Error Rate (%)",
        las = 2)

# ===== FEATURE IMPORTANCE =====

# Function to get top important features
get_top_features <- function(rf_model, n_features = 10) {
  importance_scores <- importance(rf_model)
  
  # For classification, use MeanDecreaseAccuracy
  imp_col <- "MeanDecreaseAccuracy"
  
  top_features <- head(
    importance_scores[order(importance_scores[, imp_col], decreasing = TRUE), ],
    n_features
  )
  
  return(top_features)
}

# PLOT top 5 features for each model
cat("\nTop 10 important features for each approach:\n")
for (model_name in names(rf_models)) {
  cat("\n", model_name, ":\n")
  top_feat <- get_top_features(rf_models[[model_name]])
  print(top_feat)
}

top_order_features <- row.names(get_top_features(rf_models[[3]]))
top_order_features
order_glom<- (as.data.frame(aggregated_data[3]))
order_glom$gt <- as.factor(as.data.frame(sample_data(ps_16S_rf))$gt) # check order didnt change

View(order_glom)
library(ggplot2)

for (microbe in top_order_features) { # dont forget to set counter
    counter <- counter + 1
    microbe <- paste0("Order.",microbe)
    plot <- ggplot(order_glom, aes(y = gt, x = .data[[microbe]])) +
        geom_boxplot() +
        labs(x = "GT Status", y = paste("Distribution of", microbe)) +
        theme_bw()
    
    ggsave(filename = paste0("../results/plots/microb_", counter, "_vs_gt.png"), plot = plot, width = 8, height = 6)
    if(counter == 5){
        counter <- 0
    }
}
