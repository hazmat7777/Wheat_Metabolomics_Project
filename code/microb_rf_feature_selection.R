# Save this as a new clean file

library(tidyverse)
library(phyloseq)
library(ranger)
library(caret)

ps_highdiv_relative <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")

# get a list of the genera implicated in GT suppression
# Define genera of interest
genera_implicated <- c("Pseudomonas", "Nocardioides", "Bacillus")

# Get the tax_table as a data frame
tax_df <- as.data.frame(tax_table(ps_highdiv_relative))

# Get the taxa names (OTUs/ASVs) where Genus is in your list
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]

# Now prune the phyloseq object
ps_fs <- prune_taxa(taxa_to_keep, ps_highdiv_relative)
ps_fs

# Prepare data
sample_df <- data.frame(sample_data(ps_fs))
otu_df <- data.frame(t(as.matrix(otu_table(ps_fs))))
rf_data <- cbind(sample_df, otu_df)
rf_data <- rf_data[!is.na(rf_data$gt), ]
rf_data$gt <- as.factor(rf_data$gt)

# Clean columns
rf_data$sample_name <- NULL
rf_data$og_sample_names <- NULL
if("tillage" %in% names(rf_data)) rf_data$tillage <- as.factor(rf_data$tillage)

# Cross-validation
set.seed(1250)
folds <- createFolds(rf_data$gt, k = 5, returnTrain = FALSE)
accuracies <- c()

for(i in 1:5) {
    train_data <- rf_data[-folds[[i]], ]
    test_data <- rf_data[folds[[i]], ]
    
    rf_model <- ranger(gt ~ ., data = train_data, num.trees = 3000, mtry = 5)
    predictions <- predict(rf_model, test_data)
    
    accuracy <- mean(predictions$predictions == test_data$gt)
    accuracies[i] <- accuracy
    cat("Fold", i, "accuracy:", round(accuracy, 3), "\\n")
}

cat("Mean accuracy:", round(mean(accuracies), 3), "\\n")
