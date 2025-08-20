# Save this as a new clean file

library(tidyverse)
library(phyloseq)
library(ranger)
library(caret)

ps_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")
ps_class <- tax_glom(ps_highdiv_absolute, "Class", NArm = TRUE)
taxa_names(ps_class) <- tax_table(ps_class)[, "Class"]

# Prepare data
sample_df <- data.frame(sample_data(ps_class))
otu_df <- data.frame(t(as.matrix(otu_table(ps_class))))
rf_data <- cbind(sample_df, otu_df)
rf_data <- rf_data[!is.na(rf_data$gt), ]
rf_data$gt <- as.factor(rf_data$gt)

# Clean columns
rf_data$sample_name <- NULL
rf_data$og_sample_names <- NULL
if("tillage" %in% names(rf_data)) rf_data$tillage <- as.factor(rf_data$tillage)

# Cross-validation
set.seed(123)
folds <- createFolds(rf_data$gt, k = 5, returnTrain = FALSE)
accuracies <- c()

for(i in 1:5) {
    train_data <- rf_data[-folds[[i]], ]
    test_data <- rf_data[folds[[i]], ]
    
    rf_model <- ranger(gt ~ ., data = train_data, num.trees = 100)
    predictions <- predict(rf_model, test_data)
    
    accuracy <- mean(predictions$predictions == test_data$gt)
    accuracies[i] <- accuracy
    cat("Fold", i, "accuracy:", round(accuracy, 3), "\\n")
}

cat("Mean accuracy:", round(mean(accuracies), 3), "\\n")
