# Save this as a new clean file

library(tidyverse)
library(phyloseq)
library(ranger)
library(caret)

ps_highdiv_relative <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")

ps_genus <- tax_glom(ps_highdiv_relative, "Genus", NArm = TRUE)
taxa_names(ps_genus) <- make.unique(as.character(tax_table(ps_genus)[, "Genus"])) 

cat("ESVs after tax_glom:", ntaxa(ps_genus), "\n")

# Check the taxonomy table after aggregation
tax_after <- data.frame(tax_table(ps_genus))
genus_names_after <- as.character(tax_after$Genus)
cat("Unique genera after tax_glom:", length(unique(genus_names_after)), "\n")
cat("Duplicated genera:", sum(duplicated(genus_names_after)), "\n")

# Show some examples of what's duplicated

dups <- genus_names_after[duplicated(genus_names_after)]
cat("Example duplicated genera:", head(unique(dups), 5), "\n") # only incertae sedis is duplicated.



# Prepare data
sample_df <- data.frame(sample_data(ps_genus))
otu_df <- data.frame(t(as.matrix(otu_table(ps_genus))))
rf_data <- cbind(sample_df, otu_df)
rf_data <- rf_data[!is.na(rf_data$gt), ]
rf_data$gt <- as.factor(rf_data$gt)

# Clean columns
rf_data$sample_name <- NULL
rf_data$og_sample_names <- NULL
if("tillage" %in% names(rf_data)) rf_data$tillage <- as.factor(rf_data$tillage)

# Cross-validation
set.seed(12550)
folds <- createFolds(rf_data$gt, k = 5, returnTrain = FALSE)
accuracies <- c()

for(i in 1:5) {
    train_data <- rf_data[-folds[[i]], ]
    test_data <- rf_data[folds[[i]], ]
    
    rf_model <- ranger(gt ~ ., data = train_data, num.trees = 5000, mtry = 3)
    predictions <- predict(rf_model, test_data)
    
    accuracy <- mean(predictions$predictions == test_data$gt)
    accuracies[i] <- accuracy
    cat("Fold", i, "accuracy:", round(accuracy, 3), "\\n")
}

cat("Mean accuracy:", round(mean(accuracies), 3), "\\n")
