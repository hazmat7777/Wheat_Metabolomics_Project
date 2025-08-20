library(phyloseq)
ps <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")
ps

ps_genus <- tax_glom(ps, "Genus", NArm = TRUE)
taxa_names(ps_genus) <- make.unique(as.character(tax_table(ps_genus)[, "Genus"])) 
sample_names(ps_genus) <- as.data.frame(sample_data(ps_genus))$sample_name
sample_names(ps_genus)
saveRDS(ps_genus, file = "../data/metabarcoding/ps_16S_highdiv_genus_relative.rds")

cat("ESVs after tax_glom:", ntaxa(ps_genus), "\n")


otu_tab <- otu_table(ps_genus)

prevalence_df <- data.frame(
  Prevalence = taxa_sums(otu_tab > 0),  # Number of samples where genus is present
  TotalAbundance = taxa_sums(otu_tab),  # Total abundance across all samples
  Genus = taxa_names(ps_genus)
)

# Calculate prevalence as proportion of samples
prevalence_df$PrevalenceProp <- prevalence_df$Prevalence / nsamples(ps_genus)

# Summary statistics
cat("Total samples:", nsamples(ps_genus), "\n")
cat("Total genera:", ntaxa(ps_genus), "\n")
cat("Prevalence summary:\n")
summary(prevalence_df$PrevalenceProp)

# How many genera at different prevalence thresholds?
cat("\nGenera present in at least:\n")
cat("1% of samples:", sum(prevalence_df$PrevalenceProp >= 0.01), "\n")
cat("5% of samples:", sum(prevalence_df$PrevalenceProp >= 0.05), "\n")
cat("10% of samples:", sum(prevalence_df$PrevalenceProp >= 0.10), "\n")
cat("20% of samples:", sum(prevalence_df$PrevalenceProp >= 0.20), "\n")
cat("50% of samples:", sum(prevalence_df$PrevalenceProp >= 0.50), "\n")








# running boruta algo

# Filter to genera present in at least 10% of samples
ps_genus_filtered <- prune_taxa(prevalence_df$PrevalenceProp >= 0.10, ps_genus)
cat("After 10% prevalence filtering:", ntaxa(ps_genus_filtered), "genera remain\n")

# Convert to data frame for Boruta
genus_data <- as.data.frame(t(otu_table(ps_genus_filtered)))
View(genus_data)
# Get your outcome variable (assuming it's in sample_data)
outcome_var <- as.factor(sample_data(ps_genus_filtered)$gt  )

# Combine for Boruta
boruta_data <- cbind(gt = outcome_var, genus_data)

# Run Boruta (start with fewer runs to test)
library(Boruta)
set.seed(123)
boruta_result <- Boruta(gt ~ ., data = boruta_data, 
                       maxRuns = 50,  # start smaller
                       doTrace = 2)   # show progress

print(boruta_result)
selected_genera <- getSelectedAttributes(boruta_result, withTentative = TRUE)
length(selected_genera)


# First, do train/test split
set.seed(123)  # Use same seed as metabolomics for consistency
train_indices <- sample(1:nrow(boruta_data), size = 0.7 * nrow(boruta_data))
train_data_micro <- boruta_data[train_indices, ]
test_data_micro <- boruta_data[-train_indices, ]

# Check dimensions
cat("Training set:", nrow(train_data_micro), "samples\n")
cat("Test set:", nrow(test_data_micro), "samples\n")

# Now run Boruta on training data only
boruta_result <- Boruta(gt ~ ., data = train_data_micro, 
                       maxRuns = 50,
                       doTrace = 2)


print(boruta_result)
# Evaluate selected features on test set
selected_genera <- getSelectedAttributes(boruta_result, withTentative = TRUE)

View(tax_table(ps_genus_filtered)[selected_genera, ])

# run glm
if(length(selected_genera) > 0) {
  # Wrap genus names in backticks for the formula
  genus_names_clean <- paste0("`", selected_genera, "`")
  formula <- as.formula(paste("gt ~", paste(genus_names_clean, collapse = " + ")))
  
  print(formula)  # Check what the formula looks like
  
  micro_model <- glm(formula, data = train_data_micro, family = binomial)
  
  # Test on held-out data
  test_predictions <- predict(micro_model, test_data_micro, type = "response")
  
  # Quick performance check
  library(pROC)
  test_auc <- auc(test_data_micro$gt, test_predictions)
  cat("Test AUC:", round(test_auc, 3), "\n")
  
  summary(micro_model)
}
# Simple model with just Rhodoplanes and Geomonas

simple_micro_model <- glm(gt ~ Rhodoplanes + Geomonas, data = train_data_micro, family = binomial)

# Test on held-out data
simple_test_predictions <- predict(simple_micro_model, test_data_micro, type = "response")

# Performance check
simple_test_auc <- auc(test_data_micro$gt, simple_test_predictions)
cat("Simple model Test AUC:", round(simple_test_auc, 3), "\n")

summary(simple_micro_model)

# Compare to the full tentative model
cat("Full tentative model AUC:", 0.495, "\n")
cat("Simple 2-genera model AUC:", round(simple_test_auc, 3), "\n")

# Check the actual data for these two genera
cat("Rhodoplanes summary:\n")
summary(train_data_micro$Rhodoplanes)
cat("Number of non-zero samples:", sum(train_data_micro$Rhodoplanes > 0), "\n")

cat("\nGeomonas summary:\n") 
summary(train_data_micro$Geomonas)
cat("Number of non-zero samples:", sum(train_data_micro$Geomonas > 0), "\n")

# Cross-tabulation with outcome
cat("\nRhodoplanes by outcome:\n")
table(train_data_micro$Rhodoplanes > 0, train_data_micro$gt)

cat("\nGeomonas by outcome:\n")
table(train_data_micro$Geomonas > 0, train_data_micro$gt)