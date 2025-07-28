# this script runs a randomforest model on the microobial community data (vs ITS baiting data)

library(dplyr)

# load microbial metagenomic (?) data
ps_project11 <- readRDS("../data/ps_project11.rds")

View(sample_data(ps_project11))
ps_project11 <- prune_samples(!is.na(sample_data(ps_project11)$gt), ps_project11) # remove samples without GT data

ps_project11 # 17 samples
    # way too few for a random forest model
    # want 10-20 samples PER FEATURE

# trying a simpler model first- linear model

# what data should be in the model (features)
    # taxa abundance at a certain taxonomic level (e.g. phylum, class, order, family, genus)
    # environmental factors: 
    # response variable- gt

View(tax_table(ps_project11))

View(otu_table(ps_project11))

# Convert the OTU table to a matrix
otu_matrix <- as(otu_table(ps_project11), "matrix")

# Remove taxa with zero abundance across all samples
otu_matrix <- otu_matrix[rowSums(otu_matrix) > 0, ]

# Convert back to an OTU table
otu_table(ps_project11) <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Add GT data to the OTU table
otu_gt <- t(otu_matrix)

View(otu_gt)

View(sample_data(ps_project11))

otu_gt$gt <- sample_data(ps_project11)$gt

View(otu_gt)
otu_gt <- merge(sample_data(ps_project11)[, c("sample_ID", "gt")], otu_gt, by = "row.names")

for(i in refseq(ps_project11)) {
    print(length(refseq(ps_project11)[[i]]))
}





### old one I did for mass spec data


dummy_df <- readRDS("../data/metabolomics/prepared_ms_missing_NAs.RDS")
View(dummy_df)
nrow(dummy_df) # 1239 compounds

# data prep
dummy_df_t <- t(dummy_df) # transpose
dummy_df_t <- as.data.frame(dummy_df_t) %>%
    mutate(DummyITS = sample(x = c(0,1), size = nrow(dummy_df_t), replace = TRUE)) # generate a dummy ITS value
dummy_data <- dummy_df_t[-c(1,2),] # remove first two rows (mass and RT)
View(dummy_data$DummyITS) 

# quick random forest
library(randomForest)
training <- sample(nrow(dummy_data), nrow(dummy_data)/2) # Pick some training data

model <- randomForest(DummyITS ~ ., data=dummy_data[training,], mtry=ncol(dummy_data)/2) # using all cols (.) to predict DummyITS

# Examine the model
cor.test(predict(model, dummy_data[-training,]), dummy_data$DummyITS[-training]) # correlation between predicted and actual DummyITS values crosses zero
plot(model)
text(model)
model
summary(model)

# model is rubbish because dummy data is random

## inspect data

# mr vs abundance in sample 1
plot(x = as.numeric(unlist(dummy_df[,1])), y = as.numeric(unlist(dummy_df[,3])), 
    xlab = "compound Mr", ylab = "abundance in sample 1")

apply(dummy_df, 2, function(x) sum(is.Inf(x))) # count NAs per column

str(dummy_df[,2])
table(is.finite(dummy_df[,2]))

str(dummy_df)