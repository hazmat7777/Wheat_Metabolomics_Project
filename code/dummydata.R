# this script generates dummy data for testing purposes
library(dplyr)

# load mass spec data
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