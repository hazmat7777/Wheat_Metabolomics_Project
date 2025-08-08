# load libraries


# Load necessary libraries
#install.packages("glmnet")
library(glmnet)
library(caret)

# Data preparation (as you've already done)
fk_metabolom_gt <- readRDS("../data/fk_metabolomics_gt.RDS")

## Prepare data frame
fk_metabolom_gt_t <- t(fk_metabolom_gt) # Transpose to make metabolites as columns
colnames(fk_metabolom_gt_t) <- fk_metabolom_gt_t[1, ] # Set the first row as column names
fk_metabolom_gt_t <- fk_metabolom_gt_t[-1, ] # Remove the first row (now column names)

fk_metabolom_gt_t <- as.data.frame(fk_metabolom_gt_t) # Convert to data frame
colnames(fk_metabolom_gt_t)[ncol(fk_metabolom_gt_t)] <- "gt" # Rename the last column to 'gt'
fk_metabolom_gt_t$gt <- trimws(as.character(fk_metabolom_gt_t$gt)) # Remove any extra whitespace from gt col 
fk_metabolom_gt_t$gt <- as.factor(fk_metabolom_gt_t$gt) # Convert 'gt' to a factor
levels(fk_metabolom_gt_t$gt) <- c("GT_absent", "GT_present") # Set factor levels

# Convert all metabolite columns (except gt) to numeric
fk_metabolom_gt_t[ , !colnames(fk_metabolom_gt_t) %in% "gt"] <- 
  lapply(fk_metabolom_gt_t[ , !colnames(fk_metabolom_gt_t) %in% "gt"], as.numeric)

# Re-coerce to ensure correct structure
fk_metabolom_gt_t <- as.data.frame(fk_metabolom_gt_t)

#scale data
# predictors <- fk_metabolom_gt_t[, !names(fk_metabolom_gt_t) %in% "gt"]
# predictors[] <- lapply(predictors, function(x) as.numeric(as.character(x)))
# predictors_scaled <- scale(predictors)
# fk_scaled <- data.frame(predictors_scaled, gt = fk_metabolom_gt_t$gt)
# View(fk_scaled)

#=========================================================
# LASSO REGRESSION
#=========================================================


# 1. Split data into train and test sets (70% train, 30% test)
set.seed(123)
train_idx <- createDataPartition(fk_metabolom_gt_t$gt, p = 0.7, list = FALSE)
train_data <- fk_metabolom_gt_t[train_idx, ]
test_data <- fk_metabolom_gt_t[-train_idx, ]

# 2. Prepare predictors (X) and target (y)
x_train <- as.matrix(train_data[, -which(names(train_data) == "gt")])  # predictors (excluding 'gt')
y_train <- train_data$gt  # target variable
x_test <- as.matrix(test_data[, -which(names(test_data) == "gt")])    # predictors for test
y_test <- test_data$gt    # actual target for test

# 3. Standardize the data (important for Lasso)
# glmnet standardizes by default, but let's do it manually for clarity:
# x_train_scaled <- scale(x_train)
# x_test_scaled <- scale(x_test)

# 4. Fit the Lasso model using cross-validation to find the optimal lambda
lasso_model <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial", type.measure = "class")

# 5. Plot the cross-validation results to see the best lambda
plot(lasso_model)

# 6. Best lambda
best_lambda <- lasso_model$lambda.min
cat("Best lambda:", best_lambda, "\n")

# 7. Make predictions on the test data using the best lambda
predictions <- predict(lasso_model, s = "lambda.min", newx = x_test, type = "response")

# 8. Convert predictions to class labels (GT_absent, GT_present)
predicted_class <- ifelse(predictions > 0.5, "GT_present", "GT_absent")
predicted_class <- factor(predicted_class, levels = c("GT_absent", "GT_present"))

# 9. Evaluate the model performance (confusion matrix, accuracy, etc.)
conf_matrix <- confusionMatrix(predicted_class, y_test)
print(conf_matrix)

# Optional: You can check the coefficients of the selected features
cat("Non-zero coefficients:\n")
print(coef(lasso_model, s = "lambda.min"))


#=========================================================
# TOP PREDICTOR METABOLITE VS GT
#=========================================================

rf_metab_tuned <- readRDS(file = "../results/rf_metab_tuned.rds")

# Variable importance
importance_scores <- importance(rf_final)
varImpPlot(rf_final, n.var = 20)

# rf told us that the top 5 metabolites were:
top_metabolites <- rownames(importance_scores)[order(importance_scores[,1], decreasing = TRUE)][1:5]
top_metabolites

str(fk_metabolom_gt_t)[, ncol(fk_metabolom_gt_t)-1]

class(fk_metabolom_gt_t)


top_metabolites %in% colnames(fk_metabolom_gt_t)

top_metabolites_clean <- trimws(tolower(top_metabolites))
colnames_clean <- trimws(tolower(colnames(fk_metabolom_gt_t)))

# Check for matches
top_metabolites_clean %in% colnames_clean

matches <- sapply(top_metabolites, function(metab) {
  grep(metab, colnames(fk_metabolom_gt_t), value = TRUE)
})

print(matches)

counter <- 1

# plot each of the top 5 metabolites against GT
library(ggplot2)


column_name <- top_metabolites[1]

top_metabolites

# convert the fk colnames to names
colnames(fk_metabolom_gt_t) <- make.names(colnames(fk_metabolom_gt_t))

# test
plot_df <- fk_metabolom_gt_t %>%
    select(gt, !!top_metabolites[1])

View(plot_df)
head(sort(colnames(fk_metabolom_gt_t)), n =150)

#####

for (metab in top_metabolites) { # dont forget to set counter
    if(metab == top_metabolites[1]){
        counter <- 0
    }
    counter <- counter + 1

    plot <- ggplot(fk_metabolom_gt_t, aes(y = gt, x = .data[[metab]])) +
        geom_boxplot() +
        labs(x = "GT Status", y = paste("Distribution of", actual_colname)) +
        theme_bw()
    
    ggsave(filename = paste0("../results/plots/metab_", counter, "_vs_gt.png"), plot = plot, width = 8, height = 6)
    if(counter == 5){
        counter <- 0
    }
}


#===============================================
# simple lm + boxplot of diversity vs gt?
#===============================================
