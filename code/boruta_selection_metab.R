
#install.packages("ggsignif")
# library(ggsignif)


# fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")

# all_selected_features <- c()

# for(fold in 1:5) {
#   cat("\nFold", fold, ":\n")
  
#   # Create train/test split
#   test_indices <- folds[[fold]]
#   train_data <- fk_metabolom_gt_scaled[-test_indices, ]
#   test_data <- fk_metabolom_gt_scaled[test_indices, ]
  
#   tryCatch({
#     # Run Boruta on training data only
#     set.seed(123 + fold)
#     boruta_result <- Boruta(gt ~ ., data = train_data, 
#                            maxRuns = 100, 
#                            doTrace = 0)
    
#     important_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
#     boruta_features_per_fold[[fold]] <- important_features
    
#     cat("  Boruta selected", length(important_features), "features\n")
#     if(length(important_features) > 0) {
#       cat("  Features:", paste(important_features[1:min(3, length(important_features))], collapse = ", "))
#       if(length(important_features) > 3) cat(" ...")
#       cat("\n")
#     }
    
#     if(length(important_features) > 0) {
#       # Train RF with selected features
#       train_boruta <- train_data[, c(important_features, "gt")]
#       test_boruta <- test_data[, c(important_features, "gt")]
      
#       all_selected_features <- c(all_selected_features, important_features)

#       set.seed(123 + fold)
    
#     } else {
#       # No features selected

#       cat("  No features selected\n")
#     }
    
#   }, error = function(e) {
#     cat("  Error in fold", fold, ":", e$message, "\n")

#   })
# }


# all_selected_features <- unique(all_selected_features)

# all_selected_features




#### this is its own script

# More liberal statistical threshold
all_selected_features <- c()

for(run in 1:10) {
  set.seed(100 + run)  # Different seed each time
  boruta_result <- Boruta(gt ~ ., data = train_data, 
                         maxRuns = 100, 
                         doTrace = 0)
  
  selected_features <- getSelectedAttributes(boruta_result, withTentative = TRUE)
  all_selected_features <- c(all_selected_features, selected_features)
  
  cat("Run", run, ":", length(selected_features), "features\n")
}

all_selected_features <- unique(all_selected_features)
all_selected_features

formula <- as.formula(paste("gt ~", paste(all_selected_features, collapse = " + ")))
logistic_model <- glm(formula, data = train_data, family = binomial)

summary(logistic_model)

simple_model <- glm(gt ~ Haematommic.acid + Spectral.Match.to.Pinolenic.acid.from.NIST14 + X2..4.Hydroxyphenyl.propionic.Acid, 
                   data = train_data, family = binomial)
summary(simple_model)

library(car)
vif(simple_model) # all below 5, so no multicollinearity issues

## plots


library(broom)
library(ggplot2)


# Get model coefficients with confidence intervals
model_coefs <- tidy(simple_model, conf.int = TRUE)

# # Plot coefficients
# ggplot(model_coefs[-1,], aes(x = term, y = estimate)) +  # Remove intercept
#   geom_point() +
#   geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   coord_flip() +
#   labs(title = "Logistic Regression Coefficients",
#        x = "Predictors", y = "Log Odds") +
#   theme_minimal()

# Create prediction data
train_data$predicted_prob <- predict(simple_model, type = "response")

# Plot for first predictor
p1<- ggplot(train_data, aes(x = Haematommic.acid, y = predicted_prob)) +
  geom_point(aes(color = factor(gt)), alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title ="B",
       x = "Haematommic Acid Concentration (log peak area)", y = "Predicted Disease Probability",
       color = "Disease Status") +
  theme_minimal()+
  theme(legend_position = "none")

# Plot for second predictor  
p2 <- ggplot(train_data, aes(x = Spectral.Match.to.Pinolenic.acid.from.NIST14, y = predicted_prob)) +
  geom_point(aes(color = factor(gt)), alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "C",
       x = "Pinolenic Acid concentration (log peak area)", y = "Predicted Probability",
       color = "Disease Status") +
  theme_minimal()

# Plot for second predictor  
p3 <- ggplot(train_data, aes(x = X2..4.Hydroxyphenyl.propionic.Acid, y = predicted_prob)) +
  geom_point(aes(color = factor(gt)), alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "D",
       x = "Hydroxyphenyl propionic acid concentration (log peak area)", y = "Predicted Probability",
       color = "Disease Status") +
        scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()

p3

library(cowplot)
library(gridExtra)
combined_plot <- grid.arrange(p1,p2, ncol = 2)



# First create the protective score
protective_metabolites <- c("Haematommic.acid", 
                           "Spectral.Match.to.Pinolenic.acid.from.NIST14", 
                           "X2..4.Hydroxyphenyl.propionic.Acid")

protective_score <- rowSums(scale(fk_metabolom_gt_scaled[, protective_metabolites]))

# Add it to your dataframe
gt_metab <- fk_metabolom_gt_scaled
gt_metab$gt_numeric <- as.numeric(gt_metab$gt)  # factor levels become 1, 2
gt_metab$gt_numeric <- gt_metab$gt_numeric - 1  # now it's 0/1


gt_metab$protective_score <- protective_score

simple_protective_model <- glm(gt_numeric ~ protective_score, 
                               data = gt_metab, 
                               family = binomial)

gt_metab$predicted_prob <- predict(simple_protective_model, type = "response")

summary(simple_protective_model)
gt_metab$predicted_prob <- predict(simple_protective_model, type = "response")
library(ggplot2)

# ggplot(gt_metab, aes(x = protective_score)) +
#   # Actual data points, jittered
#   geom_jitter(aes(y = gt_numeric, color = factor(gt_numeric)), 
#               height = 0.05, alpha = 0.6) +
#   # Logistic regression curve
#   geom_line(aes(y = predicted_prob), color = "black", size = 1.2) +
#   labs(
#     title = "Logistic Regression: Protective Score vs gt",
#     x = "Protective Score (Standardized Sum)",
#     y = "Actual gt (Jittered) / Predicted Probability",
#     color = "Actual gt"
#   ) +
#   theme_minimal()

# ggplot(gt_metab, aes(x = protective_score)) +
#   geom_jitter(aes(y = gt_numeric, color = factor(gt_numeric)),
#               height = 0.01, alpha = 0.9, stroke = 0) +
#   geom_line(aes(y = predicted_prob), color = "black", size = 1.5) +
#   labs(
#     title = "A",
#     x = "Protective Score (Standardized Sum)",
#     y = "Probability",
#     color = "GT"
#   ) +
#   scale_y_continuous(labels = scales::percent_format()) +
#   theme_minimal()

pA <- ggplot(gt_metab, aes(x = factor(gt), y = protective_score, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7) +
  geom_signif(comparisons = list(c("GT_absent", "GT_present")), # adjust these to match your factor levels
              annotations = "***",  # change to *, **, *** based on your p-value
              y_position = max(gt_metab$protective_score) + 0.5) +
  labs(
    title = "A",
    x = "Disease status",
    y = "Protective Score",
    fill = "Outcome"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
pA

# Fix titles for clarity
p1 <- p1 + labs(title = "B") + theme(legend.position = "none")
p2 <- p2 + labs(title = "C") + theme(legend.position = "none")
p3 <- p3 + labs(title = "D") + theme(legend.position = "none")

library(cowplot)
library(cowplot)

# Make sure your plots have legends for extraction (e.g., p2)
# If you removed legends from p2, recreate it temporarily with legend:
p2_legend <- p2 + theme(legend.position = "bottom")

# Extract the legend
shared_legend <- get_legend(p2_legend)

# Remove legends from all plots
pA_nolegend <- pA + theme(legend.position = "none")
p1_nolegend <- p1 + theme(legend.position = "none")
p2_nolegend <- p2 + theme(legend.position = "none")
p3_nolegend <- p3 + theme(legend.position = "none")

# Combine the plots without legends
combined_plots <- plot_grid(
  pA_nolegend, p1_nolegend,
  p2_nolegend, p3_nolegend,
  labels = NULL,
  ncol = 2,
  align = "hv"
)

# Add the shared legend below
final_plot <- plot_grid(
  combined_plots,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.1)  # adjust legend height
)

print(final_plot)


### fixing data snooping
# 1. Calculate protective score using ONLY training data scaling
train_metabolites <- scale(train_data[, protective_metabolites])
train_protective_score <- rowSums(train_metabolites)

# 2. Apply the SAME scaling parameters to test data
test_means <- attr(train_metabolites, "scaled:center")
test_sds <- attr(train_metabolites, "scaled:scale")
test_metabolites_scaled <- scale(test_data[, protective_metabolites], 
                                center = test_means, scale = test_sds)
test_protective_score <- rowSums(test_metabolites_scaled)

# 3. Create training dataset with protective score
train_data_protected <- train_data
train_data_protected$protective_score <- train_protective_score

# 4. Fit the protective score model on training data only
protective_model <- glm(gt ~ protective_score, 
                       data = train_data_protected, 
                       family = binomial)

summary(protective_model)

# 5. Create test dataset with protective score
test_data_protected <- test_data
test_data_protected$protective_score <- test_protective_score

# 6. Make predictions on test data
test_predictions <- predict(protective_model, 
                           newdata = test_data_protected, 
                           type = "response")

# 7. Evaluate performance on test set
library(pROC)
test_auc <- auc(test_data_protected$gt, test_predictions)
cat("Test AUC for protective score model:", round(test_auc, 3), "\n")

# 8. Create clean visualization using training data only
train_data_protected$predicted_prob <- predict(protective_model, type = "response")

# Boxplot (Panel A) using training data
pA_clean <- ggplot(train_data_protected, aes(x = factor(gt), y = protective_score, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "A",
    x = "Disease Status",
    y = "Protective Score (Standardized)",
    fill = "Outcome"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Add significance test (this is now valid since it's training data only)
library(ggsignif)
pA_clean <- pA_clean + 
  geom_signif(comparisons = list(c("GT_absent", "GT_present")), 
              test = "t.test",
              y_position = max(train_data_protected$protective_score) + 0.5)

print(pA_clean)