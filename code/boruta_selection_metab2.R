# Metabolomics Analysis: Feature Selection and Biomarker Development
# Clean script without data snooping

library(Boruta)
library(car)
library(broom)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(pROC)

# Load data (assuming train/test split already done)
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
set.seed(123)
train_indices <- sample(1:nrow(fk_metabolom_gt_scaled), size = 0.7 * nrow(fk_metabolom_gt_scaled))
train_data <- fk_metabolom_gt_scaled[train_indices, ]
test_data <- fk_metabolom_gt_scaled[-train_indices, ]

# Feature selection with Boruta (multiple runs for robustness)
all_selected_features <- c()
for(run in 1:10) {
  set.seed(100 + run)
  boruta_result <- Boruta(gt ~ ., data = train_data, 
                         maxRuns = 100, 
                         doTrace = 0)
  
  selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
  all_selected_features <- c(all_selected_features, selected_features)
  
  cat("Run", run, ":", length(selected_features), "features\n")
}

all_selected_features <- unique(all_selected_features)
cat("Total unique features selected:", length(all_selected_features), "\n")

# Build logistic regression model with selected features
formula <- as.formula(paste("gt ~", paste(all_selected_features, collapse = " + ")))
logistic_model <- glm(formula, data = train_data, family = binomial)
summary(logistic_model)

# Simple model with the 3 significant metabolites from Boruta results
# Ribonolactone, Haematomic.acid, Glycerol.Tributyrate (p < 0.05)

library(car)

# Build simple model with 3 significant metabolites
simple_model_3sig <- glm(gt ~ Ribonolactone + Haematommic.acid + Glycerol.Tributyrate, 
                        data = train_data, family = binomial)

# Model summary
summary(simple_model_3sig)

# Check VIF
cat("\nVariance Inflation Factors:\n")
vif_values <- vif(simple_model_3sig)
print(vif_values) # all below 5

# Test set performance
test_predictions_simple <- predict(simple_model_3sig, newdata = test_data, type = "response")

library(pROC)
test_auc_simple <- auc(test_data$gt, test_predictions_simple)
cat("\nTest AUC for simple 3-metabolite model:", round(test_auc_simple, 3), "\n")

# Simple 3-metabolite model (hand-selected based on Boruta results and biological relevance)
simple_model <- glm(gt ~ Haematommic.acid + Spectral.Match.to.Pinolenic.acid.from.NIST14 + X2..4.Hydroxyphenyl.propionic.Acid, 
                   data = train_data, family = binomial)
summary(simple_model)

# Check for multicollinearity
vif(simple_model)

# Create individual metabolite plots
train_data$predicted_prob <- predict(simple_model, type = "response")

# Plot B: Haematomic Acid
p1 <- ggplot(train_data, aes(x = Haematommic.acid, y = predicted_prob)) +
  geom_point(aes(color = factor(gt)), alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "B",
       x = "Haematommic Acid Concentration (log peak area)", 
       y = "Predicted Disease Probability",
       color = "Disease Status") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot C: Pinolenic Acid
p2 <- ggplot(train_data, aes(x = Spectral.Match.to.Pinolenic.acid.from.NIST14, y = predicted_prob)) +
  geom_point(aes(color = factor(gt)), alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "C",
       x = "Pinolenic Acid Concentration (log peak area)", 
       y = "Predicted Probability",
       color = "Disease Status") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot D: Hydroxyphenyl Propionic Acid
p3 <- ggplot(train_data, aes(x = X2..4.Hydroxyphenyl.propionic.Acid, y = predicted_prob)) +
  geom_point(aes(color = factor(gt)), alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "D",
       x = "Hydroxyphenyl Propionic Acid Concentration (log peak area)", 
       y = "Predicted Probability",
       color = "Disease Status") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(legend.position = "none")

# Protective Score Analysis (proper train/test split to avoid data snooping)
protective_metabolites <- c("Haematommic.acid", 
                           "Spectral.Match.to.Pinolenic.acid.from.NIST14", 
                           "X2..4.Hydroxyphenyl.propionic.Acid")

# Scale metabolites using ONLY training data
train_metabolites <- scale(train_data[, protective_metabolites])
train_protective_score <- rowSums(train_metabolites)

# Apply same scaling parameters to test data
test_means <- attr(train_metabolites, "scaled:center")
test_sds <- attr(train_metabolites, "scaled:scale")
test_metabolites_scaled <- scale(test_data[, protective_metabolites], 
                                center = test_means, scale = test_sds)
test_protective_score <- rowSums(test_metabolites_scaled)

# Create datasets with protective scores
train_data_protected <- train_data
train_data_protected$protective_score <- train_protective_score

test_data_protected <- test_data
test_data_protected$protective_score <- test_protective_score

# Fit protective score model on training data
protective_model <- glm(gt ~ protective_score, 
                       data = train_data_protected, 
                       family = binomial)
summary(protective_model)

# Evaluate on test data
test_predictions <- predict(protective_model, 
                           newdata = test_data_protected, 
                           type = "response")

test_auc <- auc(test_data_protected$gt, test_predictions)
cat("Test AUC for protective score model:", round(test_auc, 3), "\n")

# Plot A: Protective Score Boxplot
pA_clean <- ggplot(train_data_protected, aes(x = factor(gt), y = protective_score, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7) +
  geom_signif(comparisons = list(c("GT_absent", "GT_present")), 
              test = "t.test",
              y_position = max(train_data_protected$protective_score) + 0.5) +
  labs(
    title = "A",
    x = "Disease Status",
    y = "Protective Score (Standardized)",
    fill = "Outcome"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Create combined figure
# Extract legend from one plot
p2_legend <- p2 + theme(legend.position = "bottom")
shared_legend <- get_legend(p2_legend)

# Combine plots without legends
combined_plots <- plot_grid(
  pA_clean, p1,
  p2, p3,
  labels = NULL,
  ncol = 2,
  align = "hv"
)

# Add shared legend
final_plot <- plot_grid(
  combined_plots,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

print(final_plot)

# Test individual metabolite performance on test data
cat("\nIndividual metabolite test performance:\n")
for(metabolite in protective_metabolites) {
  individual_model <- glm(as.formula(paste("gt ~", metabolite)), 
                         data = train_data, family = binomial)
  individual_pred <- predict(individual_model, test_data, type = "response")
  individual_auc <- auc(test_data$gt, individual_pred)
  cat(metabolite, "Test AUC:", round(individual_auc, 3), "\n")
}