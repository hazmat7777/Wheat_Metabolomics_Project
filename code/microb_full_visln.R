# ====================================
# RANDOM FOREST BEST MODEL VISUALIZATION
# ====================================
# Run this script AFTER the feature selection comparison script

library(ranger)
library(caret)
library(ggplot2)
library(gridExtra)
library(pROC)
library(viridis)
library(phyloseq)
library(dplyr)
library(tidyr)
library(grid)  # Add this for textGrob

# ====================================
# TRAIN FINAL BEST MODEL
# ====================================
cat("\n=== TRAINING BEST MODEL FOR VISUALIZATION ===\n")

# Determine which dataset performed better from results_df
best_dataset_name <- results_df$Dataset[which.max(results_df$Mean_Accuracy)]

# Select the appropriate dataset
if(best_dataset_name == "Full_Dataset") {
  rf_data_best <- rf_data_full
  cat("Using FULL dataset for best model\n")
} else if(best_dataset_name == "Feature_Selected_2Genera") {
  rf_data_best <- rf_data_fs
  cat("Using FEATURE-SELECTED dataset for best model\n")
} else {
  # For any other feature selection method you might add
  stop("Dataset not recognized. Please update the visualization script.")
}

# Get best parameters
best_idx <- which.max(results_df$Mean_Accuracy)
best_mtry <- results_df$mtry[best_idx]
best_nodesize <- results_df$min.node.size[best_idx]
best_splitrule <- as.character(results_df$splitrule[best_idx])

cat("Best parameters: mtry =", best_mtry, ", min.node.size =", best_nodesize, 
    ", splitrule =", best_splitrule, "\n\n")

# Train final model on full dataset for visualization
set.seed(123)
final_model <- ranger(
  gt ~ ., 
  data = rf_data_best,
  num.trees = 2000,
  mtry = best_mtry,
  min.node.size = best_nodesize,
  splitrule = best_splitrule,
  importance = 'impurity',
  probability = TRUE  # Important for ROC curve
)

# Get predictions
predictions_prob <- predict(final_model, rf_data_best)$predictions
predictions_class <- ifelse(predictions_prob[,2] > 0.5, "GT_present", "GT_absent")
predictions_class <- factor(predictions_class, levels = levels(rf_data_best$gt))

# ====================================
# A) VARIABLE IMPORTANCE PLOT
# ====================================
cat("Creating variable importance plot...\n")

# Get top 8 most important variables
importance_df <- data.frame(
  Variable = names(final_model$variable.importance),
  Importance = final_model$variable.importance
) %>%
  arrange(desc(Importance)) %>%
  head(8)

# Format variable names for display
importance_df$Variable_label <- gsub("^X", "", importance_df$Variable)
# Truncate long names for better display
importance_df$Variable_label <- substr(importance_df$Variable_label, 1, 30)

# Store top 2 features for scatter plot
top_2_features <- importance_df$Variable[1:2]

p_importance <- ggplot(importance_df, aes(x = reorder(Variable_label, Importance), 
                                          y = Importance)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = round(Importance, 1)), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "A) Variable Importance",
       x = "Features",
       y = "Mean Decrease Gini") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12))

# ====================================
# B) TOP 2 FEATURES SCATTER PLOT
# ====================================
cat("Creating top 2 features scatter plot...\n")

# Create scatter plot data using top 2 most important features
scatter_data <- data.frame(
  Feature1 = rf_data_best[[top_2_features[1]]],
  Feature2 = rf_data_best[[top_2_features[2]]],
  Actual = rf_data_best$gt,
  Prob_GT_present = predictions_prob[,2]
)

# Clean feature names for axis labels
feature1_label <- gsub("^X", "", top_2_features[1])
feature1_label <- substr(feature1_label, 1, 40)
feature2_label <- gsub("^X", "", top_2_features[2])
feature2_label <- substr(feature2_label, 1, 40)

p_scatter <- ggplot(scatter_data, aes(x = Feature1, y = Feature2)) +
  geom_point(aes(color = Prob_GT_present, shape = Actual), 
             size = 3, alpha = 0.7) +
  scale_color_viridis(name = "P(GT+)", 
                      limits = c(0, 1),
                      option = "viridis") +
  scale_shape_manual(values = c(16, 17),
                     labels = c("GT absent", "GT present")) +
  labs(title = "B) Top 2 Predictive Features",
       x = feature1_label,
       y = feature2_label,
       shape = "Actual") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

# ====================================
# C) PREDICTION PROBABILITY DISTRIBUTION
# ====================================
cat("Creating probability distribution plot...\n")

prob_data <- data.frame(
  Probability = predictions_prob[,2],
  Actual = rf_data_best$gt
)

p_density <- ggplot(prob_data, aes(x = Probability, fill = Actual)) +
  geom_density(alpha = 0.6, adjust = 1.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("GT_absent" = "steelblue", 
                               "GT_present" = "coral")) +
  labs(title = "C) Prediction Probability Distribution",
       x = "P(GT present)",
       y = "Density",
       fill = "Actual") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom") +
  xlim(0, 1)

# ====================================
# D) ROC CURVE
# ====================================
cat("Creating ROC curve...\n")

# Calculate ROC
roc_obj <- roc(rf_data_best$gt, predictions_prob[,2], 
               levels = c("GT_absent", "GT_present"),
               direction = "<")
auc_value <- round(auc(roc_obj), 3)

# Create ROC data frame
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "red", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray50", size = 0.8) +
  annotate("text", x = 0.7, y = 0.3, 
           label = paste0("AUC = ", auc_value),
           size = 5, fontface = "bold") +
  labs(title = paste0("D) ROC Curve (AUC = ", auc_value, ")"),
       x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12)) +
  coord_equal() +
  xlim(0, 1) + ylim(0, 1)

# ====================================
# COMBINE ALL PLOTS
# ====================================
cat("Combining all plots...\n")

# For VSCode: explicitly open a PNG device
png("rf_best_model_visualization.png", width = 14, height = 10, units = "in", res = 300)

# Arrange plots in 2x2 grid
grid.arrange(
  p_importance, p_scatter,
  p_density, p_roc,
  ncol = 2, nrow = 2,
  top = textGrob(paste("Random Forest Model Performance -", best_dataset_name), 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Close the device
dev.off()

cat("\nVisualization saved as 'rf_best_model_visualization.png'\n")

# ====================================
# ADDITIONAL PERFORMANCE METRICS
# ====================================
cat("\n=== DETAILED PERFORMANCE METRICS ===\n")

# Confusion matrix
cm <- confusionMatrix(predictions_class, rf_data_best$gt)
cat("\nConfusion Matrix:\n")
print(cm$table)

cat("\nPerformance Metrics:\n")
cat("Accuracy:", round(cm$overall['Accuracy'], 3), "\n")
cat("Sensitivity:", round(cm$byClass['Sensitivity'], 3), "\n")
cat("Specificity:", round(cm$byClass['Specificity'], 3), "\n")
cat("Balanced Accuracy:", round(cm$byClass['Balanced Accuracy'], 3), "\n")
cat("AUC:", auc_value, "\n")

# Feature importance summary
cat("\n=== TOP 10 MOST IMPORTANT FEATURES ===\n")
top_10_importance <- data.frame(
  Rank = 1:min(10, length(importance_df$Variable_label)),
  Feature = head(importance_df$Variable_label, 10),
  Importance = round(head(importance_df$Importance, 10), 2)
)
print(top_10_importance, row.names = FALSE)

# ====================================
# SAVE INDIVIDUAL PLOTS (OPTIONAL)
# ====================================
if(TRUE) {  # Set to FALSE to skip individual saves
  cat("\nSaving individual plots...\n")
  
  # Use png() device for VSCode compatibility
  png("rf_variable_importance.png", width = 8, height = 6, units = "in", res = 300)
  print(p_importance)
  dev.off()
  
  png("rf_top_features_scatter.png", width = 8, height = 6, units = "in", res = 300)
  print(p_scatter)
  dev.off()
  
  png("rf_probability_distribution.png", width = 8, height = 6, units = "in", res = 300)
  print(p_density)
  dev.off()
  
  png("rf_roc_curve.png", width = 6, height = 6, units = "in", res = 300)
  print(p_roc)
  dev.off()
  
  cat("Individual plots saved.\n")
}

cat("\n=== VISUALIZATION COMPLETE ===\n")