# Random Forest Visualization Functions for PCA Metabolomics Data
# to be run AFTER sourcing metab_pca_then_rf.R

source("./metab_pca_then_rf2.R")  # Ensure PCA results are available

library(ggplot2)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(corrplot)

###########################################################################
# 1. VARIABLE IMPORTANCE PLOTS
###########################################################################

plot_rf_importance <- function(rf_model, title = "Random Forest Variable Importance") {
  # Get importance scores
  importance_df <- data.frame(
    PC = rownames(importance(rf_model)),
    Importance = as.vector(importance(rf_model)[,1])
  )
  
  # Order by importance
  importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE),]
  importance_df$PC <- factor(importance_df$PC, levels = importance_df$PC)
  
  # Create plot
  p <- ggplot(importance_df, aes(x = PC, y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = round(Importance, 1)), vjust = -0.3, size = 3) +
    labs(title = title,
         x = "Principal Components",
         y = "Mean Decrease in Gini") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(p)
}

###########################################################################
# 2. PARTIAL DEPENDENCE PLOTS
###########################################################################

plot_partial_dependence <- function(rf_model, train_data, target_col = "gt", 
                                   top_n = 4, grid_size = 50) {
  
  require(pdp)
  
  # Get top important variables
  importance_scores <- importance(rf_model)[,1]
  top_vars <- names(sort(importance_scores, decreasing = TRUE))[1:top_n]
  
  # Create partial dependence plots
  pd_plots <- list()
  
  for(i in 1:length(top_vars)) {
    var_name <- top_vars[i]
    
    # Calculate partial dependence
    pd_data <- partial(rf_model, pred.var = var_name, 
                      train = train_data[, !names(train_data) %in% target_col],
                      grid.size = grid_size, type = "classification", 
                      which.class = 2, prob = TRUE)
    
    # Create plot
    pd_plots[[i]] <- ggplot(pd_data, aes_string(x = var_name, y = "yhat")) +
      geom_line(color = "red", size = 1.2) +
      geom_smooth(se = TRUE, alpha = 0.3, color = "red") +
      labs(title = paste("Partial Dependence:", var_name),
           x = var_name,
           y = "Probability of GT_present") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      ylim(0, 1)
  }
  
  # Arrange plots
  do.call(grid.arrange, c(pd_plots, ncol = 2))
}

###########################################################################
# 3. PC SCATTER PLOTS WITH PREDICTIONS
###########################################################################

plot_pc_scatter <- function(pca_results, rf_model, pc_x = 1, pc_y = 2) {
  
  # Get the processed data
  data <- pca_results$data
  
  # Make predictions
  pred_prob <- predict(rf_model, data, type = "prob")[,2]
  pred_class <- predict(rf_model, data, type = "class")
  
  # Create plotting dataframe
  plot_data <- data.frame(
    PC_X = data[, paste0("PC", pc_x)],
    PC_Y = data[, paste0("PC", pc_y)],
    Actual = factor(data$gt, levels = c(0, 1), labels = c("GT_absent", "GT_present")),
    Predicted = factor(pred_class, levels = c(0, 1), labels = c("GT_absent", "GT_present")),
    Probability = pred_prob,
    Correct = ifelse(data$gt == as.numeric(pred_class) - 1, "Correct", "Incorrect")
  )
  
  # Plot 1: Actual classes with prediction probability as color
  p1 <- ggplot(plot_data, aes(x = PC_X, y = PC_Y)) +
    geom_point(aes(color = Probability, shape = Actual), size = 3, alpha = 0.8) +
    scale_color_viridis_c(name = "P(GT_present)") +
    scale_shape_manual(values = c(16, 17), name = "Actual Class") +
    labs(title = paste("PC", pc_x, "vs PC", pc_y, "- Prediction Probabilities"),
         x = paste("PC", pc_x),
         y = paste("PC", pc_y)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Plot 2: Prediction accuracy
  p2 <- ggplot(plot_data, aes(x = PC_X, y = PC_Y)) +
    geom_point(aes(color = Correct, shape = Actual), size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Correct" = "green", "Incorrect" = "red"), 
                      name = "Prediction") +
    scale_shape_manual(values = c(16, 17), name = "Actual Class") +
    labs(title = paste("PC", pc_x, "vs PC", pc_y, "- Prediction Accuracy"),
         x = paste("PC", pc_x),
         y = paste("PC", pc_y)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  grid.arrange(p1, p2, ncol = 2)
}

###########################################################################
# 4. METABOLITE CONTRIBUTION HEATMAP
###########################################################################

plot_metabolite_heatmap <- function(pca_results, rf_model, top_n_pcs = 5, top_n_metabolites = 20) {
  
  # Get loadings and importance
  loadings <- pca_results$loadings
  importance_scores <- importance(rf_model)[,1]
  
  # Get top PCs by RF importance
  top_pcs <- names(sort(importance_scores, decreasing = TRUE))[1:top_n_pcs]
  
  # For each top PC, get top contributing metabolites
  metabolite_contributions <- data.frame()
  
  for(pc in top_pcs) {
    pc_loadings <- abs(loadings[, pc])
    top_metabolites <- names(sort(pc_loadings, decreasing = TRUE))[1:top_n_metabolites]
    
    temp_df <- data.frame(
      PC = pc,
      Metabolite = top_metabolites,
      Loading = pc_loadings[top_metabolites],
      PC_Importance = importance_scores[pc]
    )
    metabolite_contributions <- rbind(metabolite_contributions, temp_df)
  }
  
  # Create heatmap data
  heatmap_data <- metabolite_contributions %>%
    select(PC, Metabolite, Loading) %>%
    tidyr::pivot_wider(names_from = PC, values_from = Loading, values_fill = 0) %>%
    column_to_rownames("Metabolite") %>%
    as.matrix()
  
  # Create heatmap
  corrplot(heatmap_data, 
           method = "color",
           col = colorRampPalette(c("white", "orange", "red"))(100),
           tl.cex = 0.8,
           tl.col = "black",
           title = "Top Metabolite Contributions to Important PCs",
           mar = c(0,0,2,0))
}

###########################################################################
# 5. PREDICTION PROBABILITY DISTRIBUTION
###########################################################################

plot_prediction_distribution <- function(rf_model, test_data, target_col = "gt") {
  
  # Get predictions
  pred_prob <- predict(rf_model, test_data, type = "prob")[,2]
  actual_class <- factor(test_data[[target_col]], levels = c(0, 1), 
                        labels = c("GT_absent", "GT_present"))
  
  # Create plotting dataframe
  plot_data <- data.frame(
    Probability = pred_prob,
    Actual = actual_class
  )
  
  # Density plot
  p1 <- ggplot(plot_data, aes(x = Probability, fill = Actual)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("GT_absent" = "blue", "GT_present" = "red")) +
    labs(title = "Distribution of Prediction Probabilities",
         x = "Probability of GT_present",
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Box plot
  p2 <- ggplot(plot_data, aes(x = Actual, y = Probability, fill = Actual)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = c("GT_absent" = "blue", "GT_present" = "red")) +
    labs(title = "Prediction Probabilities by Actual Class",
         x = "Actual Class",
         y = "Probability of GT_present") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  grid.arrange(p1, p2, ncol = 2)
}

###########################################################################
# 6. ROC CURVE AND PRECISION-RECALL CURVE
###########################################################################

# Get the data
rf_model <- pca_results$models$rf
test_data <- pca_results$test_data

# Get predictions
pred_prob <- predict(rf_model, test_data, type = "prob")[,2]
actual <- test_data$gt

cat("Predictions range:", range(pred_prob), "\n")
cat("Actual values:", table(actual), "\n")

# Create ROC curve
roc_obj <- roc(actual, pred_prob)
cat("ROC AUC:", auc(roc_obj), "\n")

# ROC plot
roc_data <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p1 <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "red", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(title = paste("ROC Curve (AUC =", round(auc(roc_obj), 3), ")"),
       x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_minimal() +
  xlim(0, 1) + ylim(0, 1)

dev.new()
print(p1)

# Manual PR curve (more reliable)
thresholds <- seq(0, 1, 0.02)
precision_vals <- numeric()
recall_vals <- numeric()

for(thresh in thresholds) {
  pred_class <- ifelse(pred_prob > thresh, 1, 0)
  
  tp <- sum(pred_class == 1 & actual == 1)
  fp <- sum(pred_class == 1 & actual == 0)
  fn <- sum(pred_class == 0 & actual == 1)
  
  if(tp + fp > 0) {
    precision <- tp / (tp + fp)
  } else {
    precision <- 1
  }
  
  if(tp + fn > 0) {
    recall <- tp / (tp + fn)
  } else {
    recall <- 0
  }
  
  precision_vals <- c(precision_vals, precision)
  recall_vals <- c(recall_vals, recall)
}

# PR plot
pr_data <- data.frame(Recall = recall_vals, Precision = precision_vals)
pr_data <- pr_data[complete.cases(pr_data), ]

p2 <- ggplot(pr_data, aes(x = Recall, y = Precision)) +
  geom_line(color = "blue", size = 1.2) +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  theme_minimal() +
  xlim(0, 1) + ylim(0, 1)

dev.new()
print(p2)

# Combine plots
combined_plot <- grid.arrange(p1, p2, ncol = 2)
dev.new()
print(combined_plot)

# Save manually
ggsave("../results/metab_rf/roc_curve.png", p1, width = 6, height = 6, dpi = 300)
ggsave("../results/metab_rf/pr_curve.png", p2, width = 6, height = 6, dpi = 300)

png("../results/metab_rf/06_roc_pr_combined.png", width = 12, height = 6, units = "in", res = 300)
grid.arrange(p1, p2, ncol = 2)
dev.off()

###########################################################################
# 7. PC BIPLOT WITH METABOLITE LOADINGS
###########################################################################

# Much simpler biplot without the complex features
plot_simple_biplot <- function(pca_results, pc_x = 1, pc_y = 2) {
  
  # Get data
  test_data <- pca_results$test_data
  pc_x_name <- paste0("PC", pc_x)
  pc_y_name <- paste0("PC", pc_y)
  
  # Simple sample plot
  plot_data <- data.frame(
    PC_X = test_data[[pc_x_name]],
    PC_Y = test_data[[pc_y_name]],
    Class = factor(test_data$gt, levels = c(0, 1), labels = c("GT_absent", "GT_present"))
  )
  
  # Create simple plot without arrows first
  p <- ggplot(plot_data, aes(x = PC_X, y = PC_Y)) +
    geom_point(aes(color = Class, shape = Class), size = 4, alpha = 0.8) +
    scale_color_manual(values = c("GT_absent" = "steelblue", "GT_present" = "darkred")) +
    scale_shape_manual(values = c("GT_absent" = 16, "GT_present" = 17)) +
    labs(title = paste("PC", pc_x, "vs PC", pc_y, "- Sample Distribution"),
         x = paste0("PC", pc_x, " (", round(pca_results$preprocessing_params$explained_variance[pc_x] * 100, 1), "%)"),
         y = paste0("PC", pc_y, " (", round(pca_results$preprocessing_params$explained_variance[pc_y] * 100, 1), "%)")) +
    theme_classic() +  # Use classic theme instead of minimal
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "white", color = "black"),
          plot.background = element_rect(fill = "white", color = NA))
  
  return(p)
}

# Test simple version first
p_simple <- plot_simple_biplot(pca_results, pc_x = 1, pc_y = 2)
dev.new()
print(p_simple)

p_simple2 <-  plot_simple_biplot(pca_results, pc_x = 10, pc_y = 7)
dev.new()
print(p_simple2)

# Save with white background explicitly
ggsave("../results/metab_rf/07_simple_biplot_pc10_pc7.png", p_simple2, 
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave("../results/metab_rf/07_simple_biplot_pc1_pc2.png", p_simple, 
       width = 10, height = 8, dpi = 300, bg = "white")

# If the simple version works, we can add metabolite arrows back
add_metabolite_arrows <- function(base_plot, pca_results, pc_x = 1, pc_y = 2, top_n = 5) {
  
  loadings <- pca_results$loadings
  
  # Get top metabolites for this PC combination
  pc_x_loadings <- abs(loadings[, pc_x])
  pc_y_loadings <- abs(loadings[, pc_y])
  
  # Combined importance for both PCs
  combined_importance <- pc_x_loadings + pc_y_loadings
  top_metabolites <- names(sort(combined_importance, decreasing = TRUE))[1:top_n]
  
  # Create arrow data
  arrow_data <- data.frame(
    x = 0,
    y = 0,
    xend = loadings[top_metabolites, pc_x] * 2,  # Scale factor
    yend = loadings[top_metabolites, pc_y] * 2,
    metabolite = substr(top_metabolites, 1, 15)  # Truncate names
  )
  
  # Add arrows to the plot
  enhanced_plot <- base_plot +
    geom_segment(data = arrow_data,
                aes(x = x, y = y, xend = xend, yend = yend),
                arrow = arrow(length = unit(0.3, "cm")),
                color = "red", alpha = 0.7, size = 1,
                inherit.aes = FALSE) +
    geom_text(data = arrow_data,
             aes(x = xend * 1.1, y = yend * 1.1, label = metabolite),
             size = 3, color = "darkred", fontface = "bold",
             inherit.aes = FALSE)
  
  return(enhanced_plot)
}

# Add arrows to the simple plot
p_with_arrows <- add_metabolite_arrows(p_simple2, pca_results, pc_x = 10, pc_y = 7, top_n = 5)
dev.new()
print(p_with_arrows)

ggsave("../results/metab_rf/09_biplot_with_arrows.png", p_with_arrows, 
       width = 12, height = 10, dpi = 300, bg = "white")

cat("Clean biplots saved!\n")


###########################################################################
# USAGE EXAMPLES
###########################################################################

# After running your main analysis, create all visualizations:

# Modified function that saves plots instead of just displaying them
create_and_save_all_rf_plots <- function(pca_results, rf_model, train_data, test_data, save_dir = "../results/metab_rf") {
  
  cat("Creating and saving Random Forest visualizations...\n")
  
  # 1. Variable importance
  cat("1. Variable Importance Plot\n")
  p1 <- plot_rf_importance(rf_model)
  ggsave(file.path(save_dir, "01_variable_importance.png"), p1, width = 10, height = 6, dpi = 300)
  
  # 2. Partial dependence plots
  cat("2. Partial Dependence Plots\n")
  if(require(pdp, quietly = TRUE)) {
    png(file.path(save_dir, "02_partial_dependence.png"), width = 12, height = 8, units = "in", res = 300)
    plot_partial_dependence(rf_model, train_data, top_n = 4)
    dev.off()
  }
  
#   # 3. PC scatter plots
#   cat("3. PC Scatter Plots\n")
#   png(file.path(save_dir, "03_pc_scatter.png"), width = 14, height = 6, units = "in", res = 300)
#   plot_pc_scatter(pca_results, rf_model, pc_x = 1, pc_y = 2)
#   dev.off()
  
#   # 4. Metabolite heatmap
#   cat("4. Metabolite Contribution Heatmap\n")
#   png(file.path(save_dir, "04_metabolite_heatmap.png"), width = 12, height = 10, units = "in", res = 300)
#   plot_metabolite_heatmap(pca_results, rf_model)
#   dev.off()
  
  # 5. Prediction distribution
  cat("5. Prediction Probability Distribution\n")
  png(file.path(save_dir, "05_prediction_distribution.png"), width = 12, height = 6, units = "in", res = 300)
  plot_prediction_distribution(rf_model, test_data)
  dev.off()
  
#   # 6. ROC and PR curves
#   cat("6. ROC and Precision-Recall Curves\n")
#   if(require(pROC, quietly = TRUE) && require(PRROC, quietly = TRUE)) {
#     png(file.path(save_dir, "06_roc_pr_curves.png"), width = 12, height = 6, units = "in", res = 300)
#     plot_roc_pr_curves(rf_model, test_data)
#     dev.off()
#   }
  
# 
}

# Call the function - you need to reconstruct train_data
# Since we don't have the original train_idx, let's use a workaround:

# SIMPLEST APPROACH - just use available data:
create_and_save_all_rf_plots(
  pca_results = pca_results,
  rf_model = pca_results$models$rf,
  train_data = pca_results$test_data,  # Not ideal but works
  test_data = pca_results$test_data,
  save_dir = "../results/metab_rf"
)

# Check what was saved
cat("Saved files:\n")
list.files("../results/metab_rf", pattern = "*.png")

# Now this will work properly!
create_and_save_all_rf_plots(
  pca_results = pca_results,
  rf_model = pca_results$models$rf,
  train_data = pca_results$train_data,  # Now properly available!
  test_data = pca_results$test_data,
  save_dir = "results"
)


















# Just run the plots individually without train_idx dependency
library(pdp)
library(pROC)
#install.packages("PRROC", dependencies = TRUE)
library(PRROC)

# 1. Variable importance
plot_rf_importance(pca_results$models$rf)

# 2. Partial dependence using ALL data (simpler)
rf_model <- pca_results$models$rf
all_features <- pca_results$data[, !names(pca_results$data) %in% "gt"]








# Get the correct PC data from your results
pc_data <- pca_results$test_data[, 1:10]  # PC1 through PC10

cat("PC data structure (correct):\n")
str(pc_data)

cat("PC data head:\n")
head(pc_data)

# Now try the partial dependence plots
top_pcs <- c("PC10", "PC5", "PC7")  # Based on your RF importance


# Create a combined dataset for all PCs
top_pcs <- c("PC10", "PC5", "PC7")
all_pd_data <- data.frame()

for(pc in top_pcs) {
  cat("Creating partial dependence data for", pc, "\n")
  
  pd_data <- partial(rf_model, pred.var = pc, train = pc_data,
                    type = "classification", which.class = 2, prob = TRUE)
  
  # Create a standardized data frame for each PC
  standardized_data <- data.frame(
    PC_name = pc,
    x_value = pd_data[[pc]],  # Get the PC values
    y_value = pd_data$yhat    # Get the predicted probabilities
  )
  
  all_pd_data <- rbind(all_pd_data, standardized_data)
}

library(ggplot2)

# Create a combined dataset for all PCs
library(ggplot2)

# Create a combined dataset for all PCs
top_pcs <- c("PC10", "PC5", "PC7")
all_pd_data <- data.frame()

for(pc in top_pcs) {
  cat("Creating partial dependence data for", pc, "\n")
  
  pd_data <- partial(rf_model, pred.var = pc, train = pc_data,
                    type = "classification", which.class = 2, prob = TRUE)
  
  # Create a standardized data frame for each PC
  standardized_data <- data.frame(
    PC_name = pc,
    x_value = pd_data[[pc]],  # Get the PC values
    y_value = pd_data$yhat    # Get the predicted probabilities
  )
  
  all_pd_data <- rbind(all_pd_data, standardized_data)
}

# Create faceted plot with 3 panels
p_faceted <- ggplot(all_pd_data, aes(x = x_value, y = y_value)) +
  geom_line(aes(color = PC_name), size = 1.2) +
  geom_smooth(se = TRUE, alpha = 0.3, color = "black") +
  facet_wrap(~ PC_name, scales = "free_x", ncol = 3) +  # 3 panels side by side
  labs(title = "Partial Dependence: Top Important PCs for GT Classification",
       x = "PC Value (Scaled)",
       y = "Probability of GT_present") +
  theme_minimal() +
  ylim(0, 1) +
  scale_color_manual(values = c("PC10" = "red", "PC5" = "blue", "PC7" = "green")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",  # Remove legend since panel titles show PC names
        strip.text = element_text(face = "bold", size = 12))

print(p_faceted)

















# Get the top 5 metabolites contributing to PC10
pc10_loadings <- abs(pca_results$loadings[, "PC10"])
top5_metabolites <- names(sort(pc10_loadings, decreasing = TRUE))[1:5]

cat("Top 5 metabolites contributing to PC10:\n")
for(i in 1:5) {
  cat(i, ".", top5_metabolites[i], "- Loading:", round(pc10_loadings[top5_metabolites[i]], 3), "\n")
}

# We need the original metabolite data (before PCA transformation)
# This should be available in your boost_data
original_metabolite_data <- boost_data[, !names(boost_data) %in% "gt"]

# Check if the metabolites exist in the original data
available_metabolites <- top5_metabolites[top5_metabolites %in% names(original_metabolite_data)]
cat("\nAvailable metabolites in original data:", length(available_metabolites), "out of", length(top5_metabolites), "\n")

if(length(available_metabolites) > 0) {
  # Create partial dependence plots for available metabolites
  all_metabolite_pd <- data.frame()
  
  for(metabolite in available_metabolites) {
    cat("Creating PDP for:", metabolite, "\n")
    
    # We need to train a simple model on the original metabolites for PDP
    # Or use the metabolite values directly with your RF model's structure
    
    # Get the metabolite values
    metabolite_values <- original_metabolite_data[[metabolite]]
    
    # Create a range of values for PDP
    metabolite_range <- seq(min(metabolite_values, na.rm = TRUE), 
                           max(metabolite_values, na.rm = TRUE), 
                           length.out = 50)
    
    # For each value, we need to transform it through PCA to get PC values
    # This is complex, so let's use a simpler approach
  }
} else {
  cat("Metabolite names don't match. Let's check the naming:\n")
  cat("First few metabolite names in loadings:\n")
  print(head(rownames(pca_results$loadings), 10))
  cat("\nFirst few names in original data:\n")
  print(head(names(original_metabolite_data), 10))
}

# Alternative approach: Show the correlation between top PC10 metabolites and GT
cat("\n=== CORRELATION ANALYSIS FOR TOP PC10 METABOLITES ===\n")

if(length(available_metabolites) > 0) {
  for(metabolite in available_metabolites[1:min(5, length(available_metabolites))]) {
    metabolite_values <- original_metabolite_data[[metabolite]]
    gt_values <- boost_data$gt
    
    # Remove missing values
    valid_idx <- !is.na(metabolite_values) & !is.na(gt_values)
    
    if(sum(valid_idx) > 0) {
      correlation <- cor(metabolite_values[valid_idx], gt_values[valid_idx])
      t_test <- t.test(metabolite_values[valid_idx] ~ gt_values[valid_idx])
      
      cat(sprintf("%s:\n", metabolite))
      cat(sprintf("  Correlation with GT: %.3f\n", correlation))
      cat(sprintf("  T-test p-value: %.4f\n", t_test$p.value))
      cat(sprintf("  Mean GT_absent: %.2f, Mean GT_present: %.2f\n\n", 
                  mean(metabolite_values[valid_idx & gt_values == 0], na.rm = TRUE),
                  mean(metabolite_values[valid_idx & gt_values == 1], na.rm = TRUE)))
    }
  }
}




# Create comprehensive plots for the top 5 PC10 metabolites
top5_metabolites <- c("X12.Hydroxyjasmonate.sulfate", "OLEOYL.GLYCEROL", "ST092794", 
                      "Massbank.PR307721.Maslinic.acid", "X2.Pentanone")

# 1. BOXPLOTS - showing distribution differences
plot_data_box <- data.frame()

for(metabolite in top5_metabolites) {
  temp_data <- data.frame(
    Metabolite = metabolite,
    Value = original_metabolite_data[[metabolite]],
    GT_Status = factor(boost_data$gt, levels = c(0, 1), labels = c("GT_absent", "GT_present"))
  )
  plot_data_box <- rbind(plot_data_box, temp_data)
}

# Create boxplots
p_metabolite_boxplots <- ggplot(plot_data_box, aes(x = GT_Status, y = Value, fill = GT_Status)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  facet_wrap(~ Metabolite, scales = "free_y", ncol = 3) +
  labs(title = "Top 5 PC10 Contributing Metabolites by GT Status",
       x = "GT Status", y = "Metabolite Abundance (Log Scale)") +
  theme_minimal() +
  scale_fill_manual(values = c("GT_absent" = "lightblue", "GT_present" = "lightcoral")) +
  scale_y_log10() +  # Log scale for better visualization
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p_metabolite_boxplots) # not convincing differences between GT_absent and GT_present
    # remember the only sig t test was for Oleoyl.glycerol

# 2. DENSITY PLOTS - showing probability distributions
p_density <- ggplot(plot_data_box, aes(x = log10(Value + 1), fill = GT_Status)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ Metabolite, scales = "free", ncol = 3) +
  labs(title = "Density Distributions of Top PC10 Metabolites",
       x = "Log10(Metabolite Abundance + 1)", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("GT_absent" = "blue", "GT_present" = "red")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom")

print(p_density)

# 3. SCATTER PLOT - Relationship with PC10 values
# Get PC10 values for each sample
pc10_values <- c(pca_results$test_data$PC10, 
                 # We need to get training PC10 values - they're in the PCA transformation
                 predict(pca_results$preprocessing_params$pca_model, 
                         as.matrix(original_metabolite_data))[setdiff(1:nrow(boost_data), 
                         as.numeric(rownames(pca_results$test_data))), "PC10"])

# Create scatter plots for each metabolite vs PC10
plot_data_scatter <- data.frame()

for(metabolite in top5_metabolites[1:3]) {  # Show top 3 to avoid overcrowding
  temp_data <- data.frame(
    Metabolite = metabolite,
    Metabolite_Value = log10(original_metabolite_data[[metabolite]] + 1),
    PC10_Value = predict(pca_results$preprocessing_params$pca_model, 
                        as.matrix(original_metabolite_data))[, "PC10"],
    GT_Status = factor(boost_data$gt, levels = c(0, 1), labels = c("GT_absent", "GT_present"))
  )
  plot_data_scatter <- rbind(plot_data_scatter, temp_data)
}

p_scatter <- ggplot(plot_data_scatter, aes(x = Metabolite_Value, y = PC10_Value, color = GT_Status)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3) +
  facet_wrap(~ Metabolite, scales = "free_x", ncol = 3) +
  labs(title = "Top PC10 Metabolites vs PC10 Scores",
       x = "Log10(Metabolite Abundance + 1)", 
       y = "PC10 Score") +
  theme_minimal() +
  scale_color_manual(values = c("GT_absent" = "blue", "GT_present" = "red")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom")

print(p_scatter)

# 4. SUMMARY TABLE of key statistics
cat("\n=== SUMMARY TABLE FOR TOP PC10 METABOLITES ===\n")
summary_table <- data.frame(
  Metabolite = top5_metabolites,
  PC10_Loading = round(pc10_loadings[top5_metabolites], 3),
  Correlation_with_GT = c(-0.115, 0.296, -0.029, -0.160, -0.206),
  P_value = c(0.3151, 0.0102, 0.8051, 0.1633, 0.0748),
  Significance = ifelse(c(0.3151, 0.0102, 0.8051, 0.1633, 0.0748) < 0.05, "***", 
                       ifelse(c(0.3151, 0.0102, 0.8051, 0.1633, 0.0748) < 0.1, "*", ""))
)

print(summary_table)









# 3. Probability distrib plots
plot_prediction_distribution(rf_model, pca_results$test_data) # good


# scatter plots- these aren't that intuitive to understand

# Get all the data we need
test_data <- pca_results$test_data
rf_model <- pca_results$models$rf

# Make predictions on test data
pred_prob <- predict(rf_model, test_data, type = "prob")[,2]  # Probability of GT_present
pred_class <- predict(rf_model, test_data, type = "class")

# Create plotting dataframe
plot_data <- data.frame(
  PC1 = test_data$PC1,
  PC2 = test_data$PC2,
  Actual = factor(test_data$gt, levels = c(0, 1), labels = c("GT_absent", "GT_present")),
  Predicted = factor(pred_class, levels = c(0, 1), labels = c("GT_absent", "GT_present")),
  Probability = pred_prob,
  Correct = ifelse(test_data$gt == as.numeric(pred_class) - 1, "Correct", "Incorrect")
)

# Plot 1: Actual classes with prediction probability as color
p1 <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Probability, shape = Actual), size = 3, alpha = 0.8) +
  scale_color_viridis_c(name = "P(GT_present)") +
  scale_shape_manual(values = c(16, 17), name = "Actual Class") +
  labs(title = "PC1 vs PC2 - Prediction Probabilities",
       x = "PC1 (36.7% variance)",
       y = "PC2 (8.6% variance)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot 2: Prediction accuracy
p2 <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Correct, shape = Actual), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Correct" = "green", "Incorrect" = "red"), 
                    name = "Prediction") +
  scale_shape_manual(values = c(16, 17), name = "Actual Class") +
  labs(title = "PC1 vs PC2 - Prediction Accuracy",
       x = "PC1 (36.7% variance)",
       y = "PC2 (8.6% variance)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Display both plots
library(gridExtra)
grid.arrange(p1, p2, ncol = 2)

# Alternative: Try with different PC combinations
# PC1 vs PC10 (since PC10 is most important)
p3 <- ggplot(plot_data, aes(x = PC1, y = test_data$PC10)) +
  geom_point(aes(color = Probability, shape = Actual), size = 3, alpha = 0.8) +
  scale_color_viridis_c(name = "P(GT_present)") +
  scale_shape_manual(values = c(16, 17), name = "Actual Class") +
  labs(title = "PC1 vs PC10 - Most Important PC",
       x = "PC1 (36.7% variance)",
       y = "PC10 (1.9% variance - Most Important)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p3)