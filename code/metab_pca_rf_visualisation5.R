# Simplified Random Forest Variable Importance Comparison
# Compares Mean Decrease Gini vs Mean Decrease Accuracy
# to be run AFTER sourcing metab_pca_then_rf.R

source("./metab_pca_then_rf3.R")  # Ensure PCA results are available

library(ggplot2)
library(gridExtra)
library(viridis)

###########################################################################
# SIMPLIFIED VARIABLE IMPORTANCE COMPARISON
###########################################################################

create_importance_comparison <- function(pca_results, rf_model, test_data, 
                                        save_path = "results/rf_importance_comparison.png") {
  
  cat("Creating variable importance comparison figure...\n")
  
  # Get importance directly from model object
  rf_importance <- rf_model$importance
  importance_gini <- rf_importance[, "MeanDecreaseGini"]
  importance_acc <- rf_importance[, "MeanDecreaseAccuracy"]
  
  # --- PANEL A: Mean Decrease Gini ---
  gini_df <- data.frame(
    PC = names(importance_gini),
    Importance = as.vector(importance_gini)
  )
  gini_df <- gini_df[order(gini_df$Importance, decreasing = TRUE),]
  gini_df <- gini_df[1:8,]  # Top 8 PCs
  gini_df$PC <- factor(gini_df$PC, levels = gini_df$PC)
  
  panel_a <- ggplot(gini_df, aes(x = PC, y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = round(Importance, 1)), vjust = -0.3, size = 3) +
    labs(title = "A) Mean Decrease Gini", 
         x = "Principal Components", 
         y = "Mean Decrease Gini") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
  
  # --- PANEL B: Mean Decrease Accuracy ---
  acc_df <- data.frame(
    PC = names(importance_acc),
    Importance = as.vector(importance_acc)
  )
  acc_df <- acc_df[order(acc_df$Importance, decreasing = TRUE),]
  acc_df <- acc_df[1:8,]  # Top 8 PCs
  acc_df$PC <- factor(acc_df$PC, levels = acc_df$PC)
  
  panel_b <- ggplot(acc_df, aes(x = PC, y = Importance)) +
    geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.8) +
    geom_text(aes(label = round(Importance, 3)), vjust = -0.3, size = 3) +
    labs(title = "B) Mean Decrease Accuracy", 
         x = "Principal Components", 
         y = "Mean Decrease Accuracy") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
  
  # --- PANEL C: PC10 vs PC5 Ordination Plot ---
  pred_prob <- predict(rf_model, test_data, type = "prob")[,2]
  pred_class <- predict(rf_model, test_data, type = "class")
  
  ordination_data <- data.frame(
    PC10 = test_data$PC10,
    PC5 = test_data$PC5,
    Actual = factor(test_data$gt, levels = c(0, 1), labels = c("GT absent", "GT present")),
    Probability = pred_prob,
    Correct = ifelse(test_data$gt == as.numeric(pred_class) - 1, "Correct", "Incorrect")
  )
  
  panel_c <- ggplot(ordination_data, aes(x = PC10, y = PC5)) +
    geom_point(aes(color = Probability, shape = Actual), size = 3, alpha = 0.8) +
    scale_color_viridis_c(name = "P(GT+)", guide = guide_colorbar(barwidth = 0.8)) +
    scale_shape_manual(values = c(16, 17), name = "Actual") +
    labs(title = "C) PC10 vs PC5 - Sample Ordination",
         x = "PC10", y = "PC5") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.key.size = unit(0.5, "cm")
    )
  
  # --- COMBINE PANELS ---
  # Top row: Both importance plots
  top_row <- grid.arrange(panel_a, panel_b, ncol = 2)
  
  # Save the combined figure (3 panels: 2 on top, 1 below)
  png(save_path, width = 12, height = 9, units = "in", res = 300, bg = "white")
  grid.arrange(top_row, panel_c, nrow = 2, heights = c(1, 1))
  dev.off()
  
  # Also display in current device
  grid.arrange(top_row, panel_c, nrow = 2, heights = c(1, 1))
  
  cat("Importance comparison figure saved to:", save_path, "\n")
  
  return(list(
    panel_a = panel_a,
    panel_b = panel_b,
    panel_c = panel_c,
    gini_top8 = gini_df,
    accuracy_top8 = acc_df,
    ordination_data = ordination_data
  ))
}

###########################################################################
# SIMPLIFIED RESULTS TABLE
###########################################################################

create_importance_table <- function(rf_model, save_path = "results/rf_importance_comparison.csv") {
  
  cat("Creating importance comparison table...\n")
  
  # Get importance directly from model object
  rf_importance <- rf_model$importance
  importance_gini <- rf_importance[, "MeanDecreaseGini"]
  importance_acc <- rf_importance[, "MeanDecreaseAccuracy"]
  
  # Create comparison table
  importance_comparison <- data.frame(
    Component = names(importance_gini),
    Mean_Decrease_Gini = round(importance_gini, 2),
    Gini_Rank = rank(-importance_gini),
    Mean_Decrease_Accuracy = round(importance_acc, 4),
    Accuracy_Rank = rank(-importance_acc)
  )
  
  # Order by Gini importance
  importance_comparison <- importance_comparison[order(importance_comparison$Mean_Decrease_Gini, decreasing = TRUE),]
  
  # Save as CSV
  write.csv(importance_comparison, save_path, row.names = FALSE)
  
  # Print summary to console
  cat("\n=== VARIABLE IMPORTANCE COMPARISON ===\n")
  cat("Top 5 Components (Mean Decrease Gini):\n")
  print(importance_comparison[1:5, c("Component", "Mean_Decrease_Gini", "Gini_Rank")])
  cat("\nTop 5 Components (Mean Decrease Accuracy):\n")
  acc_top5 <- importance_comparison[order(importance_comparison$Mean_Decrease_Accuracy, decreasing = TRUE),][1:5,]
  print(acc_top5[, c("Component", "Mean_Decrease_Accuracy", "Accuracy_Rank")])
  cat("\nFull comparison table saved to:", save_path, "\n")
  
  return(importance_comparison)
}

###########################################################################
# MAIN FUNCTION
###########################################################################

compare_rf_importance <- function(pca_results, save_dir = "results") {
  
  # Ensure save directory exists
  if(!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Extract components
  rf_model <- pca_results$models$rf
  test_data <- pca_results$test_data
  
  # Create comparison figure
  figure_path <- file.path(save_dir, "rf_importance_comparison.png")
  plots <- create_importance_comparison(pca_results, rf_model, test_data, figure_path)
  
  # Create comparison table
  table_path <- file.path(save_dir, "rf_importance_comparison.csv")
  importance_table <- create_importance_table(rf_model, table_path)
  
  cat("\n=== IMPORTANCE COMPARISON COMPLETE ===\n")
  cat("Generated files:\n")
  cat("  - Figure:", figure_path, "\n")
  cat("  - Table:", table_path, "\n")
  
  return(list(
    plots = plots,
    importance_table = importance_table,
    paths = list(figure = figure_path, table = table_path)
  ))
}

###########################################################################
# USAGE
###########################################################################

# Run the simplified comparison
# Make sure your pca_results object is loaded first!

# Execute the analysis
importance_analysis <- compare_rf_importance(pca_results, save_dir = "results")

cat("Importance comparison complete! Check the 'results' folder for outputs.\n")