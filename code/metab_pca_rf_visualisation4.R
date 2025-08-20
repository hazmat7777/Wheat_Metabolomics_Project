# Simplified Random Forest Variable Importance Comparison
# Compares Mean Decrease Gini vs Mean Decrease Accuracy
# to be run AFTER sourcing metab_pca_then_rf.R

source("./metab_pca_then_rf3.R")  # Ensure PCA results are available

library(ggplot2)
library(gridExtra)

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
  
  # --- COMBINE PANELS ---
  # Save the combined figure
  png(save_path, width = 12, height = 6, units = "in", res = 300, bg = "white")
  grid.arrange(panel_a, panel_b, ncol = 2)
  dev.off()
  
  # Also display in current device
  grid.arrange(panel_a, panel_b, ncol = 2)
  
  cat("Importance comparison figure saved to:", save_path, "\n")
  
  return(list(
    panel_a = panel_a,
    panel_b = panel_b,
    gini_top8 = gini_df,
    accuracy_top8 = acc_df
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