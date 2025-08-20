# Random Forest Visualization Functions for PCA Metabolomics Data
# to be run AFTER sourcing metab_pca_then_rf.R

source("./metab_pca_then_rf2.R")  # Ensure PCA results are available


# Combined Random Forest Analysis: Figure and Results Table
# Run this after sourcing your main analysis script

library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)
library(pROC)
library(knitr)
library(kableExtra)

###########################################################################
# COMBINED ANALYSIS FIGURE
###########################################################################

create_combined_rf_figure <- function(pca_results, rf_model, test_data, 
                                     save_path = "results/rf_combined_analysis.png") {
  
  cat("Creating combined Random Forest analysis figure...\n")
  
  # --- PANEL A: Variable Importance ---
  importance_df <- data.frame(
    PC = rownames(importance(rf_model)),
    Importance = as.vector(importance(rf_model)[,1])
  )
  importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE),]
  importance_df <- importance_df[1:8,]  # Top 8 PCs
  importance_df$PC <- factor(importance_df$PC, levels = importance_df$PC)
  
  panel_a <- ggplot(importance_df, aes(x = PC, y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = round(Importance, 1)), vjust = -0.3, size = 2.5) +
    labs(title = "A) Variable Importance", 
         x = "Principal Components", 
         y = "Mean Decrease Gini") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9)
    )
  
  # --- PANEL B: PC Scatter with Predictions ---
  pred_prob <- predict(rf_model, test_data, type = "prob")[,2]
  pred_class <- predict(rf_model, test_data, type = "class")
  
  # Use most important PCs
  top_pcs <- importance_df$PC[1:2]
  pc1_name <- as.character(top_pcs[1])
  pc2_name <- as.character(top_pcs[2])
  
  scatter_data <- data.frame(
    PC1 = test_data[[pc1_name]],
    PC2 = test_data[[pc2_name]],
    Actual = factor(test_data$gt, levels = c(0, 1), labels = c("GT absent", "GT present")),
    Probability = pred_prob,
    Correct = ifelse(test_data$gt == as.numeric(pred_class) - 1, "Correct", "Incorrect")
  )
  
  panel_b <- ggplot(scatter_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Probability, shape = Actual), size = 2.5, alpha = 0.8) +
    scale_color_viridis_c(name = "P(GT+)", guide = guide_colorbar(barwidth = 0.8)) +
    scale_shape_manual(values = c(16, 17), name = "Actual") +
    labs(title = paste0("B) ", pc1_name, " vs ", pc2_name, " - Prediction Probabilities"),
         x = pc1_name, y = pc2_name) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    )
  
  # --- PANEL C: Prediction Distribution ---
  actual_class <- factor(test_data$gt, levels = c(0, 1), labels = c("GT absent", "GT present"))
  dist_data <- data.frame(Probability = pred_prob, Actual = actual_class)
  
  panel_c <- ggplot(dist_data, aes(x = Probability, fill = Actual)) +
    geom_density(alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.8) +
    scale_fill_manual(values = c("GT absent" = "#3498db", "GT present" = "#e74c3c")) +
    labs(title = "C) Prediction Probability Distribution",
         x = "P(GT present)", y = "Density") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
  
  # --- PANEL D: ROC Curve ---
  roc_obj <- roc(test_data$gt, pred_prob, quiet = TRUE)
  roc_data <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )
  
  panel_d <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(color = "#e74c3c", size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
    labs(title = paste0("D) ROC Curve (AUC = ", round(auc(roc_obj), 3), ")"),
         x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ) +
    xlim(0, 1) + ylim(0, 1)
  
  # --- COMBINE PANELS ---
  # Create a more sophisticated layout
  top_row <- grid.arrange(panel_a, panel_b, ncol = 2, widths = c(1, 1.2))
  bottom_row <- grid.arrange(panel_c, panel_d, ncol = 2)
  
  # Save the combined figure
  png(save_path, width = 14, height = 10, units = "in", res = 300, bg = "white")
  grid.arrange(top_row, bottom_row, nrow = 2, heights = c(1, 1))
  dev.off()
  
  # Also display in current device
  dev.new(width = 14, height = 10)
  grid.arrange(top_row, bottom_row, nrow = 2, heights = c(1, 1))
  
  cat("Combined figure saved to:", save_path, "\n")
  
  return(list(
    panel_a = panel_a,
    panel_b = panel_b, 
    panel_c = panel_c,
    panel_d = panel_d,
    roc_auc = auc(roc_obj)
  ))
}

###########################################################################
# RESULTS TABLE
###########################################################################

create_results_table <- function(pca_results, rf_model, test_data, 
                                save_path = "results/rf_results_table") {
  
  cat("Creating results summary table...\n")
  
  # Get predictions
  pred_prob <- predict(rf_model, test_data, type = "prob")[,2]
  pred_class <- predict(rf_model, test_data, type = "class")
  actual <- test_data$gt
  
  # Calculate performance metrics
  cm <- table(Actual = actual, Predicted = as.numeric(pred_class) - 1)
  
  # Basic metrics
  accuracy <- sum(diag(cm)) / sum(cm)
  sensitivity <- cm[2,2] / sum(cm[2,])  # True Positive Rate
  specificity <- cm[1,1] / sum(cm[1,])  # True Negative Rate
  precision <- cm[2,2] / sum(cm[,2])    # Positive Predictive Value
  npv <- cm[1,1] / sum(cm[,1])          # Negative Predictive Value
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  
  # ROC AUC
  roc_obj <- roc(actual, pred_prob, quiet = TRUE)
  roc_auc <- as.numeric(auc(roc_obj))
  
  # Model complexity metrics
  n_features <- ncol(test_data) - 1  # Excluding target variable
  n_trees <- rf_model$ntree
  mtry <- rf_model$mtry
  
  # Variable importance summary
  importance_scores <- importance(rf_model)[,1]
  top_5_pcs <- names(sort(importance_scores, decreasing = TRUE))[1:5]
  
  # Create comprehensive results table
  results_table <- data.frame(
    Metric = c(
      "Sample Size (Test)", "Number of Features", "Number of Trees", "mtry Parameter",
      "", # Spacer
      "Accuracy", "Sensitivity (TPR)", "Specificity (TNR)", 
      "Precision (PPV)", "Negative Pred. Value", "F1-Score", "ROC AUC",
      "", # Spacer
      "True Negatives", "False Positives", "False Negatives", "True Positives"
    ),
    Value = c(
      nrow(test_data), n_features, n_trees, mtry,
      "", # Spacer
      round(accuracy, 3), round(sensitivity, 3), round(specificity, 3),
      round(precision, 3), round(npv, 3), round(f1_score, 3), round(roc_auc, 3),
      "", # Spacer
      cm[1,1], cm[1,2], cm[2,1], cm[2,2]
    ),
    Interpretation = c(
      "Test set size", "Principal components used", "Random forest trees", "Variables per split",
      "",
      "Overall classification accuracy", "Proportion of GT+ correctly identified", 
      "Proportion of GT- correctly identified", "Proportion of predicted GT+ that are true GT+",
      "Proportion of predicted GT- that are true GT-", "Harmonic mean of precision and recall",
      "Area under ROC curve",
      "",
      "Correctly classified GT absent", "Incorrectly classified as GT present",
      "Missed GT present cases", "Correctly classified GT present"
    )
  )
  
  # Create a second table for top important variables
  importance_table <- data.frame(
    Rank = 1:5,
    Component = top_5_pcs,
    Importance = round(importance_scores[top_5_pcs], 2),
    `Relative_Importance_%` = round(100 * importance_scores[top_5_pcs] / max(importance_scores), 1)
  )
  
  # Save as CSV
  write.csv(results_table, paste0(save_path, "_performance.csv"), row.names = FALSE)
  write.csv(importance_table, paste0(save_path, "_importance.csv"), row.names = FALSE)
  
  # Create formatted tables for display/reporting
  if(require(kableExtra, quietly = TRUE)) {
    
    # Performance table
    perf_table <- results_table %>%
      kable(format = "html", 
            col.names = c("Metric", "Value", "Description"),
            caption = "Random Forest Model Performance Summary") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                   full_width = FALSE) %>%
      row_spec(which(results_table$Metric == ""), background = "white", color = "white") %>%
      row_spec(1:4, background = "#f8f9fa") %>%
      row_spec(6:12, background = "#e8f4fd") %>%
      row_spec(14:17, background = "#fff3cd")
    
    # Importance table  
    imp_table <- importance_table %>%
      kable(format = "html",
            col.names = c("Rank", "Principal Component", "Gini Decrease", "Relative Importance (%)"),
            caption = "Top 5 Most Important Principal Components") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                   full_width = FALSE)
    
    # Save HTML tables
    cat(as.character(perf_table), file = paste0(save_path, "_performance.html"))
    cat(as.character(imp_table), file = paste0(save_path, "_importance.html"))
  }
  
  # Print summary to console
  cat("\n=== RANDOM FOREST ANALYSIS RESULTS ===\n")
  cat("Model Performance:\n")
  cat(sprintf("  Accuracy: %.3f\n", accuracy))
  cat(sprintf("  Sensitivity: %.3f\n", sensitivity))
  cat(sprintf("  Specificity: %.3f\n", specificity))
  cat(sprintf("  ROC AUC: %.3f\n", roc_auc))
  cat(sprintf("  F1-Score: %.3f\n", f1_score))
  cat("\nConfusion Matrix:\n")
  print(cm)
  cat("\nTop 5 Important Components:\n")
  print(importance_table)
  cat("\nFiles saved:\n")
  cat(sprintf("  %s_performance.csv\n", save_path))
  cat(sprintf("  %s_importance.csv\n", save_path))
  if(require(kableExtra, quietly = TRUE)) {
    cat(sprintf("  %s_performance.html\n", save_path))
    cat(sprintf("  %s_importance.html\n", save_path))
  }
  
  return(list(
    performance = results_table,
    importance = importance_table,
    confusion_matrix = cm,
    metrics = list(
      accuracy = accuracy, sensitivity = sensitivity, specificity = specificity,
      precision = precision, f1_score = f1_score, roc_auc = roc_auc
    )
  ))
}

###########################################################################
# MAIN EXECUTION FUNCTION
###########################################################################

create_publication_ready_analysis <- function(pca_results, save_dir = "results") {
  
  # Ensure save directory exists
  if(!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Extract components
  rf_model <- pca_results$models$rf
  test_data <- pca_results$test_data
  
  # Create combined figure
  figure_path <- file.path(save_dir, "rf_combined_analysis.png")
  plots <- create_combined_rf_figure(pca_results, rf_model, test_data, figure_path)
  
  # Create results tables
  table_path <- file.path(save_dir, "rf_results_table")
  results <- create_results_table(pca_results, rf_model, test_data, table_path)
  
  cat("\n=== PUBLICATION-READY ANALYSIS COMPLETE ===\n")
  cat("Generated files:\n")
  cat("  - Combined figure:", figure_path, "\n")
  cat("  - Results tables:", paste0(table_path, "_*.csv/html"), "\n")
  
  return(list(
    plots = plots,
    results = results,
    paths = list(figure = figure_path, tables = table_path)
  ))
}

###########################################################################
# USAGE
###########################################################################

# Run the complete analysis
# Make sure your pca_results object is loaded first!

# Execute the analysis
publication_analysis <- create_publication_ready_analysis(pca_results, save_dir = "results")

# The function will:
# 1. Create a 4-panel combined figure showing:
#    - Variable importance 
#    - PC scatter with prediction probabilities
#    - Prediction probability distributions
#    - ROC curve
# 2. Generate comprehensive results tables in CSV and HTML format
# 3. Print a summary to console

cat("Analysis complete! Check the 'results' folder for outputs.\n")