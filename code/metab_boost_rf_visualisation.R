# after you get results from metab_randforest_boost_efficient.R
#results <- evaluate_metabolomics_models(boost_data, target_col = "gt", test_prop = 0.3)

# Load additional libraries for visualization
library(ggplot2)
library(gridExtra)
#install.packages(c("kableExtra)", "corrplot", "plotly"))
library(kableExtra)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(corrplot)
library(plotly)
library(knitr)

# ===============================
# 1. RESULTS SUMMARY TABLES
# ===============================

create_model_performance_table <- function(results) {
  
  # Extract individual model results
  individual_results <- results$results$individual_results
  
  # Create performance table
  perf_table <- data.frame(
    Model = names(individual_results),
    Accuracy = sapply(individual_results, function(x) round(x$accuracy, 3)),
    Kappa = sapply(individual_results, function(x) round(x$kappa, 3)),
    stringsAsFactors = FALSE
  )
  
  # Add best threshold
  perf_table$Best_Threshold <- results$results$best_threshold
  
  # Sort by Kappa (best performance metric for classification)
  perf_table <- perf_table[order(perf_table$Kappa, decreasing = TRUE), ]
  
  return(perf_table)
}

create_feature_importance_table <- function(results, top_n = 10) {
  
  importance_data <- results$results$feature_importance
  
  # Create clean table
  importance_table <- data.frame(
    Rank = 1:min(top_n, nrow(importance_data)),
    Metabolite = importance_data$var[1:min(top_n, nrow(importance_data))],
    Relative_Importance = round(importance_data$rel.inf[1:min(top_n, nrow(importance_data))], 2),
    stringsAsFactors = FALSE
  )
  
  # Clean metabolite names for better readability
  importance_table$Metabolite_Clean <- gsub("\\.", " ", importance_table$Metabolite)
  importance_table$Metabolite_Clean <- gsub("^X", "", importance_table$Metabolite_Clean)
  
  return(importance_table)
}

create_confusion_matrix_table <- function(results) {
  
  cm <- results$results$threshold_results[[as.character(results$results$best_threshold)]]$confusion_matrix
  
  # Extract confusion matrix values
  cm_table <- as.data.frame(cm$table)
  
  # Add performance metrics
  metrics <- data.frame(
    Metric = c("Accuracy", "Kappa", "Sensitivity", "Specificity", 
               "Pos Pred Value", "Neg Pred Value"),
    Value = round(c(cm$overall['Accuracy'], cm$overall['Kappa'],
                   cm$byClass['Sensitivity'], cm$byClass['Specificity'],
                   cm$byClass['Pos Pred Value'], cm$byClass['Neg Pred Value']), 3)
  )
  
  return(list(confusion_matrix = cm_table, metrics = metrics))
}

# ===============================
# 2. PARTIAL DEPENDENCE PLOTS
# ===============================

calculate_partial_dependence <- function(model, data, feature_name, model_type = "gbm", 
                                        n_points = 50, optimal_trees = NULL) {
  
  # Get feature range
  feature_range <- seq(min(data[[feature_name]], na.rm = TRUE), 
                      max(data[[feature_name]], na.rm = TRUE), 
                      length.out = n_points)
  
  # Create grid for partial dependence
  pd_data <- data[rep(1:nrow(data), each = n_points), ]
  pd_data[[feature_name]] <- rep(feature_range, nrow(data))
  
  # Get predictions
  if(model_type == "gbm") {
    predictions <- predict(model, pd_data, n.trees = optimal_trees, type = "response")
  } else if(model_type == "rf") {
    predictions <- predict(model, pd_data, type = "prob")[,2]
  } else if(model_type == "glmnet") {
    x_pd <- as.matrix(pd_data[, !names(pd_data) %in% "gt"])
    predictions <- predict(model, x_pd, s = "lambda.min", type = "response")[,1]
  }
  
  # Calculate average prediction for each feature value
  pd_result <- data.frame(
    feature_value = feature_range,
    partial_dependence = tapply(predictions, rep(1:n_points, nrow(data)), mean)
  )
  
  return(pd_result)
}

create_partial_dependence_plots <- function(results, top_features = 6) {
  
  models <- results$models
  test_data <- results$test_data
  importance <- results$results$feature_importance
  
  # Get top features
  top_feature_names <- importance$var[1:min(top_features, nrow(importance))]
  
  plot_list <- list()
  
  for(i in 1:length(top_feature_names)) {
    feature_name <- top_feature_names[i]
    
    # Calculate PD for GBM (usually most interpretable)
    pd_data <- calculate_partial_dependence(
      models$gbm, test_data, feature_name, "gbm", optimal_trees = models$optimal_trees
    )
    
    # Clean feature name for plot
    clean_name <- gsub("\\.", " ", feature_name)
    clean_name <- gsub("^X", "", clean_name)
    
    # Create plot
    p <- ggplot(pd_data, aes(x = feature_value, y = partial_dependence)) +
      geom_line(color = "blue", size = 1.2) +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.3, color = "red") +
      labs(
        title = paste("Partial Dependence:", clean_name),
        subtitle = paste("Importance:", round(importance$rel.inf[i], 1), "%"),
        x = "Metabolite Concentration (scaled)",
        y = "Predicted Probability"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10)
      )
    
    plot_list[[i]] <- p
  }
  
  return(plot_list)
}

# ===============================
# 3. FEATURE CORRELATION HEATMAP
# ===============================

create_feature_correlation_plot <- function(results) {
  
  # Get selected features data
  feature_data <- results$test_data[, !names(results$test_data) %in% "gt"]
  
  # Calculate correlation matrix
  cor_matrix <- cor(feature_data, use = "complete.obs")
  
  # Create heatmap
  corrplot(cor_matrix, 
           method = "color",
           type = "upper",
           order = "hclust",
           tl.cex = 0.8,
           tl.col = "black",
           tl.srt = 45,
           title = "Feature Correlation Matrix",
           mar = c(0,0,1,0))
}

# ===============================
# 4. MODEL COMPARISON PLOTS
# ===============================

create_model_comparison_plot <- function(results) {
  
  individual_results <- results$results$individual_results
  
  # Create comparison data
  comparison_data <- data.frame(
    Model = names(individual_results),
    Accuracy = sapply(individual_results, function(x) x$accuracy),
    Kappa = sapply(individual_results, function(x) x$kappa)
  ) %>%
    pivot_longer(cols = c(Accuracy, Kappa), names_to = "Metric", values_to = "Value")
  
  # Create plot
  p <- ggplot(comparison_data, aes(x = reorder(Model, Value), y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(aes(label = round(Value, 3)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.25, size = 3) +
    labs(
      title = "Model Performance Comparison",
      x = "Model",
      y = "Performance Score",
      fill = "Metric"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  return(p)
}

# ===============================
# 5. PREDICTION PROBABILITY DISTRIBUTION
# ===============================

create_prediction_distribution_plot <- function(results) {
  
  test_data <- results$test_data
  predictions <- results$results$predictions$Ensemble
  
  # Create data for plotting
  plot_data <- data.frame(
    Predicted_Probability = predictions,
    True_Class = factor(test_data$gt, levels = c(0, 1), labels = c("Class 0", "Class 1"))
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Predicted_Probability, fill = True_Class)) +
    geom_histogram(alpha = 0.7, bins = 15, position = "identity") +
    geom_vline(xintercept = results$results$best_threshold, 
               linetype = "dashed", color = "red", size = 1) +
    labs(
      title = "Distribution of Prediction Probabilities",
      subtitle = paste("Red line shows optimal threshold:", results$results$best_threshold),
      x = "Predicted Probability",
      y = "Count",
      fill = "True Class"
    ) +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set1")
  
  return(p)
}

# ===============================
# 6. COMPREHENSIVE VISUALIZATION FUNCTION
# ===============================

create_metabolomics_report <- function(results, output_file = "metabolomics_report.html") {
  
  cat("Creating comprehensive metabolomics analysis report...\n")
  
  # 1. Create tables
  perf_table <- create_model_performance_table(results)
  importance_table <- create_feature_importance_table(results)
  cm_tables <- create_confusion_matrix_table(results)
  
  # 2. Create plots
  pd_plots <- create_partial_dependence_plots(results)
  comparison_plot <- create_model_comparison_plot(results)
  distribution_plot <- create_prediction_distribution_plot(results)
  
  # 3. Print summary to console
  cat("\n=== MODEL PERFORMANCE SUMMARY ===\n")
  print(kable(perf_table, caption = "Model Performance Comparison"))
  
  cat("\n=== TOP METABOLITE BIOMARKERS ===\n")
  print(kable(importance_table[,c("Rank", "Metabolite_Clean", "Relative_Importance")], 
              caption = "Feature Importance Ranking"))
  
  cat("\n=== CONFUSION MATRIX METRICS ===\n")
  print(kable(cm_tables$metrics, caption = "Classification Performance Metrics"))
  
  # 4. Display plots
  cat("\nGenerating visualizations...\n")
  
  # Arrange partial dependence plots
  if(length(pd_plots) >= 4) {
    grid.arrange(pd_plots[[1]], pd_plots[[2]], pd_plots[[3]], pd_plots[[4]], 
                 ncol = 2, top = "Partial Dependence Plots - Top 4 Metabolites")
  }
  
  # Show other plots
  print(comparison_plot)
  print(distribution_plot)
  
  # 5. Create correlation plot
  cat("Creating feature correlation heatmap...\n")
  create_feature_correlation_plot(results)
  
  return(list(
    performance_table = perf_table,
    importance_table = importance_table,
    confusion_matrices = cm_tables,
    partial_dependence_plots = pd_plots,
    comparison_plot = comparison_plot,
    distribution_plot = distribution_plot
  ))
}


# ===============================
# USAGE EXAMPLE
# ===============================


# # Generate comprehensive report
report <- create_metabolomics_report(stable_results)
report
# 
# # Create individual plots
pd_plots <- create_partial_dependence_plots(stable_results, top_features = 6)
str(pd_plots)
length(pd_plots)

# Display the pdp plots
grid.arrange(grobs = pd_plots, nrow = 2)

comparison_plot <- create_model_comparison_plot(results)
comparison_plot
# 
