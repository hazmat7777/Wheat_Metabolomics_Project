library(ggplot2)

# first load up the rf models
rf_models <- readRDS("../results/microb_rf_oob_error")

# ===== FEATURE IMPORTANCE =====

# Function to get top important features
get_top_features <- function(rf_model, n_features = 10) {
  importance_scores <- importance(rf_model)
  
  # For classification, use MeanDecreaseAccuracy
  imp_col <- "MeanDecreaseAccuracy"
  
  top_features <- head(
    importance_scores[order(importance_scores[, imp_col], decreasing = TRUE), ],
    n_features
  )
  
  return(top_features)
}

# Get top 10 features for each model
cat("\nTop 10 important features for each approach:\n")
for (model_name in names(rf_models)) {
  cat("\n", model_name, ":\n")
  top_feat <- get_top_features(rf_models[[model_name]])
  print(top_feat)
}

top_order_features <- row.names(get_top_features(rf_models[[3]]))
top_order_features
order_glom<- (as.data.frame(aggregated_data[3]))
order_glom$gt <- as.factor(as.data.frame(sample_data(ps_16S_rf))$gt) # check order didnt change

# # Plot the top 5 features for the order model (why, this one was the worst)
# counter <- 0
# for (microbe in top_order_features) { # dont forget to set counter
#     counter <- counter + 1
#     microbe <- paste0("Order.",microbe)
#     plot <- ggplot(order_glom, aes(y = gt, x = .data[[microbe]])) +
#         geom_boxplot() +
#         labs(x = "GT Status", y = paste("Distribution of", microbe)) +
#         theme_bw()
    
#     ggsave(filename = paste0("../results/plots/relative_microb_", counter, "_vs_gt.png"), plot = plot, width = 8, height = 6)
#     if(counter == 5){
#         counter <- 0
#     }
# }

## same thing for class
top_class_features <- row.names(get_top_features(rf_models[[4]]))
top_class_features
class_glom <- (as.data.frame(aggregated_data[4]))
class_glom$gt <- as.factor(as.data.frame(sample_data(ps_16S_rf))$gt) # check class didnt change

counter <- 0  # Initialize counter
for (microbe in top_class_features[1:5]) { # dont forget to set counter
    counter <- counter + 1
    microbe <- paste0("Class.",microbe)
    microbe_nicename <- gsub("\\.", " ", microbe)
    plot <- ggplot(class_glom, aes(y = gt, x = .data[[microbe]])) +
        geom_boxplot() +
        labs(x = paste("Relative Abundance of", microbe_nicename), y = "GT Status") +
        theme_bw()
    ggsave(filename = paste0("../results/plots/relative_microb_", counter, "_vs_gt.png"), plot = plot, width = 8, height = 6)
    if(counter == 5){
    counter <- 0
    }
}
