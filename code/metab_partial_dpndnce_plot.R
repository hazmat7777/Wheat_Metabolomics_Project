# need to run metab_randomforest_log.R
source("./metab_randforest_log.R") # wait a while for this to work

# Load required libraries
install.packages("pdp")
library(pdp)
library(ggplot2)
library(gridExtra)  # For arranging multiple plots

# Get the names of your top 10 features
top_10_feature_names <- names(train_subset)[1:10]  # Exclude the 'gt' column
print("Top 10 metabolites:")
print(top_10_feature_names)

# Clean names for metabolites
clean_names <- c("Toluidine", 
                 "Prothioconazole Desethio",
                 "Dihydroxybenzaldehyde", 
                 "Undecalactone", 
                 "Soyasapogenol B", 
                 "Salicylate", 
                 "Caprylate", 
                 "Syringaldehyde", 
                 "Isovalerylglycine", 
                 "Acetyl Arginine (40.0 eV)")

# ===============================
# METHOD 1: Individual PDP plots for each metabolite
# ===============================

# List to store the PDP plots
pdp_plots <- list()

# Loop to create PDP plots
for(i in 1:length(top_10_feature_names)) {
  feature_name <- top_10_feature_names[i]
  clean_name <- clean_names[i]  # Get the clean name from the list
  
  cat("Creating PDP for:", clean_name, "\n")
  
  # Create partial dependence plot
  pd_data <- partial(rf_final_10, 
                     pred.var = feature_name,
                     prob = TRUE,          # For classification probabilities
                     which.class = "GT_present")  # Probability of disease presence
  
  # Create ggplot
  p <- autoplot(pd_data) +
    labs(title = paste("PDP:", clean_name),  # Use clean name here
         x = clean_name,                     # Use clean name here
         y = "P(Take-all Disease Present)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
  
  pdp_plots[[i]] <- p
}

# Display all plots in a grid
grid.arrange(grobs = pdp_plots, ncol = 3, nrow = 4)

# ===============================
# METHOD 2: Create a multi-panel figure (publication ready)
# ===============================

# Save as high-quality figure
png("../results/plots/metabolite_topten_PDPs.png", width = 1200, height = 1600, res = 300)
grid.arrange(grobs = pdp_plots, ncol = 3, nrow = 4)
dev.off()

# ===============================
# METHOD 3: Top 6 most important metabolites only
# ===============================

# Get importance scores for the 10-feature model
importance_10 <- importance(rf_final_10)
top_6_indices <- order(importance_10[,1], decreasing = TRUE)[1:6]
top_6_names <- rownames(importance_10)[top_6_indices]

# Create PDPs for top 6 only
pdp_top6 <- list()
for(i in 1:6) {
  feature_name <- top_6_names[i]
  
  pd_data <- partial(rf_final_10, 
                     pred.var = feature_name,
                     prob = TRUE,
                     which.class = "GT_present")
  
  p <- autoplot(pd_data) +
    labs(title = paste("PDP:", feature_name),
         x = feature_name,
         y = "P(Take-all Disease)") +
    theme_minimal()
  
  pdp_top6[[i]] <- p
}

# Display top 6 in a 2x3 grid
grid.arrange(grobs = pdp_top6, ncol = 3, nrow = 2)

# ===============================
# METHOD 4: Alternative using plotPartial for base R plots
# ===============================

# If you prefer base R plots instead of ggplot
par(mfrow = c(3, 4))  # 3 rows, 4 columns for 10 plots

for(i in 1:length(top_10_feature_names)) {
  feature_name <- top_10_feature_names[i]
  
  pd_data <- partial(rf_final_10, 
                     pred.var = feature_name,
                     prob = TRUE,
                     which.class = "GT_present")
  
  plotPartial(pd_data, 
              main = paste("PDP:", feature_name),
              xlab = feature_name,
              ylab = "P(Disease Present)")
}
par(mfrow = c(1, 1))  # Reset layout












##################################

# diagnosing what's up with salicylate

feature_name <- "Salicylate"  # Use exact name from your data

# 1. Check basic statistics
cat("Salicylate summary:\n")
print(summary(train_subset[[feature_name]]))
