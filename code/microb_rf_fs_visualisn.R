# boxplots on the top 2 features from the Random Forest feature selection
# boxplots didnt work due to all the zeroes.

# Load required libraries
library(phyloseq)
library(ggplot2)

# Extract data from phyloseq object
sample_df <- data.frame(sample_data(ps_fs))
otu_df <- data.frame(t(as.matrix(otu_table(ps_fs))))

# Combine sample data with OTU data
combined_data <- cbind(sample_df, otu_df)

# Get the top 20 feature names from your Random Forest results
# (assuming you have top_features from your RF analysis)
top_20_names <- names(top_features)

# Perform t-tests for each of the top 20 features
ttest_results <- data.frame(
  feature = character(0),
  t_statistic = numeric(0),
  p_value = numeric(0),
  mean_GT_present = numeric(0),
  mean_GT_absent = numeric(0),
  mean_diff = numeric(0),
  abs_mean_diff = numeric(0),
  stringsAsFactors = FALSE
)

cat("Performing t-tests for top 20 features...\n")

for(feature in top_20_names) {
  # Get data for each group
  present_group <- combined_data[combined_data$gt == "GT_present", feature]
  absent_group <- combined_data[combined_data$gt == "GT_absent", feature]
  
  # Remove any NAs
  present_group <- present_group[!is.na(present_group)]
  absent_group <- absent_group[!is.na(absent_group)]
  
  # Perform t-test
  if(length(present_group) > 1 && length(absent_group) > 1) {
    ttest <- t.test(present_group, absent_group)
    
    # Calculate means and difference
    mean_present <- mean(present_group)
    mean_absent <- mean(absent_group)
    mean_diff <- mean_present - mean_absent
    
    # Store results
    ttest_results <- rbind(ttest_results, data.frame(
      feature = feature,
      t_statistic = as.numeric(ttest$statistic),
      p_value = ttest$p.value,
      mean_GT_present = mean_present,
      mean_GT_absent = mean_absent,
      mean_diff = mean_diff,
      abs_mean_diff = abs(mean_diff),
      stringsAsFactors = FALSE
    ))
  }
}

# Sort by absolute mean difference (biggest difference between groups)
ttest_results <- ttest_results[order(ttest_results$abs_mean_diff, decreasing = TRUE), ]

# Display results
cat("\nT-test results for top 20 features (sorted by absolute mean difference):\n")
print(ttest_results)

# Get the top 2 features with biggest difference
top_2_features <- ttest_results$feature[1:2]
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("TOP 2 FEATURES WITH BIGGEST DIFFERENCE BETWEEN GROUPS:\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("1.", top_2_features[1], "\n")
cat("   Mean difference:", round(ttest_results$abs_mean_diff[1], 6), "\n")
cat("   p-value:", format.pval(ttest_results$p_value[1], digits = 3), "\n\n")
cat("2.", top_2_features[2], "\n")
cat("   Mean difference:", round(ttest_results$abs_mean_diff[2], 6), "\n")
cat("   p-value:", format.pval(ttest_results$p_value[2], digits = 3), "\n")

# Get taxonomy information for the top 2 features
tax_df <- as.data.frame(tax_table(ps_fs))

# Function to get taxonomic info for a feature
get_tax_info <- function(esv_name) {
  if(esv_name %in% rownames(tax_df)) {
    genus <- tax_df[esv_name, "Genus"]
    species <- tax_df[esv_name, "Species"]
    
    # Handle NAs or empty values
    genus <- ifelse(is.na(genus) || genus == "", "Unknown_Genus", genus)
    species <- ifelse(is.na(species) || species == "", "Unknown_Species", species)
    
    return(paste(genus, species))
  } else {
    return("Unknown_Taxonomy")
  }
}

# Get taxonomy for top 2 features
tax_info_1 <- get_tax_info(top_2_features[1])
tax_info_2 <- get_tax_info(top_2_features[2])

# Create boxplots for the top 2 features
# Feature 1
plot_data_1 <- data.frame(
  gt = combined_data$gt,
  abundance = combined_data[[top_2_features[1]]]
)

p1 <- ggplot(plot_data_1, aes(x = gt, y = abundance, fill = gt)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("GT_present" = "#E74C3C", "GT_absent" = "#3498DB")) +
  labs(
    title = paste(tax_info_1, "\n(", top_2_features[1], ")", sep=""),
    subtitle = paste("p-value:", format.pval(ttest_results$p_value[1], digits = 3),
                    "| Mean diff:", round(ttest_results$abs_mean_diff[1], 6)),
    x = "Ground Truth",
    y = "Relative Abundance"
  ) +
  theme_bw() +  # White background theme
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )

# Feature 2
plot_data_2 <- data.frame(
  gt = combined_data$gt,
  abundance = combined_data[[top_2_features[2]]]
)

p2 <- ggplot(plot_data_2, aes(x = gt, y = abundance, fill = gt)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("GT_present" = "#E74C3C", "GT_absent" = "#3498DB")) +
  labs(
    title = paste(tax_info_2, "\n(", top_2_features[2], ")", sep=""),
    subtitle = paste("p-value:", format.pval(ttest_results$p_value[2], digits = 3),
                    "| Mean diff:", round(ttest_results$abs_mean_diff[2], 6)),
    x = "Ground Truth",
    y = "Relative Abundance"
  ) +
  theme_bw() +  # White background theme
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )

# Save plots to files (this should always work)
ggsave("../results/microb/top_feature_1_boxplot.png", p1, width = 8, height = 6, dpi = 300)
ggsave("../results/microb/top_feature_2_boxplot.png", p2, width = 8, height = 6, dpi = 300)

cat("\nPlots saved as:\n")
cat("- top_feature_1_boxplot.png (", top_2_features[1], ")\n")
cat("- top_feature_2_boxplot.png (", top_2_features[2], ")\n")
cat("Check your working directory for the image files.\n")

# Save the t-test results
#write.csv(ttest_results, "ttest_results_top20_features.csv", row.names = FALSE)
cat("\nT-test results saved to: ttest_results_top20_features.csv\n")