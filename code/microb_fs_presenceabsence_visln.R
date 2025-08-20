# Load required libraries
library(phyloseq)
library(ggplot2)

# Extract data from phyloseq object (assuming you have ps_fs or use your main phyloseq object)
# If you want to analyze ALL features, use your main phyloseq object instead of ps_fs
# Change ps_fs to your main phyloseq object name if needed
sample_df <- data.frame(sample_data(ps_fs))  # Change ps_fs if using different object
otu_df <- data.frame(t(as.matrix(otu_table(ps_fs))))  # Change ps_fs if using different object

# Combine sample data with OTU data
combined_data <- cbind(sample_df, otu_df)

# Get all feature names (excluding sample metadata columns)
feature_names <- colnames(otu_df)

# Calculate detection frequency differences for ALL features
detection_results <- data.frame(
  feature = character(0),
  prev_GT_present = numeric(0),
  prev_GT_absent = numeric(0),
  prev_diff = numeric(0),
  abs_prev_diff = numeric(0),
  fisher_p = numeric(0),
  n_detected_present = numeric(0),
  n_detected_absent = numeric(0),
  n_total_present = numeric(0),
  n_total_absent = numeric(0),
  stringsAsFactors = FALSE
)

cat("Analyzing detection frequencies for", length(feature_names), "features...\n")

for(feature in feature_names) {
  # Get data for each group
  present_group <- combined_data[combined_data$gt == "GT_present", feature]
  absent_group <- combined_data[combined_data$gt == "GT_absent", feature]
  
  # Remove any NAs
  present_group <- present_group[!is.na(present_group)]
  absent_group <- absent_group[!is.na(absent_group)]
  
  if(length(present_group) > 0 && length(absent_group) > 0) {
    # Calculate detection (presence/absence)
    detected_present <- sum(present_group > 0)
    detected_absent <- sum(absent_group > 0)
    total_present <- length(present_group)
    total_absent <- length(absent_group)
    
    # Calculate prevalence (detection frequency)
    prev_present <- detected_present / total_present
    prev_absent <- detected_absent / total_absent
    prev_diff <- prev_present - prev_absent
    
    # Fisher's exact test for presence/absence
    fisher_test <- NA
    if(total_present > 0 && total_absent > 0) {
      contingency_table <- matrix(c(
        detected_present, total_present - detected_present,
        detected_absent, total_absent - detected_absent
      ), nrow = 2, byrow = TRUE)
      
      fisher_test <- fisher.test(contingency_table)$p.value
    }
    
    # Store results
    detection_results <- rbind(detection_results, data.frame(
      feature = feature,
      prev_GT_present = prev_present,
      prev_GT_absent = prev_absent,
      prev_diff = prev_diff,
      abs_prev_diff = abs(prev_diff),
      fisher_p = fisher_test,
      n_detected_present = detected_present,
      n_detected_absent = detected_absent,
      n_total_present = total_present,
      n_total_absent = total_absent,
      stringsAsFactors = FALSE
    ))
  }
}

# Sort by absolute prevalence difference (biggest detection frequency difference)
detection_results <- detection_results[order(detection_results$abs_prev_diff, decreasing = TRUE), ]

# Display top results
cat("\nTop 20 features with biggest detection frequency differences:\n")
print(head(detection_results, 20))

# Get the top 2 features with biggest detection frequency difference
top_2_features <- detection_results$feature[1:2]

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("TOP 2 FEATURES WITH BIGGEST DETECTION FREQUENCY DIFFERENCES:\n")
cat(paste(rep("=", 80), collapse=""), "\n")

for(i in 1:2) {
  feature <- detection_results$feature[i]
  cat(i, ".", feature, "\n")
  cat("   Detection in GT_present:", round(detection_results$prev_GT_present[i] * 100, 1), "% (", 
      detection_results$n_detected_present[i], "/", detection_results$n_total_present[i], ")\n")
  cat("   Detection in GT_absent: ", round(detection_results$prev_GT_absent[i] * 100, 1), "% (", 
      detection_results$n_detected_absent[i], "/", detection_results$n_total_absent[i], ")\n")
  cat("   Difference:", round(detection_results$prev_diff[i] * 100, 1), "percentage points\n")
  cat("   Fisher's p-value:", format.pval(detection_results$fisher_p[i], digits = 3), "\n\n")
}

# Get taxonomy information for the top 2 features
tax_df <- as.data.frame(tax_table(ps_fs))  # Change ps_fs if using different object

get_tax_info <- function(esv_name) {
  if(esv_name %in% rownames(tax_df)) {
    genus <- tax_df[esv_name, "Genus"]
    species <- tax_df[esv_name, "Species"]
    
    genus <- ifelse(is.na(genus) || genus == "", "Unknown_Genus", genus)
    species <- ifelse(is.na(species) || species == "", "Unknown_Species", species)
    
    return(paste(genus, species))
  } else {
    return("Unknown_Taxonomy")
  }
}

# Create detection frequency bar plots
plot_data_combined <- data.frame(
  Feature = rep(c(paste(get_tax_info(top_2_features[1]), "\n(", top_2_features[1], ")", sep=""),
                  paste(get_tax_info(top_2_features[2]), "\n(", top_2_features[2], ")", sep="")), each = 2),
  Group = rep(c("GT_absent", "GT_present"), 2),
  Detection_Frequency = c(
    detection_results$prev_GT_absent[1] * 100,
    detection_results$prev_GT_present[1] * 100,
    detection_results$prev_GT_absent[2] * 100,
    detection_results$prev_GT_present[2] * 100
  ),
  P_value = rep(c(
    format.pval(detection_results$fisher_p[1], digits = 3),
    format.pval(detection_results$fisher_p[2], digits = 3)
  ), each = 2)
)

# Create the plot
p <- ggplot(plot_data_combined, aes(x = Group, y = Detection_Frequency, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = paste0(round(Detection_Frequency, 1), "%")), 
            vjust = -0.3, size = 4, fontface = "bold") +
  facet_wrap(~ Feature, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("GT_absent" = "#3498DB", "GT_present" = "#E74C3C")) +
  labs(
    title = "Detection Frequency: Top 2 Features with Largest Differences",
    subtitle = "Percentage of samples where each feature is detected (abundance > 0)",
    x = "Ground Truth Group",
    y = "Detection Frequency (%)",
    fill = "Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylim(0, max(plot_data_combined$Detection_Frequency) * 1.1)

# Add p-values as strip labels
p <- p + geom_text(data = data.frame(
  Feature = unique(plot_data_combined$Feature),
  P_value = c(format.pval(detection_results$fisher_p[1], digits = 3),
             format.pval(detection_results$fisher_p[2], digits = 3)),
  x = 1.5, y = rep(max(plot_data_combined$Detection_Frequency) * 1.05, 2)
), aes(x = x, y = y, label = paste("p =", P_value)), 
inherit.aes = FALSE, size = 3.5, fontface = "italic")

# Save the plot
ggsave("../results/microb/detection_frequency_top2_features.png", p, width = 12, height = 8, dpi = 300)

cat("Plot saved as: results/microb/detection_frequency_top2_features.png\n")

# Save the results
write.csv(detection_results, "../results/microb/detection_frequency_all_features.csv", row.names = FALSE)
cat("Full results saved to: detection_frequency_all_features.csv\n")

    library(ggplot2)
# Create the plot
p <- ggplot(detection_results, aes(x = prev_GT_absent * 100, y = prev_GT_present * 100)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "ESV Prevalence: GT_present vs GT_absent",
       x = "Prevalence in GT_absent (%)",
       y = "Prevalence in GT_present (%)",
       subtitle = paste("Total ESVs analyzed:", nrow(detection_results))) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

# Save directly to file
ggsave("../results/microb/prevalence_comparison.png", p, width = 8, height = 6, dpi = 300)

cat("Prevalence comparison plot saved as: prevalence_comparison.png\n")

