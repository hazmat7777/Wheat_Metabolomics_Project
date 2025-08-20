# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

# Extract data from phyloseq object
sample_df <- data.frame(sample_data(ps_fs))
otu_df <- data.frame(t(as.matrix(otu_table(ps_fs))))

# Get the top discriminatory ESVs (from your previous analysis)
# Let's take top 30 features with biggest detection frequency differences
top_n <- 30
top_esvs <- head(detection_results$feature, top_n)

# Subset to top ESVs only
otu_subset <- otu_df[, top_esvs]

# Add sample metadata
heatmap_data <- cbind(sample_df[, "gt", drop = FALSE], otu_subset)

# Transform abundances - log10 transform with pseudocount for better visualization
# Add small pseudocount to avoid log(0)
otu_subset_log <- log10(otu_subset + 1e-6)

# Prepare data for heatmap
heatmap_data_log <- cbind(sample_df[, "gt", drop = FALSE], otu_subset_log)

# Order samples by group for better visualization
heatmap_data_log <- heatmap_data_log[order(heatmap_data_log$gt), ]

# Create sample names with group info
heatmap_data_log$sample_id <- paste0(rownames(heatmap_data_log), "_", heatmap_data_log$gt)

# Melt data for ggplot
heatmap_melted <- melt(heatmap_data_log, 
                       id.vars = c("gt", "sample_id"),
                       variable.name = "ESV", 
                       value.name = "log_abundance")

# Get taxonomy info for ESV labels
tax_df <- as.data.frame(tax_table(ps_fs))
get_tax_label <- function(esv_name) {
  if(esv_name %in% rownames(tax_df)) {
    genus <- tax_df[esv_name, "Genus"]
    species <- tax_df[esv_name, "Species"]
    
    genus <- ifelse(is.na(genus) || genus == "", "Unknown", genus)
    species <- ifelse(is.na(species) || species == "", "sp.", species)
    
    # Truncate long names
    if(nchar(paste(genus, species)) > 25) {
      return(paste0(substr(genus, 1, 10), "... (", substr(esv_name, 5, 10), ")"))
    } else {
      return(paste0(genus, " ", species, " (", substr(esv_name, 5, 10), ")"))
    }
  } else {
    return(esv_name)
  }
}

# Create ESV labels with taxonomy
esv_labels <- sapply(top_esvs, get_tax_label)
names(esv_labels) <- top_esvs

# Update factor levels for better ordering
heatmap_melted$ESV <- factor(heatmap_melted$ESV, levels = top_esvs)

# Create the heatmap
p1 <- ggplot(heatmap_melted, aes(x = ESV, y = sample_id, fill = log_abundance)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_viridis_c(name = "Log10\nAbundance", 
                       trans = "identity",
                       na.value = "grey90",
                       option = "plasma") +
  facet_grid(gt ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0, size = 12, face = "bold"),
    panel.spacing = unit(0.1, "cm"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = paste("Top", top_n, "Discriminatory ESVs: Abundance Heatmap"),
    subtitle = "Ordered by detection frequency difference between groups",
    x = "ESV (with taxonomy)",
    y = "Samples"
  ) +
  scale_x_discrete(labels = esv_labels)

# Alternative version: Presence/Absence heatmap (cleaner for sparse data)
# Convert to presence/absence (0/1)
otu_subset_binary <- otu_subset
otu_subset_binary[otu_subset_binary > 0] <- 1

heatmap_data_binary <- cbind(sample_df[, "gt", drop = FALSE], otu_subset_binary)
heatmap_data_binary <- heatmap_data_binary[order(heatmap_data_binary$gt), ]
heatmap_data_binary$sample_id <- paste0(rownames(heatmap_data_binary), "_", heatmap_data_binary$gt)

# Melt binary data
heatmap_melted_binary <- melt(heatmap_data_binary, 
                              id.vars = c("gt", "sample_id"),
                              variable.name = "ESV", 
                              value.name = "presence")

heatmap_melted_binary$ESV <- factor(heatmap_melted_binary$ESV, levels = top_esvs)
heatmap_melted_binary$presence <- factor(heatmap_melted_binary$presence, levels = c(0, 1))

# Create presence/absence heatmap
p2 <- ggplot(heatmap_melted_binary, aes(x = ESV, y = sample_id, fill = presence)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_manual(values = c("0" = "white", "1" = "darkblue"), 
                    name = "Detected",
                    labels = c("Absent", "Present")) +
  facet_grid(gt ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0, size = 12, face = "bold"),
    panel.spacing = unit(0.1, "cm"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = paste("Top", top_n, "Discriminatory ESVs: Presence/Absence Heatmap"),
    subtitle = "Blue = detected, White = not detected",
    x = "ESV (with taxonomy)",
    y = "Samples"
  ) +
  scale_x_discrete(labels = esv_labels)

# Save both versions
ggsave("abundance_heatmap_log.png", p1, width = 16, height = 10, dpi = 300)
ggsave("presence_absence_heatmap.png", p2, width = 16, height = 10, dpi = 300)

cat("Heatmaps saved as:\n")
cat("- abundance_heatmap_log.png (log10 transformed abundances)\n")
cat("- presence_absence_heatmap.png (binary presence/absence)\n")

# Print summary
cat("\nHeatmap Summary:\n")
cat("- Samples are grouped by GT status (present/absent)\n")
cat("- ESVs are ordered by detection frequency difference\n")
cat("- Each row represents one sample\n")
cat("- Each column represents one ESV\n")
cat("- Top panel: GT_absent samples\n")
cat("- Bottom panel: GT_present samples\n")

# Also create a summary table of the top ESVs
summary_table <- detection_results[1:top_n, c("feature", "prev_GT_present", "prev_GT_absent", "prev_diff", "fisher_p")]
summary_table$prev_GT_present_pct <- round(summary_table$prev_GT_present * 100, 1)
summary_table$prev_GT_absent_pct <- round(summary_table$prev_GT_absent * 100, 1)
summary_table$taxonomy <- sapply(summary_table$feature, get_tax_label)

write.csv(summary_table, "top_discriminatory_esvs_summary.csv", row.names = FALSE)
cat("- top_discriminatory_esvs_summary.csv (summary table)\n")