library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Load data
ps_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
fk_metabolom_gt_scaled$sample_name <- rownames(fk_metabolom_gt_scaled)

# ===============================
# PANEL A: MICROBIAL ORDINATION BY PHYLUM
# ===============================

# # Agglomerate taxa to phylum level
# ps_phylum <- tax_glom(ps_highdiv_absolute, taxrank = "Phylum")

# # Transform to relative abundance
# ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

# # Get OTU table for ordination
# otu_matrix <- as.matrix(otu_table(ps_phylum))
# if(taxa_are_rows(ps_phylum)) {
#   otu_matrix <- t(otu_matrix)
# }

# # Perform PCoA with Bray-Curtis distance
# microbe_dist <- vegdist(otu_matrix, method = "bray")
# microbe_pcoa <- cmdscale(microbe_dist, eig = TRUE)

# # Calculate percentage of variance explained
# microbe_var_explained <- round(microbe_pcoa$eig[1:2] / sum(microbe_pcoa$eig) * 100, 1)

# # Create dataframe for plotting
# microbe_ord_df <- data.frame(
#   PC1 = microbe_pcoa$points[, 1],
#   PC2 = microbe_pcoa$points[, 2],
#   gt = sample_data_microbe$gt,
#   sample_name = rownames(microbe_pcoa$points)
# )


# # Create microbial ordination plot
# p1 <- ggplot(microbe_ord_df, aes(x = PC1, y = PC2, color = factor(gt))) +
#   geom_point(size = 3, alpha = 0.7) +
#   labs(
#     title = "A",
#     x = paste0("PC1 (", microbe_var_explained[1], "%)"),
#     y = paste0("PC2 (", microbe_var_explained[2], "%)"),
#     color = "Disease Status"
#   ) +
#   scale_color_manual(
#     values = c("GT_absent" = "#E69F00", "GT_present" = "#56B4E9"),
#     labels = c("GT_absent" = "GT Absent", "GT_present" = "GT Present")
#   ) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.05, size = 14, face = "bold"),
#     legend.position = "bottom",
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 9)
#   ) +
#   stat_ellipse(aes(group = factor(gt)), type = "norm", level = 0.68, alpha = 0.3)

ps_phylum <- tax_glom(ps_16S_highdiv_absolute, taxrank = "Phylum")

# ordinate the data using NMDS
ps_phylum_ord <- ordinate(ps_phylum, "NMDS", "bray")
p_phylum <- plot_ordination(ps_phylum, ps_phylum_ord, type = "sample", color = "gt", title = "A")

plot(p_phylum)

# ===============================
# PANEL B: METABOLITE ORDINATION
# ===============================

# Prepare metabolite data
metabolite_cols <- names(fk_metabolom_gt_scaled)[names(fk_metabolom_gt_scaled) != "gt" & 
                                                names(fk_metabolom_gt_scaled) != "sample_name"]

# Create metabolite matrix (samples as rows, metabolites as columns)
metabolite_data <- fk_metabolom_gt_scaled %>%
  filter(!is.na(gt)) %>%
  dplyr::select(all_of(metabolite_cols))

# Convert to numeric matrix
metabolite_matrix <- as.matrix(metabolite_data)
metabolite_matrix[is.na(metabolite_matrix)] <- 0

# Perform PCA on metabolite data
metabolite_pca <- prcomp(metabolite_matrix, scale. = TRUE, center = TRUE)

# Calculate percentage of variance explained
metabolite_var_explained <- round(summary(metabolite_pca)$importance[2, 1:2] * 100, 1)

# Create dataframe for plotting
metabolite_ord_df <- data.frame(
  PC1 = metabolite_pca$x[, 1],
  PC2 = metabolite_pca$x[, 2],
  gt = fk_metabolom_gt_scaled$gt[!is.na(fk_metabolom_gt_scaled$gt)],
  sample_name = rownames(metabolite_matrix)
)

# Create metabolite ordination plot
p2 <- ggplot(metabolite_ord_df, aes(x = PC1, y = PC2, color = factor(gt))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "B",
    x = paste0("PC1 (", metabolite_var_explained[1], "%)"),
    y = paste0("PC2 (", metabolite_var_explained[2], "%)"),
    color = "Disease Status"
  ) +
  scale_color_manual(
    values = c("GT_absent" = "#E69F00", "GT_present" = "#56B4E9"),
    labels = c("GT_absent" = "GT Absent", "GT_present" = "GT Present")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.05, size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  stat_ellipse(aes(group = factor(gt)), type = "norm", level = 0.68, alpha = 0.3)

# ===============================
# COMBINE PANELS
# ===============================

# Combine both panels
combined_ordination <- grid.arrange(p1, p2, ncol = 2)

# Save the figure
ggsave("../results/plots/ordination_microbe_metabolite_composition.png", 
       combined_ordination, width = 14, height = 6, dpi = 300)

# Print the combined plot
combined_ordination

# ===============================
# OPTIONAL: STATISTICAL TESTS
# ===============================

# PERMANOVA for microbial communities
microbe_permanova <- adonis2(microbe_dist ~ gt, data = sample_data_microbe, permutations = 999)
cat("Microbial community PERMANOVA results:\n")
print(microbe_permanova)

# PERMANOVA for metabolite communities
metabolite_dist <- vegdist(metabolite_matrix, method = "euclidean")
metabolite_permanova <- adonis2(metabolite_dist ~ gt, data = metabolite_ord_df, permutations = 999)
cat("\nMetabolite composition PERMANOVA results:\n")
print(metabolite_permanova)