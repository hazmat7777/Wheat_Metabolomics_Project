library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

# Load data
ps_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
fk_metabolom_gt_scaled$sample_name <- rownames(fk_metabolom_gt_scaled)

# ===============================
# PANEL A: MICROBIAL ORDINATION BY PHYLUM
# ===============================

# Agglomerate taxa to phylum level
ps_phylum <- tax_glom(ps_highdiv_absolute, taxrank = "Phylum")

# Ordinate the data using NMDS
ps_phylum_ord <- ordinate(ps_phylum, "NMDS", "bray")

# Extract stress value for NMDS
nmds_stress <- round(ps_phylum_ord$stress, 3)

# Create microbial ordination plot
p1 <- plot_ordination(ps_phylum, ps_phylum_ord, type = "sample", color = "gt", title = "A") +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    x = paste0("NMDS1"),
    y = paste0("NMDS2"),
    subtitle = paste0("Stress = ", nmds_stress),
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

# Perform NMDS on metabolite data
metabolite_nmds <- metaMDS(metabolite_matrix, distance = "euclidean", k = 2, trymax = 100)

# Extract stress value for NMDS
metabolite_stress <- round(metabolite_nmds$stress, 3)

# Create dataframe for plotting
metabolite_ord_df <- data.frame(
  NMDS1 = metabolite_nmds$points[, 1],
  NMDS2 = metabolite_nmds$points[, 2],
  gt = fk_metabolom_gt_scaled$gt[!is.na(fk_metabolom_gt_scaled$gt)],
  sample_name = rownames(metabolite_matrix)
)

# Create metabolite ordination plot
p2 <- ggplot(metabolite_ord_df, aes(x = NMDS1, y = NMDS2, color = factor(gt))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "B",
    x = "NMDS1",
    y = "NMDS2",
    subtitle = paste0("Stress = ", metabolite_stress),
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

# Extract legend from p2
legend <- get_legend(p2)

# Remove legend from both plots
p1_no_legend <- p1 + theme(legend.position = "none")
p2_no_legend <- p2 + theme(legend.position = "none")

# Arrange plots and legend
combined_ordination <- grid.arrange(
  arrangeGrob(p1_no_legend, p2_no_legend, ncol = 2),
  legend,
  nrow = 2,
  heights = c(10, 1)
)

# Save the figure
ggsave("../results/plots/ordination_microbe_metabolite_composition3.png", 
       combined_ordination, width = 14, height = 6, dpi = 300)

# Print the combined plot
combined_ordination

# ===============================
# OPTIONAL: STATISTICAL TESTS
# ===============================

# PERMANOVA for microbial communities (using phyloseq object)
microbe_dist <- phyloseq::distance(ps_phylum, method = "bray")
sample_df <- data.frame(sample_data(ps_phylum))
microbe_permanova <- adonis2(microbe_dist ~ gt, data = sample_df, permutations = 999)
cat("Microbial community PERMANOVA results:\n")
print(microbe_permanova)

# PERMANOVA for metabolite communities
metabolite_dist <- vegdist(metabolite_matrix, method = "euclidean")
metabolite_permanova <- adonis2(metabolite_dist ~ gt, data = metabolite_ord_df, permutations = 999)
cat("\nMetabolite composition PERMANOVA results:\n")
print(metabolite_permanova)