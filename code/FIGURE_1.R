# Load all required packages
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(terra)
library(sf)
library(units)
library(geodata)
library(raster)
library(tidyverse)
library(readxl)
library(readODS)

sf_use_s2(FALSE)

# ===============================
# LOAD DATA FOR ALL PANELS
# ===============================

# Load ordination data
ps_phylum <- readRDS("../data/metabarcoding/ps_16S_highdiv_genus_relative.rds")
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
fk_metabolom_gt_scaled$sample_name <- rownames(fk_metabolom_gt_scaled)

# Load mapping data
env_data <- read_xlsx("../data/env_data/sample_field_specific_environmental_data.xlsx")
env_data_filtered <- env_data %>%
  filter(!is.na(gps_latitude) & !is.na(gps_longitude)) %>%
  filter(project == "farmer_kit") %>%
  dplyr::select(gps_latitude, gps_longitude, sample_name)

# Add GT data for mapping
gt <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
gt$sample_name <- rownames(gt)

gt <- gt %>%
  dplyr::select(sample_name, gt) %>%
  mutate(gt = as.character(gt)) %>%
  filter(!is.na(gt)) %>%
  filter(gt != "")
colnames(gt) <- c("sample_name", "gt")

# Merge mapping data
env_data_filtered <- env_data_filtered %>%
  left_join(gt, by = "sample_name") %>%
  mutate(gt = as.factor(gt)) %>%
  dplyr::select(sample_name, gps_latitude, gps_longitude, gt) %>%
  filter(!is.na(gt)) %>%
  filter(gt != "")

# Convert to sf object
env_data_sf <- st_as_sf(env_data_filtered, coords = c("gps_longitude", "gps_latitude"), crs = 4326)

# Load UK map
ne_10 <- st_read('../../CMEECourseWork/week09/data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp')
st_agr(ne_10) <- 'constant'
ne_10_uk <- ne_10 %>%
  filter(ADMIN == "United Kingdom") %>%
  st_transform(crs = 4326) %>%
  dplyr::select(ADMIN) %>%
  st_crop(xmin = -4, ymin = 50, xmax = 2, ymax = 55)

# Add jitter to the sf object
env_data_sf_jittered <- env_data_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(
    X_jittered = X + runif(nrow(.), -0.1, 0.1),
    Y_jittered = Y + runif(nrow(.), -0.1, 0.1)
  ) %>%
  st_as_sf(coords = c("X_jittered", "Y_jittered"), crs = st_crs(env_data_sf)) %>%
  bind_cols(st_drop_geometry(env_data_sf))

# ===============================
# PANEL A: SAMPLE MAP
# ===============================

pA <- ggplot() +
  geom_sf(data = ne_10_uk, fill = "honeydew", color = "black") +
  geom_sf(data = env_data_sf_jittered, aes(color = gt), size = 3.5, shape = 17, alpha = 0.9) +
  labs(title = "A") +
  scale_color_manual(
    values = c("GT_absent" = "#E69F00", "GT_present" = "#56B4E9"),
    labels = c("GT_absent" = "GT Absent", "GT_present" = "GT Present"),
    name = "Disease Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.05, size = 28, face = "bold"),
    legend.position = "none"
  )
pA
# ===============================
# PANEL B: MICROBIAL ORDINATION BY PHYLUM
# ===============================

# Ordinate the data using NMDS
ps_phylum_ord <- ordinate(ps_phylum, "NMDS", "bray")

# Extract stress value for NMDS
nmds_stress <- round(ps_phylum_ord$stress, 3)

# Create microbial ordination plot
pB <- plot_ordination(ps_phylum, ps_phylum_ord, type = "sample", color = "gt", title = "B") +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    x = paste0("NMDS1"),
    y = paste0("NMDS2"),
    subtitle = paste0("Stress = ", nmds_stress),
    color = "Disease Status"
  ) +
  scale_color_manual(
    values = c("GT_absent" = "#E69F00", "GT_present" = "#56B4E9"),
    labels = c("GT_absent" = "GT Absent", "GT_present" = "GT Present"),
    name = "Disease Status"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.05, size = 28, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    legend.position = "none"
  ) +
  stat_ellipse(aes(group = factor(gt)), type = "norm", level = 0.68, alpha = 0.3)

# ===============================
# PANEL C: METABOLITE ORDINATION
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
pC <- ggplot(metabolite_ord_df, aes(x = PC1, y = PC2, color = factor(gt))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "C",
    x = paste0("PC1 (", metabolite_var_explained[1], "%)"),
    y = paste0("PC2 (", metabolite_var_explained[2], "%)"),
    color = "Disease Status"
  ) +
  scale_color_manual(
    values = c("GT_absent" = "#E69F00", "GT_present" = "#56B4E9"),
    labels = c("GT_absent" = "GT Absent", "GT_present" = "GT Present"),
    name = "Disease Status"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.05, size = 28, face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    legend.position = "bottom",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  ) +
  stat_ellipse(aes(group = factor(gt)), type = "norm", level = 0.68, alpha = 0.3)

# ===============================
# COMBINE ALL PANELS
# ===============================

# Extract legend from pC
legend <- get_legend(pC)

# Remove legends from all plots
pA_no_legend <- pA + theme(legend.position = "none")
pB_no_legend <- pB + theme(legend.position = "none")
pC_no_legend <- pC + theme(legend.position = "none")

# Create left panel (B and C stacked)
left_panel <- arrangeGrob(pB_no_legend, pC_no_legend, nrow = 2)

# Combine with map on right and legend at bottom
combined_figure <- grid.arrange(
  arrangeGrob(pA_no_legend, left_panel, ncol = 2, widths = c(1, 1)),
  legend,
  nrow = 2,
  heights = c(10, 1)
)

combined_figure

# Save the figure
ggsave("../figures_tables/figure1.png", 
       combined_figure, width = 16, height = 10, dpi = 600)

# Print the combined plot
combined_figure

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