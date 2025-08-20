# Load all required libraries
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tibble)
library(readxl)
library(readODS)

# ===============================
# DATA LOADING AND PREPARATION
# ===============================

ps_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
fk_metabolom_gt_scaled$sample_name <- rownames(fk_metabolom_gt_scaled)

# ===============================
# LINEAR MODEL DATA PREPARATION
# ===============================




# load it 
merged_div_df <- readRDS("../data/merged_metab_microb_diversity_data.rds")


# ===============================
# PANEL 1: LINEAR MODEL (BIGGEST)
# ===============================

# Plot the diversity metrics
p_linear <- ggplot(merged_div_df, aes(x = Shannon, y = metab_diversity)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title ="A",
        y = "Metabolomics Diversity (Shannon)",
        x = "16S Diversity (Shannon)") +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# run a simple linear regression
lm_result <- lm(metab_diversity ~ Shannon, data = merged_div_df)

p_linear

# ===============================
# PANEL 2: UNFILTERED MICROBIAL DIVERSITY
# ===============================
# Calculate diversity for all taxa (unfiltered)

View(merged_div_df)

p1 <- ggplot(merged_div_df, aes(x = factor(gt), y = Shannon, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(title = "B",
    x = "",
    y = "Shannon Diversity (All Taxa)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )


# ===============================
# PANEL 3: UNFILTERED METABOLOMIC DIVERSITY
# ===============================
# Calculate total metabolite abundance (unfiltered - all metabolites)
metabolite_cols <- names(fk_metabolom_gt_scaled)[names(fk_metabolom_gt_scaled) != "gt" & names(fk_metabolom_gt_scaled) != "sample_name"]

unfiltered_metabolite_diversity <- fk_metabolom_gt_scaled %>%
  rowwise() %>%
  mutate(
    total_metabolite_abundance = sum(c_across(all_of(metabolite_cols)), na.rm = TRUE)
  ) %>%
  dplyr::select(gt, total_metabolite_abundance) %>%
  filter(!is.na(gt))

p2 <- ggplot(unfiltered_metabolite_diversity, aes(x = factor(gt), y = total_metabolite_abundance, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(title = "D",
    x = "",
    y = "Total Metabolite Abundance"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# ===============================
# PANEL 4: FILTERED METABOLOMIC DIVERSITY
# ===============================
interesting_compounds <- c("Salicylate", "ISOLEUCINE", "Ornithine", "THREONINE", 
"ASPARTATE", "NICOTINAMIDE", "TRYPTOPHAN", "VALINE", "succinic acid", 
"LACTIC ACID", "Linalool", "Geranyl acetate", "trans-Nerolidol", "Nerol", 
"FERULATE", "CAFFEATE", "caffeine", "Cinnamic acid - 40.0 eV", 
"Vanillic acid", "BETAINE", "CHOLINE", "SARCOSINE", "N,N-Dimethylglycine", 
"3-HYDROXYBENZOATE", "4-HYDROXYBENZOATE", "Sinapic acid", 
"Hydroquinone", "TYROSINE", "trans-3-Hydroxycinnamic acid", 
"3,4-DIHYDROXYBENZOATE", "COUMARIN", "L-Arginine", "L-Proline", 
"PHENYLALANINE", "berberine", "VANILLIN", "TAURINE", "THYMIDINE", 
"QUINOLINE", "Usnic.acid", "Chrysin", "trans.3.Hydroxycinnamic.acid")

fk_metabolom_gt_t <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
available_antifungal <- intersect(interesting_compounds, colnames(fk_metabolom_gt_t))
antifungal_data <- fk_metabolom_gt_t[, c(available_antifungal, "gt")]

predictor_cols <- names(antifungal_data)[names(antifungal_data) != "gt"]
antifungal_data[predictor_cols] <- lapply(antifungal_data[predictor_cols], function(x) {
as.numeric(as.character(x))
})

metabolite_diversity <- antifungal_data %>%
  rowwise() %>%
  mutate(
    metabolite_abundance = sum(c_across(all_of(available_antifungal)), na.rm = TRUE)
  ) %>%
  dplyr::select(gt, metabolite_abundance)

p3 <- ggplot(metabolite_diversity, aes(x = factor(gt), y = metabolite_abundance, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(title = "E",
    x = "",
    y = "Feature-selected Metabolite Abundance"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# ===============================
# PANEL 5: FILTERED MICROBIAL DIVERSITY
# ===============================
genera_implicated <- c("Pseudomonas", "Bacillus")
ps_16S_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")
sd <- as.data.frame(sample_data(ps_highdiv_absolute))
sample_names(ps_highdiv_absolute) <- sd$sample_name

tax_df <- as.data.frame(tax_table(ps_highdiv_absolute))
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]
ps_fs <- prune_taxa(taxa_to_keep, ps_highdiv_absolute)

ps_fs_div_df <- data.frame(estimate_richness(ps_fs, measures = c("Shannon", "Observed")))
ps_fs_div_df$sample_name <- rownames(ps_fs_div_df)
ps_fs_div_df$gt <- fk_metabolom_gt_scaled$gt[match(ps_fs_div_df$sample_name, fk_metabolom_gt_scaled$sample_name)]
ps_fs_div_df <- ps_fs_div_df %>%
    dplyr::filter(!is.na(gt))

View(ps_fs_div_df)

p4 <- ggplot(ps_fs_div_df, aes(x = factor(gt), y = Observed, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(title = "C",
    x = "",
    y = "Pseudomonas and Bacillus Richness"
  ) +
#   scale_fill_manual(values = c("0" = "#E69F00", "1" = "#56B4E9")) +
#   scale_x_discrete(labels = c("0" = "GT_absent", "1" = "GT_present")) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# p4

# ##
# genera_implicated <- c("Pseudomonas", "Bacillus")
# tax_df <- as.data.frame(tax_table(ps_highdiv_absolute))
# taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]
# ps_fs <- prune_taxa(taxa_to_keep, ps_highdiv_absolute)

# ps_fs_div_df <- data.frame(estimate_richness(ps_fs, measures = c("Shannon", "Observed")))
# ps_fs_div_df$sample_name <- rownames(ps_fs_div_df)
# ps_fs_div_df$gt <- fk_metabolom_gt_scaled$gt[match(ps_fs_div_df$sample_name, fk_metabolom_gt_scaled$sample_name)]
# ps_fs_div_df <- ps_fs_div_df %>%
#     filter(!is.na(gt))

# p4 <- ggplot(ps_fs_div_df, aes(x = factor(gt), y = Observed, fill = factor(gt))) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
#   labs(
#     x = "GT Status",
#     y = "Observed Pseudomonas and Bacillus Richness"
#   ) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
#   )

# ===============================
# COMBINE ALL FIVE PANELS
# ===============================

# # Create the bottom four panel grid
# bottom_four_panels <- grid.arrange(p1, p4, p2, p3, ncol = 2, nrow = 2)
# # Create the bottom four panel grid with shared x-axis label
# bottom_four_panels <- arrangeGrob(p1, p4, p2, p3, ncol = 2, nrow = 2,
#                                  bottom = textGrob("GT Status", gp = gpar(fontsize = 12)))

bottom_four_panels <- grid.arrange(p1, p4, p2, p3, ncol = 2, nrow = 2,
                                  bottom = textGrob("GT Status", gp = gpar(fontsize = 12)))


# Combine the linear model (top, bigger) with the four panels (bottom)
combined_five_panel <- grid.arrange(
  p_linear,
  bottom_four_panels,
  nrow = 2,
  heights = c(2, 3)  # Make linear model panel bigger
)

# Save the figure
ggsave("../results/plots/five_panel_diversity_figure.png", combined_five_panel, 
       width = 10, height = 15, dpi = 600)

# Print the combined plot
combined_five_panel

# ===============================
# STATISTICAL RESULTS
# ===============================

cat("Linear regression results:\n")
print(summary(lm_result))