library(phyloseq)
library(vegan)
library(dplyr)

ps_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
View(fk_metabolom_gt_scaled)
fk_metabolom_gt_scaled$sample_name <- rownames(fk_metabolom_gt_scaled)

# feature selection- only important genera.
genera_implicated <- c("Pseudomonas", "Bacillus") # removed "Bacillus"

tax_df <- as.data.frame(tax_table(ps_highdiv_absolute))
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]
ps_fs <- prune_taxa(taxa_to_keep, ps_highdiv_absolute)

ps_fs

ps_fs_div_df <- data.frame(estimate_richness(ps_fs, measures = c("Shannon", "Observed")))
ps_fs_div_df$sample_name <- rownames(ps_fs_div_df)
#View(ps_fs_div_df)


ps_fs_div_df$gt <- fk_metabolom_gt_scaled$gt[match(ps_fs_div_df$sample_name, fk_metabolom_gt_scaled$sample_name)]
ps_fs_div_df <- ps_fs_div_df %>%
    filter(!is.na(gt)) # remove rows with NA in gt
ps_fs_div_df

t_test_result <- t.test(Observed ~ gt, data = ps_fs_div_df)
cat("T-test result:\n")
print(t_test_result)

# make these into boxplots.

interesting_compounds <- c("Salicylate", "ISOLEUCINE", "Ornithine", "THREONINE", 
"ASPARTATE", "NICOTINAMIDE", "TRYPTOPHAN", "VALINE", "succinic acid", 
"LACTIC ACID", "Linalool", "Geranyl acetate", "trans-Nerolidol", "Nerol", 
"FERULATE", "CAFFEATE", "caffeine", "Cinnamic acid - 40.0 eV", 
"Vanillic acid", "BETAINE", "CHOLINE", "SARCOSINE", "N,N-Dimethylglycine", 
"3-HYDROXYBENZOATE", "4-HYDROXYBENZOATE", "Sinapic acid", 
"Hydroquinone", "TYROSINE", "trans-3-Hydroxycinnamic acid", 
"3,4-DIHYDROXYBENZOATE", "COUMARIN", "L-Arginine", "L-Proline", 
"PHENYLALANINE",   "Salicylate",                    # [1] - keratolytic antifungal
  "berberine",                     # [630] - broad spectrum antifungal
  "VANILLIN",                     # [339] - demonstrated antifungal activity
  "caffeine",                     # [485] - direct antifungal vs Candida
  "COUMARIN",                     # [674] - antifungal vs dermatophytes
  "TAURINE",                      # [619] - antifungal properties
  "NICOTINAMIDE",                 # [81] - B-vitamin with antifungal activity
  "THYMIDINE",                    # [214] - nucleoside with antifungal properties
  "QUINOLINE",                    # [498] - scaffold in antifungal medications
  "Usnic.acid",                   # [573] - lichen metabolite, antifungal
  "Chrysin",                      # [477] - flavonoid with antifungal activity
  "CAFFEATE",                     # [426] - phenolic with antifungal properties
  "trans.3.Hydroxycinnamic.acid"  # [267] - related to caffeic acid
  )

# load data
fk_metabolom_gt_t <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
View(fk_metabolom_gt_t)
# save it as a csv
write.csv(fk_metabolom_gt_t, "../data/metabolomics/fk_metabolomics_gt_noNA.csv", row.names = FALSE)


# Check which antifungal metabolites are present in your dataset
available_antifungal <- intersect(interesting_compounds, colnames(fk_metabolom_gt_t))
available_antifungal

# ===============================
# DATA PREPARATION WITH ANTIFUNGAL SUBSET
# ===============================
head(colnames(fk_metabolom_gt_t))
# Create subset with only antifungal metabolites + target variable
antifungal_data <- fk_metabolom_gt_t[, c(available_antifungal, "gt")]

# Scale the predictors
predictor_cols <- names(antifungal_data)[names(antifungal_data) != "gt"]
antifungal_data[predictor_cols] <- lapply(antifungal_data[predictor_cols], function(x) {
as.numeric(as.character(x))
})

##################################################


# CHANGE THE LEVELS
look^ at that

library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)

ps_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")
fk_metabolom_gt_scaled <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")
fk_metabolom_gt_scaled$sample_name <- rownames(fk_metabolom_gt_scaled)

# ===============================
# PANEL 1: UNFILTERED MICROBIAL DIVERSITY
# ===============================
# Calculate diversity for all taxa (unfiltered)
ps_all_div_df <- data.frame(estimate_richness(ps_highdiv_absolute, measures = c("Shannon", "Observed")))
ps_all_div_df$sample_name <- rownames(ps_all_div_df)
ps_all_div_df$gt <- fk_metabolom_gt_scaled$gt[match(ps_all_div_df$sample_name, fk_metabolom_gt_scaled$sample_name)]
ps_all_div_df <- ps_all_div_df %>%
    filter(!is.na(gt))

p1 <- ggplot(ps_all_div_df, aes(x = factor(gt), y = Observed, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(
    x = "GT Status",
    y = "Observed Richness (All Taxa)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# ===============================
# PANEL 2: UNFILTERED METABOLOMIC DIVERSITY
# ===============================
# Calculate total metabolite abundance (unfiltered - all metabolites)
metabolite_cols <- names(fk_metabolom_gt_scaled)[names(fk_metabolom_gt_scaled) != "gt" & names(fk_metabolom_gt_scaled) != "sample_name"]

unfiltered_metabolite_diversity <- fk_metabolom_gt_scaled %>%
  rowwise() %>%
  mutate(
    total_metabolite_abundance = sum(c_across(all_of(metabolite_cols)), na.rm = TRUE)
  ) %>%
  select(gt, total_metabolite_abundance) %>%
  filter(!is.na(gt))

p2 <- ggplot(unfiltered_metabolite_diversity, aes(x = factor(gt), y = total_metabolite_abundance, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(
    x = "GT Status",
    y = "Total Metabolite Abundance"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# ===============================
# PANEL 3: FILTERED METABOLOMIC DIVERSITY (Your existing code)
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
  select(gt, metabolite_abundance)

p3 <- ggplot(metabolite_diversity, aes(x = factor(gt), y = metabolite_abundance, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(
    x = "GT Status",
    y = "Feature-selected Metabolite Abundance"
  ) +
#   scale_fill_manual(values = c("0" = "#E69F00", "1" = "#56B4E9")) +
#   scale_x_discrete(labels = c("0" = "GT_absent", "1" = "GT_present")) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# ===============================
# PANEL 4: FILTERED MICROBIAL DIVERSITY (Your existing code)
# ===============================
genera_implicated <- c("Pseudomonas", "Bacillus")
tax_df <- as.data.frame(tax_table(ps_highdiv_absolute))
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]
ps_fs <- prune_taxa(taxa_to_keep, ps_highdiv_absolute)

ps_fs_div_df <- data.frame(estimate_richness(ps_fs, measures = c("Shannon", "Observed")))
ps_fs_div_df$sample_name <- rownames(ps_fs_div_df)
ps_fs_div_df$gt <- fk_metabolom_gt_scaled$gt[match(ps_fs_div_df$sample_name, fk_metabolom_gt_scaled$sample_name)]
ps_fs_div_df <- ps_fs_div_df %>%
    filter(!is.na(gt))

p4 <- ggplot(ps_fs_div_df, aes(x = factor(gt), y = Observed, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(
    x = "GT Status",
    y = "Observed Pseudomonas and Bacillus Richness"
  ) +
#   scale_fill_manual(values = c("0" = "#E69F00", "1" = "#56B4E9")) +
#   scale_x_discrete(labels = c("0" = "GT_absent", "1" = "GT_present")) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# ===============================
# COMBINE ALL FOUR PANELS
# ===============================
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Save the figure
ggsave("../results/plots/four_panel_unfiltered_filtered_gt.png", combined_plot, 
       width = 12, height = 10, dpi = 300)

# Print the combined plot
combined_plot
library(ggplot2)
# library(gridExtra)



