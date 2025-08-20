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

# Load the metabolomics data for linear model
fk_metabolom_gt <- readRDS("../data/fk_metabolomics_gt.RDS")

## Prepare data frame
fk_metabolom_gt_t <- t(fk_metabolom_gt) # Transpose to make metabolites as columns
colnames(fk_metabolom_gt_t) <- fk_metabolom_gt_t[1, ] # Set the first row as column names
fk_metabolom_gt_t <- fk_metabolom_gt_t[-1, ] # Remove the first row (now column names)

fk_metabolom_gt_t <- as.data.frame(fk_metabolom_gt_t) # Convert to data frame
colnames(fk_metabolom_gt_t)[ncol(fk_metabolom_gt_t)] <- "gt" # Rename the last column to 'gt'
fk_metabolom_gt_t$gt <- trimws(as.character(fk_metabolom_gt_t$gt)) # Remove any extra whitespace from gt col 
fk_metabolom_gt_t$gt <- as.factor(fk_metabolom_gt_t$gt) # Convert 'gt' to a factor
levels(fk_metabolom_gt_t$gt) <- c("GT_absent", "GT_present") # Set factor levels

# make sure all metabolomics data is numeric
fk_metabolom_gt_t[] <- lapply(fk_metabolom_gt_t, function(x) {as.numeric(x)} )

# get a measure of diversity in the metabolomics data
metab_diversity <- diversity(fk_metabolom_gt_t[, -ncol(fk_metabolom_gt_t)], index = "shannon")
met_div_df <- as.data.frame(metab_diversity)
met_div_df$sample_name <- rownames(met_div_df)

# Load the phyloseq object
ps_16S <- readRDS(file = "../data/metabarcoding/ps_16S_highdiv_absolute.rds")

# Get the sample names from the phyloseq object
sample_df <- as.data.frame(sample_data(ps_16S))
sample_names_16S<- sample_df$sample_name
og_sample_names <- sample_names_16S

# select last 3 characters of sample names
sample_names_16S <- substr(sample_names_16S, nchar(sample_names_16S) - 2, nchar(sample_names_16S))
# trim whitespace
sample_names_16S <- trimws(sample_names_16S)
# remove the "." prefix where it exists
sample_names_16S <- gsub("^\\.", "", sample_names_16S)
# remove the "K" prefix where it exists
sample_names_16S <- gsub("K", "", sample_names_16S)
# add "fk" prefix to match metabolomics data
sample_names_16S <- paste0("fk", sample_names_16S)
sample_df$sample_name<- sample_names_16S

sample_names_16S[6] <- "fk161" # manually change one which went wrong
sample_names_16S[12] <- "fk140" # manually change one which went wrong
sample_names_16S[75] <- "fk88" # manually change one which went wrong

# add changes
sample_df$sample_name<- sample_names_16S

# remove the duplicate sample name
sample_to_remove<- rownames(sample_df[sample_df$sample_name == "fk88", ][2]) # check the sample names
ps_16S <- prune_samples(sample_names(ps_16S) != sample_to_remove, ps_16S)
sample_df <- sample_df[rownames(sample_df) != sample_to_remove, ] # remove the sample from the sample data frame

# update the sample names in the phyloseq object
sample_names(ps_16S) <- sample_df$sample_name

## make a table of diversity metrics
# get diversity metrics for the 16S data
seq_counts <- sample_sums(ps_16S)

div_16S_df <- data.frame(estimate_richness(ps_16S, measures = c("Shannon", "Observed")), seq_count = seq_counts)
div_16S_df$sample_name <- rownames(div_16S_df)

# Merge metabolomics and 16S data frames by sample name
merged_div_df <- left_join(met_div_df, div_16S_df, by = "sample_name")
merged_div_df <- merged_div_df %>%
    filter(!is.na(Shannon) & !is.na(metab_diversity)) # remove rows with NA in either diversity metric

# add gt column from metabolomics data
merged_div_df$gt <- fk_metabolom_gt_scaled$gt[match(merged_div_df$sample_name, fk_metabolom_gt_scaled$sample_name)]
merged_div_df <- merged_div_df %>%
    filter(!is.na(gt)) # remove rows with NA in gt

# is there are relationship between diversity and gt?
merged_div_df$gt <- as.factor(merged_div_df$gt) # convert gt to factor
levels(merged_div_df$gt) <- c("GT_absent", "GT_present") # set factor levels  

# ===============================
# PANEL 1: LINEAR MODEL (BIGGEST)
# ===============================

# Plot the diversity metrics
p_linear <- ggplot(merged_div_df, aes(x = Shannon, y = metab_diversity)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(y = "Metabolomics Diversity (Shannon)",
       x = "16S Diversity (Shannon)") +
  theme_bw()

# run a simple linear regression
lm_result <- lm(metab_diversity ~ Shannon, data = merged_div_df)

# ===============================
# PANEL 2: UNFILTERED MICROBIAL DIVERSITY
# ===============================
# Calculate diversity for all taxa (unfiltered) - use the processed data from merged_div_df
ps_all_div_df <- merged_div_df %>%
  select(sample_name, Shannon, Observed, gt) %>%
  rename(Shannon_16S = Shannon, Observed_16S = Observed)

p1 <- ggplot(ps_all_div_df, aes(x = factor(gt), y = Observed_16S, fill = factor(gt))) +
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
# PANEL 3: UNFILTERED METABOLOMIC DIVERSITY
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
  select(gt, metabolite_abundance)

p3 <- ggplot(metabolite_diversity, aes(x = factor(gt), y = metabolite_abundance, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(
    x = "GT Status",
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
# For the filtered microbial diversity, we need to calculate it from the processed phyloseq
# but use only samples that are in our merged_div_df

# Get sample names that are in our processed dataset
valid_samples <- merged_div_df$sample_name

# Load the processed phyloseq object (this should match the sample naming)
ps_16S <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")

# Apply the same sample name processing as in the original script
sample_df <- as.data.frame(sample_data(ps_16S))
sample_names_16S <- sample_df$sample_name

# Apply the same transformations
sample_names_16S <- substr(sample_names_16S, nchar(sample_names_16S) - 2, nchar(sample_names_16S))
sample_names_16S <- trimws(sample_names_16S)
sample_names_16S <- gsub("^\\.", "", sample_names_16S)
sample_names_16S <- gsub("K", "", sample_names_16S)
sample_names_16S <- paste0("fk", sample_names_16S)

sample_names_16S[6] <- "fk161"
sample_names_16S[12] <- "fk140"
sample_names_16S[75] <- "fk88"

sample_df$sample_name <- sample_names_16S

# Remove duplicate
sample_to_remove <- rownames(sample_df[sample_df$sample_name == "fk88", ][2])
ps_16S <- prune_samples(sample_names(ps_16S) != sample_to_remove, ps_16S)
sample_df <- sample_df[rownames(sample_df) != sample_to_remove, ]

# Update sample names
sample_names(ps_16S) <- sample_df$sample_name

# Now filter for genera
genera_implicated <- c("Pseudomonas", "Bacillus")
tax_df <- as.data.frame(tax_table(ps_16S))
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% genera_implicated]
ps_fs <- prune_taxa(taxa_to_keep, ps_16S)

# Calculate diversity and filter to only samples in our merged dataset
ps_fs_div_df <- data.frame(estimate_richness(ps_fs, measures = c("Shannon", "Observed")))
ps_fs_div_df$sample_name <- rownames(ps_fs_div_df)

# Filter to only samples that are in our merged dataset and add GT status
ps_fs_div_df <- ps_fs_div_df %>%
  filter(sample_name %in% valid_samples) %>%
  left_join(merged_div_df %>% select(sample_name, gt), by = "sample_name") %>%
  filter(!is.na(gt))

p4 <- ggplot(ps_fs_div_df, aes(x = factor(gt), y = Observed, fill = factor(gt))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(
    x = "GT Status",
    y = "Observed Pseudomonas and Bacillus Richness"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# ===============================
# COMBINE ALL FIVE PANELS
# ===============================

# Create the bottom four panel grid
bottom_four_panels <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Combine the linear model (top, bigger) with the four panels (bottom)
combined_five_panel <- grid.arrange(
  p_linear,
  bottom_four_panels,
  nrow = 2,
  heights = c(2, 3)  # Make linear model panel bigger
)

# Save the figure
ggsave("../results/plots/five_panel_combined_figure.png", combined_five_panel, 
       width = 14, height = 16, dpi = 300)

# Print the combined plot
combined_five_panel

# ===============================
# STATISTICAL RESULTS
# ===============================

cat("Linear regression results:\n")
print(summary(lm_result))