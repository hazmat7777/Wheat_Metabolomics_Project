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

# update the sample names in the phyloseq object
sample_df <- as.data.frame(sample_data(ps_16S))
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
 dplyr::select(sample_name, Shannon, Observed, gt) %>%
 rename(Shannon_16S = Shannon, Observed_16S = Observed)

p1 <- ggplot(ps_all_div_df, aes(x = factor(gt), y = Observed_16S, fill = factor(gt))) +
 geom_boxplot(alpha = 0.7, outlier.shape = NA) +
 geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
 labs(
 x = NULL,  # Remove individual x-axis label
 y = "Observed Richness (All Taxa)"
 ) +
 theme_classic() +
 theme(
 legend.position = "none",
 plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
 axis.title.x = element_blank()  # Ensure x-axis title is removed
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
 labs(
 x = NULL,  # Remove individual x-axis label
 y = "Total Metabolite Abundance"
 ) +
 theme_classic() +
 theme(
 legend.position = "none",
 plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
 axis.title.x = element_blank()  # Ensure x-axis title is removed
 )

# ===============================
# COMBINE ALL THREE PANELS WITH SHARED X-AXIS LABEL
# ===============================
# Create the bottom panel grid with shared x-axis label
bottom_TWO_panels <- arrangeGrob(
  p1, p2, 
  ncol = 2, 
  nrow = 1,
  bottom = textGrob("GT Status", gp = gpar(fontsize = 12))  # Add shared x-axis label
)

# Combine the linear model (top, bigger) with the two panels (bottom)
combined_three_panel <- grid.arrange(
 p_linear,
 bottom_TWO_panels,
 nrow = 2,
 heights = c(2, 3) # Make linear model panel bigger
)

# Save the figure
ggsave("../results/plots/diversity_3panel_combined_figure2.png", combined_three_panel,
 width = 7, height = 7, dpi = 300)

# Print the combined plot
plot(combined_three_panel)
# ===============================
# STATISTICAL RESULTS
# ===============================

cat("Linear regression results:\n")
print(summary(lm_result))