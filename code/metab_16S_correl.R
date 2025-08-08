# comparing metabolomics and 16S data
# Load necessary libraries
library(phyloseq)
library(dplyr)
library(tibble)
library(readxl)
library(readODS)
library(vegan)

#=================================
# METABOLOMICS DATA
#=================================

# Load the metabolomics data
fk_metabolom_gt <- readRDS("../data/fk_metabolomics_gt.RDS")

View(fk_metabolom_gt)


## Prepare data frame
fk_metabolom_gt_t <- t(fk_metabolom_gt) # Transpose to make metabolites as columns
colnames(fk_metabolom_gt_t) <- fk_metabolom_gt_t[1, ] # Set the first row as column names
fk_metabolom_gt_t <- fk_metabolom_gt_t[-1, ] # Remove the first row (now column names)

fk_metabolom_gt_t <- as.data.frame(fk_metabolom_gt_t) # Convert to data frame
colnames(fk_metabolom_gt_t)[ncol(fk_metabolom_gt_t)] <- "gt" # Rename the last column to 'gt'
fk_metabolom_gt_t$gt <- trimws(as.character(fk_metabolom_gt_t$gt)) # Remove any extra whitespace from gt col 
fk_metabolom_gt_t$gt <- as.factor(fk_metabolom_gt_t$gt) # Convert 'gt' to a factor
levels(fk_metabolom_gt_t$gt) <- c("GT_absent", "GT_present") # Set factor levels

View(fk_metabolom_gt_t)
str(fk_metabolom_gt_t)
# make sure all metabolomics data is numeric
fk_metabolom_gt_t[] <- lapply(fk_metabolom_gt_t, function(x) {as.numeric(x)} )

# get a measure of diversity in the metabolomics data
metab_diversity <- diversity(fk_metabolom_gt_t[, -ncol(fk_metabolom_gt_t)], index = "shannon")
met_div_df <- as.data.frame(metab_diversity)
met_div_df$sample_name <- rownames(met_div_df)

#==================================
# 16S DATA
#==================================

# Load the phyloseq object
ps_16S <- readRDS(file = "../data/metabarcoding/ps_16S_filtered_allsamples.rds")
ps_16S
# View(otu_table(ps_16S))
# view(sample_data(ps_16S))

# Get the sample names from the phyloseq object
sample_df <- as.data.frame(sample_data(ps_16S))
sample_names_16S<- sample_df$sample_name
og_sample_names <- sample_names_16S
str(og_sample_names)
# select last 3 characters of sample names
sample_names_16S <- substr(sample_names_16S, nchar(sample_names_16S) - 2, nchar(sample_names_16S))
sample_names_16S
# trim whitespace
sample_names_16S <- trimws(sample_names_16S)
# remove the "." prefix where it exists
sample_names_16S <- gsub("^\\.", "", sample_names_16S)
# remove the "K" prefix where it exists
sample_names_16S <- gsub("K", "", sample_names_16S)
# add "fk" prefix to match metabolomics data
sample_names_16S <- paste0("fk", sample_names_16S)
sample_df$sample_name<- sample_names_16S
View(sample_df) # check the sample names

sample_names_16S[6] <- "fk161" # manually change one which went wrong
sample_names_16S[12] <- "fk140" # manually change one which went wrong
sample_names_16S[75] <- "fk88" # manually change one which went wrong

# add changes
sample_df$sample_name<- sample_names_16S

# remove the duplicate sample name
#View((sample_df[sample_df$sample_name == "fk88", ])) # check the sample names
sample_to_remove<- rownames(sample_df[sample_df$sample_name == "fk88", ][2]) # check the sample names
sample_to_remove
ps_16S <- prune_samples(sample_names(ps_16S) != sample_to_remove, ps_16S)
sample_df <- sample_df[rownames(sample_df) != sample_to_remove, ] # remove the sample from the sample data frame

View(sample_df)

# update the sample names in the phyloseq object
sample_names(ps_16S) <- sample_df$sample_name

## make a table of diversity metrics
# get diversity metrics for the 16S data
seq_counts <- sample_sums(ps_16S)
View(seq_counts)

div_16S_df <- data.frame(estimate_richness(ps_16S, measures = c("Shannon", "Observed")), seq_count = seq_counts)
div_16S_df$sample_name <- rownames(div_16S_df)

View(div_16S_df) # check the diversity metrics

write.csv(div_16S_df, file = "../results/16S_div_seqcounts.csv", row.names = FALSE)

#==================================
# MERGE THE DATA
#==================================
# Merge metabolomics and 16S data frames by sample name
merged_div_df <- left_join(met_div_df, div_16S_df, by = "sample_name")

merged_div_df <- merged_div_df %>%
    filter(!is.na(Shannon))

View(merged_div_df) # check the merged data frame

library(ggplot2)
# Plot the diversity metrics
ggplot(merged_div_df, aes(x = metab_diversity, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Diversity Metrics Comparison",
       x = "Metabolomics Diversity (Shannon)",
       y = "16S Diversity (Shannon)") +
  theme_minimal()

# run a simple linear regression
lm_result <- lm(Shannon ~ metab_diversity, data = merged_div_df)
summary(lm_result) # significant slope indicates a relationship between the two diversity metrics


str(merged_div_df)

# remove zero seq_depth from merged_div_df
merged_div_df_someseqs <- merged_div_df %>%
    filter(seq_count > 0)

# sequencing depth effect
plot(log(merged_div_df$seq_count), merged_div_df$Shannon, xlab = "Log sequencing Depth", ylab = "16S Diversity (Shannon)", main = "Sequencing Depth vs 16S Diversity")
lm_log_depth_vs_16S_div <- lm(Shannon ~ log(seq_count), data = merged_div_df_someseqs)
summary(lm_log_depth_vs_16S_div) # significant slope indicates a relationship between sequencing depth

# adding sequencing depth as a covariate
lm_result_depth <- lm(metab_diversity ~ Shannon + seq_count, data = merged_div_df)
lm_result_depth2 <- lm(Shannon ~ metab_diversity + seq_count, data = merged_div_df)


lm_result_logdepth <- lm(metab_diversity ~ Shannon + log(seq_count), data = merged_div_df_someseqs)
lm_result_logdepth2 <- lm(Shannon ~ metab_diversity + log(seq_count), data = merged_div_df_someseqs)

summary(lm_result_depth) # significant slope indicates a relationship between the two diversity
summary(lm_result_depth2) # significant slope indicates a relationship between the two diversity metrics

