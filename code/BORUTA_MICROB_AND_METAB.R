library(tidyverse)
library(phyloseq)

# load data
fk_metabolom_gt_logged <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")

ps <- readRDS("../data/metabarcoding/ps_16S_highdiv_genus_relative.rds")

# load important metabolites
important_metabs <- c("Haematommic.acid","X3.4.Dihydroxybenzaldehyde", "Ribonolactone")
fk_metabs <- fk_metabolom_gt_logged %>%
  mutate(sample_name = rownames(.)) %>%
  select(sample_name, important_metabs, gt)

View(fk_metabs)

# load important genera
important_genera <- c("Acidiferrimicrobium","Incertae Sedis.317","Mariniblastus")


taxa_names(ps)
ps_genus <- ps %>%
  subset_taxa(taxa_names(ps) %in% important_genera)

otu_df <- as.data.frame(otu_table(ps_genus) )

sample_df <- as.data.frame(sample_data(ps_genus))
View(sample_df)
View(otu_df)

#merge sample and otu data
otu_gt_df <- otu_df %>%
  t() %>%
  as.data.frame() %>%
  mutate(sample_name = rownames(.)) %>%
  left_join(sample_df, by = "sample_name") %>%
  select(sample_name, important_genera, gt)

View(otu_gt_df)

## merge metabolomics and metabarcoding data
merged_df <- otu_gt_df %>%
  left_join(fk_metabs, by = "sample_name") %>%
  filter(!is.na("gt"))

merged_df_protection <- merged_df %>%
  mutate(
    metab_protection = rowSums(select(., all_of(important_metabs))),
    microb_protection = rowSums(select(., all_of(c("Incertae Sedis.317", "Mariniblastus")))) - 
                        select(., "Acidiferrimicrobium") %>% pull()
  ) %>%
  filter(!is.na(metab_protection) & !is.na(microb_protection))


str(merged_df_protection)

cor(merged_df_protection$metab_protection, merged_df_protection$microb_protection)

library(tidyverse)
library(corrplot)

# Rename Incertae Sedis.317 to Unclassified Paracoccaceae
merged_df_protection <- merged_df_protection %>%
  rename(`Unclassified Paracoccaceae` = `Incertae Sedis.317`)

# Extract just the genera and metabolites for correlation
cor_data <- merged_df_protection %>%
  select(Acidiferrimicrobium, `Unclassified Paracoccaceae`, Mariniblastus,
         Haematommic.acid, X3.4.Dihydroxybenzaldehyde, Ribonolactone)

# Calculate correlation matrix
cor_matrix <- cor(cor_data, use = "complete.obs")

# Extract just the genera-metabolite correlations
genera_names <- c("Acidiferrimicrobium", "Unclassified Paracoccaceae", "Mariniblastus")
metabolite_names <- c("Haematommic.acid", "X3.4.Dihydroxybenzaldehyde", "Ribonolactone")

# Subset correlation matrix (genera vs metabolites)
genera_metab_cor <- cor_matrix[genera_names, metabolite_names]

# Option 1: Using corrplot
corrplot(genera_metab_cor, 
         method = "color",
         type = "full",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         title = "Genera-Metabolite Correlations",
         mar = c(0,0,2,0))

# Option 3: Using ggplot2 for more control
library(reshape2)
cor_melted <- melt(genera_metab_cor)
colnames(cor_melted) <- c("Genus", "Metabolite", "Correlation")

ggplot(cor_melted, aes(x = Metabolite, y = Genus, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 3)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Genera-Metabolite Correlations",
       x = "Metabolites", y = "Genera")

# Print the correlation values
print("Correlation matrix:")
print(round(genera_metab_cor, 3))


##########
#GLM



# Create combined protection score and standardize variables
merged_df_combinedprotection <- merged_df_protection %>%
  mutate(
    # Standardize individual scores
    metab_protection_scaled = scale(metab_protection)[,1],
    microb_protection_scaled = scale(microb_protection)[,1],
    # Combined protection score (equal weighting)
    combined_protection = metab_protection_scaled + microb_protection_scaled,
    # gt
    gt.x = as.factor(gt.x)
  ) %>%
  select(sample_name, metab_protection_scaled, microb_protection_scaled, combined_protection, gt.x)

View(merged_df_combinedprotection)

cor(merged_df_combinedprotection$metab_protection_scaled, merged_df_combinedprotection$microb_protection_scaled)

# Run GLM with combined protection score
glm_combined <- glm(gt.x ~ combined_protection, 
                    data = merged_df_combinedprotection, 
                    family = binomial)

# Summary of the model
summary(glm_combined)

# Optional: Compare with individual models
glm_metab <- glm(gt.x ~ metab_protection_scaled, 
                 data = merged_df_combinedprotection, 
                 family = binomial)

glm_microb <- glm(gt.x ~ microb_protection_scaled, 
                  data = merged_df_combinedprotection, 
                  family = binomial)

glm_both <- glm(gt.x ~ metab_protection_scaled + microb_protection_scaled, 
                data = merged_df_combinedprotection, 
                family = binomial)

# Compare AIC values
cat("AIC values:\n")
cat("Combined model:", AIC(glm_combined), "\n")
cat("Metabolite only:", AIC(glm_metab), "\n")
cat("Microbe only:", AIC(glm_microb), "\n")
cat("Both separate:", AIC(glm_both), "\n")