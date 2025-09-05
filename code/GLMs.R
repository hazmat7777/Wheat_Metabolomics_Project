## load data. NOTE i need to add the additive GLMs for microb and metabs.
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

# load important genera
important_genera <- c("Acidiferrimicrobium","Incertae Sedis.317","Mariniblastus")
taxa_names(ps)
ps_genus <- ps %>%
  subset_taxa(taxa_names(ps) %in% important_genera)
otu_df <- as.data.frame(otu_table(ps_genus) )
sample_df <- as.data.frame(sample_data(ps_genus))

#merge sample and otu data
otu_gt_df <- otu_df %>%
  t() %>%
  as.data.frame() %>%
  mutate(sample_name = rownames(.)) %>%
  left_join(sample_df, by = "sample_name") %>%
  select(sample_name, important_genera, gt)

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



library(corrplot)

# Rename Incertae Sedis.317 to Unclassified Paracoccaceae
merged_df_protection <- merged_df_protection %>%
  rename(`Paracoccaceae` = `Incertae Sedis.317`)

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


merged_df_newnames <- merged_df  %>%
  rename(`Paracoccaceae` = `Incertae Sedis.317`) %>%
  rename(`gt` = `gt.x`)

str(merged_df_newnames)


## GLMs (so far have only combined one)

# all metabolites
glm_metab <- glm(as.factor(gt) ~ Haematommic.acid + X3.4.Dihydroxybenzaldehyde + Ribonolactone,
                 data = merged_df_newnames,
                 family = binomial)

summary(glm_metab)

# check variable inflation factors
library(car)
vif_values <- vif(glm_metab)
vif_values

# just the most important metabolite
glm_34 <- glm(as.factor(gt) ~ X3.4.Dihydroxybenzaldehyde,
                 data = merged_df_newnames,
                 family = binomial)

summary(glm_34) # higher AIC

# all genera
glm_microb <- glm(as.factor(gt) ~  Acidiferrimicrobium + Paracoccaceae + Mariniblastus,
                  data = merged_df_newnames,
                  family = binomial)

summary(glm_microb)

glm_acid <- glm(as.factor(gt) ~  Acidiferrimicrobium,
                  data = merged_df_newnames,
                  family = binomial)

summary(glm_acid)

# metabs combined
str(merged_df_combinedprotection)
glm_metab_combined <- glm(gt.x ~ metab_protection_scaled,
                    data = merged_df_combinedprotection, 
                    family = binomial)

summary(glm_metab_combined)

#microbs combined
glm_microb_combined <- glm(gt.x ~ microb_protection_scaled,
                    data = merged_df_combinedprotection, 
                    family = binomial)

summary(glm_microb_combined)

# combining metabs and microbs
glm_combined_adding <- glm(gt.x ~ microb_protection_scaled + metab_protection_scaled, 
                    data = merged_df_combinedprotection, 
                    family = binomial)

summary(glm_combined_adding)
vif(glm_combined_adding)

glm_combined <- glm(gt.x ~ combined_protection, 
                    data = merged_df_combinedprotection, 
                    family = binomial)

summary(glm_combined)

## everything added
glm_everything <- glm(as.factor(gt) ~  Acidiferrimicrobium + Paracoccaceae + Mariniblastus + Haematommic.acid + X3.4.Dihydroxybenzaldehyde + Ribonolactone,
                  data = merged_df_newnames,
                  family = binomial)


summary(glm_everything)

# pseudo-r2
nagelkerke_all_r2 <- PseudoR2(glm_everything, which = "Nagelkerke")
nagelkerke_all_r2
nagelkerke_metab_r2 <- PseudoR2(glm_metab, which = "Nagelkerke")
nagelkerke_metab_r2
nagelkerke_microb_r2 <- PseudoR2(glm_microb, which = "Nagelkerke")
nagelkerke_microb_r2


# Calculate pseudo-R² values
#install.packages("DescTools")
library(DescTools)

# McFadden's pseudo-R²
mcfadden_r2 <- PseudoR2(glm_everything, which = "McFadden")
mcfadden_r2
# Nagelkerke's pseudo-R²  
nagelkerke_r2 <- PseudoR2(glm_combined, which = "Nagelkerke")
nagelkerke_r2
# Or calculate McFadden's manually:
null_deviance <- glm_combined$null.deviance
residual_deviance <- glm_combined$deviance
mcfadden_manual <- 1 - (residual_deviance / null_deviance)

cat("McFadden's pseudo-R²:", round(mcfadden_manual, 3), "\n")


# Create summary dataframe of model comparisons
model_summary <- data.frame(
 Model = c("Combined Protection", "Metabolite Only", "Microbe Only", "Both Separate"),
 AIC = c(AIC(glm_combined), AIC(glm_metab), AIC(glm_microb), AIC(glm_both)),
 Delta_AIC = c(0, NA, NA, NA)
)

# Calculate delta AIC (difference from best model)
best_aic <- min(model_summary$AIC)
model_summary$Delta_AIC <- model_summary$AIC - best_aic

# Add pseudo-R² (McFadden's)
model_summary$Pseudo_R2 <- c(
 1 - (glm_combined$deviance / glm_combined$null.deviance),
 1 - (glm_metab$deviance / glm_metab$null.deviance),
 1 - (glm_microb$deviance / glm_microb$null.deviance),
 1 - (glm_both$deviance / glm_both$null.deviance)
)

# Round values
model_summary$AIC <- round(model_summary$AIC, 2)
model_summary$Delta_AIC <- round(model_summary$Delta_AIC, 2)
model_summary$Pseudo_R2 <- round(model_summary$Pseudo_R2, 3)

# Display table
print(model_summary)

write.csv(model_summary, "../supplementary/borutamodel_AIC_comparison")




