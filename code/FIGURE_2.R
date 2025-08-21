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
    metab_protection = Haematommic.acid + X3.4.Dihydroxybenzaldehyde + Ribonolactone,
    microb_protection = `Incertae Sedis.317` + Mariniblastus - Acidiferrimicrobium
  ) %>%
  filter(!is.na(metab_protection) & !is.na(microb_protection))


str(merged_df_protection)
######################### Microbial abundance plots

library(ggsignif)
library(scales)
library(cowplot)

# Prepare data for microbial plots
# Rename for clarity
merged_df_protection <- merged_df_protection %>%
  rename(`Unclassified Paracoccaceae` = `Incertae Sedis.317`)

# Get the three genera
important_genera_renamed <- c("Acidiferrimicrobium", "Unclassified Paracoccaceae", "Mariniblastus")

# Create GLM model for predictions
microbial_glm <- glm(factor(gt.x) ~ Acidiferrimicrobium + `Unclassified Paracoccaceae` + Mariniblastus, 
                     data = merged_df_protection, 
                     family = binomial)

# Add predictions to data
merged_df_protection$predicted_prob <- predict(microbial_glm, type = "response")

# Plot D: Protective score boxplot
pA_microb <- ggplot(merged_df_protection, aes(x = factor(gt.x), y = microb_protection, fill = factor(gt.x))) +
  geom_boxplot(alpha = 0.7) +
  geom_signif(comparisons = list(c("GT_absent", "GT_present")), 
              test = "t.test",
              map_signif_level = TRUE,
              y_position = max(merged_df_protection$microb_protection, na.rm = TRUE) + 0.01) +
  labs(
    title = "D",
    x = "Disease Status",
    y = "Microbial Protective Score",
    fill = "Disease Status"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 19),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20))

# Plot A: Acidiferrimicrobium
p1_microb <- ggplot(merged_df_protection, aes(x = Acidiferrimicrobium, y = predicted_prob)) +
  geom_point(aes(color = factor(gt.x)), alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "A",
       x = "Acidiferrimicrobium",
       y = "Predicted Disease Probability",
       color = "Disease Status") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 19),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20))

# Plot B: Unclassified Paracoccaceae
p2_microb <- ggplot(merged_df_protection, aes(x = `Unclassified Paracoccaceae`, y = predicted_prob)) +
  geom_point(aes(color = factor(gt.x)), alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "B",
       x = "Unclassified Paracoccaceae",
       y = "Predicted Probability",
       color = "Disease Status") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 19),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20))

# Plot C: Mariniblastus
p3_microb <- ggplot(merged_df_protection, aes(x = Mariniblastus, y = predicted_prob)) +
  geom_point(aes(color = factor(gt.x)), alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "C",
       x = "Mariniblastus",
       y = "Predicted Probability",
       color = "Disease Status") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 19),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20))

# Create legend from one plot
p2_legend_microb <- p2_microb + theme(legend.position = "bottom")
shared_legend_microb <- get_legend(p2_legend_microb)

# Combine plots
combined_plots_microb <- plot_grid(
  p1_microb, p2_microb,
  p3_microb, pA_microb,
  labels = NULL,
  ncol = 2,
  align = "hv"
)

# Add shared legend
final_plot_microb <- plot_grid(
  combined_plots_microb,
  shared_legend_microb,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

print(final_plot_microb)

# Save the plot
ggsave("../figures_tables/figure2.png",
       final_plot_microb, width = 8, height = 10, dpi = 300)