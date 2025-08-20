###########################
library(phyloseq)
library(dplyr)

# load ps
ps <- readRDS("../data/metabarcoding/ps.rds")
print(ps)

## PREPROCESSING

# what to filter
tax_data <- as.data.frame(tax_table(ps))

unique(tax_data$Kingdom)
unique(tax_data$Phylum)
unique(tax_data$Class)
sort(unique(tax_data$Order))
sort(unique(tax_data$Family))

# filter out unclassified taxa, mitochondria and chloroplasts

ps_project2 <- ps %>%
    phyloseq::subset_taxa(Order != "Incertae Sedis") %>% 
    phyloseq::subset_taxa(Class != "Incertae Sedis") %>%  
    phyloseq::subset_taxa(Phylum != "unclassified_Bacteria") %>%
    phyloseq::subset_taxa(Kingdom != "Eukaryota") %>% # only proks
    phyloseq::subset_taxa(Order != "Chloroplast") %>% # only interested in bacterial microbiome here
    phyloseq::subset_taxa(Family != "Mitochondria") # not organelles

ps
ps_project2 # still got plenny taxa

# transform to relative abundance using a custom fn
ps_project3 <- transform_sample_counts(ps_project2,
    function(x) x / sum(x)) 

sum(is.na(otu_table(ps_project3)))
sum(is.nan(otu_table(ps_project3)))
sum(is.infinite(otu_table(ps_project3)))
    # no division by NA or zero

#View(otu_table(ps_project3))

# filter rarest taxa (whose abundance is less than 5e-4 across samples)
ps_project4 <- filter_taxa(ps_project3, function(x) mean (x) > 5e-7, TRUE) 

ps_project4 # heavily filtered- from 5k to 324 to 255

## merging the rare taxa into one 'other' group

# get the names of the rare taxa
all_taxa <- taxa_names(ps_project3)
abundant_taxa <- taxa_names(ps_project4)

rare_taxa <- setdiff(all_taxa, abundant_taxa)

rare_taxa
length(rare_taxa) # 59k lol

ps_project4

# merge and rename the rare taxa
ps_project5 <- merge_taxa(ps_project2, rare_taxa)

merged_group_name <- intersect(taxa_names(ps_project5), rare_taxa) # the name given to the rare_taxa merged group

merged_index <- which(taxa_names(ps_project5) == merged_group_name)
taxa_names(ps_project5)[merged_index] <- "Rare_taxa"

"Rare_taxa" %in% taxa_names(ps_project5)  # TRUE, means rename was successful

# change the taxonomic table entries for the merged group
tax_table(ps_project5)["Rare_taxa", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- "Other"

# change to relative rather than absolute abundance
ps_project6 <- transform_sample_counts(ps_project5,
    function(x) x / sum(x)) # first transform to relative abundance using a custom fn








#############################









# filtering out samples with too low diversity
richness_df <- estimate_richness(ps_project2, measures = c("Observed", "Shannon", "Simpson"))
#View(richness_df)
saveRDS(richness_df, file = "../data/all_16s_richness_df_new.rds")
write.csv(richness_df, file = "../results/all_16s_richness_df_new.csv")

ps_16S_highdiv_absolute <- prune_samples(richness_df$Observed > 300, ps_project2) # keep samples with more than 6 observed taxa
ps_16S_highdiv_relative <- prune_samples(richness_df$Observed > 300, ps_project3) # keep samples with more than 6 observed taxa

# filter out samples whose
ps_16S_highdiv_absolute <- prune_samples(!is.na(sample_data(ps_16S_highdiv_absolute)$gt), ps_16S_highdiv_absolute)
ps_16S_highdiv_relative <- prune_samples(!is.na(sample_data(ps_16S_highdiv_relative)$gt), ps_16S_highdiv_relative)

# should have done this ages ago- filter out OTUs with zero abundance across all samples
ps_16S_highdiv_absolute <- prune_taxa(taxa_sums(ps_16S_highdiv_absolute) > 0, ps_16S_highdiv_absolute)
ps_16S_highdiv_relative <- prune_taxa(taxa_sums(ps_16S_highdiv_relative) > 0, ps_16S_highdiv_relative)

# save filtered objects
saveRDS(ps_16S_highdiv_absolute, file = "../data/metabarcoding/ps_16S_highdiv_absolute.rds")
saveRDS(ps_16S_highdiv_relative, file = "../data/metabarcoding/ps_16S_highdiv_relative.rds")

# load the filtered objects
ps_16S_highdiv_absolute <- readRDS("../data/metabarcoding/ps_16S_highdiv_absolute.rds")
ps_16S_highdiv_relative <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")

ps_16S_highdiv_absolute
ps_16S_highdiv_relative









# DO ORDINATION NOW

# Summarize the data at the phylum level
ps_phylum <- tax_glom(ps_16S_highdiv_absolute, taxrank = "Phylum")
taxa_names(ps_phylum) <- tax_table(ps_phylum)[, "Phylum"] #rename taxa to phylum names

# ordinate the data using NMDS
ps_phylum_ord <- ordinate(ps_phylum, "NMDS", "bray")
p_phylum <- plot_ordination(ps_phylum, ps_phylum_ord, type = "sample", color = "gt", title = "samples")
  
plot(p_phylum)

saveRDS(p_phylum, file = "../results/plots/ps_phylum_ordination.rds")
ggsave("../results/plots/ordination_highdiv_phylum_gt.png", plot = p_phylum, width = 8, height = 6, dpi = 300)

View(sample_data(ps_phylum))

# ordination at ASV level


# Option A: Keep top N most abundant taxa
top_n <- 100000  # Adjust based on your needs
top_taxa <- names(sort(taxa_sums(ps_16S_highdiv_absolute), decreasing = TRUE))[1:top_n]
ps_filtered_abs <- prune_taxa(top_taxa, ps_16S_highdiv_absolute)
ps_filtered_rel <- prune_taxa(top_taxa, ps_16S_highdiv_relative)


# ordinate at esv level
ps_asv_ord <- ordinate(ps_filtered_rel, "NMDS", "bray")
p_asv_ord <- ordinate(ps_16S_highdiv_absolute, "NMDS", "bray")
p_asv <- plot_ordination(ps_filtered_rel, ps_asv_ord, type = "samples", color = "gt", title = "samples")
p_asv


#ggsave("../results/plots/ordination_highdiv_gt.png", plot = p_asv, width = 8, height = 6, dpi = 300)

###############################

## heatmap?
plot_heatmap(ps_phylum, low = "blue", high = "red")


#############################

# barplots

# really hard to do these with my very high diversity

plot_bar(ps_filtered_rel, fill = "Phylum") # top N taxa, showing absolute abundance
# plot_bar(ps_project4, fill = "Phylum") # no rare taxa, showing absolute abundance
#     # a bit weird that most of the rare taxa are found in one plate of samples
# plot_bar(ps_project5, fill = "Phylum") # rare taxa merged (other), showing absolute abundance
# plot_bar(ps_project6, fill = "Phylum") # rare taxa merged (other), showing relative abundance
plot_bar(ps_filtered_abs, fill = "Class") # rare taxa merged (other), showing relative abundance
# plot_bar(ps_project6, fill = "Family") # rare taxa merged (other), showing relative abundance
# plot_bar(ps_project2, fill = "Class") # all taxa, showing absolute abundance
# plot_bar(ps_project2, fill = "Family") # all taxa, showing absolute abundance




















#############################

# filtering out samples with too low diversity
nrow(tax_table(ps_project2)) # 1545 taxa
richness_df <- estimate_richness(ps_project2, measures = c("Observed", "Shannon", "Simpson"))
saveRDS(richness_df, file = "../data/all_16s_richness_df_new.rds")
write.csv(richness_df, file = "../results/all_16s_richness_df_new.csv")

# filter out low richness samples- do this once I've done the ordination
ps_16S_highdiv_absolute <- prune_samples(richness_df$Observed > 6, ps_project2) # keep samples with more than 6 observed taxa
    # when I filtered >5,  one point with richness 6 messed it up
    # will need to justify this filtering step in the report

ps_16S_highdiv_relative <- prune_samples(richness_df$Observed > 6, ps_project3)
ps_16S_highdiv_absolute
ps_16S_highdiv_relative

########### CHECKPOINT ###########

# save filtered objects
saveRDS(ps_16S_highdiv_absolute, file = "../data/ps_16S_highdiv_absolute.rds")
saveRDS(ps_16S_highdiv_relative, file = "../data/ps_16S_highdiv_relative.rds")


#################################

## how many sequences are there in each sample? (number of sequences found in all esvs per sample)
View(sample_sums(ps_16S_highdiv_absolute)) # quite high variation, also very low generally- we asked for 10k

# view community composition
plot_bar(ps_16S_highdiv_absolute, fill = "Phylum")
plot_bar(ps_16S_highdiv_relative, fill = "Phylum")
plot_bar(ps_16S_highdiv_absolute, fill = "Class")
plot_bar(ps_16S_highdiv_relative, fill = "Class")
plot_bar(ps_16S_highdiv_absolute, fill = "Family")
plot_bar(ps_16S_highdiv_relative, fill = "Family")

View(tax_table(ps_16S_highdiv_absolute))
View(otu_table(ps_16S_highdiv_absolute))

sample_names(ps_16S_highdiv_absolute)



## Ordination plots


##################################################




