# phyloseq object from simple meta pipeline
# this time from all 16S sequences
# here- filtering taxa and samples, and visualising composition

# load dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

library(BiocManager)
library(phyloseq)

library(dplyr)
library(tibble)

#load data
ps_16S <- readRDS("../data/metabarcoding/all_farmkits_16S_ps_project.rds")

ps_16S

View(tax_table(ps_16S))

# load tax table
tax_table(ps_16S) <- read.csv("../results/ESV_tax_table.csv", row.names = 1, header = TRUE) %>%
    as.matrix() %>%
    tax_table()

View(tax_table(ps_16S))

## datawrang to merge forward/reverse reads

# Extract row names (taxa names) from the tax_table
rownames_taxa <- rownames(tax_table(ps_16S))

# Remove the "ESV_" prefix
rownames_taxa_no_prefix <- gsub("ESV_", "", rownames_taxa)

# Convert the row names to numeric for sorting
rownames_numeric <- as.numeric(rownames_taxa_no_prefix)

# Sort the row names in ascending numeric order
sorted_indices <- order(rownames_numeric)

# Apply the sorted order to the tax_table
tax_table_sorted <- tax_table(ps_16S)[sorted_indices, ]
View(tax_table_sorted)
write.csv(tax_table_sorted, file = "../results/ESV_tax_table_sorted.csv")

# Get the indices of odd-numbered taxa
odd_indices <- which((sort(rownames_numeric) %%2 !=0))
tax_table_filtered <- tax_table_sorted[odd_indices, ]
View(tax_table_filtered)

# Update the phyloseq object with the new tax_table
ps_16S@tax_table <- tax_table_filtered

# View the updated tax_table
View(tax_table(ps_16S))

# trim the esvs to only the ones in the tax table
ps_16S_merged <- prune_taxa(taxa_names(ps_16S) %in% rownames(tax_table(ps_16S)), ps_16S) # think this will fuck up my data only slightly
    # a small number of esvs didnt have both forward and reverse strands- may have missed a couple
View(tax_table(ps_16S_merged))

## PREPROCESSING

# what to filter
tax_data <- as.data.frame(tax_table(ps_16S_merged))

unique(tax_data$Kingdom)
unique(tax_data$Phylum)
unique(tax_data$Class)
sort(unique(tax_data$Order))
sort(unique(tax_data$Family))

# filter out unclassified taxa, mitochondria and chloroplasts

ps_project2 <- ps_16S_merged %>%
    phyloseq::subset_taxa(Order != "Incertae Sedis") %>% 
    phyloseq::subset_taxa(Class != "Incertae Sedis") %>%  
    phyloseq::subset_taxa(Phylum != "unclassified_Bacteria") %>%
    phyloseq::subset_taxa(Kingdom != "Eukaryota") %>% # only proks
    phyloseq::subset_taxa(Order != "Chloroplast") %>% # only interested in bacterial microbiome here
    phyloseq::subset_taxa(Family != "Mitochondria") # not organelles

# check the taxonomic table
View(tax_table(ps_project2))

ps_16S
ps_16S_merged
ps_project2

# filter out OTUs with no hits on any sample
ps_project2 <- prune_taxa(taxa_sums(ps_project2) > 0, ps_project2) #removes missing taxa

ps_project2 <- prune_samples(sample_sums(ps_project2) > 0, ps_project2) # removes empty samples

# Check the counts of taxa remaining after filtering
dim(tax_table(ps_project2))
dim(otu_table(ps_project2))

ps_project2 # 2 samples got filtered (not sure what went wrong with them!)

empty_samples <- setdiff(sample_names(ps_16S), sample_names(ps_project2)) # check which samples were removed
empty_samples # 3 empties

# transform to relative abundance using a custom fn
ps_project3 <- transform_sample_counts(ps_project2,
    function(x) x / sum(x)) 

sum(is.na(otu_table(ps_project3)))
sum(is.nan(otu_table(ps_project3)))
sum(is.infinite(otu_table(ps_project3)))
    # no division by NA or zero

# filter rarest taxa (whose abundance is less than 5e-4 across samples)
ps_project4 <- filter_taxa(ps_project3, function(x) mean (x) > 5e-4, TRUE) 

ps_project4 # heavily filtered- from 5k to 324 to 255

## merging the rare taxa into one 'other' group

# get the names of the rare taxa
all_taxa <- taxa_names(ps_project3)
abundant_taxa <- taxa_names(ps_project4)

rare_taxa <- setdiff(all_taxa, abundant_taxa)

rare_taxa
length(rare_taxa) # 69

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


# graphics
plot_bar(ps_project2, fill = "Phylum") # all taxa, showing absolute abundance
plot_bar(ps_project4, fill = "Phylum") # no rare taxa, showing absolute abundance
    # a bit weird that most of the rare taxa are found in one plate of samples
plot_bar(ps_project5, fill = "Phylum") # rare taxa merged (other), showing absolute abundance
plot_bar(ps_project6, fill = "Phylum") # rare taxa merged (other), showing relative abundance
plot_bar(ps_project6, fill = "Class") # rare taxa merged (other), showing relative abundance
plot_bar(ps_project6, fill = "Family") # rare taxa merged (other), showing relative abundance
plot_bar(ps_project2, fill = "Class") # all taxa, showing absolute abundance
plot_bar(ps_project2, fill = "Family") # all taxa, showing absolute abundance

# filtering out samples with too low diversity
View(tax_table(ps_project2)) # 255 taxa
richness_df <- estimate_richness(ps_project2, measures = c("Observed", "Shannon", "Simpson"))
saveRDS(richness_df, file = "../data/all_16s_richness_df.rds")
write.csv(richness_df, file = "../results/all_16s_richness_df.csv")

View(otu_table(ps_project2))
otu_df <- as.data.frame(otu_table(ps_project2))

View(richness_df)

# filter out low richness samples
ps_16S_highdiv_absolute <- prune_samples(richness_df$Observed > 6, ps_project2) # keep samples with more than 6 observed taxa
    # when I filtered >5,  one point with richness 6 messed it up
    # will need to justify this filtering step in the report

ps_16S_highdiv_relative <- prune_samples(richness_df$Observed > 6, ps_project3)
ps_16S_highdiv_absolute
ps_16S_highdiv_relative

########### CHECKPOINT ###########

# save filtered objects
saveRDS(ps_16S_highdiv_absolute, file = "../data/ps_16S_highdiv_absolute_odd.rds")
saveRDS(ps_16S_highdiv_relative, file = "../data/ps_16S_highdiv_relative_odd.rds")

# load em
# ps_project7 <- readRDS("../data/ps_project7.rds")
# ps_project8 <- readRDS("../data/ps_project8.rds")
# ps_project2 <- readRDS("../data/ps_project2.rds")

#################################

## how many sequences are there in each sample?
View(sample_sums(ps_16S_highdiv_absolute)) # quite high variation

# view community composition
plot_bar(ps_16S_highdiv_absolute, fill = "Phylum")
plot_bar(ps_16S_highdiv_relative, fill = "Phylum")
plot_bar(ps_16S_highdiv_absolute, fill = "Class")
plot_bar(ps_16S_highdiv_relative, fill = "Class")
plot_bar(ps_16S_highdiv_absolute, fill = "Family")
plot_bar(ps_16S_highdiv_relative, fill = "Family")

View(tax_table(ps_16S_highdiv_absolute))

## Ordination plots

# Summarize the data at the phylum level
ps_phylum <- tax_glom(ps_16S_highdiv_relative, taxrank = "Phylum")
taxa_names(ps_phylum) <- tax_table(ps_phylum)[, "Phylum"] #rename taxa to phylum names

# ordinate the data using NMDS
ps_phylum_ord <- ordinate(ps_phylum, "NMDS", "bray")
p_phylum <- plot_ordination(ps_phylum, ps_phylum_ord, type = "sample", color = "sample_name", title = "samples")
  
plot(p_phylum)

View(sample_data(ps_phylum))

# ordination at ASV level (is that what this is doing?)
ps_asv_ord <- ordinate(ps_16S_highdiv_absolute, "NMDS", "bray")
p_asv <- plot_ordination(ps_16S_highdiv_absolute, ps_asv_ord, type = "samples", color = "sample_name", title = "samples")
p_asv

View(richness_df) # a had 13 ESVs, b had 7 ESVs

ps_asv_relative_ord <- ordinate(ps_16S_highdiv_relative, "NMDS", "bray")
p_asv_relative_tillage <- plot_ordination(ps_16S_highdiv_relative, ps_asv_relative_ord, type = "samples", color = "tillage", title = "samples")
p_asv_relative_samples <- plot_ordination(ps_16S_highdiv_relative, ps_asv_relative_ord, type = "samples", color = "sample_name", title = "samples")

p_asv_relative_tillage # maybe some clustering by tillage?
p_asv_relative_samples

###############################

## heatmap?
plot_heatmap(ps_phylum, low = "blue", high = "red")

##################################################




