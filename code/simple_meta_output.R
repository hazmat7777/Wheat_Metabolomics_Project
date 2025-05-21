# load dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

library(BiocManager)
library(phyloseq)

library(dplyr)
library(tibble)

#load data
ps_project <- readRDS("../data/ps_farmkits.rds")

ps_project 

## PREPROCESSING

# filter out unclassified taxa, mitochondria and chloroplasts
ps_project2 <- ps_project2 %>%
    phyloseq::subset_taxa(order != "unclassified_Cyanobacteriia") %>% 
    phyloseq::subset_taxa(class != "unclassified_Cyanobacteria") %>%  
    phyloseq::subset_taxa(phylum != "unclassified_Bacteria") %>% # havent removed NAs, useful to see where they are and how many
    phyloseq::subset_taxa(domain != "unclassified_Root") %>% 
    phyloseq::subset_taxa(domain != "Eukaryota") %>%
    phyloseq::subset_taxa(order != "Chloroplast") %>% # only interested in bacterial microbiome here
    phyloseq::subset_taxa(family != "Mitochondria") # not organelles

# filter out OTUs with no hits on any sample
ps_project2 <- prune_taxa(taxa_sums(ps_project) > 0, ps_project) 

ps_project
ps_project2 # can see how many taxa got filtered

# first transform to relative abundance using a custom fn
ps_project3 <- transform_sample_counts(ps_project2,
    function(x) x / sum(x)) 

dim(otu_table(ps_project3))
rowSums(otu_table(ps_project3))

# filter rarest taxa (whose abundance is less than 5e-4 across samples)
ps_project4 <- filter_taxa(ps_project3, function(x) mean (x) > 5e-4, TRUE) 

ps_project4 # heavily filtered- from 5k to 469 to 307

## merging the rare taxa into one 'other' group

# get the names of the rare taxa
all_taxa <- taxa_names(ps_project2)
abundant_taxa <- taxa_names(ps_project4)

rare_taxa <- setdiff(all_taxa, abundant_taxa)

rare_taxa
length(rare_taxa) # 162

# merge and rename the rare taxa
ps_project5 <- merge_taxa(ps_project2, rare_taxa)

merged_group_name <- intersect(taxa_names(ps_project5), rare_taxa) # the name given to the rare_taxa merged group

merged_index <- which(taxa_names(ps_project5) == merged_group_name)
taxa_names(ps_project5)[merged_index] <- "Rare_taxa"

"Rare_taxa" %in% taxa_names(ps_project5)  # Should return TRUE

# change the taxonomic table entries for the merged group
tax_table(ps_project5)["Rare_taxa", c("domain", "phylum", "class", "order", "family", "genus")] <- "Other"

# change to relative rather than absolute abundance
ps_project6 <- transform_sample_counts(ps_project5,
    function(x) x / sum(x)) # first transform to relative abundance using a custom fn


# # MAKE A NEW BINARY VARIABLE- HUMAN
# sample_data(GP)$human = factor( get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue") )

## TAXA TABLE
taxa_names(ps_project6)

rank_names(ps_project4) # same as:  colnames(tax_table(ps_project2))

View(tax_table(ps_project6))


# unique(taxdata_df$rootrank) # all that col is root
unique(taxdata_df$phylum)

## ENVIRONMENTAL SAMPLE DATA
sampledata <-(sample_data(ps_project4)) # looks like the sample data were already there

sample_names(ps_project4)
view(sampledata)

# ****** sample data doesnt have sampleID or plate ********
# ****** sample data are messy from farmers- clean on env_datawrang *********
    # leave it for now

## TREE
phy_tree(ps_project4) # empty

# SEQUENCES
refseq(ps_project4)
    # actual nucleotide sequence of each OTU/ASV

### GRAPHICS

## Bar plot
plot_bar(ps_project2, fill = "class") 

plot_bar(ps_project4, fill = "class")
    # useful to see that some samples lost more species than others. NEED TO MAKE AN OTHER GROUP!
    # note that abundance is now a proportion- is that good?
    # 1024:1033 is ALL NAs!

plot_bar(ps_project4, fill = "phylum") # maybe omit taxa which are fewer than a threshold
plot_bar(ps_project5, fill = "phylum") # maybe omit taxa which are fewer than a threshold
plot_bar(ps_project6, fill = "phylum") # maybe omit taxa which are fewer than a threshold
plot_bar(ps_project6, fill = "class") # maybe omit taxa which are fewer than a threshold
plot_bar(ps_project6, fill = "family") # maybe omit taxa which are fewer than a threshold

# samples to remove based on these:
    # 1024-1033 all NAs
    # 1008- 1044
    # 1005- 1044 small abundance one class
    #1020-1096 all one class
    #1005-1110 and 1007-1045 all one class

sample_names(ps_project6)[15] #dont get the order but heres one to remove

## Ordination plot
ps_ord <- ordinate(ps_project2, "NMDS", "bray")
# what is ordination
    # makes a distance matrix
    # based on how similar two communities are
    # two communities with similar amounts of the same ASV = low distance
    # so why (below) are many points overlapping?

p1 <- plot_ordination(ps_project2, ps_ord, type = "sample", color = "sample_name", title = "samples")
    # warning because may be too few data
print(p1)
    # p1 is very different when you lump together the rare taxa (it's worse, distances get massive)
# why are points overlapping
head(ps_ord$points)

library(ggplot2)
print(p1 + geom_jitter(width = 0.005, height = 0.005))



p2 <- plot_ordination(ps_project2, ps_ord, type = "taxa", color = "phylum", title = "taxa")

print(p2) # points overlap again
print(p2+ geom_jitter(width = 0.005, height = 0.005))
#

ordination_coords <- ps_ord$points # or ps_ord@points
nrow(unique(ordination_coords))


ntaxa(ps_project6)
head(sample_data(ps_project6)[,1:5])

## heatmap?
plot_heatmap(ps_project4) # dont think this worked


