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

## PREPROCESSING

# what to filter
tax_data <- as.data.frame(tax_table(ps_16S))

unique(tax_data$Kingdom)
unique(tax_data$Phylum)
unique(tax_data$Class)
sort(unique(tax_data$Order))
sort(unique(tax_data$Family))

# filter out unclassified taxa, mitochondria and chloroplasts

ps_project2 <- ps_16S %>%
    phyloseq::subset_taxa(Order != "Incertae Sedis") %>% 
    phyloseq::subset_taxa(Class != "Incertae Sedis") %>%  
    phyloseq::subset_taxa(Phylum != "unclassified_Bacteria") %>%
    phyloseq::subset_taxa(Kingdom != "Eukaryota") %>% # only proks
    phyloseq::subset_taxa(Order != "Chloroplast") %>% # only interested in bacterial microbiome here
    phyloseq::subset_taxa(Family != "Mitochondria") # not organelles

ps_16S
ps_project2

# filter out OTUs with no hits on any sample
ps_project2 <- prune_taxa(taxa_sums(ps_project2) > 0, ps_project2) #removes missing taxa

ps_project2 <- prune_samples(sample_sums(ps_project2) > 0, ps_project2) # removes empty samples

# Check the counts of taxa remaining after filtering
dim(tax_table(ps_project2))
dim(otu_table(ps_project2))

ps_project2 # 2 samples got filtered (not sure what went wrong with them!)

empty_samples <- setdiff(sample_names(ps_16S), sample_names(ps_project2)) # check which samples were removed

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

richness_df <- estimate_richness(ps_project2, measures = c("Observed", "Shannon", "Simpson"))
saveRDS(richness_df, file = "../data/all_16s_richness_df.rds")

View(otu_table(ps_project2))

View(richness_df)

##################################################


# filter out low richness samples
ps_16S_highdiv_absolute <- prune_samples(richness_df$Observed > 6, ps_project2) # keep samples with more than 6 observed taxa
    # when I filtered >5,  one point with richness 6 messed it up
    # will need to justify this filtering step in the report

ps_16S_highdiv_relative <- prune_samples(richness_df$Observed > 6, ps_project3)
ps_16S_highdiv_absolute
ps_16S_highdiv_relative


########### CHECKPOINT ###########

# save filtered objects
# saveRDS(ps_project7, file = "../data/ps_project7.rds")
# saveRDS(ps_project8, file = "../data/ps_project8.rds")
# saveRDS(ps_project2, file = "../data/ps_project2.rds")

# load em
# ps_project7 <- readRDS("../data/ps_project7.rds")
# ps_project8 <- readRDS("../data/ps_project8.rds")
# ps_project2 <- readRDS("../data/ps_project2.rds")


#################################

# view community composition
plot_bar(ps_project7, fill = "phylum")
plot_bar(ps_project7, fill = "class") 
plot_bar(ps_project7, fill = "family") 
plot_bar(ps_project8, fill = "phylum")
plot_bar(ps_project8, fill = "class") 
plot_bar(ps_project8, fill = "family") 

View(tax_table(ps_project8))

## Ordination plot

# Summarize the data at the phylum level
ps_phylum <- tax_glom(ps_project8, taxrank = "phylum")
taxa_names(ps_phylum) <- tax_table(ps_phylum)[, "phylum"] #rename taxa to phylum names

# ordinate the data using NMDS
ps_phylum_ord <- ordinate(ps_phylum, "NMDS", "bray")
p3 <- plot_ordination(ps_phylum, ps_phylum_ord, type = "sample", color = "sample_name", title = "samples")
  
plot(p3)

# ordination at ASV level
ps_asv_ord <- ordinate(ps_project2, "NMDS", "bray")
p4 <- plot_ordination(ps_project2, ps_asv_ord, type = "samples", color = "sample_name", title = "samples") # warning because may be too few data
p4
# same but with rare taxa merged
ps_asv_ord_2 <- ordinate(ps_project7, "NMDS", "bray")
p5 <- plot_ordination(ps_project7, ps_asv_ord_2, type = "samples", color = "sample_name", title = "samples") # warning because may be too few data
p5 # even fewer points, more overlap

head(ps_asv_ord$points)

# 1_demultiplex.bc1005--bc1101.hifi_reads -0.03200422 -0.01219933
# 1_demultiplex.bc1007--bc1056.hifi_reads -0.03207408 -0.01222850
# 1_demultiplex.bc1007--bc1096.hifi_reads -0.03209634 -0.01225365
# 1_demultiplex.bc1008--bc1082.hifi_reads -0.03223436 -0.01236374
# 1_demultiplex.bc1012--bc1065.hifi_reads -0.03200789 -0.01224109
# 1_demultiplex.bc1015--bc1044.hifi_reads -0.03218156 -0.01220171

# what reads are in these samples?

View(data.frame(OTU = rownames(otu_table(ps_project7)[, "1_demultiplex.bc1005--bc1101.hifi_reads"])[otu_table(ps_project7)[, "1_demultiplex.bc1005--bc1101.hifi_reads"] > 0], 
           Count = otu_table(ps_project7)[, "1_demultiplex.bc1005--bc1101.hifi_reads"][otu_table(ps_project7)[, "1_demultiplex.bc1005--bc1101.hifi_reads"] > 0]))

View(data.frame(OTU = rownames(otu_table(ps_project7)[, "1_demultiplex.bc1007--bc1056.hifi_reads"])[otu_table(ps_project7)[, "1_demultiplex.bc1007--bc1056.hifi_reads"] > 0], 
           Count = otu_table(ps_project7)[, "1_demultiplex.bc1007--bc1056.hifi_reads"][otu_table(ps_project7)[, "1_demultiplex.bc1007--bc1056.hifi_reads"] > 0]))

# look at the otu table
View(otu_table(ps_project7))[, c("1_demultiplex.bc1005--bc1101.hifi_reads", "1_demultiplex.bc1007--bc1056.hifi_reads", "1_demultiplex.bc1007--bc1096.hifi_reads", "1_demultiplex.bc1008--bc1082.hifi_reads", "1_demultiplex.bc1012--bc1065.hifi_reads", "1_demultiplex.bc1015--bc1044.hifi_reads")]

# they seem to have different ASVs, so why are they so close together in the ordination space?

# why so overlapping?



################################################

### the old thing i did

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

# what is this actually showing

# why are points overlapping
head(ps_ord$points)

str(ps_ord$points)

# Sort the ordination points by the first column
sorted_points <- as.data.frame(ps_ord$points) %>%
    arrange(.[, 2]) # Sort by the first column

# View the sorted data
View(sorted_points)

library(ggplot2)
print(p1 + geom_jitter(width = 0.005, height = 0.005))

p2 <- plot_ordination(ps_project2, ps_ord, type = "taxa", color = "phylum", title = "taxa")

print(p2) # points overlap again
print(p2+ geom_jitter(width = 0.005, height = 0.005))

# what is this actually showing

ordination_coords <- ps_ord$points # or ps_ord@points
nrow(unique(ordination_coords))


ntaxa(ps_project6)
head(sample_data(ps_project6)[,1:5])

###############################



## heatmap?
plot_heatmap(ps_phylum, low = "blue", high = "red")

colnames(tax_table(ps_project8))



## how many sequences are there in each sample?
View(sample_sums(ps_project7))

# how many unique taxa are there per sample?
View(estimate_richness(ps_project7, measures = c("Observed", "Shannon", "Simpson")))
# doesnt work with ps_project8 because it got changed to relative abundance

richness_df <- estimate_richness(ps_project7, measures = c("Observed", "Shannon", "Simpson"))
saveRDS(richness_df, file = "../data/richness_df.rds")

View(otu_table(ps_project7))

View(richness_df)
View(sample_sums(ps_project7))

##################################################

# the ordination plot had overlap because the low richness samples were outliers, and caused the high richness samples to squeeze together
library(ggplot2)

# removing low richness samples
# load em
ps_project7 <- readRDS("../data/ps_project7.rds")
ps_project8 <- readRDS("../data/ps_project8.rds")

# merge richness data into sample data of ps_project7 and ps_project8
sample_data(ps_project7)$Observed_richness <- richness_df$Observed
sample_data(ps_project8)$Observed_richness <- richness_df$Observed

# filter out low richness samples
ps_project9 <- prune_samples(richness_df$Observed > 6, ps_project7) # keep samples with more than 10 observed taxa
ps_project10 <- prune_samples(richness_df$Observed > 6, ps_project8) # keep samples with more than 10 observed taxa
    # when I filtered >5,  one point with richness 6 messed it up
    # will need to justify this filtering step in the report

ps_project9 # 22 samples

# save the filtered object
saveRDS(ps_project9, file = "../data/ps_project9.rds")
saveRDS(ps_project10, file = "../data/ps_project10.rds")

# plot ordination again
ps_phylum_ord2 <- ordinate(ps_project9, "NMDS", "bray")

# Create the ordination plot
p6 <- plot_ordination(ps_project9, ps_phylum_ord2, type = "sample", color = "sample_name", title = "Samples After Removing Low Richness Samples") +
    geom_text(aes(label = Observed_richness), size = 3, vjust = -0.5)  # Add richness labels

# Plot the updated ordination
plot(p6)
plot(p6)



