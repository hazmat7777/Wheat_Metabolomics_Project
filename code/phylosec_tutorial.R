
library("phyloseq")

packageVersion("phyloseq")

library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

# tools for constructing phyloseq component data
?phyloseq

# create a pretend otu table that you read from a file called otumat
otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
otumat # shows abundances of microbial taxa (OTUs) across samples

# it needs sample and OTU names
rownames(otumat) <- paste0("OTU",1:nrow(otumat))
colnames(otumat) <- paste0("Sample",1:nrow(otumat))
otumat

# now a pretend taxonomy table
taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat

class(otumat)
class(taxmat)

# combine these vanilla r matrices into a phyloseq object
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
TAX

physeq = phyloseq(OTU,TAX)
physeq

# visualise family level composition
plot_bar(physeq, fill = "Family")

# create random sample data and add to the combined dataset
sampledata <- sample_data(data.frame(
    Location = sample(LETTERS[1:4], size = nsamples(physeq), replace = TRUE),
    Depth = sample(50:1000, size = nsamples(physeq), replace = TRUE),
    row.names=sample_names(physeq),
    stringsAsFactors = FALSE
))
sampledata

# create random phylo tree with ape package and add it to dataset
library("ape")
random_tree = rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
plot(random_tree)

# combine new data components to the phyloseq object

    #using merge_phyloseq
physeq1 <- merge_phyloseq(physeq, sampledata, random_tree)
physeq1

    #fresh call to phyloseq() to build it again from scratch
physeq2 <- phyloseq(OTU, TAX, sampledata, random_tree)
physeq2

#theyre identical
identical(physeq1, physeq2)

# buidlign a tree plot
plot_tree(physeq1,color = "Location", label.tips = "taxa_names", ladderize = "left", plot.margin = 0.3)
plot_tree(physeq1, color="Depth", shape="Location", label.tips="taxa_names", ladderize="right", plot.margin=0.3)

plot_heatmap(physeq1)
plot_heatmap(physeq1, taxa.label="Phylum")

### DATA IMPORT

## importing .biom files

# first define the file paths
rich_dense_biom  = system.file("extdata", "rich_dense_otu_table.biom",  package="phyloseq")
rich_sparse_biom = system.file("extdata", "rich_sparse_otu_table.biom", package="phyloseq")
min_dense_biom   = system.file("extdata", "min_dense_otu_table.biom",   package="phyloseq")
min_sparse_biom  = system.file("extdata", "min_sparse_otu_table.biom",  package="phyloseq")
treefilename = system.file("extdata", "biom-tree.phy",  package="phyloseq")
refseqfilename = system.file("extdata", "biom-refseq.fasta",  package="phyloseq")

# now can use filepaths as an argument to the import_biom fn
import_biom(rich_dense_biom, treefilename, refseqfilename, parseFunction = parse_taxonomy_greengenes)

import_biom(rich_sparse_biom, treefilename, refseqfilename, parseFunction=parse_taxonomy_greengenes)

import_biom(min_dense_biom, treefilename, refseqfilename, parseFunction=parse_taxonomy_greengenes)

import_biom(min_sparse_biom, treefilename, refseqfilename, parseFunction=parse_taxonomy_greengenes)

myData = import_biom(rich_dense_biom, treefilename, refseqfilename, parseFunction=parse_taxonomy_greengenes)
myData # saving it for later

# making plots
plot_tree(myData, color = "Genus", shape ="BODY_SITE", size= "abundance")

# why doesnt this work, its from the tutorial
plot_richness(myData, x="BODY_SITE", color="Description")

plot_bar(myData, fill = "Genus")

refseq(myData)

## import_mothur
# software package to process barcoded amplicon sequences and perform OTU clustering
mothlist  = system.file("extdata", "esophagus.fn.list.gz", package="phyloseq")
mothgroup = system.file("extdata", "esophagus.good.groups.gz", package="phyloseq")
mothtree  = system.file("extdata", "esophagus.tree.gz", package="phyloseq")
show_mothur_cutoffs(mothlist)

cutoff  ="0.10"
x = import_mothur(mothlist, mothgroup, mothtree, cutoff)
x

plot_tree(x, color = "samples")

SDF = data.frame(samples=sample_names(x), row.names=sample_names(x))
sample_data(x) <- sample_data(SDF)
plot_richness(x)

#### ACCESSING BITS OF DATA

data("GlobalPatterns")

class(GlobalPatterns) # the data comes from the phyloseq package

GlobalPatterns

ls("package:phyloseq") # get a list of all the commands/fns

ntaxa(GlobalPatterns) # all the (OTUs)

nsamples(GlobalPatterns)

sample_names(GlobalPatterns)[1:5] # first 5 samples

rank_names(GlobalPatterns) # instead of doing colnames in taxa_table

sample_variables(GlobalPatterns) # instead of colnames pf sample_data

otu_table(GlobalPatterns)[1:5, 1:5]
tax_table(GlobalPatterns)[1:5, 1:4]
phy_tree(GlobalPatterns)

taxa_names(GlobalPatterns)[1:10]

##### PREPROCESSING DATA

# use prune_taxa()/samples() when the complete seubset of OTUs/samples is available
myTaxa = names(sort(taxa_sums(GlobalPatterns),
    decreasing = TRUE)[1:10]) # the lowest 10
ex1 = prune_taxa(myTaxa, GlobalPatterns) # narrow GP to my named taxa
plot(phy_tree(ex1), show.node.label = TRUE)

plot_tree(ex1, color = "SampleType",
        label.tips = "Phylum", ladderize = "left",
        justify = "left" , size = "Abundance")

# use filter_taxa()/ filter_samples to subset based on data in that table
GPr = transform_sample_counts(GlobalPatterns,
    function(x) x / sum(x)) # first transform to relative abundance using a custom fn
GPfr = filter_taxa(GPr, function(x) mean (x) > 1e-5, TRUE) #highly-subsetted object

GPfr # contains just 4264 of original 19216

# subset then prune
GP.chl = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
GP.chl = prune_samples(sample_sums(GP.chl)>=20, GP.chl)

GP.chl

# merge taxa
GP.chl.merged = merge_taxa(GP.chl, taxa_names(GP.chl)[1:5]) # merges the first 5 OTUs in the gpchl dataset

GP.chl.merged # but why would it be useful to merge them?
    # maybe could merge the rare species into an "other"?

# merge (agglomerate) based on phylo/tax thresh
gpsfb = subset_taxa(GlobalPatterns, Phylum == "Bacteroidetes")
gpsfbg = tax_glom(gpsfb, "Family") # agglom this bd dataset at family rank

plot_tree(gpsfbg, color = "SampleType", shape = "Class", size = "abundance") # visualise

# transform abundance values
transform_sample_counts(GP.chl, function(OTU) OTU/sum(OTU) ) # convert to fractional abundance
    # I can't do this but could I make a new column of fractional abundance?

# remove taxa not seen more than 3 times in at least 20% of samples
GP = filter_taxa(GlobalPatterns, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
    # protects vs an OTU with small mean and large CV (?)

# MAKE A NEW BINARY VARIABLE- HUMAN
sample_data(GP)$human = factor( get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue") )

# standardize abundances to the median sequencing depth
total = median(sample_sums(GP))
standf = function(x, t=total) {
  round(t * (x / sum(x)))
}
gps = transform_sample_counts(GP, standf)

# filter the taxa using a cutoff of 3 for the coefficient of variation
gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)

# subset the data to bacteroidetes
gpsfb = subset_taxa(gpsf, Phylum=="Bacteriodetes")

# graphics
title = "plot_bar; Bacteroidetes-only"

plot_bar(gpsfb, "SampleType", "Abundance", title=title)
plot_bar(gpsfb, "SampleType", "Abundance", "Family", title=title)
plot_bar(gpsfb, "Family", "Abundance", "Family", 
         title=title, facet_grid="SampleType~.")

##### MAKING GRAPHICS 

library("plyr"); packageVersion("plyr")

theme_set(theme_bw())

# preprocessing like above
GP = GlobalPatterns
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP)) # remove otus that dont appear more than 5 times in more than half the samples
GP1 = prune_taxa(wh0, GP)

GP1

# tranform to even sampling depth
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

# only 5 most abundant phyla
phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

GP1 # leaves 204 taxa

# define a human-associated categorical variable
human = get_variable(GP1, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
sample_data(GP1)$human <- factor(human)

GP1 # Added a new sample variable

View(sample_data(GP1))

# make the ordination plot
GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)
    # apparently there is overplotting/occlusion
    # (too many points -> cant understand data visdually)
    # soln:

p1 + facet_wrap(~Phylum, 3) # 3 is the number of rows

