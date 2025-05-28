# the goal of this script
    # generate the three files needed for picrust pipeline
        # metadata
        # seqs.fna
        # table.biom

library(phyloseq)
library(Biostrings)

### getting seqs.fna

# load filtered data
ps <- readRDS("../data/ps_project7.rds")

# remove the rare_taxa bit
taxa_to_keep <- setdiff(taxa_names(ps), "Rare_taxa")

ps <- prune_taxa(taxa_to_keep, ps)

"Rare_taxa" %in% taxa_to_keep # removed successfully

# extract reference sequences
seqs <- refseq(ps)  # get sequences from ps (DNAStringSet object)

seqs <- seqs[taxa_names(ps)]  # Reorder and drop any extras

# Name the sequences with the OTU or ASV IDs
names(seqs) <- taxa_names(ps) # names are ESV numbers

# Write to FASTA file
writeXStringSet(seqs, filepath = "../data/farmkitseqs.fna")

### getting table.biom

library(biomformat)

otu_table <- otu_table(ps)
otu_mat <- as(otu_table, "matrix")
otu_biom <- make_biom(data = otu_mat)

write_biom(otu_biom, biom_file = "../data/farmkitstable.biom")
