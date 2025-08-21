# phyloseq object from simple meta pipeline
# this time from all 16S sequences
# here- collating taxa data, gt status, esv abundance.

#load data

# Load required libraries
library(phyloseq)
library(readr)
library(dplyr)
library(readxl)
library(readODS)

# Read the ESV abundance table
# Assuming your file is saved as "esv_abundance_table.txt" or similar
abundance_data <- read_delim("../data/metabarcoding/Harry_filtered_farmkit_ESVs.tsv", delim = "\t")
str(abundance_data)
# Read the taxonomy table 
# Assuming your taxonomy file has columns like: ESV, Kingdom, Phylum, Class, Order, Family, Genus, Species
taxonomy_data <- read_xlsx("../data/metabarcoding/ESV_tax_table_new.xlsx")

# Prepare the abundance matrix (OTU table)
# Set ESV as rownames and remove the Sequence column
abundance_matrix <- abundance_data %>%
  select(-Sequence) %>%  # Remove sequence column
  column_to_rownames("ESV") %>%  # Set ESV as rownames
  as.matrix()

# Prepare the taxonomy matrix
# Ensure ESV names match between abundance and taxonomy tables
taxonomy_matrix <- taxonomy_data %>%
  column_to_rownames("ESV") %>%  # Set ESV as rownames
  as.matrix()

############################################ SAMPLE DATA of GT STATUS
ps_project2 <- readRDS("../data/metabarcoding/all_farmkits_16S_ps_project.rds")
ps_project2

sample_names(ps_project2)
sample_names(ps)

clean_sample_names <- sample_names(ps) %>%
  gsub("^Sample_", "", .) %>%
  gsub("_R1_001__Run[0-9]+$", "", .)

clean_sample_names %in% sample_names(ps_project2) # missing the last four
sample_names(ps) <- clean_sample_names

sample_df<- as_tibble(sample_data(ps_project2))
sample_df <- sample_df %>%
    mutate(sample_name = sample_names(ps_project2)) %>%
    select(sample_name, tillage)

View(sample_df)
gt <- read_ods("../data/baiting/baiting_farmkits.ods", sheet = 1) %>%
    as.data.frame() %>%
    select(sample_name, gt)
View(gt)

View(sample_data(ps_project2))
sample_gt_df <- left_join(sample_df, gt, by = "sample_name")
sample_gt_df <- as.data.frame(sample_gt_df)
rownames(sample_gt_df) <- sample_gt_df$sample_name

View(sample_gt_df)
sample_data(ps_project2) <- sample_gt_df
View(sample_data(ps_project2))

########################

# Create phyloseq components
OTU <- otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy_matrix)
#SAMPLE <- sample_data(sample_data)

# Create the phyloseq object
ps <- phyloseq(OTU, TAX)
ps
View(sample_sums(ps))


# Optional: Filter out ESVs with very low abundance
# Remove ESVs present in fewer than 2 samples or with total abundance < 10
ps_filtered <- filter_taxa(ps, function(x) sum(x > 0) >= 2 & sum(x) >= 10, TRUE)

cat("\nAfter filtering - Number of taxa:", ntaxa(ps_filtered))



############################# change the sample_data
ps_project2 <- readRDS("../data/metabarcoding/all_farmkits_16S_ps_project.rds")
sample_names(ps_project2)

# Get the sample names from the phyloseq object
sample_df <- as.data.frame(sample_data(ps_project2))
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
View((sample_df[sample_df$sample_name == "fk88", ])) # check the sample names
sample_to_remove<- rownames(sample_df[sample_df$sample_name == "fk88", ][2]) # check the sample names
sample_to_remove
ps_project2 <- prune_samples(sample_names(ps_project2) != sample_to_remove, ps_project2)
sample_df <- sample_df[rownames(sample_df) != sample_to_remove, ] # remove the sample from the sample data frame

View(sample_df)

# make the old sample names a column
sample_df$og_sample_names <- sample_names(ps_project2)

sample_df <- as_tibble(sample_df)
sample_df <- sample_df %>%
    select(sample_name, tillage, og_sample_names)

# NOW MERGE WITH GT

View(sample_df)
gt <- read_ods("../data/baiting/baiting_farmkits.ods", sheet = 1) %>%
    as.data.frame() %>%
    select(sample_name, gt)
View(gt)


sample_gt_df <- left_join(sample_df, gt, by = "sample_name")
sample_gt_df <- as.data.frame(sample_gt_df)
rownames(sample_gt_df) <- sample_gt_df$og_sample_names

View(sample_gt_df)

sample_data(ps) <- sample_gt_df

saveRDS(ps, file = "../data/metabarcoding/ps.rds")

sample_names(ps)
View(otu_table(ps))
View(tax_table(ps))
View(sample_sums(ps)) # check the sample sums

View(sample_data(ps)) # check the sample data

sample_names(ps) <- sample_data(ps)$sample_name # set the sample names to the new ones
sample_names(ps)

# filter out OTUs with zero abundance across all samples
ps_filtered <- prune_taxa(taxa_sums(ps) > 0, ps)

