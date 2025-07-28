### do the samples cluster by geographic location?

# load, prep data
ps_project7 <- readRDS("../data/metabarcoding/ps_project7.rds")
ps_project8 <- readRDS("../data/metabarcoding/ps_project8.rds")
ps_phylum <- tax_glom(ps_project8, taxrank = "phylum")
taxa_names(ps_phylum) <- tax_table(ps_phylum)[, "phylum"] #rename taxa to phylum names

# Remove samples without GPS coordinates
sample_data_df <- as.data.frame(sample_data(ps_phylum)) # Convert sample_data to a data frame
ps_phylum_gps <- prune_samples(!is.na(sample_data_df$gps_latitude), ps_phylum)
nrow(sample_data(ps_phylum_gps)) # 40 samples with GPS
View(sample_data(ps_phylum_gps))

## making a distance decay plot
library(vegan)
library(ggplot2)
install.packages("geosphere") # for geographic distance calculations
library(geosphere)

# geographic distance matrix
gps_coords <- sample_data(ps_phylum_gps)[, c("gps_latitude", "gps_longitude")]
View(gps_coords)
geo_dist_m <- distm(x = gps_coords) # from geosphere
rownames(geo_dist_m) <- sample_names(ps_phylum_gps)
colnames(geo_dist_m) <- sample_names(ps_phylum_gps)
View(geo_dist_m) # values in metres

# phylum composition distance matrix
phylum_dist_m <- phyloseq::distance(ps_phylum_gps, method = "bray", type = "samples")
View(as.matrix(phylum_dist_m))

# Convert distance matrices to vectors of pairwise distances
geo_dist_vec <- as.vector(geo_dist_m[lower.tri(geo_dist_m)]) # extracts lower triangle, no zeroes or duplicates
phylum_dist_vec <- as.vector(as.matrix(phylum_dist_m)[lower.tri(as.matrix(phylum_dist_m))])

# Plot the relationship
plot(geo_dist_vec, phylum_dist_vec,
     xlab = "Geographic Distance (meters)",
     ylab = "Phylum Composition Distance (Bray-Curtis)",
     main = "Distance Decay Relationship")

## distance decay plot for class level

# Summarize the data at the class level
ps_class <- tax_glom(ps_project8, taxrank = "class")

# filter to samples with GPS coordinates
ps_class_gps <- prune_samples(!is.na(sample_data_df$gps_latitude), ps_class)

# phylum composition distance matrix
class_dist_m <- phyloseq::distance(ps_class_gps, method = "bray", type = "samples")
View(as.matrix(phylum_dist_m))

# Convert distance matrices to vectors of pairwise distances
class_dist_vec <- as.vector(as.matrix(class_dist_m)[lower.tri(as.matrix(phylum_dist_m))])

# Plot the relationship
plot(geo_dist_vec, class_dist_vec,
     xlab = "Geographic Distance (meters)",
     ylab = "Phylum Composition Distance (Bray-Curtis)",
     main = "Distance Decay Relationship")

# same thing for family level

# Summarize the data at the family level
ps_family <- tax_glom(ps_project8, taxrank = "family")

# filter to samples with GPS coordinates
ps_family_gps <- prune_samples(!is.na(sample_data_df$gps_latitude), ps_family)
# family composition distance matrix
family_dist_m <- phyloseq::distance(ps_family_gps, method = "bray", type = "samples")
View(as.matrix(family_dist_m))
# Convert distance matrices to vectors of pairwise distances
family_dist_vec <- as.vector(as.matrix(family_dist_m)[lower.tri(as.matrix(family_dist_m))])
# Plot the relationship
plot(geo_dist_vec, family_dist_vec,
     xlab = "Geographic Distance (meters)",
     ylab = "Family Composition Distance (Bray-Curtis)",
     main = "Distance Decay Relationship")

# same thing for genus level
# Summarize the data at the genus level
ps_genus <- tax_glom(ps_project8, taxrank = "genus")
# filter to samples with GPS coordinates
ps_genus_gps <- prune_samples(!is.na(sample_data_df$gps_latitude), ps_genus)
# genus composition distance matrix
genus_dist_m <- phyloseq::distance(ps_genus_gps, method = "bray", type = "samples")
# Convert distance matrices to vectors of pairwise distances
genus_dist_vec <- as.vector(as.matrix(genus_dist_m)[lower.tri(as.matrix(genus_dist_m))])
# Plot the relationship
plot(geo_dist_vec, genus_dist_vec,
     xlab = "Geographic Distance (meters)",
     ylab = "Genus Composition Distance (Bray-Curtis)",
     main = "Distance Decay Relationship", )

# Fit a linear model
lobf_model <- lm(genus_dist_vec ~ geo_dist_vec)

# Add the line of best fit to the plot
abline(lobf_model, col = "blue", lwd = 2)

### is anything significant?
cor.test(geo_dist_vec, phylum_dist_vec, method = "spearman")
cor.test(geo_dist_vec, class_dist_vec, method = "spearman")
cor.test(geo_dist_vec, family_dist_vec, method = "spearman")
cor.test(geo_dist_vec, genus_dist_vec, method = "spearman")
summary(lobf_model) # slope is positive, but not significant
    # no sig correlation
    # expect a positive correlation- the further apart the samples, the more different they are
