# Observing data, processing data, prepping it for analysis

# Install MSprep
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("MSPrep")

# load data
mass_spec_df <- read.csv("../data/metabolomics/20250327_Soil_Bact_OB.csv")

nrow(mass_spec_df) #37187
head(colnames(mass_spec_df))

# checking how many have formulae

library(dplyr)

ms_formulae_df <- filter(mass_spec_df, Formula != "")
(nrow(mass_spec_df) - nrow(ms_formulae_df))/ nrow(mass_spec_df) * 100
    # getting rid of 8359 rows (22%)
    # let's just go with it for now

# checking how many have names
ms_named_df <- filter(ms_formulae_df, Name != "")

# some seem to have names but no formulae
ms_named_no_formulae <- filter(mass_spec_df, Name != "" & Formula == "")
nrow(ms_named_no_formulae)
View(ms_named_no_formulae)
    # chuck these 87 ambiguous ones


nrow(ms_named_df)
    # only 2525 have names

# how many different compounds are represented by these?
length(unique(ms_named_df$Name))
    # only 1428 diff named compounds

# how many diff formulae are unnamed?
nrow(ms_formulae_df)
length(unique(ms_formulae_df$Formula))
    # also cuts it down a bit but still
    # 19293 different formulae
    # maximum 28828 diff compounds (9000 duplicate formulae could be isomers)

# filter out the unrelated the soybean study
soil_columns <- grep("Soil", colnames(ms_formulae_df), value = TRUE)
soil_columns
# Subset the data frame to include only soil-related columns
ms_formulae_df <- ms_formulae_df[, c("Tags", "Checked", "Name", "Formula", "Annot..Source..Predicted.Compositions", 
                              "Annot..Source..mzCloud.Search", "Annot..Source..mzVault.Search", "Annot..DeltaMass..ppm.", 
                              "Calc..MW", "m.z", "RT..min.", "Area..Max..", "X..mzCloud.Results", 
                              "X..mzVault.Results", "mzCloud.Best.Match", "mzCloud.Best.Match.Confidence", 
                              "mzCloud.Best.Sim..Match", "mzCloud.Best.Tree.Match", "mzVault.Best.Match", 
                              "mzCloud.Library.Match..Autoprocessed", "mzCloud.Library.Match..Reference", 
                              "mzVault.Library.Match..20220905_ALL_GNPS_v2c_with_CompoundMetadata_20240813", 
                              "mzVault.Library.Match..20240730_IMPACT_ALL_v2_with_CompoundMetadata", 
                              "MS2", "Reference.Ion", soil_columns)]

soil_columns

## USE MSPREP
library(MSPrep)

?msPrepare # fn that does it all

# df needs only the peak data for the soil
peak_columns <- grep("Peak.Rating", colnames(ms_formulae_df), value = TRUE)
ms_formulae_prep <- ms_formulae_df[, c("Calc..MW", "RT..min.", peak_columns)]

?sub()

# getting colnames in the right format
df <- ms_formulae_prep
colnames <- colnames(df)[3:length(colnames(df))] # just want to change peak area names
colnames <- substr(colnames, 1, nchar(colnames) - 4) # remove the last 4 letters
split_names <- strsplit(colnames, split = "_") # split to 3 substrings
split_names <- lapply(split_names, function(x) { # Modify the second part (x[2]) of the split
  x[2] <- sub("\\.", "-", x[2]) # had to escape \\ from . (which means any character)
  return(x)
})
new_colnames <- sapply(split_names, function(x) paste(x[3], x[1], x[2], sep = "_")) # rearrange
new_colnames # it worked
colnames(ms_formulae_prep)[3:length(colnames(df))] <- new_colnames # apply changes
head(colnames(ms_formulae_prep)) # done

?msSummarize()

# running msprep
summarizedDF <- msSummarize(ms_formulae_prep,
                            cvMax = 0.50,
                            minPropPresent = 2/3,
                            compVars = c("Calc..MW", "RT..min."),
                            sampleVars = c("subject_id","replicate"),
                            colExtraText = "Soil.raw.._Peak.Rating..20250318_",
                            separator = "-",
                            missingValue = 1 # placeholder. I have loads though (how many? I should check)
                            )

    # i notice that there are lots of NAs, can filter them with minproppresent and change them with missingValue

sum(ms_formulae_prep == 0, na.rm = TRUE)  # Count of zeros
sum(is.na(ms_formulae_prep))  

# full pipeline
prepared_ms <- msPrepare(ms_formulae_prep,
                            cvMax = 0.50,
                            minPropPresent = 2/3,
                            compVars = c("Calc..MW", "RT..min."),
                            sampleVars = c("subject_id","replicate"),
                            colExtraText = "Soil.raw.._Peak.Rating..20250318_",
                            separator = "-",
                            missingValue = NA, # treat NAs as missing. I have MANY NAs
                            imputeMethod = "halfmin", # replacing missing values with half the minimum observed value in that variable
                            nPcs = 3,
                            normalizeMethod = "median",
                            transform = "none"
                            )
View(prepared_ms)

nrow(prepared_ms) # 1239 compounds. But maybe my filtering too strong.. think about imputing?

## other args I could use
    # prepared_ms <- msPrepare(ms_formulae_df,
    #                         minPropPresent = 1/3, # below this the compound will be set to 0
    #                         filterPercent = 0.8,
    #                         missingValue = 0,
    #                         imputeMethod = "halfmin",
    #                         nPcs = 3,
    #                         normalizeMethod = "median",
    #                         transform = "none",
    #                         compVars = c("Mass", "Retention.Time",
    #                                      "Compound.Name"),
    #                         sampleVars = c("subject_id", "replicate"),
    #                         colExtraText = "X",
    #                         separator = "_",
    #                         returnToSE = TRUE)

saveRDS(prepared_ms, "../data/metabolomics/prepared_ms_missing_NAs.RDS")
prepared_ms <- readRDS("../data/metabolomics/prepared_ms_missing_NAs.RDS")

View(prepared_ms)


# bcpa explained:

# Bayesian approach to PCA: BPCA combines an EM approach for PCA with a Bayesian model bpca: Bayesian PCA missing value estimation in pcaMethods: A collection of PCA methods
# dont want to use this if missing values are not random:

# Calculate the proportion of missing values per metabolite
missing_per_metabolite <- apply(ms_formulae_prep[-c(1,2)], 1, function(x) sum(is.na(x))/length(x))

View(missing_per_metabolite)

length(missing_per_metabolite)

# Calculate mean intensity per metabolite (ignoring NAs)
mean_intensity_per_metabolite <- apply(ms_formulae_prep[-c(1,2)], 1, function(x) mean(x, na.rm = TRUE))
length(mean_intensity_per_metabolite)
# Test correlation
cor.test(missing_per_metabolite, mean_intensity_per_metabolite, method = "spearman")
    # barely any correlation, so missing values are not just low intensity
    # means we can use bcpa

# If you know your detection limits, check if missing values occur 
# primarily in the low-intensity range

# Find the minimum observed value for each metabolite
min_observed <- apply(ms_formulae_prep, 1, function(x) min(x, na.rm = TRUE))

# Compare this to the proportion of missing values
plot(log(min_observed + 1), missing_per_metabolite,
     xlab = "Log(Minimum Observed Value + 1)", 
     ylab = "Proportion Missing",
     main = "Missing Pattern vs Detection Threshold")

# If there's a strong negative correlation, it suggests values are 
# missing because they're below detection limit (MNAR)
cor.test(log(min_observed + 1), missing_per_metabolite, method = "spearman")

# Check if there are many zeros that might have been converted to NA
sum(your_data == 0, na.rm = TRUE)

# Look at the distribution of minimum values
hist(min_observed, main = "Distribution of Minimum Observed Values")
summary(min_observed)

## inspect the data

dummy_df <- readRDS("../data/metabolomics/prepared_ms")

# (plotted by accident) Mr vs retention time
plot(x = dummy_df[,1], y = dummy_df[,2], 
    xlab = "compound Mr", ylab = "retention time",
    ylim = c(0,16))
    # why are there so few compounds with medium RT? it's hydrophobic vs hydrophilic, to do with elution method

# mr vs abundance in sample 1
plot(x = dummy_df[,1], y = dummy_df[,3], 
    xlab = "compound Mr", ylab = "abundance in sample 1")

View(dummy_df)