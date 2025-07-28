# load packages
rm(list = ls()) # clear workspace
library(readxl) # for reading xlsx files
library(openxlsx) # for reading xlsx files
library(dplyr)    # Data manipulation
library(readODS)

# load data
df <- read.xlsx("../data/metabolomics/farmkit_peak_areas.xlsx") # this is the full metabolomics dataset I collected from each farmkit using msdial
gt <- read_ods("../data/baiting/baiting_farmkits.ods", sheet = 1) %>%
    as_tibble() %>%
    select(sample_name, gt)

# cleaning baiting data
View(gt)
length(gt$sample_name) # 70 rows in df
length(unique(gt$sample_name)) # 69 uniques... so one set of duplicates (error in collection/ data entry)
duplicates <- gt[duplicated(gt$sample_name) | duplicated(gt$sample_name, fromLast = TRUE), ]
View(duplicates)
    # fk93 had 1 and 0- this is an error in collection. Remove it:
gt <- gt %>%
    filter(sample_name != "93") # remove the fk93 sample

# which gt sample names are not in fk_colnames
fk_colnames <- grep("fk", colnames(df), value = TRUE) #  82 farmkits had mass spec data collected
fk_colnames <- gsub("fk", "", fk_colnames) # remove the fk_ prefix so same as gt$sample_name
gt$sample_name[!gt$sample_name %in% fk_colnames] # 207, 243 are missing
    # 243 is my fault
    # 3 was meant to be 103
    # 207- apparently they baited it but it wasnt sent for mass spec analysis. 

# which farmkits have baiting data
baited_fks <- fk_colnames[fk_colnames %in% gt$sample_name]
baited_fks 

# df of only the sample metabolomics for farmkits with GT data
baited_fk_colnames <- paste0("fk", baited_fks) # re-add the fk_ prefix
df_fk <- df[, c("Metabolite.name", fk_colnames)]
#transpose the GT data frame to have samples as rows and metabolites as columns
gt_t <- gt %>%
    mutate(sample_name = paste0("fk", sample_name)) %>% # re-add the fk_ prefix to match df_fk 
    t() %>%
    as.data.frame()

colnames(gt_t) <- gt_t[1, ] # set the first row as column names
gt_t <- gt_t[-1, ] # remove the first row (now column names)

# add a 'Metabolite.name' column to gt_t
gt_t <- gt_t %>% 
  mutate(Metabolite.name = NA) %>%  # Add a column with the metabolite names
  select(Metabolite.name, everything()) %>%  # Move it to the front
  select(-`fk243`, -`fk207`)  # Remove the 2 farmkits that don't have GT data

# convert character columns to numeric in gt_t
gt_t[] <- lapply(gt_t, function(x) if(is.character(x)) as.numeric(x) else x)

# add the GT data to the bottom of the metabolomics data
fk_metabolom_gt <- bind_rows(df_fk, as.data.frame(gt_t)) # bind the GT data to the metabolomics data

View(tail(fk_metabolom_gt))

saveRDS(fk_metabolom_gt, "../data/fk_metabolomics_gt.RDS") # save the combined data














?order()

library(dplyr)
glimpse(df)

length(unique(df$Formula))

dim(df2) # 225 columns, 37187 samples
length(unique(df$Calc..MW)) # all unique masses

df2<- df[with(df, order(Calc..MW)), ] # ordered by MW

df2$Calc..MW[1]
min(df$Calc..MW) # it worked

