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
    select(sample_name, gt_yn)

# cleaning baiting data
View(gt)
length(gt$sample_name) # 77 rows in df
length(unique(gt$sample_name)) # 77 farmkits
duplicates <- gt[duplicated(gt$sample_name) | duplicated(gt$sample_name, fromLast = TRUE), ]
duplicates # no duplicates

colnames(df)
# which gt sample names are not in fk_colnames
fk_colnames <- grep("fk", colnames(df), value = TRUE) 
length(fk_colnames) #  82 farmkits had mass spec data collected
gt$sample_name[!gt$sample_name %in% fk_colnames] # 94, 206, 207, 243 are missing
    # 243 is my fault
    # 207- apparently they baited it but it wasnt sent for mass spec analysis. 

# which farmkits have baiting data
baited_fks <- fk_colnames[fk_colnames %in% gt$sample_name]
baited_fks # 73 

# df of only the sample metabolomics for farmkits with GT data
df_fk <- df[, c("Metabolite.name", fk_colnames)]
#transpose the GT data frame to have samples as rows and metabolites as columns

gt_t <- gt %>%
    t() %>%
    as.data.frame()

colnames(gt_t) <- gt_t[1, ] # set the first row as column names
gt_t <- gt_t[-1, ] # remove the first row (now column names)

# add a 'Metabolite.name' column to gt_t
gt_t <- gt_t %>% 
  mutate(Metabolite.name = NA) %>%  # Add a column with the metabolite names
  select(Metabolite.name, everything()) %>%  # Move it to the front
  select(-`fk243`, -`fk207`)  # Remove the 2 farmkits that don't have GT data

View(gt_t)

# convert character columns to numeric in gt_t
gt_t[] <- lapply(gt_t, function(x) if(is.character(x)) as.numeric(x) else x)

# add the GT data to the bottom of the metabolomics data
fk_metabolom_gt <- bind_rows(df_fk, as.data.frame(gt_t)) # bind the GT data to the metabolomics data

# filter columns for which there are NA values
fk_metabolom_gt_noNA <- fk_metabolom_gt[, c(TRUE, colSums(is.na(fk_metabolom_gt[, -1])) == 0)]

View(fk_metabolom_gt_noNA)
View(fk_metabolom_gt[nrow(fk_metabolom_gt),])


saveRDS(fk_metabolom_gt_noNA, "../data/fk_metabolomics_gt_noNA.RDS") # save the combined data














?order()

library(dplyr)
glimpse(df)

length(unique(df$Formula))

dim(df2) # 225 columns, 37187 samples
length(unique(df$Calc..MW)) # all unique masses

df2<- df[with(df, order(Calc..MW)), ] # ordered by MW

df2$Calc..MW[1]
min(df$Calc..MW) # it worked

