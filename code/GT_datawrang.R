# ADDING GT VARIABLE TO THE SAMPLE DATA OF PHYLOSEQ OBJECT

# install.packages("openxlsx")
library(openxlsx)  # Reading data from Excel files
library(dplyr)    # Data manipulation
# install.packages("readODS")
library(readODS)




# load farmkits data
farmkits = read.xlsx("../data/GreenMicrobiome_FarmKits_EnvData_Harry.xlsx", sheet = 1)

farmkits <- as_tibble(farmkits) %>%
    dplyr::filter(!is.na(gps_latitude) & !is.na(gps_longitude) & !is.na(sample_name)) #%>%
    #dplyr::filter("sample_ID" %in% sample_names(ps_project9))

View(farmkits)

# load the GT data
gt <- read_ods("../data/baiting_farmkits.ods", sheet = 1) %>%
    as_tibble() %>%
    select(sample_name, gt_yn)

gt$sample_name <- as.character(gt$sample_name)

View(gt)

# merge the GT data with farmkits to get the sample names
gt_farmkits <- left_join(farmkits, gt, by = "sample_name")

View(gt_farmkits)
sum(is.na(gt_farmkits$gt)) # 3 farmkits have no GT data

# merge the GT dataframe with the sample data from the phyloseq object
ps_project9 <- readRDS("../data/ps_project9.rds")
view(sample_data(ps_project9))

sample_data(ps_project9)$sample_ID <- rownames(sample_data(ps_project9))

length(sample_data(ps_project9)$sample_ID)  # Should equal nsamples(ps_project9)
nsamples(ps_project9)

# pruning samples to only those with GT data
ps_project11 <- prune_samples(sample_data(ps_project9)$sample_ID %in% gt_farmkits$sample_ID, ps_project9)
ps_project11

View(sample_data(ps_project11))

# Filter gt_farmkits to include only matching sample_IDs
filtered_gt_farmkits <- gt_farmkits %>%
    filter(sample_ID %in% sample_data(ps_project11)$sample_ID)

# Convert sample_data to a data frame
sample_data_df <- as.data.frame(sample_data(ps_project11))

# Perform the left join
sample_data_joined <- left_join(sample_data_df, 
                                filtered_gt_farmkits[, c("sample_ID", "gt")], 
                                by = "sample_ID")

str(sample_data_joined)

# Update the sample_data in the phyloseq object
sample_data(ps_project11)$gt <- sample_data_joined$gt

View(sample_data(ps_project11))

saveRDS(ps_project11, "../data/ps_project11.rds")







# # load GT baiting data
# gt_1 <- read.xlsx("../data/baIting30_04_25_Harry.xlsx", sheet = 1)[, 1:4] %>% # one half of the baiting data
#     as_tibble() %>%
#     dplyr::select(SAMPLE, gt) %>%
#     dplyr::filter(!grepl("[A-Za-z]", SAMPLE))

# # 1 for presence, 0 for absence
# gt_1$gt[!is.na(gt_1$gt)] <- 1
# gt_1$gt[is.na(gt_1$gt)] <- 0 # previously NA meant absence

# View(gt_1)

# # Join the GT baiting data with farmkits data
# gt_farmkits <- farmkits[, c(1,4,6,7)] %>%
#     left_join(gt_1, by = c("sample_name" = "SAMPLE"))

# # Check the joined data
# View(gt_farmkits)

# sum(is.na(gt_farmkits$gt)) # 14 farmkits have no GT data. I only have 30 and some will be filtered bc saequencing didnt work
# View(gt_farmkits[is.na(gt_farmkits$gt),])


### Mass spectrometry random forest

### Community composition random forest

### do the sequenced samples match the env samples?
ps_project2 <- readRDS("../data/ps_project2.rds")

# Extract sample names from the phyloseq object
samples <- sample_names(ps_project2)
print(samples)

gt_farmkits$sample_ID %in% samples # one sample in farmkits not in phyloseq object
samples %in% gt_farmkits$sample_ID # some samples in phyloseq not in farmkits





















samples <- as.vector(farmkits$sample_name)

class(samples)
print(samples)

cleaned_names <- gsub("[A-Za-z]", "", samples) # remove letters
print(cleaned_names)

cleaned_names <- cleaned_names[!grepl("\\.", samples)]
print(cleaned_names) # no dots

length(cleaned_names) # 11

baiting <- read.xlsx("../data/baIting30_04_25.xlsx", sheet = 1)

View(baiting[, 1:4])

baiting_samples <- as.vector(baiting$SAMPLE)

cleaned_bait <- gsub("[A-Za-z]", "", baiting_samples) # remove letters

cleaned_bait <- cleaned_bait[!grepl("\\.", baiting_samples)] # remove dots
print(cleaned_bait)
print(cleaned_names)

cleaned_names %in% cleaned_bait # check if cleaned names are in baiting samples