# install packages
# install.packages("terra")    # Core raster GIS data package
# install.packages('sf')       # Core vector GIS data package
# install.packages('raster')   # Older raster GIS package required by some packages
# install.packages('geodata')  # Data downloader

# install.packages('sp')        # Older vector GIS package - replaced by sf in most cases
# install.packages('rgdal')     # Interface to the Geospatial Data Abstraction Library
# install.packages('lwgeom')    # Lightweight geometry engine

# install.packages('ggplot2')    # Plotting package
# install.packages('gridExtra')  # Extensions to ggplot

# load packages
library(terra)     # core raster GIS package
library(sf)        # core vector GIS package
library(units)     # used for precise unit conversion

library(geodata)   # Download and load functions for core datasets
library(raster)    # older raster GIS package

library(tidyverse)  # data manipulation and visualization
sf_use_s2(FALSE)
library(readxl)    # for reading Excel files
library(readODS)   # for reading ODS files
# load the data
env_data <- read_xlsx("../data/env_data/GreenMicrobiome_FarmKits_EnvData_Harry.xlsx")
env_data <- read_xlsx("../data/env_data/sample_field_specific_environmental_data.xlsx")
# check the data

# filter out the rows with no GPS coordinates
env_data_filtered <- env_data %>%
    filter(!is.na(gps_latitude) & !is.na(gps_longitude)) %>%
    filter(project == "farmer_kit") %>%
    select(gps_latitude, gps_longitude, sample_name)

View(env_data_filtered)

# add in GT data here
# gt <- read_ods("../data/env_data/gt_oldnames.ods")

# View(gt)
# View(env_data_filtered)

# colnames(gt) <- c("sample_name", "gt")
gt <- readRDS("../data/metabolomics/fk_metabolomics_gt_logged.RDS")

gt$sample_name <- rownames(gt)
tail(colnames(gt))
gt <- gt %>%
    select(Farmkit, gt) %>%
    mutate(gt = as.character(gt)) %>%
    filter(!is.na(gt)) %>%
    filter(gt != "")

View(gt)
colnames(gt) <- c("sample_name", "gt")

# merge them
env_data_filtered <- env_data_filtered %>%
    left_join(gt, by = "sample_name") %>%
    mutate(gt = as.factor(gt)) %>%
    select(sample_name, gps_latitude, gps_longitude, gt)

env_data_filtered <- env_data_filtered %>%
    filter(!is.na(gt)) %>%
    filter(gt != "")
    
View(env_data_filtered)

write.csv(env_data_filtered, "../data/env_data/filtered_env_data.xlsx")



# convert to sf object
env_data_sf <- st_as_sf(env_data_filtered, coords = c("gps_longitude", "gps_latitude"), crs = 4326)

plot(env_data_sf,
     xlab = "Longitude", ylab = "Latitude",
     main = "Sample Locations")


# # plot background map cropped to the UK
# world_map <- geodata::world(country = "United Kingdom", resolution = "low") %>%
#     st_transform(crs = 4326)

ne_10 <- st_read('../../CMEECourseWork/week09/data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp')

st_agr(ne_10) <- 'constant'

ne_10_uk <- ne_10 %>%
    filter(ADMIN == "United Kingdom") %>%
    st_transform(crs = 4326) %>%
    dplyr::select(ADMIN) %>%
    st_crop(xmin = -4, ymin = 50, xmax = 2, ymax = 55)

glimpse(ne_10_uk)

str(env_data_sf)

# Add jitter to the sf object
env_data_sf_jittered <- env_data_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(
    X_jittered = X + runif(nrow(.), -0.1, 0.1),  # Adjust jitter amount as needed
    Y_jittered = Y + runif(nrow(.), -0.1, 0.1)
  ) %>%
  st_as_sf(coords = c("X_jittered", "Y_jittered"), crs = st_crs(env_data_sf)) %>%
  bind_cols(st_drop_geometry(env_data_sf))

# Plot with ggplot2
p1 <- ggplot() +
  geom_sf(data = ne_10_uk, fill = "khaki", color = "black") +  # UK map
  geom_sf(data = env_data_sf_jittered, aes(color = gt), size = 3, shape = 17, alpha = 0.7) +  # Sample points
  theme_minimal()

pdf("../results/sampling_locations2.pdf", width = 8, height = 6)
print(p1)

dev.off()

