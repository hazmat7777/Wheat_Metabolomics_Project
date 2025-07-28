# install packages
install.packages("terra")    # Core raster GIS data package
install.packages('sf')       # Core vector GIS data package
install.packages('raster')   # Older raster GIS package required by some packages
install.packages('geodata')  # Data downloader

install.packages('sp')        # Older vector GIS package - replaced by sf in most cases
install.packages('rgdal')     # Interface to the Geospatial Data Abstraction Library
install.packages('lwgeom')    # Lightweight geometry engine

install.packages('ggplot2')    # Plotting package
install.packages('gridExtra')  # Extensions to ggplot

# load packages
library(terra)     # core raster GIS package
library(sf)        # core vector GIS package
library(units)     # used for precise unit conversion

library(geodata)   # Download and load functions for core datasets
library(raster)    # older raster GIS package

library(tidyverse)  # data manipulation and visualization
sf_use_s2(FALSE)

# load the data
env_data <- read.csv("../data/GreenMicrobiome_FarmKits_EnvData.csv")

# check the data
glimpse(env_data)

# filter out the rows with no GPS coordinates
env_data_filtered <- env_data %>%
    filter(!is.na(gps_latitude) & !is.na(gps_longitude)) %>%
    select(gps_latitude, gps_longitude, sample_ID)

# convert to sf object
env_data_sf <- st_as_sf(env_data_filtered, coords = c("gps_longitude", "gps_latitude"), crs = 4326)

plot(env_data_sf,
     xlab = "Longitude", ylab = "Latitude",
     main = "Sample Locations")


# plot background map cropped to the UK
world_map <- geodata::world(country = "United Kingdom", resolution = "low") %>%
    st_transform(crs = 4326)

ne_10 <- st_read('../../CMEECourseWork/week09/data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp')

st_agr(ne_10) <- 'constant'

ne_10_uk <- ne_10 %>%
    filter(ADMIN == "United Kingdom") %>%
    st_transform(crs = 4326) %>%
    select(ADMIN) %>%
    st_crop(xmin = -4, ymin = 50, xmax = 2, ymax = 55)

glimpse(ne_10_uk)

# Plot with ggplot2
p1 <- ggplot() +
  geom_sf(data = ne_10_uk, fill = "khaki", color = "black") +  # UK map
  geom_sf(data = env_data_sf, color = "red", size = 2, pch = 4) +  # Sample points
  theme_minimal()

pdf("../results/sampling_locations.pdf", width = 8, height = 6)
print(p1)
dev.off()
