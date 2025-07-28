library(raster)

## clim data

# Load the DTR variable from NetCDF file
dtr_brick <- brick("../data/env_data/clim_data/cru_ts4.09.2021.2024.dtr.dat.nc", varname = "dtr")
    # a brick- each layer is a timestep

# Check the names of the layers (timesteps)
names(dtr_brick) # 4 years 

uk_extent <- extent(-10, 3, 49, 61) 
dtr_uk <- crop(dtr_brick, uk_extent) # crop to just UK

library(sp)

list.files("../data/env_data/clim_data/") # check the files in the directory

# Load and crop the RasterBrick (climate data)
dtr_brick <- brick("../data/env_data/clim_data/cru_ts4.09.2021.2024.dtr.dat.nc", varname = "dtr")
uk_extent <- extent(-10, 3, 49, 61) 
dtr_brick <- crop(dtr_brick, uk_extent) # crop to just UK

# load and filter the GPS data
gps_df <- read.csv("../data/GreenMicrobiome_FarmKits_EnvData.csv")
gps_df <- gps_df %>%
    filter(!is.na(gps_latitude) & !is.na(gps_longitude)) %>%
    as_tibble() %>%
    dplyr::select(gps_latitude, gps_longitude, sample_ID)

# Convert to SpatialPoints
coordinates(gps_df) <- ~gps_longitude + gps_latitude
crs(gps_df) <- crs(dtr_brick)  # Set the coordinate reference system

# Assuming dtr_brick is already loaded as the raster data
timestep <- dtr_brick[[37:42]] # first half of 2024

# Extract climate data at each GPS point for the last timestep
climate_data <- extract(mean(timestep), gps_df)

# Add climate data to the GPS DataFrame
gps_df$diurnal_temprange <- climate_data

# View the updated dataframe
View(gps_df)

summary(timestep)
summary(values(timestep))  # Check for NA values in the raster layer

is.na(timestep)

plot(timestep)
















# Load and crop the RasterBrick (climate data)
dtr_brick <- brick("../data/env_data/clim_data/cru_ts4.09.2021.2024.dtr.dat.nc", varname = "dtr")
uk_extent <- extent(-10, 3, 49, 61)
dtr_brick <- crop(dtr_brick, uk_extent) # crop to just UK

# Load additional climate variables
pre_brick <- brick("../data/env_data/clim_data/cru_ts4.09.2021.2024.pre.dat.nc", varname = "pre")
pre_brick <- crop(pre_brick, uk_extent)

tmn_brick <- brick("../data/env_data/clim_data/cru_ts4.09.2021.2024.tmn.dat.nc", varname = "tmn")
tmn_brick <- crop(tmn_brick, uk_extent)

tmx_brick <- brick("../data/env_data/clim_data/cru_ts4.09.2021.2024.tmx.dat.nc", varname = "tmx")
tmx_brick <- crop(tmx_brick, uk_extent)

wet_brick <- brick("../data/env_data/clim_data/cru_ts4.09.2021.2024.wet.dat.nc", varname = "wet")
wet_brick <- crop(wet_brick, uk_extent)

# load and filter the GPS data
gps_df <- read.csv("../data/GreenMicrobiome_FarmKits_EnvData.csv")
gps_df <- gps_df %>%
 filter(!is.na(gps_latitude) & !is.na(gps_longitude)) %>%
 as_tibble() %>%
 dplyr::select(gps_latitude, gps_longitude, sample_ID)

# Convert to SpatialPoints
coordinates(gps_df) <- ~gps_longitude + gps_latitude
crs(gps_df) <- crs(dtr_brick) # Set the coordinate reference system

# Extract data for first half of 2024 (timesteps 37:42)
timestep <- 37:42

# Extract climate data at each GPS point and calculate means
dtr_data <- extract(mean(dtr_brick[[timestep]]), gps_df)
pre_data <- extract(mean(pre_brick[[timestep]]), gps_df)
tmn_data <- extract(mean(tmn_brick[[timestep]]), gps_df)
tmx_data <- extract(mean(tmx_brick[[timestep]]), gps_df)
wet_data <- extract(mean(wet_brick[[timestep]]), gps_df)

# Add all climate data to the GPS DataFrame
gps_df$diurnal_temp_range <- dtr_data
gps_df$precipitation <- pre_data
gps_df$min_temp <- tmn_data
gps_df$max_temp <- tmx_data
gps_df$monthly_wet_days <- wet_data

# Convert back to regular dataframe if needed
climate_df <- as.data.frame(gps_df)

View(climate_df)

# notes on data
    # dtr: Diurnal Temperature Range (DTR) is the difference between the maximum and minimum temperature in a day, averaged over the first 6 months of 2024.
    # pre: Precipitation (PRE) is the mean monthly rainfall/snowfall (mm) in the first half of 2024
    # tmn: Mean Minimum Temperature (TMN) is the mean of the monthly minimum temperatures in a month.
    # tmx: Mean Maximum Temperature (TMX) is the mean of the monthly maximum temperatures in a month.
    # wet: Monthly Wet Days (WET) is the number of days in a month with precipitation greater than 1mm, averaged over the first half of 2024.







### plotting

# Convert the raster (timestep) to a data frame
raster_df <- as.data.frame(rasterToPoints(timestep), xy = TRUE)
colnames(raster_df) <- c("longitude", "latitude", "value")

# Create a SpatialPointsDataFrame from gps_df (if it's not already)
gps_df <- as.data.frame(gps_df)  
ggplot() +
  # Raster as background
  geom_tile(data = raster_df, aes(x = longitude, y = latitude, fill = value)) +
  scale_fill_viridis_c() +  # You can use any color scale you prefer
  # GPS points on top of the raster
  geom_point(data = gps_df, aes(x = gps_longitude, y = gps_latitude), color = "red", size = 2) +
  coord_equal() +  # Equal scaling for both axes
  labs(title = "GPS Points Overlaid on Raster", x = "Longitude", y = "Latitude") +
  theme_minimal()






## soil geochemistry data

# load the soil data
mg <- raster("../data/env_data/soil_properties/UK_topsoil_data_geotiff_MgO_v1/MgO_v1.tif")
zn <- raster("../data/env_data/soil_properties/UK_topsoil_data_geotiff_Zn_v1/Zn_v1.tif")

# Plot the raster data
plot(mg, main = "Magnesium Oxide (MgO) in UK Topsoil",
    xlab = "Longitude", ylab = "Latitude")
plot(zn, main = "Zinc (Zn) in UK Topsoil",
    xlab = "Longitude", ylab = "Latitude")
# Check the CRS of the rasters

# can stack em
chem_stack <- stack(mg, zn)

# not much point doing this because farms will fertilize the soil