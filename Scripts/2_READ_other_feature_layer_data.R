# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Read in all feature layer data (apart from agricultural
# land classification)
#
# Script Description: 
# - This script reads, and converts to GB hectads, previously downloaded data
# - Species distributions are not read in here as this data has already been
# processed (distribution data not publicly available)

# LOAD PACKAGES -----------------------------------------------------------

# Source library
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(data.table)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)
library(caTools)
library(abind)

# LOAD FUNCTIONS AND SOURCE DATA -------------------------------------------

# Create GB base maps
source("Functions/SNIPPET_create_GB_shapes.R") 

# Overwrite GB.rast with 'trimmed' GB raster
GB.rast <- raster("../Data/Spatial_datasets/GB.raster.tif")

# READ IN MULTI LAYER 2015 EU DATASET ---------------------------------------
# N.B. This data isn't used within analysis

# This dataset is explained in the report
# (https://op.europa.eu/en/publication-detail/-/publication/e6480d73-2a48-4249-8a3e-f94ea635475b/language-en)
# and available here
# (https://data.jrc.ec.europa.eu/dataset/7e3f0681-5967-41f7-ae9b-87f1c3cfac4f)

# Set up file locations
source.file.loc <- "R:/rsrch/cb751/lab/Charles/PhD/Chapter4/Data/Spatial_datasets/2015_ES_services"
shp.file.name <- "10kmgrideu28"
csv.file.name <- "DownscaledData.csv"
save.name <- paste0("GB_res_", shp.file.name)

# Read in shapefile
EU.hectad.poly <- readOGR(source.file.loc,
                          shp.file.name)
# Change projection...
GB.hectad.poly <- spTransform(EU.hectad.poly,
                              bng) %>%
  crop(., extent(base.grid)) %>% # then crop the extent
  .[GB.poly,] # then mask to GB

# Read in table
EU.hectad.data <- read.table(paste0(source.file.loc,"/", csv.file.name),
                             sep = ",",
                             skip = 1,
                             header = TRUE)

# Join data
GB.data <- merge(x = GB.hectad.poly, y = EU.hectad.data,
                 by.x = "CellCode", by.y =  "CELLCODE")

# Save layer  
writeOGR(GB.data,
         layer =  save.name,
         dsn = source.file.loc,
         driver="ESRI Shapefile",
         overwrite_layer = TRUE)

# WILDERNESS ------------------------------------------------------------------

# Set up file locations
source.file.loc <- "R:/rsrch/cb751/lab/Charles/PhD/Chapter4/Data/Spatial_datasets/Wilderness"
file.name <- "BritishIsles_equal_osgb.tif"
save.file.loc <- paste0(source.file.loc, "/", "processed_", file.name)

# Read in raster
wild.data <- raster(paste0(source.file.loc, "/", file.name)) %>%
  projectRaster(.,
                disaggregate(base.grid,10),
                res = 1000,
                crs = bng) %>%
  mask(., GB.poly)
plot(wild.data)

# Change to GB.rast resolution
wild.data <- resample(wild.data, GB.rast) 

# Remove 'NA' hectads, and mask it using GB.rast
wild.data[is.na(wild.data)] <- 0
wild.data <- mask(wild.data, GB.rast)

# Create normalised values
normalised = (values(wild.data) - min(na.omit(values(wild.data))))/
  (max(na.omit(values(wild.data))) - min(na.omit(values(wild.data))))

#Assign to raster
wild.data[] <- normalised

# Save layer
writeRaster(wild.data,
            filename = save.file.loc,
            overwrite=TRUE,
            format = "GTiff")


# POTENTIAL VISITS PER YEAR ---------------------------------------------------

# This is modelled visits per km2 in 2016. Layer available here
# (https://www.sciencedirect.com/science/article/pii/S1617138116300152#fig0015)

# Set up file locations
source.file.loc <- "R:/rsrch/cb751/lab/Charles/PhD/Chapter4/Data/Spatial_datasets/Recreation_visits_per_ha/Visit_ha_NP_Euro/pr_kr_bt_nc"
file.name <- "w001001.adf"
save.file.loc <- paste0(source.file.loc, "/", "GB_res_", file.name)

# Read in raster (EU projection and resolution)
EU.rec.data <- raster(paste0(source.file.loc, "/", file.name))
crs(EU.rec.data) <- crs(EU.hectad.poly)
plot(EU.rec.data)

# Convert to points
EU.rec.points <- as(EU.rec.data, "SpatialPointsDataFrame")

# Change projection
GB.rec.points <- spTransform(EU.rec.points,
                              bng) %>%
  crop(., extent(base.grid)) # then crop the extent

# Convert back to raster
GB.data <- rasterize(GB.rec.points, base.grid, field  = "w001001",
                         fun = mean)

# Remove 'NA' hectads, and mask it using GB.rast
GB.data[is.na(GB.data)] <- 0
GB.data <- mask(GB.data, GB.rast)

# Create normalised values
normalised = (values(GB.data) - min(na.omit(values(GB.data))))/
  (max(na.omit(values(GB.data))) - min(na.omit(values(GB.data))))

#Assign to raster
GB.data[] <- normalised

# Save layer
writeRaster(GB.data,
            filename = save.file.loc,
            overwrite=TRUE,
            format = "GTiff")

# LANDSCAPE VALUE --------------------------------------------------------------------

# This is an proxy of landscape value by taking the average value of several social media platforms
# Available here (https://www.pnas.org/content/early/2016/10/25/1614158113.abstract)
# Calculation needed on raw data.

### Read in data

# Set up file locations
source.file.loc <- "R:/rsrch/cb751/lab/Charles/PhD/Chapter4/Data/Spatial_datasets/Landscape_value"
file.name <- "flickr_5k_m/w001001.adf"
save.file.loc <- paste0(source.file.loc, "/", "aesthetic_value")

# Read in Flickr
flickr.r <- raster(paste0(source.file.loc, "/", file.name)) %>%
  projectRaster(.,
                disaggregate(base.grid,2),
                res = 1000,
                crs = bng) %>%
  mask(., disaggregate(GB.rast,2))
plot(flickr.r)

# Read in Instagram
file.name <- "insta_5k_m/w001001.adf"

# Read in raster, change projection...
insta.r <- raster(paste0(source.file.loc, "/", file.name)) %>%
  projectRaster(.,
                disaggregate(base.grid,2),
                res = 1000,
                crs = bng) %>%
  mask(., disaggregate(GB.rast,2))
plot(insta.r)

# Read in Panoramio
file.name <- "pano_5k_m/w001001.adf"

# Read in raster, change projection...
panom.r <- raster(paste0(source.file.loc, "/", file.name)) %>%
  projectRaster(.,
                disaggregate(base.grid,2),
                res = 1000,
                crs = bng) %>%
  mask(., disaggregate(GB.rast,2))
plot(panom.r)

# Rank values
flickr.r[!is.na(flickr.r)] <- rank(na.omit(flickr.r[]))
insta.r[!is.na(insta.r)] <- rank(na.omit(insta.r[]))
panom.r[!is.na(panom.r)] <- rank(na.omit(panom.r[]))

# Mean rank

lands.val.r <- calc(stack(flickr.r, insta.r, panom.r), mean)
plot(lands.val.r)

### Merge and aggregate

# Map to base.grid
lands.val.big <- resample(lands.val.r, base.grid)

# Remove 'NA' hectads, and mask it using GB.rast
lands.val.big[is.na(lands.val.big)] <- 0
lands.val.big <- mask(lands.val.big, GB.rast)

# Create normalised values
normalised = (values(lands.val.big) - min(na.omit(values(lands.val.big))))/
  (max(na.omit(values(lands.val.big))) - min(na.omit(values(lands.val.big))))

#bAssign to raster
lands.val.big[] <- normalised
plot(lands.val.big)
hist(values(lands.val.big))

# Save layer
writeRaster(lands.val.big,
            filename = save.file.loc,
            overwrite=TRUE,
            format = "GTiff")

# FLOOD REGULATION ES -------------------------------------------------------------------

# This is an estimation of flow of flood regulation ES by weighting the supply by the catchment demand.
# Calculation needed on raw data.

### Read in data

# Read in supply
source.file.loc <- "R:/rsrch/cb751/lab/Charles/PhD/Chapter4/Data/Spatial_datasets/Flood_regulation/FRsupply_data"
file.name <- "supply.asc"
save.file.loc <- paste0(source.file.loc, "/", "flood_reg_flow")

# Read in raster, change projection...
F.supply.r <- raster(paste0(source.file.loc, "/", file.name)) %>%
  projectRaster(.,
                disaggregate(base.grid,10),
                res = 1000,
                crs = bng) %>%
  mask(., GB.poly)

# Read in demand
source.file.loc <- "R:/rsrch/cb751/lab/Charles/PhD/Chapter4/Data/Spatial_datasets/Flood_regulation/FRdemand_data"
file.name <- "demand.asc"

# Read in raster, change projection...
F.demand.r <- raster(paste0(source.file.loc, "/", file.name)) %>%
  projectRaster(.,
                disaggregate(base.grid,10),
                res = 1000,
                crs = bng) %>%
  mask(., GB.poly)

# Resample supply and demand
F.supply.r <- aggregate(F.supply.r, fact = 10, fun = sum) %>%
  resample(., base.grid) %>%
  mask(., GB.rast)
F.demand.r <- aggregate(F.demand.r, fact = 10, fun = sum) %>%
  resample(., base.grid) %>%
  mask(., GB.rast) 

plot(F.supply.r)
plot(F.demand.r)

# Rank supply and demand
F.supply.r[!is.na(F.supply.r)] <- rank(na.omit(F.supply.r[]))
F.demand.r[!is.na(F.demand.r)] <- rank(na.omit(F.demand.r[]))

# Find min
F.flow.r <- min(F.supply.r,F.demand.r )
plot(F.flow.r)

# Remove 'NA' hectads, and mask it using GB.rast
F.flow.r[is.na(F.flow.r)] <- 0
F.flow.r <- mask(F.flow.r, GB.rast)

# Create normalised values
normalised = (values(F.flow.r) - min(na.omit(values(F.flow.r))))/
  (max(na.omit(values(F.flow.r))) - min(na.omit(values(F.flow.r))))

#Assign to raster
F.flow.r[] <- normalised
plot(F.flow.r)
hist(values(F.flow.r))

# Save layer
writeRaster(F.flow.r,
            filename = save.file.loc,
            overwrite=TRUE,
            format = "GTiff")


# CARBON ------------------------------------------------------------------------

### Read in above-ground carbon (tonnes per hectare)

# Set up file locations
source.file.loc <- "R:/rsrch/cb751/lab/Charles/PhD/Chapter4/Data/Spatial_datasets/Carbon/Above-ground_carbon"
file.name <- "Carbon"
save.file.loc <- paste0(source.file.loc, "/", "processed_", file.name)

# Read in shapefile
ag.carbon.data <- readOGR(source.file.loc,
                          file.name)

# Convert to raster (1km resolution)
ag.carbon.r <- rasterFromXYZ(data.frame(x = coordinates(ag.carbon.data)[,1],
                                        y = coordinates(ag.carbon.data)[,2],
                                        Carbon = ag.carbon.data$Carbon),
                             crs = bng)
plot(ag.carbon.r)

### Read in below-ground carbon (ktC km-2 kilotons per km)

# Set up file locations
source.file.loc <- "R:/rsrch/cb751/lab/Charles/PhD/Chapter4/Data/Spatial_datasets/Carbon/Below-ground_carbon"
file.name <- "gbstock100"

# Read in shapefile
bg.carbon.data <- readOGR(source.file.loc,
                          file.name)

# Calculation to sum 0-30cm, and 30-100cm, to find 0-100cm carbon
bg.carbon.data$totalCarbon <- bg.carbon.data$MAP30 + bg.carbon.data$MAP100

# Calculation convert from kilotons per km^2 to tonnes per hectare
bg.carbon.data$totalCarbon <- (bg.carbon.data$totalCarbon * 1000) / 100

# Convert from km to m scale, and get centroid
bg.carbon.data$EASTING <- (bg.carbon.data$EASTING * 1000) + 500
bg.carbon.data$NORTHING <- (bg.carbon.data$NORTHING * 1000) + 500

# Convert to raster (1km resolution)
bg.carbon.r <- rasterFromXYZ(data.frame(x = bg.carbon.data$EASTING,
                                        y = bg.carbon.data$NORTHING,
                                        totalCarbon = bg.carbon.data$totalCarbon),
                             crs = bng)
plot(bg.carbon.r)

### Merge and aggregate

# Sum above and below ground carbon
total.carbon.r <- bg.carbon.r + ag.carbon.r

# Map to base.grid
total.carbon.r <- resample(total.carbon.r, base.grid)

# Remove 'NA' hectads, and mask it using GB.rast
total.carbon.r[is.na(total.carbon.r)] <- 0
total.carbon.r <- mask(total.carbon.r, GB.rast)

# Create normalised values
normalised = (values(total.carbon.r) - min(na.omit(values(total.carbon.r))))/
  (max(na.omit(values(total.carbon.r))) - min(na.omit(values(total.carbon.r))))

#Assign to raster
total.carbon.r[] <- normalised
plot(total.carbon.r)

# Save layer
writeRaster(total.carbon.r,
            filename = save.file.loc,
            overwrite=TRUE,
            format = "GTiff")

# ES FLOW FOR CROP POLLINATION -----------------------------------------------------------

# This is an estimation of flow of pollination ES by weighting the supply by the demand.
# Available here (https://www.sciencedirect.com/science/article/pii/S1470160X13002768?via%3Dihub)
# Calculation needed on raw data.

### Read in data

# Read in supply ( area percentage pollinator habitat per km2 grid cell (%) )
source.file.loc <- "R:/rsrch/cb751/lab/Charles/PhD/Chapter4/Data/Spatial_datasets/Pollination/Pollin_Dryad"
file.name <- "3b_visitprob.tif"
save.file.loc <- paste0(source.file.loc, "/", "poll_flow")

# Read in raster, change projection...
P.supply.r <- raster(paste0(source.file.loc, "/", file.name)) %>%
  projectRaster(.,
                disaggregate(base.grid,10),
                res = 1000,
                crs = bng) %>%
  mask(., GB.poly)

# Read in demand (Coverage pollinator dependent crops 
# weighted by pollinator dependency level, normalized to 0-10,000 scale.)
file.name <- "2b_poldephotspot.tif"

# Read in raster, change projection...
P.demand.r <- raster(paste0(source.file.loc, "/", file.name)) %>%
  projectRaster(.,
                disaggregate(base.grid,10),
                res = 1000,
                crs = bng) %>%
  mask(., GB.poly)

# Resample supply and demand
P.supply.r <- aggregate(P.supply.r, fact = 10, fun = sum) %>%
  resample(., base.grid)
P.demand.r <- aggregate(P.demand.r, fact = 10, fun = sum) %>%
  resample(., base.grid)

plot(P.supply.r)
plot(P.demand.r)

# Rank supply and demand
P.supply.r[!is.na(P.supply.r)] <- rank(na.omit(P.supply.r[]))
P.demand.r[!is.na(P.demand.r)] <- rank(na.omit(P.demand.r[]))

# Find min
P.flow.r <- min(P.supply.r,P.demand.r )
plot(P.flow.r)

# Remove 'NA' hectads, and mask it using GB.rast
P.flow.r[is.na(P.flow.r)] <- 0
P.flow.r <- mask(P.flow.r, GB.rast)

# Create normalised values
normalised = (values(P.flow.r) - min(na.omit(values(P.flow.r))))/
  (max(na.omit(values(P.flow.r))) - min(na.omit(values(P.flow.r))))

#Assign to raster
P.flow.r[] <- normalised

plot(P.flow.r)
hist(values(P.flow.r))

# Save layer
writeRaster(P.flow.r,
            filename = save.file.loc,
            overwrite=TRUE,
            format = "GTiff")

# AGRICULTURAL VALUE -----------------------------------------------------------

# Already calculated in previous script
