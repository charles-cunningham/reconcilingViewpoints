# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Read in agricultural classification data
#

# LOAD PACKAGES -----------------------------------------------------------

# Source library
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(magrittr)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)
library(dplyr)
library(gdalUtils)

# LOAD FUNCTIONS AND SET PARAMETERS -------------------------------------

# Create GB base maps
source("Functions/SNIPPET_create_GB_shapes.R") 

# Disaggregation factor
dis.factor <- 200

# Land area mask proportion
land.cutoff <- 0.5

# CREATE TRIMMED RASTER FOR REST OF ANALYSIS -----------------------------

#Create land area raster
land.area <- rasterize(GB.poly, GB.rast, getCover = TRUE)

# Create land area mask (> land.cutoff)
land.mask <- land.area
land.mask[land.mask < land.cutoff] <- NA
land.mask[land.mask >= land.cutoff] <- 1

# Write GB raster
writeRaster(land.mask, "../Data/Spatial_datasets/GB.raster.tif",
            overwrite=TRUE)

# Overwrite GB.rast with 'trimmed' GB raster
GB.rast <- raster("../Data/Spatial_datasets/GB.raster.tif")

# READ IN ALCs -------------------------------------------------------------

# England: https://data.gov.uk/dataset/952421ec-da63-4569-817d-4d6399df40a1/provisional-agricultural-land-classification-alc
# Wales: http://lle.gov.wales/catalogue/item/PredictiveAgriculturalLandClassificationALCMap2/?lang=en
# Scotland: https://www.hutton.ac.uk/learning/natural-resource-datasets/landcover/land-capability-agriculture

# Read in agricultural land classification shapefiles
ALC.e <- readOGR("../Data/Spatial_datasets/Land_classes/England",
                  "Agricultural_Land_Classification_Provisional_England",
                  p4s = projection(bng) )
ALC.w <- readOGR("../Data/Spatial_datasets/Land_classes/Wales",
                  "wg_alc_predictive_agricultural_land_classification_2Polygon",
                  p4s = projection(bng) )
ALC.s <- readOGR("../Data/Spatial_datasets/Land_classes/Scotland",
                  "LCA_250K",
                  p4s = projection(bng) )

# CREATE INTEROPERABLE CLASS COLUMN ------------------------------------------

# Read in .csv with interoperable codes
# (edit this to change classes rather than script

ALC.lookup <- read.csv(file = "../Data/Spatial_datasets/Land_classes/Land_class_conversion.csv",
                       header = TRUE, na.strings = "")

# England

ALC.e.comp <- merge(ALC.e, unique(ALC.lookup[,c("England.Code","Interoperable.Code")]),
                    by.x = "alc_grade", by.y = "England.Code")
ALC.e.comp <- ALC.e.comp["Interoperable.Code"]

# Wales

ALC.w.comp <- merge(ALC.w, unique(ALC.lookup[,c("Wales.Code","Interoperable.Code")]),
                    by.x = "predictive", by.y = "Wales.Code")
ALC.w.comp <- ALC.w.comp["Interoperable.Code"]

# Scotland

ALC.s.comp <- merge(ALC.s, unique(ALC.lookup[,c("Scotland.Code","Interoperable.Code")]),
                    by.x = "LCCODE", by.y = "Scotland.Code")
ALC.s.comp <- ALC.s.comp["Interoperable.Code"]

# COMBINE NATION DATA --------------------------------------------------------

# Fix geometry errors 
ALC.e.comp <- gBuffer(ALC.e.comp, byid=TRUE, width=0)
ALC.w.comp <- gBuffer(ALC.w.comp, byid=TRUE, width=0)
ALC.s.comp <- gBuffer(ALC.s.comp, byid=TRUE, width=0)

# Remove overlapping areas
ALC.w.comp <- ALC.w.comp - ALC.e.comp
ALC.s.comp <- ALC.s.comp - ALC.e.comp

# List nation spdfs
ALC.list <- list(ALC.e.comp,ALC.w.comp,ALC.s.comp)

# Combine
ALC.GB <- Reduce(bind, ALC.list)

# Rename value name (must be <10 characters)
names(ALC.GB) <- "Code"

# Clip to PA Terrestrial .SHP
ALC.GB <- raster::intersect( ALC.GB, GB.poly )

# N.B. The Wales spdf  causes local crash when using writeOGR() function after above code
# This can be avoided by saving to an .Rdata file then loading it back after a crash

# Save to Rdata file
save(ALC.GB, file = "../Data/Spatial_datasets/Land_classes/ALC_GB_polygon.Rdata")

##### CHECKPOINT #####

# WRITE TO SHAPEFILE -----------------------------------------------------------

# Load .Rdata file
load(file = "../Data/Spatial_datasets/Land_classes/ALC_GB_polygon.Rdata")

# Save GB polygon to shapefile
writeOGR(ALC.GB,
         layer = "GB.Agri.Land.Classes",
         dsn = "../Data/Spatial_datasets/Land_classes/GB",
         driver="ESRI Shapefile",
         overwrite_layer = TRUE)

##### CHECKPOINT #####

# CREATE FINE-SCALE RASTER -------------------------------------------------------

# Use 200 for disaggregation factor
# i.e. 50m resolution (Wales is a polygonised 50m raster)

gdal_rasterize(src_datasource = paste0("../Data/Spatial_datasets/Land_classes/GB/GB.Agri.Land.Classes.shp"), # Input
               dst_filename = paste0("../Data/Spatial_datasets/Land_classes/GB/GB.ALC.small.tif"), # Output
               a = "Code", # Value column
               a_nodata = 999, # Need this parameter to separate 0s and NAs
               tr = c(10000/dis.factor, 10000/dis.factor), # Resolution
               tap = TRUE) # Keep alignments tidy

##### CHECKPOINT #####

# PROCESS RASTER -------------------------------------------------------------------

### Read in
ALC.base.small <- raster(paste0("../Data/Spatial_datasets/Land_classes/GB/GB.ALC.small.tif"))

# Align with other rasters
ALC.base.small <- resample(ALC.base.small, disaggregate(base.grid, dis.factor))

# Check
plot(ALC.base.small)

ALC.scaled <- ALC.base.small

# Create normalised values
normalised = (values(ALC.scaled) - min(na.omit(values(ALC.scaled))))/
   (max(na.omit(values(ALC.scaled)))-min(na.omit(values(ALC.scaled))))

# Assign to raster
ALC.scaled[] <- normalised

# Invert so higher values are greater agri value
ALC.scaled[] <- 1- ALC.scaled[]

# Check
plot(ALC.scaled)

### Aggregate 
ALC.base.big <- aggregate(ALC.scaled, fun = mean, fact = dis.factor)

# Mask using GB.rast
# (Note '1' is equivalent to '0' for positive-weighted rasters)
ALC.base.big <- mask(ALC.base.big, GB.rast)

plot(ALC.base.big)

# Write rasters
writeRaster(ALC.base.big,
            "../Data/Spatial_datasets/Land_classes/GB/ALC.rast.normd.tif", 
            format = "GTiff", overwrite=TRUE)
