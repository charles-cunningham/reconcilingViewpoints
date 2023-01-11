# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Read in protected area data
#

# LOAD PACKAGES AND SET PARAMETERS---------------------------------------

# Source library
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(plyr)
library(dplyr)
library(magrittr)
library(reshape2)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)
library(parallel)
library(gdalUtils)

# SET GB SPATIAL OBJECTS AND PARAMETERS --------------------------

source("Functions/SNIPPET_create_GB_shapes.R") 

# Set Data folder
dataDir <- "../Data/"

# Save GB.poly
GB.poly.toWrite <- as(GB.poly, "SpatialPolygonsDataFrame")
names(GB.poly.toWrite) <- "GB.poly"

writeOGR(GB.poly.toWrite,
         layer =  "GB.poly",
         dsn = "Spatial_datasets",
         driver="ESRI Shapefile",
         overwrite_layer = TRUE)

#Set disaggregation factor (10k / 400 = 25m resolution)
dis.factor <- 400

# READ IN PA .SHP --------------------------------------------------------
# N.B. Ignore warning message: doesn't change projection, just p4s string

# SSSI
SSSI.e <- readOGR("Spatial_datasets/PAs/PA_shps/SSSI_england",
                  "Sites_of_Special_Scientific_Interest__England____Natural_England",
                  p4s = as.character(bng))
SSSI.s <- readOGR("Spatial_datasets/PAs/PA_shps/SSSI_scotland",
                  "SSSI_SCOTLAND",
                  p4s = as.character(bng))
SSSI.w <- readOGR("Spatial_datasets/PAs/PA_shps/SSSI_wales",
                  "SSSI_June2020",
                  p4s = as.character(bng))

# NNR
NNR.e <- readOGR("Spatial_datasets/PAs/PA_shps/NNR_england",
                 "National_Nature_Reserves___Natural_England",
                 p4s = as.character(bng))
NNR.s <- readOGR("Spatial_datasets/PAs/PA_shps/NNR_scotland",
                 "NNR_SCOTLAND",
                 p4s = as.character(bng))
NNR.w <- readOGR("Spatial_datasets/PAs/PA_shps/NNR_wales",
                 "NRW_NNRPolygon",
                 p4s = as.character(bng))

# MERGE -------------------------------------------------------

# Join
PA.list <- list(SSSI.e,SSSI.s,SSSI.w,
                NNR.e, NNR.s, NNR.w)

# Merge all the PAs supplied (this is not a union - polygons will overlap)
PA.merge <- Reduce(bind, PA.list)
PA.merge <- gBuffer( PA.merge, width=0, byid=TRUE )

### Need to clip to PA Terrestrial .SHP
PA.merge <- raster::crop(PA.merge, base.grid)
PA.merge <- raster::intersect( PA.merge, GB.poly )

# Write to .shp (get a SHAPE_AREA warning message, measurement is accurate after check) 
writeOGR(PA.merge,
         layer =  "PA_merge",
         dsn = "Spatial_datasets/PAs/Merged_PAs",
         driver="ESRI Shapefile",
         overwrite_layer = TRUE)

# RASTERIZE ---------------------------------------------------------------------

# Make fine-scale raster, need to rename
gdal_rasterize(src_datasource = paste0("Spatial_datasets/PAs/Merged_PAs/PA_merge.shp"),
               dst_filename = paste0("Spatial_datasets/PAs/Merged_PAs/PA_merge_small.tif"),
               burn = 1,
               a_nodata = 0,
               tr = c(10000 / dis.factor, 10000 / dis.factor),
               tap = TRUE,
               output_Raster = TRUE
)

# PROCESS RASTER ----------------------------------------------------------------
  
# Disaggregate and format objects
base.grid.small <- disaggregate(base.grid, fact = dis.factor)
GB.rast.small <- disaggregate(GB.rast, fact = dis.factor)
GB.rast.small <- extend(GB.rast.small, extent(base.grid))

# PROCESS SMALL RASTER ----------------------------------------------------------

PA.small <- raster(paste0("Spatial_datasets/PAs/Merged_PAs/PA_merge_small.tif"))

#Make sure no outlying coastal PAs
PA.small <- crop(PA.small, extent(base.grid))

# Make sure all extents the same
PA.small <- extend(PA.small, extent(base.grid))

#Take out marine PAs by cropping to GB
PA.small <- mask(PA.small, GB.rast.small)

# Write
writeRaster(
  PA.small,
  paste0("Spatial_datasets/PAs/PA_merge_small.tif"),
  format = "GTiff",
  overwrite = TRUE
)

# CREATE HECTAD RASTER ------------------------------------------------------

PA.small <- raster(paste0("Spatial_datasets/PAs/PA_merge_small.tif"))

# Aggregate
PA.big <- aggregate(PA.small, fun = sum, fact = dis.factor)

# to get proportion need to divide by dis.factor^2
PA.big <- PA.big / (dis.factor ^ 2)

# Mask it using base.grid
PA.big <- mask(PA.big, GB.rast)

#Write rasters
PA.raster <- list(PA.big, GB.rast) %>%
  do.call(merge, .)

writeRaster(
  PA.raster,
  paste0("Spatial_datasets/PAs/PA_merge_big.tif"),
  format = "GTiff",
  overwrite = TRUE
) 