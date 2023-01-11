# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Calculate PA performance under different viewpoints 
# and integration approaches
#

# LOAD PACKAGES AND SET PARAMETERS---------------------------------------

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
library(tibble)

### GB
source("Functions/SNIPPET_create_GB_shapes.R") 

# Set Data folder
dataDir <- "../Data/"

# Overwrite GB.rast with 'trimmed' GB raster
GB.rast <- raster(paste0(dataDir, "Spatial_datasets/GB.raster.tif"))

# Viewpoints
runs <- c("TRAD", "NEW", "ECON","SOC") 

### READ IN PA AND VIEWPOINTS ----------------------------------------------------

# Read in PA
PA.rast <- raster(paste0(dataDir, "Spatial_datasets/PAs/PA_merge_big.tif"))

# Read in prioritisation runs

# Raster
SP.r.stack <- lapply(runs, function(x) {
  
  paste0(dataDir, "Zonation/", x, "/CAZ/outputfile.CAZ_E.rank.compressed.tif") %>% 
    raster(.)
  
  }) %>% stack(.)

names(SP.r.stack) <- runs

SP.stack = list()

for (j in 1:length(runs)) {
  
  SP.stack_j <- readGDAL(paste0(
    dataDir, "Zonation/",
    runs[j],
    "/CAZ/outputfile.CAZ_E.rank.compressed.tif")) 
  SP.stack[[j]] <- SP.stack_j}

# SIMILARITY OF PA NETWORK AND APPROACHES ----------------------------------------

# Make data frame
PA.compare <- data.frame(values(PA.rast),
                           values(SP.r.stack[[1]]),
                           values(SP.r.stack[[2]]),
                           values(SP.r.stack[[3]]),
                           values(SP.r.stack[[4]]))
   
# Set column names            
colnames(PA.compare) <- c( "PA", runs)

# Remove NA values
PA.compare <- na.omit(PA.compare)

# Correlation of GB PAs with prioritisation runs

# Correlation matrix
cor(PA.compare, method = "spearman") 

# Hence correlation between PA network distribution and different approaches is
PA.cor <- cor(PA.compare, method = "spearman")[1, -1]
PA.cor %>% round(., 3)
