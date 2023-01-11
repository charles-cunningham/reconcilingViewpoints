# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Create GB vectors required for analysis
#

### LOAD LIBRARIES & INSTALL PACKAGES ---------------------------------

library(raster)
library(rgdal)
library(rgeos)
library(magrittr)

### SNIPPET ------------------------------------------------------------

# Load in CRS
bng = CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 
          +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")

# Set base grid (10x10km resolution raster)
base.grid <- raster(ncol=70, nrow=140, 
                    xmn=0, xmx=700000, 
                    ymn=-100000, ymx=1300000,
                    crs = bng)

projection(base.grid) <- bng

# Read in GB polygon from GADM, removing small islands
GBR1 <- raster::getData("GADM", country="GBR", level=1) # Read in
GBR2 <- GBR1[GBR1$NAME_1 != "Northern Ireland",] # Remove NI
GBR3 <- disaggregate(GBR2) # Separate into different polygons
GBR3$area_sqkm <- area(GBR3) / 1000000 # Calculate area
GBR4 <- GBR3[GBR3$area_sqkm > 20,] # Only include islands > 20km
GB.poly <- aggregate(GBR4) # Aggregate back

# Simplify slightly to make rasterisation faster
GB.poly <- gSimplify(GB.poly,  0.005) 
GB.poly <- spTransform(GB.poly, proj4string(base.grid))

# Create raster of hectads with >0 GB land within them
GB1 <- GB.poly %>%
  aggregate(.) %>%
  rasterize(.,base.grid)
GB2 <- GB.poly %>%
  aggregate(.) %>%
  as(., "SpatialLines") %>%
  rasterize(.,base.grid)
GB.rast <- (!is.na(GB1)) + (!is.na(GB2))
GB.rast[values(GB.rast) > 0] <- 1
GB.rast[values(GB.rast) == 0] <- NA
values(GB.rast)[!is.na(values(GB.rast))] <- 0

# Tidy environment
rm(GB1, GB2, GBR1, GBR2, GBR3, GBR4)

