# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Run Zonation for different viewpoints
#

# LOAD PACKAGES -----------------------------------------------------------

# Source library
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)
library(abind)

# Set working directory for Zonation
setwd("../Data")

#Create Zonation folder
dir.create(paste0("Zonation"), showWarnings = TRUE)

# DEFINE LAYER WEIGHTS -------------------------------------------------------

### Choose variables for Zonation run, only need to change here 

# Non-bio layers
non.bio.layers <- c( "carbon", "rec.visits", "flood.ES.flow", "poll.ES.flow",
                     "wilderness", "lands.val", "agri.val")

carbon <- "Spatial_datasets\\Carbon\\Above-ground_carbon\\processed_Carbon.tif"
rec.visits <- "Spatial_datasets\\Recreation_visits_per_ha\\Visit_ha_NP_Euro\\pr_kr_bt_nc\\GB_res_w001001.tif"
flood.ES.flow <- "Spatial_datasets\\Flood_regulation\\FRsupply_data\\flood_reg_flow.tif"
poll.ES.flow <- "Spatial_datasets\\Pollination\\Pollin_Dryad\\poll_flow.tif"
wilderness <- "Spatial_datasets\\Wilderness\\processed_BritishIsles_equal_osgb.tif"
lands.val <- "Spatial_datasets\\Landscape_value\\aesthetic_value.tif"
agri.val <- "Spatial_datasets\\Land_classes\\GB\\ALC.rast.normd.tif"

# Check out non bio layers
to.check <- agri.val
raster(gsub("\\\\", "/", to.check)) %>% plot(.)

###!!!### Assign weights here - need to think about these 
layer.weights <- data.frame(layers = c("Species_dists", non.bio.layers),
                              TRAD = c(1  , 0, 0, 0, 0  , 0.25,0 ,0 ),
                              NEW  = c(1  , 1, 1, 1, 0.5, 0   ,1 ,-0.5 ), 
                              ECON = c(1  , 1, 0, 0, 0.5, 0   ,0 ,-1 ),
                              SOC  = c(0.5, 0, 1, 1, 0  , 0   ,1 ,0 ))

# Add average column for one of the compromise options
layer.weights$MEAN <- rowMeans(layer.weights[,-1])

# Save layer weights for other scripts
save(layer.weights, file = "Zonation/layer_weights.Rda")

# ZONATION RUN FOR ALL -------------------------------------------------------

# Objects for loop
all.approaches <- names(layer.weights)[-1]
rem.rule <- c("CAZ", "ABF") # 1 = CAZ, 2 = ABF , 3 = TBF, 4 = GBF, 5 = Random

# Run loop
for (i in all.approaches) {
  for (j in 1:length(rem.rule)) {

# Create folders
dir.create(paste0("Zonation/", i, "/", rem.rule[j] ), recursive = TRUE)

### SETTINGS FILE

# Create settings file. Load with vmat in settings file - speeds up zonation run
settings.text <- paste0(
  
"[Settings]

removal rule = ", j, # 1 = CAZ, 2 = ABF , 3 = TBF, 4 = GBF, 5 = Random
"\nwarp factor = 1
edge removal = 1
add edge points = 0
use SSI = 0
use planning unit layer = 0
use cost = 0 # 1 
use mask = 0
use groups = 0

mask missing areas = 1     #1=use mask for missing areas (Scotland, islands)
area mask file = Spatial_datasets/GB.raster.tif  #name of mask file, hash out if not needed

use boundary quality penalty = 0
BQP mode = 2
BLP = 0
use tree connectivity = 0
use interactions = 0

annotate name = 1
logit space = 0
treat zero-areas as missing data = 0
z = 0.25
resample species = 0

[Info-gap settings]
Info-gap proportional = 0
use info-gap weights  = 0

[vmat]
save vmat = _save_vmat  ")

# Write settings to file
settings <- write(settings.text, file = paste0("Zonation/",
                                      i,
                                      "/",
                                      rem.rule[j],
                                      "/settings.dat"))
    
### SPECIES AND FEATURE LISTS AND FILES     

# Get species files      
sp.list <- paste0("Species_datasets\\SDMforZonation\\",
                  dir("Species_datasets/SDMforZonation"))
      
# Number of species      
nsp <- NROW(sp.list)

# Add in non-biodiversity feature layers, have to 'mget' on non.bio.layers    
feature.list <- c(sp.list, mget(non.bio.layers) %>%
                    unname(.) %>%
                    unlist(.))
                
nTotal <- NROW(feature.list)

# Create feature list file with associated weights(.spp)  
feature.list.df <- data.frame(weights = c(rep((1/nsp) * layer.weights[1,i], nsp), 
                                          layer.weights[-1,i]),
                              alpha = rep(0, nTotal),
                              BQP = rep(1, nTotal),
                              NQP = rep(1, nTotal),
                              funct = rep(1, nTotal),
                              spFile = feature.list)

# Write table
write.table(feature.list.df,
            file = paste0("Zonation/",
                          i,
                          "/",
                          rem.rule[j],
                          "/FileLst.spp"),
            row.names = FALSE,
            col.names = FALSE)
 
# Output file     
temp.str <- paste0("\\Zonation\\", i, "\\", rem.rule[j])

# Run Zonation from R
system(paste0("c:\\progra~1\\zonation\\bin\\zig4 -r ",
              getwd(), paste0(temp.str, "\\settings.dat "),
              getwd(), paste0(temp.str, "\\FileLst.spp "),
              getwd(), paste0(temp.str, "\\outputfile.txt 0.0 0 1.0 0 ")))

  }
}

# Switch working directory back
setwd("../ViewpointsAnalysis")
