# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Evaluate approaches to viewpoint integration
#

# LOAD PACKAGES AND SET PARAMETERS---------------------------------------

# Source library
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(ggplot2)
library(magrittr)
library(rgdal)
library(rgeos)
library(raster)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggpubr)

### GB
source("Functions/SNIPPET_create_GB_shapes.R") 

# Set Data folder
dataDir <- "../Data/"

# Overwrite GB.rast with 'trimmed' GB raster
GB.rast <- raster(paste0(dataDir, "Spatial_datasets/GB.raster.tif"))

#List the quantiles we're going to extract accumulation data for
top.percent <- c(5,10,17,30)
quants <- 1 - (top.percent/100)
quants <- sort(quants, decreasing = TRUE)

# Define specifics
my.theme <- theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.title = element_text(size=8),
                  legend.text = element_text(size=8),
                  legend.key = element_blank(),
                  strip.text.x = element_blank(),
                  strip.background = element_blank(),
                  panel.background = element_blank(),
                  legend.justification=c(1,0.8),
                  legend.position=c(1,0.8),
                  plot.margin= unit(rep(0, 4), "lines"))

cutoff <- "mod.10.ModandFiltRaw"
taxa <- c("Merge_all")
runs <- c("TRAD", "NEW", "ECON", "SOC")

# Removal rule
rem.rule <- c("CAZ", "ABF") # 1 = CAZ, 2 = ABF , 3 = TBF, 4 = GBF, 5 = Random

# Colour schemes
pal <- brewer.pal(n = length(quants), name = "RdBu")
rev.pal <- rev(pal)

# New colour scheme
reds <- c("#fdc0aa", "#FB6A4A", "#DE2D26", "#8e0d12")
reds.rev <- rev(reds)

# DEFINE VARIABLES HERE -----------------------------------------------------

# Bio layer
bio <- raster( paste0(dataDir, "Species_datasets/wrscr.tif"))

# Non-bio layers
non.bio.layers <- c( "carbon", "rec.visits", "flood.ES.flow", "poll.ES.flow",
                     "wilderness", "lands.val", "agri.val")

carbon <- raster(paste0(dataDir, "Spatial_datasets/Carbon/Above-ground_carbon/processed_Carbon.tif"))
rec.visits <- raster(paste0(dataDir, "Spatial_datasets/Recreation_visits_per_ha/Visit_ha_NP_Euro/pr_kr_bt_nc/GB_res_w001001.tif"))
flood.ES.flow <- raster(paste0(dataDir, "Spatial_datasets/Flood_regulation/FRsupply_data/flood_reg_flow.tif"))
poll.ES.flow <- raster(paste0(dataDir, "Spatial_datasets/Pollination/Pollin_Dryad/poll_flow.tif"))
wilderness <- raster(paste0(dataDir, "Spatial_datasets/Wilderness/processed_BritishIsles_equal_osgb.tif"))
lands.val <- raster(paste0(dataDir, "Spatial_datasets/Landscape_value/aesthetic_value.tif"))
agri.val <- raster(paste0(dataDir, "Spatial_datasets/Land_classes/GB/ALC.rast.normd.tif"))

# List the layers of interest
layers <- c("bio", non.bio.layers)
layer.shorts <- c("B", "C", "R", "F", "P", "W", "L", "A*")

### CREATE SPDF OF RUNS (VIEWPOINTS)

for( i in rem.rule) {
  
  SP.stack = list()
  
  for (j in 1:length(runs)) {
    
    SP.stack_j <- readGDAL(paste0(
      dataDir, "Zonation/",
      runs[j], "/", i, 
      "/outputfile.", i, "_E.rank.compressed.tif")) 
    SP.stack[[j]] <- SP.stack_j}
  
  assign(x = paste0(i, "_stack"), value = SP.stack)
  
}

#### LOOP THROUGH CAZ AND ABF HERE
for ( x in 1:length(rem.rule)) { # Runs use the same variable names and they are reassigned to rule specific names at the end of the loop

# Assign which removal rule we will look at
  SP.stack <- get(paste0(rem.rule[x], "_stack"))

# INCLUSIVE ----------------------------------------------------------------
# N.B. 'vote-counting' method
 
### Calculate values
  
# Sum viewpoints across cells
SP.agg <- c(raster(SP.stack[[1]]),
            raster(SP.stack[[2]]),
            raster(SP.stack[[3]]), 
            raster(SP.stack[[4]])) %>% 
  stack(.) %>% 
  calc(., sum)

### Convert from aggregate values to quantile values

### Reset quants based off aggregate sgdf
agg.quants <- vector(length = length(quants))
for(i in 1:length(quants))
  {
  agg.quants[i] <- quantile(na.omit(SP.agg[]), probs = quants[i])
}

ranks <- SP.agg[]

# Change highest quantile
ranks[SP.agg[] >= agg.quants[1] & 
        !is.na(SP.agg[])] <- quants[1]

# Change other quantiles
for (j in 2:length(agg.quants)){
  
  ranks[SP.agg[] >= agg.quants[j] &
          SP.agg[] < agg.quants[j-1] &
          !is.na(SP.agg[])] <- quants[j] }

# Change other values to NA
ranks[SP.agg[] < min(agg.quants)] <- NA 

SP.agg[] <- ranks

# Convert to data frame
plotData <- data.frame(SP.agg[])
plotData[,c("x","y")] <- coordinates(SP.agg)

### Reshape the data for ggplot
plotData = reshape2::melt(plotData, id.vars = c('x','y'))
levels(plotData$variable) <- runs

# Crop to non-NA values only
plotData <- subset(plotData, x>= 75000 &
                     x <= 645000 &
                     y>= 25000 &
                     y <= 1205000)

# Plot
INC.map <- ggplot(aes(x = x, y = y), data = plotData) +
  geom_tile(aes(fill = as.factor(value))) + facet_wrap(~ variable) +
  scale_fill_manual(na.value = "white",
                    values = rev.pal) +
  coord_equal() +
  ggtitle("") +
  my.theme + 
  theme(legend.position = "none" )
INC.map <- INC.map + geom_polygon(data = GB.poly, aes(long,lat,group=group), color = "black", fill = NA)
INC.map

# PLURALIST --------------------------------------------------------
# N.B. method addressing equity

### Calculate values

# Read in layer weightings for each run (remove MEAN as not needed here)
load(paste0(dataDir, "Zonation/layer_weights.Rda"))
layer.weights$MEAN <- NULL

# Carry out PCA
layers.pca <- prcomp(layer.weights[,-1])
layers.pca
summary(layers.pca) # 1st PC is 61% variance explained

# How many PCs to include?
# Ensure all viewpoints have their largest PC component included

lapply(runs, function(y) { # Loop through runs
  
  #For each run...
  layers.pca$rotation[y,] %>% # Extract PC components
    abs(.) %>% # Take absolute value
    which.max(.) %>% # Find which PC is the maximum
    as.numeric(.) # Convert to number
  
}) %>% unlist(.) %>% # Unlist
  max(.) # Find max of all runs, i.e. how many PCs you should include

### Weight ranks with PC eigenvectors, and sum (dot product)

SP.agg.PC1 <- c(raster(SP.stack[[1]]) * layers.pca$rotation[,"PC1"][[1]],
                raster(SP.stack[[2]]) * layers.pca$rotation[,"PC1"][[2]],
                raster(SP.stack[[3]]) * layers.pca$rotation[,"PC1"][[3]],
                raster(SP.stack[[4]]) * layers.pca$rotation[,"PC1"][[4]]) %>%
  stack(.) %>% 
  calc(., sum)

SP.agg.PC2 <- c(raster(SP.stack[[1]]) * layers.pca$rotation[,"PC2"][[1]],
                raster(SP.stack[[2]]) * layers.pca$rotation[,"PC2"][[2]],
                raster(SP.stack[[3]]) * layers.pca$rotation[,"PC2"][[3]],
                raster(SP.stack[[4]]) * layers.pca$rotation[,"PC2"][[4]]) %>%
  stack(.) %>% 
  calc(., sum)

SP.agg.PC3 <- c(raster(SP.stack[[1]]) * layers.pca$rotation[,"PC3"][[1]],
                raster(SP.stack[[2]]) * layers.pca$rotation[,"PC3"][[2]],
                raster(SP.stack[[3]]) * layers.pca$rotation[,"PC3"][[3]],
                raster(SP.stack[[4]]) * layers.pca$rotation[,"PC3"][[4]]) %>%
  stack(.) %>% 
  calc(., sum)

# Sum absolute dot product values
SP.agg.PC <- abs(SP.agg.PC1) + abs(SP.agg.PC2) + abs(SP.agg.PC3)

### Convert from aggregate values to quantile values

### Reset quants based off aggregate sgdf
agg.quants <- vector(length = length(quants))
for(i in 1:length(quants))
{
  agg.quants[i] <- quantile(na.omit(SP.agg.PC[]), probs = quants[i])
}

### Change from rank values to quantile values

ranks <- SP.agg.PC[]

# Change highest quantile
ranks[SP.agg.PC[] >= agg.quants[1] & 
        !is.na(SP.agg.PC[])] <- quants[1]

# Change other quantiles
for (j in 2:length(agg.quants)){
  
  ranks[SP.agg.PC[] >= agg.quants[j] &
          SP.agg.PC[] < agg.quants[j-1] &
          !is.na(SP.agg.PC[])] <- quants[j] }

# Change other values to NA
ranks[SP.agg.PC[] < min(agg.quants)] <- NA 

SP.agg.PC[] <- ranks

# Convert to data frame
plotData <- data.frame(SP.agg.PC[])
plotData[,c("x","y")] <- coordinates(SP.agg.PC)

### Reshape the data for ggplot
plotData = reshape2::melt(plotData, id.vars = c('x','y'))
levels(plotData$variable) <- runs

# Crop to non-NA values only
plotData <- subset(plotData, x>= 75000 &
                     x <= 645000 &
                     y>= 25000 &
                     y <= 1205000)

# Plot
PLUR.map <- ggplot(aes(x = x, y = y), data = plotData) +
  geom_tile(aes(fill = as.factor(value))) + facet_wrap(~ variable) +
  scale_fill_manual(na.value = "white",
                    values = rev.pal) +
  coord_equal() +
  ggtitle("") +
  my.theme +
  theme(legend.position = "none" )
PLUR.map <- PLUR.map + geom_polygon(data = GB.poly, aes(long,lat,group=group), color = "black", fill = NA)
PLUR.map

# PRE-PRIORITISATION MEAN (MEAN) --------------------------------------------
# N.B. Aggregating viewpoints prior to viewpoint prioritisation

### Calculate values

SP.agg.mean <- raster(paste0(
  dataDir, "Zonation/", "MEAN/", rem.rule[x],
  "/outputfile.", rem.rule[x], "_E.rank.compressed.tif"))

### Convert from aggregate values to quantile values

# Reset quants based off aggregate sgdf
agg.quants <- vector(length = length(quants))
for(i in 1:length(quants))
{
  agg.quants[i] <- quantile(na.omit(SP.agg.mean[]), probs = quants[i])
}

ranks <- SP.agg.mean[]

# Change highest quantile
ranks[SP.agg.mean[] >= agg.quants[1] & 
        !is.na(SP.agg.mean[])] <- quants[1]

# Change other quantiles
for (j in 2:length(agg.quants)){
  
  ranks[SP.agg.mean[] >= agg.quants[j] &
          SP.agg.mean[] < agg.quants[j-1] &
          !is.na(SP.agg.mean[])] <- quants[j] }

# Change other values to NA
ranks[SP.agg.mean[] < min(agg.quants)] <- NA 

SP.agg.mean[] <- ranks

# Convert to data frame
plotData <- data.frame(SP.agg.mean[])
plotData[,c("x","y")] <- coordinates(SP.agg.mean)

### Reshape the data for ggplot
plotData = reshape2::melt(plotData, id.vars = c('x','y'))
levels(plotData$variable) <- runs

# Crop to non-NA values only
plotData <- subset(plotData, x>= 75000 &
                     x <= 645000 &
                     y>= 25000 &
                     y <= 1205000)

# Plot
mean.map <- ggplot(aes(x = x, y = y), data = plotData) +
  geom_tile(aes(fill = as.factor(value))) + facet_wrap(~ variable) +
  scale_fill_manual(na.value = "white",
                    values = rev.pal) +
  coord_equal() +
  ggtitle("") +
  my.theme + 
  theme(legend.position = "none" )
mean.map <- mean.map + geom_polygon(data = GB.poly, aes(long,lat,group=group), color = "black", fill = NA)
mean.map

# RANK PRIORITISATION (RANK) -------------------------------------------
# N.B. Carrying out additional prioritisation on the viewpoint rankings

### Calculate values

# Set working directory for Zonation
setwd(dataDir)

# Create folder
dir.create(paste0("Zonation/PRIORITIES_PRIORITISATION/", rem.rule[x]), recursive = TRUE)

feature.list = c(paste0("Zonation/", runs, "/", rem.rule[x],
    "/outputfile.", rem.rule[x], "_E.rank.compressed.tif")) %>%
  gsub("/", "\\\\", .)

### SETTINGS FILE

# Create settings file. Load with vmat in settings file - speeds up zonation run massively
settings.text <- paste0(
  
  "[Settings]

removal rule = ", x, # 1 = CAZ, 2 = ABF , 3 = TBF, 4 = GBF, 5 = Random
"\nwarp factor = 1
edge removal = 1
add edge points = 0
use SSI = 0
use planning unit layer = 0
use cost = 0 
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
settings <- write(settings.text, file = paste0("Zonation/PRIORITIES_PRIORITISATION/", rem.rule[x], "/settings.dat"))

### SPECIES AND FEATURE LISTS AND FILES     

# Add in non-biodiversity feature layers, have to 'mget' on non.bio.layers    
nTotal <- NROW(feature.list)

# To write Biodiversity feature list file (.spp)  
#-- Change weights to add different criteria
feature.list.df <- data.frame(weights = rep(1, nTotal),
                              alpha = rep(0, nTotal),
                              BQP = rep(1, nTotal),
                              NQP = rep(1, nTotal),
                              funct = rep(1, nTotal),
                              spFile = feature.list)

# Write table
write.table(feature.list.df,
            file = paste0("Zonation/PRIORITIES_PRIORITISATION/", rem.rule[x], "/FileLst.spp"),
            row.names = FALSE,
            col.names = FALSE)

# Set output folder
temp.str <- paste0("\\Zonation\\PRIORITIES_PRIORITISATION\\", rem.rule[x])

# Call Zonation from R
system(paste0("c:\\progra~1\\zonation\\bin\\zig4 -r ",
              getwd(), paste0(temp.str, "\\settings.dat "),
              getwd(), paste0(temp.str, "\\FileLst.spp "),
              getwd(), paste0(temp.str, "\\outputfile.txt 0.0 0 1.0 0 ")))

# Load Zonation output
SP.agg.rank <- raster(paste0(
  "Zonation/", "PRIORITIES_PRIORITISATION/", rem.rule[x],
  "/outputfile.", rem.rule[x], "_E.rank.compressed.tif"))

### Convert from aggregate values to quantile values

# Reset quants based off aggregate sgdf
agg.quants <- vector(length = length(quants))
for(i in 1:length(quants))
{
  agg.quants[i] <- quantile(na.omit(SP.agg.rank[]), probs = quants[i])
}

ranks <- SP.agg.rank[]

# Change highest quantile
ranks[SP.agg.rank[] >= agg.quants[1] & 
        !is.na(SP.agg.rank[])] <- quants[1]

# Change other quantiles
for (j in 2:length(agg.quants)){
  
  ranks[SP.agg.rank[] >= agg.quants[j] &
          SP.agg.rank[] < agg.quants[j-1] &
          !is.na(SP.agg.rank[])] <- quants[j] }

# Change other values to NA
ranks[SP.agg.rank[] < min(agg.quants)] <- NA 

SP.agg.rank[] <- ranks

# Convert to data frame
plotData <- data.frame(SP.agg.rank[])
plotData[,c("x","y")] <- coordinates(SP.agg.rank)

### Reshape the data for ggplot
plotData = reshape2::melt(plotData, id.vars = c('x','y'))
levels(plotData$variable) <- runs

# Crop to non-NA values only
plotData <- subset(plotData, x>= 75000 &
                     x <= 645000 &
                     y>= 25000 &
                     y <= 1205000)

# Plot
rank.map <- ggplot(aes(x = x, y = y), data = plotData) +
  geom_tile(aes(fill = as.factor(value))) + facet_wrap(~ variable) +
  scale_fill_manual(na.value = "white",
                    values = rev.pal) +
  coord_equal() +
  ggtitle("") +
  my.theme + 
  theme(legend.position = "none" )
rank.map <- rank.map + geom_polygon(data = GB.poly, aes(long,lat,group=group), color = "black", fill = NA)
rank.map

# Change working directory back
setwd("../ViewpointsAnalysis")

# CREATE QUANTILE ACCUMULATION TABLE -------------------------------

# Stack rasters
aggs <- c("INCLUSIVE", "PLURALIST", "MEAN", "RANK")
SP.stack <- list(SP.agg, SP.agg.PC, SP.agg.mean, SP.agg.rank)

# Create layers plus blank vector
plot.layers <- c( "bio", "carbon", "rec.visits", "flood.ES.flow", "poll.ES.flow",
                  "wilderness", "lands.val", "blank", "agri.val")
plot.layer.shorts <- c("B", "C", "R", "F", "P", "W", "L", " ", "A*")

# Create the table which we are going to populate

num.rows <-  length(plot.layers) * length(aggs) * length(quants)

acc.agg <- data.frame("aggs" = rep(aggs, each = length(plot.layers) * length(quants)),
                      "layers" = rep(plot.layers, times = length(aggs), each = length(quants)),
                      "quants" = rep(quants, num.rows/length(quants)),
                      "total.value" = rep(NA,num.rows),
                      "added.value" = rep(NA,num.rows),
                      "pos.or.neg" = rep(NA,num.rows))

# State levels
acc.agg$aggs <- factor(acc.agg$aggs,
                        levels = aggs)
acc.agg$layers <- factor(acc.agg$layers,
                          levels = plot.layers)
acc.agg$quants <- as.factor(acc.agg$quants)

# Create a loop to populate the table (apart from added.value)

iter = 1

for (i in 1 : length(aggs)) {
  for (j in 1 : length(plot.layers)) {
    for (k in 1 : length(quants)) {
      
      ### Read in quantile mask
 
      # Get raster from run i
      r <- SP.stack[[i]]
     
      # Extract only data above quantile
      r[r<quants[k]] <- NA
      r[!is.na(r)] <- 1
      
      ### Positive or negative resource
      
      if(plot.layers[j] != "agri.val") {
        acc.agg$pos.or.neg[iter] <- "Positive"
      }
      
      if(plot.layers[j] == "agri.val") {
        acc.agg$pos.or.neg[iter] <- "Negative"
      }
      
      # Only carry out if not the blank layer
      if(plot.layers[j] != "blank") {
      
      ### Mask layer raster and sum
      
      # If layer is positive resource, mask to included amount
      if(plot.layers[j] != "agri.val") {
        acc.value <- raster::mask(get(plot.layers[j]), mask = r, inverse = FALSE) %>% cellStats(., sum)
      }
      # If layer is negative resource, mask to excluded amount
      if(plot.layers[j] == "agri.val") {
        acc.value <- raster::mask(get(plot.layers[j]), mask = r, inverse = TRUE) %>% cellStats(., sum)
      }
      
      ### Find total quantile proportion covered
      
      # Divide by total layer 'resource' to get proportion
      acc.prop <- acc.value / cellStats(get(plot.layers[j]), sum)
      
      # Add to total prop data table
      acc.agg$total.value[iter] <- acc.prop
      }
      iter = iter + 1
    }
  }
}

### Add "added.value" column ( N.B. just for plot as "agri.val" quants are switched to plot)

iter = 1

for (i in 1 : length(aggs)) {
  for (j in 1 : length(plot.layers)) {
    for (k in 1 : length(quants)) {
      
      if (plot.layers[j] != "agri.val") {
        
        if (k == 1) { acc.added.prop <- acc.agg$total.value[iter]}
        
        if (k>1) { # if not the first quantile...
          
          acc.added.prop <- acc.agg$total.value[iter] - #take away the sum of all other higher quantiles
            acc.agg$total.value[acc.agg$aggs == aggs[i] &
                                   acc.agg$layers == plot.layers[j] &
                                   acc.agg$quants == quants[k-1]]
        }}
      
      if (plot.layers[j] == "agri.val"){
        
        if (k == 1) {
          
          acc.added.prop <- acc.agg$total.value[acc.agg$aggs == aggs[i] &
                                                  acc.agg$layers == plot.layers[j] &
                                                  acc.agg$quants == quants[length(quants)]]
        }
        
        if (k  > 1) { 
          
          acc.added.prop <- acc.agg$total.value[acc.agg$aggs == aggs[i] &
                                                  acc.agg$layers == plot.layers[j] &
                                                  acc.agg$quants == quants[length(quants) - k + 1]] - 
            acc.agg$total.value[acc.agg$aggs == aggs[i] &
                                  acc.agg$layers == plot.layers[j] &
                                  acc.agg$quants == quants[length(quants) - k + 2]]
        }}
      
      # Assign to acc.table
      acc.agg$added.value[iter] <- acc.added.prop
      
      iter = iter + 1
    }
  }
}

# Save for later efficiency analysis
save(acc.agg, file = paste0(dataDir, "Spatial_datasets/", rem.rule[x], "_aggregate_resource_coverage.Rdata"))

# MAPS --------------------------------------------------------------

for (i in 1:length(aggs)) {
 
### Create plot data
plotData <- data.frame(SP.stack[[i]][])
plotData[,c("x","y")] <- coordinates(SP.stack[[i]])
### Reshape the data for ggplot
plotData = reshape2::melt(plotData, id.vars = c('x','y'))
levels(plotData$variable) <- aggs[i]

# Crop to non-NA values only
plotData <- subset(plotData, x>= 75000 &
                     x <= 645000 &
                     y>= 25000 &
                     y <= 1205000)

### Map
Zon.map <- ggplot(aes(x = x, y = y), data = plotData) +
  geom_tile(aes(fill = as.factor(value))) + facet_wrap(~ variable) +
  scale_fill_manual(na.value = "white",
                    values = reds) +
  coord_equal() +
  my.theme  + theme(legend.position = "none",
                    plot.background = element_blank(),
                    panel.grid = element_blank(),
                    plot.title = element_blank(),
                    strip.text.x = element_blank())

Zon.map <- Zon.map + geom_polygon(data = GB.poly, aes(long,lat,group=group), color = "black", fill = NA)

# BARPLOT --------------------------------------------------------------

ES.barplot <- ggplot(acc.agg[acc.agg$aggs == aggs[i],], aes( y= added.value, 
                                                         x= factor(layers, labels = plot.layer.shorts), 
                                                         fill = interaction(quants, pos.or.neg))) + ## change here
  geom_bar(position="stack", stat="identity", col = "black") +
  scale_fill_manual( values = c(reds.rev, reds)) +
  facet_wrap(~aggs) +
  theme_minimal() +
  xlab("Resource") +
  ylab("Proportion resource included") +
  scale_y_continuous( breaks = seq(0,1,0.1),
                      minor_breaks = seq(0,1,0.05),
                      limits = c(0,1)) +
  theme(axis.title.x = element_text(size=30,margin = margin(t=-10,b=10)),
        axis.title.y = element_text(size=30),
        axis.text.x = element_text(size=30, margin = margin(t=-15,b=10), colour = "black"),
        axis.text.y = element_text(size=26, colour = "black"),
        plot.title = element_blank(),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "black", size= 0.2, linetype="dashed"),
        panel.grid.minor.y = element_line(colour = "black", size= 0.2, linetype="dashed"),
        panel.ontop = TRUE,
        legend.position = "none")

# COMBINE MAP AND BARPLOT ------------------------------------------

ES.inset <- ggplotGrob(ES.barplot)
main.map <- ggplotGrob(Zon.map)

# Create empty grob to fill
empty_base <- ggplot() + geom_rect(mapping=aes(xmin = extent(GB.poly)@xmin,
                                               xmax = 800000,
                                               ymin = extent(GB.poly)@ymin,
                                               ymax = extent(GB.poly)@ymax),
                                   fill = NA) + 
  my.theme +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank()) + 
  theme(plot.margin = margin(-6,-5,-6,-5, "lines"))

# Fill grob
inlaid.plot <- empty_base + 
  annotation_custom(grob = main.map,
                    xmin = extent(GB.poly)@xmin,
                    xmax = extent(GB.poly)@xmax,
                    ymin = extent(GB.poly)@ymin,
                    ymax = extent(GB.poly)@ymax) +
  annotation_custom(grob = ES.inset,
                    xmin = 440000,
                    xmax = 800000,
                    ymin = 500000,
                    ymax = 1100000)


assign(paste0("inlaid.plot.", aggs[i]), inlaid.plot, envir = .GlobalEnv)

} # warning messages okay, just removing NA values from blank which is what we want!

# Dummy plot for legend.grob
dummy.legend <- ggplot(acc.agg, aes( y= added.value, x= layers, fill = interaction(quants, pos.or.neg))) + 
  geom_bar(position="stack", stat="identity", col = "black") +
  scale_fill_manual( values = c(reds.rev, reds),
                     breaks = c("0.95.Positive","0.9.Positive", "0.83.Positive", "0.7.Positive" ),
                     labels = c("5%", "10%", "17%", "30%"),
                     name = "Total land coverage ") +
  facet_wrap(~aggs) + 
  theme(legend.position="bottom",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40))

### Plot 2 main integration approaches with 'ggarrange' for CAZ only
if (rem.rule[x] == "CAZ") {

main.inlaid.plots <- ggarrange(inlaid.plot.INCLUSIVE, inlaid.plot.PLURALIST,
                              ncol=2, nrow=1, 
                              labels = aggs[1:2], label.x = 0.55, label.y = 0.4, 
                              font.label = c(size = 40), hjust = c(-0.49,-0.45),
                              common.legend = TRUE, legend= "bottom", legend.grob = get_legend(dummy.legend)) + 
  theme(plot.margin = margin(0,4,1,1, "lines"),
        plot.background = element_rect(fill = "white", colour = "white"))

main.inlaid.plots

# Save
ggsave(main.inlaid.plots, filename = paste0("../Writing/Figures/", rem.rule[x], "_integrated_inlaid_main.png"),
       width = 56, height = 40, dpi = 400, units = "cm")
ggsave(main.inlaid.plots, filename = paste0("../Writing/Figures/", rem.rule[x], "_integrated_inlaid_main.pdf"),
       width = 56, height = 40, dpi = 600, units = "cm")
}

### Plot all 4 integration approaches with 'ggarrange' for supplement
all.inlaid.plots <- ggarrange(inlaid.plot.INCLUSIVE, inlaid.plot.PLURALIST,inlaid.plot.MEAN, inlaid.plot.RANK,
                              ncol=2, nrow=2,
                              labels = aggs, label.x = 0.5, label.y = 0.4,
                              font.label = c(size = 40), hjust = c(-0.58,-0.52, -1.41, -1.46),
                              common.legend = TRUE, legend= "bottom", legend.grob = get_legend(dummy.legend))+ 
  theme(plot.margin = margin(2,5,2,1, "lines"),
        plot.background = element_rect(fill = "white", colour = "white"))

all.inlaid.plots

# Save
ggsave(all.inlaid.plots, filename = paste0("../Writing/Figures/", rem.rule[x], "_integrated_inlaid_supplemental.png"),
       width = 50, height = 65, dpi = 400, units = "cm")
ggsave(all.inlaid.plots, filename = paste0("../Writing/Figures/", rem.rule[x], "_integrated_inlaid_supplemental.pdf"),
       width = 50, height = 65, dpi = 600, units = "cm")

### Reassign all variables to removal rule specific names to look through if needed
assign(x = paste0(rem.rule[x], ".SP.agg"), value = SP.agg)
assign(x = paste0(rem.rule[x], "INC.map"), value = INC.map)
assign(x = paste0(rem.rule[x], "SP.agg.PC"), value = SP.agg.PC)
assign(x = paste0(rem.rule[x], "PLUR.map"), value = PLUR.map)
assign(x = paste0(rem.rule[x], "SP.agg.mean"), value = SP.agg.mean)
assign(x = paste0(rem.rule[x], "mean.map"), value = mean.map)
assign(x = paste0(rem.rule[x], "SP.agg.rank"), value = SP.agg.rank)
assign(x = paste0(rem.rule[x], "rank.map"), value = rank.map)
assign(x = paste0(rem.rule[x], "acc.agg"), value = acc.agg)
assign(x = paste0(rem.rule[x], "all.inlaid.plots"), value = all.inlaid.plots)

}
