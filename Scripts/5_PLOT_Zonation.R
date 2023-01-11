# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Plot feature layers and Zonation viewpoint 
# prioritisation outputs
#

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
library(gridExtra)
library(grid)
library(RColorBrewer)
library(tibble)
library(ggpubr)

### GB
source("Functions/SNIPPET_create_GB_shapes.R") 

# Set Data folder
dataDir <- "../Data/"

# Overwrite GB.rast with >50% land cover GB raster
GB.rast <- raster(paste0(dataDir, "Spatial_datasets/GB.raster.tif"))

# List the quantiles we're going to extract accumulation data for
top.percent <- c(5,10,17,30)
quants <- 1 - (top.percent/100)
quants <- sort(quants, decreasing = TRUE)

# Plot specifics
my.theme <- theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.grid = element_blank(),
                  title = element_blank(),
                  legend.key.width = unit(1,"cm"),
                  legend.title = element_text(size=42),
                  legend.text = element_text(size=42),
                  strip.text.x = element_blank(),
                  strip.background = element_blank(),
                  panel.background = element_blank())

pal <- brewer.pal(n = length(quants), name = "RdBu")

rev.pal <- brewer.pal(n = length(quants), name = "RdBu") %>% 
  rev(.)

runs <- c("TRAD", "NEW", "ECON", "SOC")

# CALCULATE WRSCR -------------------------------------------------
# N.B. WRSCR provided as require permissions to use raw distributions

# ### Calculate 'biological value' (mean species distribution protected)
# 
# # Make brick of distributions
# sp.brick <- list.files("../Data/Species_datasets/SDMforZonation", full.names = TRUE) %>%
#   stack(.) %>%
#   brick(.)
# sp.wrscr <- mask(sp.brick, GB.rast)
# 
# # Divide distribution by total extent
# for(i in 1:nlayers(sp.wrscr)) {
# 
#   sp.sum <- cellStats(sp.wrscr[[i]], sum)
# 
#   sp.wrscr[[i]] <- sp.wrscr[[i]] / sp.sum
# 
# }
# 
# # Remove weird NA values
# sp.wrscr[is.na(sp.wrscr)] <- 0
# 
# # Find average
# bio <- calc(sp.wrscr, mean)
# 
# #Normalise bio for plots
# normalised = (values(bio) - min(na.omit(values(bio))))/
#  (max(na.omit(values(bio))) - min(na.omit(values(bio))))
# 
# # Assign to raster
# bio[] <- normalised
# bio <- mask(bio, GB.rast)
# 
# # Save
# writeRaster(bio,
#             filename = "../Data/Species_datasets/wrscr.tif",
#             overwrite=TRUE,
#             format = "GTiff")

# DEFINE VARIABLES HERE ----------------------------------------------

# Bio layer (WRSCR)
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

# List the all feature layers
layers <- c("bio", non.bio.layers)
layer.shorts <- c("B", "C", "R", "F", "P", "W", "L", "A*")

# Removal rule
rem.rule <- c("CAZ", "ABF") # 1 = CAZ, 2 = ABF , 3 = TBF, 4 = GBF, 5 = Random

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

# PLOT FEATURE LAYERS -----------------------------------------------------

# ES data to plot from rasters
ES.df <- data.frame(as(bio,"SpatialGridDataFrame")@data,
                   as(carbon,"SpatialGridDataFrame")@data,
                   as(rec.visits,"SpatialGridDataFrame")@data,
                   as(flood.ES.flow,"SpatialGridDataFrame")@data,
                   as(poll.ES.flow,"SpatialGridDataFrame")@data,
                   as(wilderness,"SpatialGridDataFrame")@data,
                   as(lands.val,"SpatialGridDataFrame")@data,
                   as(agri.val,"SpatialGridDataFrame")@data)
names(ES.df) <- layers
ES.df[,c("x","y")] <- coordinates(as(wilderness,"SpatialGridDataFrame"))

# Crop to non-NA values only
ES.df <- subset(ES.df, x>= 75000 &
                  x <= 645000 &
                  y>= 25000 &
                  y <= 1205000)

### Reshape the data for ggplot
plotData = reshape2::melt(ES.df, id.vars = c('x','y'))

# Plot
ES.maps <- ggplot(aes(x = x, y = y), data = plotData) +
  geom_tile(aes(fill = value)) + facet_wrap(~ variable) +
  scale_fill_gradient2("Rescaled\nvalue",
                       na.value = "white",
                       low = "white",
                       high = "#5C821A", # base green "forestgreen"
                       guide = guide_colourbar(ticks = FALSE,
                                               title.position = "top",
                                               title.vjust = 0.2,
                                               title.hjust = 0.2,
                                               scale_colour = "black",
                                               barwidth = unit(1, "cm"),
                                               barheight = unit(8, "cm"),
                                               label.hjust=0),
                       labels=c("0", "0.5", "1"),
                       breaks=c(0.05,0.5,0.95)) +
  coord_equal() +
  my.theme +
  theme(legend.position=c(0.86,0.18),
        panel.spacing.x= unit(-1.5, "lines"),
        panel.spacing.y= unit(-3, "lines"))

# Add GB polygon
ES.maps <- ES.maps + geom_polygon(data = GB.poly, aes(long,lat,group=group), color = "black", fill = NA)

# Add in layer.shorts text
dat_text <- data.frame(label = layer.shorts,
                       variable   = layers)
dat_text$variable = factor(dat_text$variable, levels = layers)

ES.maps <-  ES.maps + geom_text(data = dat_text,
                                mapping = aes(x = 200000,
                                              y = 1100000,
                                              label = label),
                                inherit.aes = TRUE,
                                size = 16)

# Save (N.B. for all plots: one low res .png and one high res .pdf)
ggsave(plot = ES.maps,
       filename = "../Writing/Figures/ES_maps.png",
       width = 25, height = 50, dpi = 300, units = "cm")
ggsave(plot = ES.maps,
       filename = "../Writing/Figures/ES_maps.pdf",
       width = 25, height = 50, dpi = 600, units = "cm")

# CREATE QUANTILE ACCUMULATION TABLE ---------------------------------------------

# Loop through CAZ and ABF
for (i in rem.rule) {

SP.stack <- get(paste0(i, "_stack"))

# Create layers plus blank vector
plot.layers <- c( "bio", "carbon", "rec.visits", "flood.ES.flow", "poll.ES.flow",
                  "wilderness", "lands.val", "blank", "agri.val")
plot.layer.shorts <- c("B", "C", "R", "F", "P", "W", "L", " ", "A*")

# Create the table which we are going to populate
num.rows <-  length(plot.layers)*length(runs)*length(quants)

# Create data frame (add in blank column to separate A* from the rest)
acc.data <- data.frame("runs" = rep(runs, each = length(plot.layers) * length(quants)),
                       "layers" = rep(plot.layers, times = length(runs), each = length(quants)),
                       "quants" = rep(quants, num.rows/length(quants)),
                       "total.value" = rep(NA,num.rows),
                       "added.value" = rep(NA,num.rows),
                       "pos.or.neg" = rep(NA,num.rows))

# State levels
acc.data$runs <- factor(acc.data$runs,
                        levels = runs)
acc.data$layers <- factor(acc.data$layers,
                          levels = plot.layers)
acc.data$quants <- as.factor(acc.data$quants)

# Create a loop to populate the table (apart from added.value)
iter = 1

for (k in 1 : length(runs)) {
  for (l in 1 : length(plot.layers)) {
    for (m in 1 : length(quants)) {
      
      ### Read in quantile mask
      
      # Get raster from run k
      r <- raster(SP.stack[[k]]) 
      
      # Extract only data above quantile
      r[r<quants[m]] <- NA
      r[!is.na(r)] <- 1
      
      ### Positive or negative resource
      
      if(plot.layers[l] != "agri.val") {
        acc.data$pos.or.neg[iter] <- "Positive"
      }
      
      if(plot.layers[l] == "agri.val") {
        acc.data$pos.or.neg[iter] <- "Negative"
        }
      
      # Only carry out if not the blank layer
      if(plot.layers[l] != "blank") {
        
      ### Mask layer raster and sum

      # If layer is positive resource, mask to included amount
      if(plot.layers[l] != "agri.val") {
        acc.value <- mask(get(plot.layers[l]), mask = r, inverse = FALSE) %>% cellStats(., sum)
      }
      # If layer is negative resource, mask to excluded amount
      if(plot.layers[l] == "agri.val") {
        acc.value <- mask(get(plot.layers[l]), mask = r, inverse = TRUE) %>% cellStats(., sum)
      }
      
      ### Find cumulative quantile proportion covered
      
      # Divide by total layer 'resource' to get proportion
      acc.prop <- acc.value / cellStats(get(plot.layers[l]), sum)
      
      # Add to total prop data table
      acc.data$total.value[iter] <- acc.prop
      
      }
      iter = iter + 1
    }
  }
}

### Add "added.value" column 
# ( N.B. just for plot as "agri.val" quants are switched on plot)

iter = 1

for (k in 1 : length(runs)) {
  for (l in 1 : length(plot.layers)) {
    for (m in 1 : length(quants)) {
      
      if (plot.layers[l] != "agri.val") {
        
        if (m == 1) { acc.added.prop <- acc.data$total.value[iter]}
        
        if (m > 1) { # if not the first quantile...
          
          acc.added.prop <- acc.data$total.value[iter] - #take away the sum of all other higher quantiles
            acc.data$total.value[acc.data$runs == runs[k] &
                                   acc.data$layers == plot.layers[l] &
                                    acc.data$quants == quants[m-1]]
        }}
      
      if (plot.layers[l] == "agri.val") {
        
        if (m == 1) {
          
          acc.added.prop <- acc.data$total.value[acc.data$runs == runs[k] &
                                                   acc.data$layers == plot.layers[l] &
                                                   acc.data$quants == quants[length(quants)]]
          }
      
      if (m > 1) { 
          
          acc.added.prop <- acc.data$total.value[acc.data$runs == runs[k] &
                                                    acc.data$layers == plot.layers[l] &
                                                    acc.data$quants == quants[length(quants) - m + 1]] - 
            acc.data$total.value[acc.data$runs == runs[k] &
                                   acc.data$layers == plot.layers[l] &
                                   acc.data$quants == quants[length(quants) - m + 2]]
        }}
        

      # Assign to acc.table
      acc.data$added.value[iter] <- acc.added.prop
  
      iter = iter + 1
      
    }}}

# Save for later efficiency analysis for each removal rule
save(acc.data, file = paste0(dataDir, "Spatial_datasets/", i, "_viewpoint_resource_coverage.Rdata"))

}

# PLOT INLAID RESOURCE COVERAGE BARPLOT ------------------------------------------

### MAPS

for (i in rem.rule){

  SP.stack <- get(paste0(i, "_stack"))
  load(file = paste0(dataDir, "Spatial_datasets/", i, "_viewpoint_resource_coverage.Rdata"))
 
 # Loop to convert continuous ranks to discrete quantiles
  for (k in 1:length(SP.stack)) {
    
    # Reassign ranks to quantiles
    ranks <- SP.stack[[k]]@data
    
    # Change highest quantile
    ranks[ranks >= quants[1] & !is.na(ranks)] <- quants[1]
    
    # Change other quantiles
    for (l in 2:length(quants)) {
      ranks[ranks >= quants[l] &
              ranks < quants[l - 1] &
              !is.na(ranks)] <- quants[l]
    }
    
    # Change other values to NA
    ranks[ranks < min(quants)] <- NA
    
    # Convert back
    SP.stack[[k]]@data <- ranks
  }
 
  for (m in 1:length(runs)) {
    
    ### Create plot data
    plotData <- data.frame(SP.stack[[m]]@data)
    plotData[, c("x", "y")] <- coordinates(SP.stack[[m]])
    ### Reshape the data for ggplot
    plotData <- reshape2::melt(plotData, id.vars = c('x', 'y'))
    levels(plotData$variable) <- runs[m]
    
    # Crop to non-NA values only
    plotData <- subset(plotData, x >= 75000 &
                         x <= 645000 &
                         y >= 25000 &
                         y <= 1205000)
    
    ### Map
    Zon.map <- ggplot(aes(x = x, y = y), data = plotData) +
      geom_tile(aes(fill = as.factor(value))) + facet_wrap( ~ variable) +
      scale_fill_manual(na.value = "white",
                        values = rev.pal) +
      coord_equal() +
      my.theme  + theme(
        legend.position = "none",
        plot.background = element_blank(),
        plot.title = element_blank(),
        strip.text.x = element_blank() )
    
    Zon.map <- Zon.map + geom_polygon(
      data = GB.poly,
      aes(long, lat, group = group),
      color = "black",
      fill = NA)

### INLAID BARPLOT

ES.barplot <- ggplot(acc.data[acc.data$runs == runs[m],],
                     aes( y = added.value,
                          x = factor(layers, labels = plot.layer.shorts),
                          fill = interaction(quants, pos.or.neg))) + 
  geom_bar(position="stack", stat="identity", col = "black") +
  scale_fill_manual( values = c(pal, rev.pal) ) +
  facet_wrap(~runs) +
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

### COMBINE MAP AND BARPLOT

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

assign(paste0("inlaid.plot.", runs[m]), inlaid.plot, envir = .GlobalEnv)

} # warning messages okay, just removing NA values from blank rows

# Dummy plot for legend.grob
 dummy.legend <- ggplot(acc.data, aes( y=added.value, x=layers, fill = interaction(quants, pos.or.neg))) + 
   geom_bar(position="stack", stat="identity", col = "black") +
   scale_fill_manual( values = c(pal, rev.pal),
                      breaks = c("0.95.Positive","0.9.Positive", "0.83.Positive", "0.7.Positive" ),
                      labels = c("5%", "10%", "17%", "30%"),
                      name = "Total land coverage ") +
   facet_wrap(~runs) + 
   theme(legend.position="bottom",
         legend.key.width = unit(1,"cm"),
         legend.title = element_text(size = 40),
         legend.text = element_text(size = 40),)
 
# Plot all 4 viewpoints with 'ggarrange'
all.inlaid.plots <- ggarrange(inlaid.plot.TRAD, inlaid.plot.NEW, inlaid.plot.ECON, inlaid.plot.SOC,
           ncol=2, nrow=2, 
           labels = runs,
           label.x = 0.68, label.y = 0.4, 
           font.label = c(size = 40), hjust = c(-0.48,-0.65, -0.46, -0.74),
           common.legend = TRUE, legend= "bottom", legend.grob = get_legend(dummy.legend))+ 
  theme(plot.margin = margin(2,5,2,1, "lines"),
        plot.background = element_rect(fill = "white", colour = "white"))

all.inlaid.plots

# Save
ggsave(all.inlaid.plots, filename = paste0("../Writing/Figures/", i,"_viewpoints_inlaid.png"),
       width = 50, height = 65, dpi = 400, units = "cm")
ggsave(all.inlaid.plots, filename = paste0("../Writing/Figures/", i,"_viewpoints_inlaid.pdf"),
       width = 50, height = 65, dpi = 600, units = "cm")
}
