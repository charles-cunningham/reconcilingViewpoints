# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Evaluate viewpoint and integration approach efficiency
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
library(gtable)
library(RColorBrewer)
library(ggpubr)

### GB
source("Functions/SNIPPET_create_GB_shapes.R") 

# Set Data folder
dataDir <- "../Data/"

# Overwrite GB.rast with 'trimmed' GB raster
GB.rast <- raster(paste0(dataDir, "Spatial_datasets/GB.raster.tif"))

#####
### DEFINE SPECIFICS
#####

#List the quantiles we're going to extract accumulation data for
top.percent <- c(5,10,17,30)
quants <- 1 - (top.percent/100)
quants <- sort(quants, decreasing = TRUE)

# Layers and layer names
cutoff <- "mod.10.ModandFiltRaw"
taxa <- c("Merge_all")
runs <- c("TRAD", "NEW", "ECON", "SOC")
aggs <- c("INCLUSIVE", "PLURALIST", "MEAN", "RANK") #c("INCLUSIVE", "PLURALIST")
rem.rule <- c("CAZ", "ABF")

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

# Create layers plus blank vector
plot.layers <- c( "bio", "carbon", "rec.visits", "flood.ES.flow", "poll.ES.flow",
                  "wilderness", "lands.val", "blank", "agri.val")
plot.layer.shorts <- c("B", "C", "R", "F", "P", "W", "L", " ", "A*")

# Define theme
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
                  panel.background = element_rect(fill = NA, color = "white"),
                  legend.justification=c(1,0.8),
                  legend.position=c(1,0.8),
                  plot.margin= unit(rep(0, 4), "lines"))

# Colour schemes
pal <- brewer.pal(n = length(quants), name = "RdBu")

rev.pal <- brewer.pal(n = length(quants), name = "RdBu") %>% 
  rev(.)

reds <- brewer.pal(n = length(quants), name = "Reds")
reds.rev <- brewer.pal(n = length(quants), name = "Reds") %>% 
  rev(.)

# CALCULATE EFFICIENCY TABLES -------------------------------------------------

# Loop through removal rules
for (i in rem.rule) {

### Load in resource coverage tables from previous script

# Load viewpoint coverage data from PLOT_Zonation.R script
load(paste0(dataDir, "Spatial_datasets/", i, "_viewpoint_resource_coverage.Rdata"))
# Load aggregate coverage data from ANALYSIS_trade_offs.R script
load(paste0(dataDir, "Spatial_datasets/", i, "_aggregate_resource_coverage.Rdata"))

# Change "agg" column name of acc.agg to be able to join to acc.data
colnames(acc.agg)[1] <- "runs.and.aggs"
colnames(acc.data)[1] <- "runs.and.aggs"

# Merge viewpoints and aggregates
acc.all <- rbind(acc.data, acc.agg)

# Remove "added.value" column to avoid confusion, this was just for plot
acc.all$added.value <- NULL

### Maximum proportion resource coverage possible

# Create column
acc.all$max <- NA

# Populate new column
for (j in 1:NROW(acc.all)) {
  
  # Only carry out if not the blank layer
  if (acc.all$layers[j] != "blank") {
    
    # Extract relevant quantile
    quant <- acc.all$quants[j] %>%
      as.character(.) %>%
      as.numeric(.)
    
    # Extract ES raster
    r <- as.character(acc.all$layers[j]) %>% get(.)
    
    if (acc.all$pos.or.neg[j] == "Positive") {
      
      # if positive ES...
      
      # Get top values possible at quantile coverage
      top.values <- values(r) %>%
        .[!is.na(.)] %>%
        .[. >= quantile(., prob = quant)] %>%
        sum(.)
      
      # Assign top proportion to df
      acc.all$max[j] <- top.values / cellStats(r, sum)
    }
    
    else if (acc.all$pos.or.neg[j] == "Negative") {
      
      # if negative ES, i.e. agri.val...
      
      # Get bottom values possible at quantile coverage
      bottom.values <- values(r) %>%
        .[!is.na(.)] %>%
        .[. < quantile(., prob = (1 - quant))] %>%
        sum(.)
      
      # Assign 1 minus minimum possible proportion (i.e. max proportion excluded) to df
      acc.all$max[j] <- 1 - (bottom.values / cellStats(r, sum))
      
    }
  }
}

# Create new column of viewpoint/aggregate "efficiency"
acc.all$eff <- acc.all$total.value/ acc.all$max

# Reverse quant factor levels for plot
acc.all$quants <- factor(acc.agg$quants, levels = quants)

### Create data frame of summary values 

# Create empty data frame
agg.stats <- data.frame("runs.and.aggs" =  levels(acc.all$runs.and.aggs) %>% rep(., each = length(quants)),
                        "quants" = rep(quants, times = length(levels(acc.all$runs.and.aggs))),
                        "Mean" = NA,
                        "Median" = NA,
                        "Minimum" = NA,
                        "grouping" = rep(c("A","A","A","A","B","B","B","B") %>% rep(., each = length(quants))))

# Populate data frame with mean, median, and min
for(j in 1:NROW(agg.stats)) {
  
  agg.stats$Mean[j] <- acc.all[acc.all$runs.and.aggs == agg.stats$runs.and.aggs[j] &
                                 acc.all$quants == agg.stats$quants[j],]$eff %>% 
    na.omit(.) %>% mean(.)
  
  agg.stats$Median[j] <- acc.all[acc.all$runs.and.aggs == agg.stats$runs.and.aggs[j] &
                                   acc.all$quants == agg.stats$quants[j],]$eff %>% 
    na.omit(.) %>% median(.)
  
  agg.stats$Minimum[j] <- acc.all[acc.all$runs.and.aggs == agg.stats$runs.and.aggs[j] &
                                    acc.all$quants == agg.stats$quants[j],]$eff %>% 
    na.omit(.) %>% min(.)
}

# Sort levels and melt
agg.stats$runs.and.aggs <- factor(agg.stats$runs.and.aggs,
                                  levels = c(runs, aggs))
agg.stats$quants <- factor(agg.stats$quants, levels = quants)

agg.stats <- reshape2::melt(agg.stats, id = c("runs.and.aggs", "quants", "grouping"))

# Assign
assign(paste0(i, ".acc.all"), acc.all, envir = globalenv())
assign(paste0(i, ".agg.stats"), agg.stats, envir = globalenv())

}

# Save data frames
save(CAZ.acc.all, file = paste0(dataDir, "Spatial_datasets/CAZ_efficiency.Rdata"))
save(ABF.acc.all, file = paste0(dataDir, "Spatial_datasets/ABF_efficiency.Rdata"))
save(CAZ.agg.stats, file = paste0(dataDir, "Spatial_datasets/CAZ_efficiency_summary.Rdata"))
save(ABF.agg.stats, file = paste0(dataDir, "Spatial_datasets/ABF_efficiency_summary.Rdata"))

# PLOTS -------------------------------------------------------------------------

# Load data frames (also load these for extracting values)
load(file = paste0(dataDir, "Spatial_datasets/CAZ_efficiency.Rdata"))
load(file = paste0(dataDir, "Spatial_datasets/ABF_efficiency.Rdata"))
load(file = paste0(dataDir, "Spatial_datasets/CAZ_efficiency_summary.Rdata"))
load(file = paste0(dataDir, "Spatial_datasets/ABF_efficiency_summary.Rdata"))

for (i in rem.rule) {

# Load in relevant 'acc.all' file
acc.all <- get(paste0(i, ".acc.all"))
agg.stats <- get(paste0(i, ".agg.stats"))

### DUMMY LEGEND

# Dummy plot for legend.grob
dummy.legend <- ggplot(data = agg.stats[agg.stats$variable == "Minimum",], 
                       aes(y = value,
                           x = factor(runs.and.aggs, labels = c(runs, aggs)),
                           fill = interaction(quants,runs.and.aggs))) +
  geom_bar(position="stack", stat="identity", col = "black") +
  scale_fill_manual( values = rep(pal, times = length(c(runs, aggs))),
                     breaks = c("0.95.TRAD", "0.9.TRAD", "0.83.TRAD", "0.7.TRAD"),
                     labels = c("5%", "10%", "17%", "30%"),
                     name = "Total land coverage ") +
  theme(legend.position="bottom",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40))

### Plot efficiency summary stats

mean.plot <- ggplot(data = agg.stats[agg.stats$variable == "Mean",], 
                    aes(y = value,
                        x = factor(runs.and.aggs, labels = c(runs, "INC", "PLUR", "MEAN", "RANK")),
                        fill = interaction(quants, runs.and.aggs))) +
  geom_bar(position="dodge", stat="identity", col = "black") +
  facet_wrap(~grouping,  scales = "free_x") + # facet_grid(~grouping,  scales = "free_x", space = "free_x") for equal bar widths
  scale_fill_manual( values = c(pal, pal, pal, pal, pal, pal, pal, pal)) +
  scale_y_continuous( breaks = seq(0,1,0.2),
                      minor_breaks = seq(0,0.8,0.05),
                      limits = c(0,0.8)) + 
  scale_x_discrete(expand = expansion(add = 0.75)) +
  ylab("Mean coverage efficiency") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40 , vjust = 2),
        axis.text.x = element_text(size = 24, colour = "black", vjust = 5),
        axis.text.y = element_text(size = 30, colour = "black"),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.grid.minor.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.ontop = TRUE,
        plot.margin = margin(0,1,0,0, "lines"),
        panel.spacing = unit(0, "lines"))

min.plot <- ggplot(data = agg.stats[agg.stats$variable == "Minimum",], 
                    aes(y = value,
                        x = factor(runs.and.aggs, labels = c(runs, "INC", "PLUR", "MEAN", "RANK")),
                        fill = interaction(quants, runs.and.aggs))) +
  geom_bar(position="dodge", stat="identity", col = "black") +
  facet_wrap(~grouping,  scales = "free_x") + # facet_grid(~grouping,  scales = "free_x", space = "free_x") for equal bar widths
  scale_fill_manual( values = c(pal, pal, pal, pal, pal, pal, pal, pal)) +
  scale_y_continuous( breaks = seq(0,1,0.2),
                      minor_breaks = seq(0,0.8,0.05),
                      limits = c(0,0.8)) + 
  scale_x_discrete(expand = expansion(add = 0.75)) +
  ylab("Minimum coverage efficiency") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40 , vjust = 2),
        axis.text.x = element_text(size = 24, colour = "black", vjust = 5),
        axis.text.y = element_text(size = 30, colour = "black"),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.grid.minor.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.ontop = TRUE,
        plot.margin = margin(0,0,0,1, "lines"),
        panel.spacing = unit(0, "lines"))

ggarrange(mean.plot, min.plot, ncol = 2,
          common.legend = TRUE, legend= "top", legend.grob = get_legend(dummy.legend)) + 
  theme(plot.margin = margin(0,1,0,1, "lines"))

ggsave(last_plot(), filename = paste0("../Writing/Figures/", i, "_efficiency_summary_supplemental.png"),
       width = 60, height = 30, dpi = 300, units = "cm")#
ggsave(last_plot(), filename = paste0("../Writing/Figures/", i, "_efficiency_summary_supplemental.pdf"),
       width = 60, height = 30, dpi = 600, units = "cm")

### Plot efficiency summary stats for main figure (CAZ only)
if (i == "CAZ") {
  
mean.plot <- ggplot(data = agg.stats[agg.stats$variable == "Mean" &
                                       !(agg.stats$runs.and.aggs %in% c("MEAN", "RANK")),], 
                    aes(y = value,
                        x = factor(runs.and.aggs, labels = c(runs, aggs[1:2])),
                        fill = interaction(quants, runs.and.aggs))) +
  geom_bar(position="dodge", stat="identity", col = "black") +
  facet_wrap(~grouping,  scales = "free_x") + # facet_grid(~grouping,  scales = "free_x", space = "free_x") for equal bar widths
  scale_fill_manual( values = c(pal, pal, pal, pal, pal, pal)) +
  scale_y_continuous( breaks = seq(0,1,0.2),
                      minor_breaks = seq(0,0.8,0.05),
                      limits = c(0,0.8)) + 
  ylab("Mean coverage efficiency") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40 , vjust = 2),
        axis.text.x = element_text(size = 28, colour = "black", vjust = 5),
        axis.text.y = element_text(size = 30, colour = "black"),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.grid.minor.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.ontop = TRUE,
        plot.margin = margin(0,1,0,0, "lines"),
        panel.spacing = unit(0, "lines"))

min.plot <- ggplot(data = agg.stats[agg.stats$variable == "Minimum" &
                                      !(agg.stats$runs.and.aggs %in% c("MEAN", "RANK")),], 
                   aes(y = value,
                       x = factor(runs.and.aggs, labels = c(runs, aggs[1:2])),
                       fill = interaction(quants, runs.and.aggs))) +
  geom_bar(position="dodge", stat="identity", col = "black") +
  facet_wrap(~grouping,  scales = "free_x") + # facet_grid(~grouping,  scales = "free_x", space = "free_x") for equal bar widths
  scale_fill_manual( values = c(pal, pal, pal, pal, pal, pal)) +
  scale_y_continuous( breaks = seq(0,1,0.2),
                      minor_breaks = seq(0,0.8,0.05),
                      limits = c(0,0.8)) + 
  ylab("Minimum coverage efficiency") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40 , vjust = 2),
        axis.text.x = element_text(size = 26, colour = "black", vjust = 5),
        axis.text.y = element_text(size = 30, colour = "black"),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.grid.minor.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.ontop = TRUE,
        plot.margin = margin(0,0,0,1, "lines"),
        panel.spacing = unit(0, "lines"))

ggarrange(mean.plot, min.plot,ncol = 2,
          common.legend = TRUE, legend= "top", legend.grob = get_legend(dummy.legend)) + 
  theme(plot.margin = margin(0,1,0,1, "lines"))

ggsave(last_plot(), filename = paste0("../Writing/Figures/", i, "_efficiency_summary_main.png"),
       width = 60, height = 30, dpi = 300, units = "cm")
ggsave(last_plot(), filename = paste0("../Writing/Figures/", i, "_efficiency_summary_main.pdf"),
       width = 60, height = 30, dpi = 600, units = "cm")
}

### Plot all resource efficiency
eff.plot <- ggplot(data = acc.all, aes( y = eff,
                                        x = factor(layers, labels = plot.layer.shorts),
                                        fill = interaction(quants,runs.and.aggs))) + ## change here
  geom_bar(position="dodge", stat="identity", col = "black") +
  facet_wrap(~runs.and.aggs, ncol = 2, scales='free') +
  scale_fill_manual( values = rep(pal, times = length(c(runs, aggs)))) +
  scale_y_continuous( breaks = seq(0,1,0.2),
                      minor_breaks = seq(0,1,0.1),
                      limits = c(0,1)) +
  theme_minimal() +
  xlab(element_blank()) +
  ylab("Coverage Efficiency") +
  theme(axis.title.x = element_text(size=30,margin=margin(t=-10,b=10)),
        axis.title.y = element_text(size=40, ),
        axis.text.x = element_text(size=30, margin=margin(t=-5,b=10), colour = "black"),
        axis.text.y = element_text(size=24, colour = "black"),
        strip.text = element_text(size=30, colour = "black", margin=margin(t=40,b=10)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.grid.minor.y = element_line(colour = "black", size=0.1, linetype="dashed"),
        panel.ontop = TRUE)

eff.plot <- ggarrange(eff.plot, ncol = 1, common.legend = TRUE, legend= "top", 
                         legend.grob = get_legend(dummy.legend)) + 
  theme(plot.margin = margin(0,1,0,1, "lines"))

ggsave(eff.plot, filename = paste0("../Writing/Figures/", i, "_efficiency_supplemental.png"),
       width = 35, height = 50, dpi = 300, units = "cm")
ggsave(eff.plot, filename = paste0("../Writing/Figures/", i, "_efficiency_supplemental.pdf"),
       width = 35, height = 50, dpi = 600, units = "cm")

}

# ANALYSIS ----------------------------------------------------------------

# Load dataframes (also load these for extracting values)
load(file = paste0(dataDir, "Spatial_datasets/CAZ_efficiency.Rdata"))
load(file = paste0(dataDir, "Spatial_datasets/ABF_efficiency.Rdata"))
load(file = paste0(dataDir, "Spatial_datasets/CAZ_efficiency_summary.Rdata"))
load(file = paste0(dataDir, "Spatial_datasets/ABF_efficiency_summary.Rdata"))

CAZ.agg.stats$percent <- (CAZ.agg.stats$value * 100) %>% round(., 1) 

# Mean feature coverage
CAZ.agg.stats[CAZ.agg.stats$variable == "Mean" &
                CAZ.agg.stats$quants == 0.83,]

# Minimum feature coverage
CAZ.agg.stats[CAZ.agg.stats$variable == "Minimum" &
                CAZ.agg.stats$quants == 0.83,]


# Range within ABF mean viewpoints compared to CAZ

ABF.agg.stats[ABF.agg.stats$variable == "Mean" &
                ABF.agg.stats$quants == 0.83 &
                ABF.agg.stats$grouping == "B",] %>%
  .["value"] %>%
  range(.) %>%
  round(., 3)


CAZ.agg.stats[CAZ.agg.stats$variable == "Mean" &
                CAZ.agg.stats$quants == 0.83 &
                CAZ.agg.stats$grouping == "B",] %>%
  .["value"] %>%
  range(.) %>%
  round(., 3)

# Range within ABF min viewpoints compared to CAZ

ABF.agg.stats[ABF.agg.stats$variable == "Minimum" &
                ABF.agg.stats$quants == 0.7 &
                ABF.agg.stats$grouping == "B",] %>%
  .["value"] %>%
  range(.) %>%
  round(., 3)


CAZ.agg.stats[CAZ.agg.stats$variable == "Minimum" &
                CAZ.agg.stats$quants == 0.7 &
                CAZ.agg.stats$grouping == "B",] %>%
  .["value"] %>%
  range(.) %>%
  round(., 3)

# Loss in coverage of biodiversity between TRAD and compromise approaches
CAZ.acc.all[CAZ.acc.all$layers == "bio",] # total value column

