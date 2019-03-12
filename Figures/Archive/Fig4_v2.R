###################################################################################
#Name: Wetlandscape Hydrology - Figure 4 Seasonality in depth and fluxes

#Date: 14 JAN 2019
#Purpose: Create figure 4 to plot observed/simulated 

##################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
rm(list=ls(all=TRUE))

library(reshape2)
library(ggplot2)
library(dplyr)

# var.name = 'y_w'
var.name = 'runoff_in'

## For Delmarva --------------------------------------------------------------------
# 1. Load data ---------------------------------------------------------------------

rawdata <- read.table('/nfs/WHC-data/Figure Generation/delmarva_output.csv', sep=",", header=TRUE)
rawdata$value <- as.numeric(as.character(rawdata$value))

levels <- subset(rawdata, scale == 'weltand' & var == var.name)   # need to fix typo!

#level_stat[,2:4] <- t(apply(levels, 1, quantile, probs=c(0.25,0.5, 0.75), na.rm=TRUE))

stat <- levels %>%
  group_by(day) %>%
  summarise(x25th = quantile(value, probs = 0.25, na.rm = T), 
            median = quantile(value, probs = 0.5, na.rm = T),
            x75th = quantile(value, probs = 0.75, na.rm = T))

stat2 <- melt(stat, by_id = 'day')
stat2$day <- as.numeric(as.character(stat$day))
stat$day <- as.numeric(as.character(stat$day))

# 2. Plot Data ---------------------------------------------------------------------
ggplot() + 
  geom_line(data = stat2, aes(x = day, y = value, group = variable),
              size = 1) +
  geom_ribbon(data = stat, aes(x = day, ymin = x25th, ymax = x75th),
              fill = 'blue', alpha = 0.1) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold')) +
  labs(x = 'Day of Year', y = 'Normalized Water Level', title = 'Figure 4a - Delmarva') 



## For PPR----- --------------------------------------------------------------------
# 1. Load data ---------------------------------------------------------------------
rawdata <- read.table('/nfs/WHC-data/Figure Generation/ppr_2.csv', sep=",", header=TRUE)
# rawdata <- read.table('/nfs/WHC-data/Figure Generation/ppr_output.csv', sep=",", header=TRUE)
rawdata$value <- as.numeric(as.character(rawdata$value))

levels <- subset(rawdata, scale == 'weltand' & var == var.name)   # need to fix typo!

#level_stat[,2:4] <- t(apply(levels, 1, quantile, probs=c(0.25,0.5, 0.75), na.rm=TRUE))

stat <- levels %>%
  group_by(day) %>%
  summarise(x25th = quantile(value, probs = 0.25, na.rm = T), 
            median = quantile(value, probs = 0.5, na.rm = T),
            x75th = quantile(value, probs = 0.75, na.rm = T))

stat2 <- melt(stat, by_id = 'day')
stat2$day <- as.numeric(as.character(stat$day))
stat$day <- as.numeric(as.character(stat$day))

# 2. Plot Data ---------------------------------------------------------------------
ggplot() + 
  geom_line(data = stat2, aes(x = day, y = value, group = variable, 
                              linetype=ifelse(variable=='median', 'longdash', 'solid')),
            size = 1) +
  geom_ribbon(data = stat, aes(x = day, ymin = x25th, ymax = x75th),
              fill = 'blue', alpha = 0.1) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.position = "none") +
  labs(x = 'Day of Year', y = 'Runoff_in', title = 'Figure 4b - PPR') 








# =================================================================================
# Catchment Scale Values
# =================================================================================

rm(list=ls(all=TRUE))

## For Delmarva --------------------------------------------------------------------
# 1. Load data ---------------------------------------------------------------------

rawdata <- read.table('/nfs/WHC-data/Figure Generation/delmarva_output.csv', sep=",", header=TRUE)
rawdata$value <- as.numeric(as.character(rawdata$value))

levels <- subset(rawdata, scale == 'catchment' & var == 'y_wt')   

#level_stat[,2:4] <- t(apply(levels, 1, quantile, probs=c(0.25,0.5, 0.75), na.rm=TRUE))

stat <- levels %>%
  group_by(day) %>%
  summarise(x25th = quantile(value, probs = 0.05, na.rm = T), 
            median = quantile(value, probs = 0.5, na.rm = T),
            x75th = quantile(value, probs = 0.95, na.rm = T))

stat2 <- melt(stat, by_id = 'day')
stat2$day <- as.numeric(as.character(stat$day))
stat$day <- as.numeric(as.character(stat$day))

# 2. Plot Data ---------------------------------------------------------------------
ggplot() + 
  geom_line(data = stat2, aes(x = day, y = value, group = variable),
            size = 1) +
  geom_ribbon(data = stat, aes(x = day, ymin = x25th, ymax = x75th),
              fill = 'blue', alpha = 0.1) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold')) +
  labs(x = 'Day of Year', y = 'Normalized Wetland Stage', title = 'Figure 4b - Delmarva Catchment Scale') 


rm(list=ls(all=TRUE))

## For PPR ---- --------------------------------------------------------------------
# 1. Load data ---------------------------------------------------------------------
rawdata <- read.table('~/Wetland_Hydrologic_Capacitance_Model/data/ppr_2.csv', sep=",", header=TRUE)
# rawdata <- read.table('/nfs/WHC-data/Figure Generation/ppr_output.csv', sep=",", header=TRUE)
rawdata$value <- as.numeric(as.character(rawdata$value))

levels <- subset(rawdata, scale == 'catchment' & var == 'y_wt')   

#level_stat[,2:4] <- t(apply(levels, 1, quantile, probs=c(0.25,0.5, 0.75), na.rm=TRUE))

stat <- levels %>%
  group_by(day) %>%
  summarise(x25th = quantile(value, probs = 0.25, na.rm = T), 
            median = quantile(value, probs = 0.5, na.rm = T),
            x75th = quantile(value, probs = 0.75, na.rm = T))

stat2 <- melt(stat, by_id = 'day')
stat2$day <- as.numeric(as.character(stat$day))
stat$day <- as.numeric(as.character(stat$day))

# 2. Plot Data ---------------------------------------------------------------------
ggplot() + 
  geom_line(data = stat2, aes(x = day, y = value, group = variable, 
                              linetype=ifelse(variable=='median', 'longdash', 'solid')),
            size = 1) +
  geom_ribbon(data = stat, aes(x = day, ymin = x25th, ymax = x75th),
              fill = 'blue', alpha = 0.1) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.position = "none") +
  labs(x = 'Day of Year', y = 'Normalized Wetland Stage', title = 'Figure 4b - PPR Catchment Scale') 










# plot climate data
load('/nfs/WHC-data/Regional_Analysis/PPR/inputs/climate.Rdata')
a <- rep(1:365, times = 1000)
b <- data.frame(a, pet.VAR, precip.VAR)







