###################################################################################
#Name: Wetlandscape Hydrology - Figure 4b Seasonality in GW_local

#Date: 14 JAN 2019
#Purpose: Create figure 4b

##################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
rm(list=ls(all=TRUE))

library(reshape2)
library(ggplot2)
library(dplyr)

variable_name <- 'runoff_in'
label_x <- 'Runoff In'

# 1. Load data ---------------------------------------------------------------------

rawdata <- read.table('/nfs/WHC-data/Figure Generation/delmarva_output.csv', sep=",", header=TRUE)
rawdata$value <- as.numeric(as.character(rawdata$value))

levels <- subset(rawdata, scale == 'weltand' & var == variable_name )   # need to fix typo!

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
  labs(x = 'Day of Year', y = label_x, title = 'Figure 4c - Delmarva') 



# 1. Load data ---------------------------------------------------------------------

rawdata <- read.table('/nfs/WHC-data/Figure Generation/ppr_output.csv', sep=",", header=TRUE)
rawdata$value <- as.numeric(as.character(rawdata$value))

levels <- subset(rawdata, scale == 'weltand' & var == variable_name)   # need to fix typo!

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
  labs(x = 'Day of Year', y = label_x, title = 'Figure 4d - PPR') 
