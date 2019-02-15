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

# Load data; using PPR as test case
load("/nfs/WHC-data/Regional_Analysis/Delmarva/backup/results.RData")

level <- as.data.frame(sapply(results[,41:405], as.numeric)) #<- sapply is here
level <- data.frame(t(level))

# are results already normalized? if not...
#scaled_level <- data.frame(scale(level))
scaled_level <- level
scaled_level$DOY <- 1:365

level_stat <- data.frame(rowMeans(scaled_level[1:ncol(level)],na.rm=TRUE))
level_stat[,2:4] <- t(apply(scaled_level, 1, quantile, probs=c(0.25,0.5, 0.75), na.rm=TRUE))

colnames(level_stat) <- c('mean','25th', 'median','75th')
level_stat$DOY <- 1:365
mean_level$DOY <- 1:365

plot_data <- melt(scaled_level, id.vars = "DOY")
plot_data2 <- melt(level_stat, id.vars = "DOY")

# Plot Data ------------------------------------------------------------------------
ggplot() + 
  geom_line(data = subset(plot_data2, variable == 'median'), 
            aes(x = DOY, y = value),
                size = 1)  + 
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black', size = 0.5) +
  geom_ribbon(aes(x = level_stat$DOY, ymin = level_stat[,2],
                  ymax = level_stat[,4]) ,
              fill="grey", alpha=0.6) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold')) +
  labs(x = 'Day of Year', y = 'Normalized Water Level', title = 'Figure 4a - Delmarva') 


# Seasonal Water Balance Plots -------------------------------------------------------
# need daily SW/GW flux (median or mean across the wetlands) 














