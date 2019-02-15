###################################################################################
#Name: Wetlandscape Hydrology - Figure 5 Wetland & Inter-Catchment Comparison

#Date: 14 JAN 2019
#Purpose: Create figure 5 to plot observed/simulated 

##################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
rm(list=ls(all=TRUE))

library(reshape2)
library(ggplot2)

# Load data; using PPR as test case
load("/nfs/WHC-data/Regional_Analysis/Delmarva/backup/results.RData")



# Seasonal Water Balance Plots -------------------------------------------------------
fluxes <- as.data.frame(sapply(results[,25:28], as.numeric)) #<- sapply is here
colnames(fluxes) <- c('SW_in', 'SW_out', 'GW_out', 'GW_in')

plot_data_fluxes <- melt(fluxes)

ggplot(plot_data_fluxes, aes( x = value, fill = variable)) + 
  geom_histogram(color = 'black', alpha = 0.7, breaks=seq(-501,500,25)) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold')) +
  labs(x = 'Day of Year', y = 'Normalized Water Level', title = 'Figure 5a - Delmarva Wetlands') 

summary(fluxes)


plot_data_fluxes$value = abs(plot_data_fluxes$value)
ggplot(plot_data_fluxes, aes( x = variable, y = value, fill = variable)) + 
  geom_boxplot(outlier.shape = NA) +
  ylim(-50, 500) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold')) +
  labs(x = 'Day of Year', y = 'Normalized Water Level', title = 'Figure 5a - Delmarva Wetlands') 
