###################################################################################
#Name: Wetlandscape Hydrology - Figure 5 Regional Comparisons

#Date: 14 JAN 2019
#Purpose: Create figure 5

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
rawdata$region <- 'DELMARVA'

rawdata2 <- read.table('/nfs/WHC-data/Figure Generation/ppr_output.csv', sep=",", header=TRUE)
rawdata2$value <- as.numeric(as.character(rawdata2$value))
rawdata2$region <- 'PPR'

data <- rbind(rawdata, rawdata2)

levels <- subset(data, scale == 'weltand' & day == 0)   # need to fix typo!
levels <- subset(data, var == 'qf_in' | var == 'qf_out' | 
                          var == 'SW_in' | var == 'SW_out' |
                          var == 'GW_out')


ggplot(levels, aes( x = var, y = value, fill = region)) + 
  geom_boxplot(color = 'black', alpha = 0.7, outlier.shape = NA) +
  scale_y_continuous(limits = quantile(levels$value, c(0.1, 0.8)))+
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold')) +
  labs(x = '', y = 'Annual Flux', title = 'Figure 5a - Wetland Scale') 



# Catchment scale -----------------------------------------------------------------
# 1. Load data ---------------------------------------------------------------------

levels <- subset(data, scale == 'catchment' & day == 0)   # need to fix typo!
levels <- subset(data, var == 'sw_out' | var == 'gw_out' | 
                   var == 'sw_out_day' | var == 'gw_out_day')


ggplot(levels, aes( x = var, y = value, fill = region)) + 
  geom_boxplot(color = 'black', alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~ var, ncol = 2, scales='free') +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold')) +
  labs(x = '', y = 'Flux', title = 'Figure 5b - Catchment Scale') 

