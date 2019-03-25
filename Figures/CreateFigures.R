###################################################################################
#Name: Wetlandscape Hydrology - Main Script to Create Figures

#Date: 14 JAN 2019
#Purpose: Call scripts to create individual figures
##################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
rm(list=ls(all=TRUE))

library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

# Create Figure 4 - Seasonality Plots ------------------------------------------------
# Note: Output plots are saved in nfs/WHC-data/Figure Generation/Output
setwd("~/Wetland_Hydrologic_Capacitance_Model/Figures")
source('Fig4_FUN.R')

Fig4_FUN('ppr')
Fig4_FUN('delmarva')
Fig4_FUN('florida')

# Create Figure 5 - Seasonality Plots ------------------------------------------------
# Note: Output plots are saved in nfs/WHC-data/Figure Generation/Output

# 1. Load data (to be used for figure 6 as well) -------------------------------------
rawdata <- read.table('/nfs/WHC-data/Figure Generation/delmarva.csv', sep=",", header=TRUE)
rawdata$value <- as.numeric(as.character(rawdata$value))
rawdata$region <- 'DELMARVA'

rawdata2 <- read.table('/nfs/WHC-data/Figure Generation/ppr.csv', sep=",", header=TRUE)
rawdata2$value <- as.numeric(as.character(rawdata2$value))
rawdata2$region <- 'PPR'

rawdata3 <- read.table('/nfs/WHC-data/Figure Generation/florida.csv', sep=",", header=TRUE)
rawdata3$value <- as.numeric(as.character(rawdata3$value))
rawdata3$region <- 'FLORIDA'

data <- do.call("rbind", list(rawdata, rawdata2, rawdata3))
# remove(rawdata, rawdata2, rawdata3) #don't actually need

setwd("~/Wetland_Hydrologic_Capacitance_Model/Figures")
source('Fig5_FUN.R')

Fig5_FUN(data)


# BROKEN ----NEED TO FIX!!!!! ):     ): --------------------------------------
# Create Figure 6 - PCA -------------------------------------------------------
# Note: Output plots are saved in nfs/WHC-data/Figure Generation/Output
setwd("~/Wetland_Hydrologic_Capacitance_Model/Figures")
source('Fig6_FUN.R')

Fig6_FUN(data)

