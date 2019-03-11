###################################################################################
#Name: Wetlandscape Hydrology - Figure 3 Region Validation

#Date: 14 JAN 2019
#Purpose: Create figure 3 to plot observed/simulated 

##################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
rm(list=ls(all=TRUE))
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd('/nfs/WHC-data/Figure Generation/')
# 1. Load data ---------------------------------------------------------------------
rawdata1 <- read.table('/nfs/WHC-data/Figure Generation/delmarva_output.csv', sep=",", header=TRUE)
rawdata1$value <- as.numeric(as.character(rawdata1$value))
rawdata1$region <- 'DEL'


rawdata2 <- read.table('/nfs/WHC-data/Figure Generation/ppr_output.csv', sep=",", header=TRUE)
rawdata2$value <- as.numeric(as.character(rawdata2$value))
rawdata2$region <- 'PPR'

data <- rbind(rawdata1, rawdata2)


df <- data %>%
  filter(day == 0) %>%
  filter(!var %in% c("GW_in")) %>%
  filter(scale %in% c("weltand")) %>%
  select(-one_of('day', 'X'))

df.wide <- spread(df, var, value)
regional.pca <- prcomp(df.wide[,c(6, 8, 12, 14)], center = TRUE,scale. = TRUE)

#setwd("~/Wetland_Hydrologic_Capacitance_Model/data")
#save(regional.pca, df.wide,file =  'PCA2.rdata')


#install.packages('ggfortify')
#library(ggfortify)

autoplot(regional.pca, data = df.wide, colour = 'region',
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)

