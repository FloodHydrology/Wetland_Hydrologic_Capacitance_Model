####################################################################################
# Name: HUC12 Analysis
# Coder: C. Nathan Jones & Fred Cheng
# Date: 10 Jan 2019
# Purpose: Complete Regional Analysis  
####################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
# 1a. Clear Memory  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls(all=TRUE))

# 1b. Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sp)         # for spatial analysis
library(raster)     # for spatial analysis
library(rgdal)      # for spatial analysis
library(rgeos)      # for spatial analysis
library(geosphere)  # for spatial analysis
library(tidyverse)  # for data wrangling
library(MASS)       #distribrutional analysis
library(markovchain)# for MCMC modeling
library(Evapotranspiration)
library(lubridate)  # dealing with dates
library(parallel)   # parallel computing
library(rslurm)     # parallel computing

####################################################################################
# Step 2: Regional simulations------------------------------------------------------
####################################################################################
#2.1 Delmarva~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1.a run script to prep input data
source("analysis/HUC12_Model_Input_Delmarva.R")

#2.1.b run the model
#grab function
source("R/regional_analysis.R")
source("R/WHC_2.R")
# regional_analysis(WetID=1,n.years =10, pet.VAR,precip.VAR,wetlands.shp,HUC12.shp, 
#                   catchments.shp, flowlines.shp,fac.grd, soils.shp, dem.grd, 
#                   nfw_centroid.shp, rootdepth.grd)



sopts <- list(partition = "sesynctest", time = "1:00:00" )
params<-data.frame(WetID=wetlands.shp$WetID[1:16],
                   n.years=10)
delmarva<- slurm_apply(regional_analysis, params,
                       add_objects = c(
                         #Functions
                         "wetland.hydrology",
                         #Spatial data
                         "fac.grd","catchments.shp","flowlines.shp",
                         "soils.shp","wetlands.shp","dem.grd", "rootdepth.grd",
                         #Climate data
                         "precip.VAR", "pet.VAR"),
                       nodes = 2, cpus_per_node=8,
                       pkgs=c('sp','raster','rgdal','rgeos','dplyr'),
                       slurm_options = sopts)



#2.2 Run model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~









#####################################################################################
# Step 3:PPR Analysis--------------------------------------------------------------
#####################################################################################





