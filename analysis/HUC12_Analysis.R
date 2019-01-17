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
#2.0 Global Options
n.years<-10
n.nodes<-16
n.cpus<-8

#2.1 Delmarva~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1.a run script to prep input data
source("analysis/HUC12_Model_Input_Delmarva.R")

#temporary backup for testing
save.image("backup_temp.RData")
#-------------------------
remove(list=ls())
load("backup_temp.RData")

#2.1.b run the model
#define functions from file 
source("R/regional_analysis.R")
source("R/WHC_2.R")

#Create wrapper function 
fun<-function(ID){
   regional_analysis(WetID=ID,n.years, pet.VAR,precip.VAR,wetlands.shp,HUC12.shp, 
                     catchments.shp, flowlines.shp,fac.grd, soils.shp, dem.grd, 
                     nfw_centroid.shp, rootdepth.grd)
}

#run using SLURM
sopts <- list(partition = "sesync", time = "12:00:00" )
params<-data.frame(ID=wetlands.shp$WetID[1:16])
t0<-Sys.time()
delmarva<- slurm_apply(fun, 
                       params,
                       add_objects = c(
                         #Functions
                         "wetland.hydrology", "regional_analysis",
                         #Spatial data
                         "fac.grd","catchments.shp","flowlines.shp",
                         "soils.shp","wetlands.shp","dem.grd", "rootdepth.grd", 'n.years',
                         #Climate data
                         "precip.VAR", "pet.VAR"),
                       nodes = n.nodes, cpus_per_node=n.cpus,
                       pkgs=c('sp','raster','rgdal','rgeos','dplyr'),
                       slurm_options = sopts)

# 3.4 Retrieve results
print_job_status(delmarva)
results <- get_slurm_out(delmarva, outtype = "table")
cleanup_files(delmarva)
tf<-Sys.time()
tf-t0

#2.2 Run model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~









#####################################################################################
# Step 3:PPR Analysis----------------------------------------------------------------
#####################################################################################





