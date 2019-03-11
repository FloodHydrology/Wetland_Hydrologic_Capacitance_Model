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
#2.1 Define global simulation options------------------------------------
n.years<-10
n.nodes<-16
n.cpus<-8

#define functions from file 
source("R/regional_analysis.R")
source("R/WHC_2.R")

#2.2 Delmarva----------------------------------------------------------------------
#a run script to prep input data
source("analysis/HUC12_Model_Input_Delmarva.R")

#b Create wrapper function 
dmv_fun<-function(ID){
   regional_analysis(WetID=ID,n.years, pet.VAR,precip.VAR,wetlands.shp,HUC12.shp, 
                     catchments.shp, flowlines.shp,fac.grd, soils.shp, dem.grd, 
                     nfw_centroid.shp, rootdepth.grd)
}

#c run using SLURM (this will take ~2.5 hrs)
sopts <- list(partition = "sesync", time = "12:00:00")
params<-data.frame(ID=wetlands.shp$WetID)
delmarva<- slurm_apply(dmv_fun, 
                       params,
                       add_objects = c(
                         #Functions
                         "wetland.hydrology", "regional_analysis",
                         #Spatial data
                         "fac.grd","catchments.shp","flowlines.shp","HUC12.shp",
                         "soils.shp","wetlands.shp","dem.grd", "rootdepth.grd", 'n.years',
                         #Climate data
                         "precip.VAR", "pet.VAR"),
                       nodes = n.nodes, cpus_per_node=n.cpus,
                       pkgs=c('sp','raster','rgdal','rgeos','tidyverse'),
                       slurm_options = sopts)

#e Retrieve results
print_job_status(delmarva)
results <- get_slurm_out(delmarva, outtype = "table")
cleanup_files(delmarva)
tf<-Sys.time()
tf-t0

#f Write output to "data" folder
write.csv(results,"data/delmarva.csv")

#2.3 PPR----------------------------------------------------------------------
#a remove previous data
remove(list=ls()[ls()!= 'n.years' &
                 ls()!= 'n.nodes' &
                 ls()!= 'n.cpus'  &
                 ls()!= 'regional_analysis' &
                 ls()!= 'wetland.hydrology'])

#b run script to prep input data
source("analysis/HUC12_Model_Input_PPR.R")

#c Create wrapper function 
ppr_fun<-function(ID){
  regional_analysis(WetID=ID,n.years, pet.VAR,precip.VAR,wetlands.shp,HUC12.shp, 
                    catchments.shp, flowlines.shp,fac.grd, soils.shp, dem.grd, 
                    nfw_centroid.shp, rootdepth.grd)}

#d run using SLURM (this will take ~2.5 hrs)
sopts  <- list(partition = "sesync", time = "12:00:00")
params <- data.frame(ID=wetlands.shp$WetID)
ppr    <- slurm_apply(ppr_fun, 
                      params,
                      add_objects = c(
                         #Functions
                         "wetland.hydrology", "regional_analysis",
                         #Spatial data
                         "fac.grd","catchments.shp","flowlines.shp","HUC12.shp",
                         "soils.shp","wetlands.shp","dem.grd", "rootdepth.grd", 'n.years',
                         #Climate data
                         "precip.VAR", "pet.VAR"),
                      nodes = n.nodes, cpus_per_node=n.cpus,
                      pkgs=c('sp','raster','rgdal','rgeos','tidyverse'),
                      slurm_options = sopts)

#e Retrieve results
print_job_status(ppr)
results <- get_slurm_out(ppr, outtype = "table")
cleanup_files(ppr)

#f Write output to "data" folder
write.csv(results,"data/ppr.csv")

#2.3 PPR----------------------------------------------------------------------
#a remove previous data
remove(list=ls()[ls()!= 'n.years' &
                   ls()!= 'n.nodes' &
                   ls()!= 'n.cpus'  &
                   ls()!= 'regional_analysis' &
                   ls()!= 'wetland.hydrology'])

#b run script to prep input data
source("analysis/HUC12_Model_Input_Florida.R")

#c Create wrapper function 
florida_fun<-function(ID){
  regional_analysis(WetID=ID,n.years, pet.VAR,precip.VAR,wetlands.shp,HUC12.shp, 
                    catchments.shp, flowlines.shp,fac.grd, soils.shp, dem.grd, 
                    nfw_centroid.shp, rootdepth.grd)}

#d run using SLURM (this will take ~2.5 hrs)
sopts   <- list(partition = "sesync", time = "12:00:00")
params  <- data.frame(ID=wetlands.shp$WetID)
florida <- slurm_apply(florida_fun, 
                       params,
                       add_objects = c(
                        #Functions
                        "wetland.hydrology", "regional_analysis",
                        #Spatial data
                        "fac.grd","catchments.shp","flowlines.shp","HUC12.shp",
                        "soils.shp","wetlands.shp","dem.grd", "rootdepth.grd", 'n.years',
                        #Climate data
                        "precip.VAR", "pet.VAR"),
                       nodes = n.nodes, cpus_per_node=n.cpus,
                       pkgs=c('sp','raster','rgdal','rgeos','tidyverse'),
                       slurm_options = sopts)

#e Retrieve results
print_job_status(florida)
results <- get_slurm_out(florida, outtype = "table")
cleanup_files(florida)

#f Write output to "data" folder
write.csv(results,"data/florida.csv")