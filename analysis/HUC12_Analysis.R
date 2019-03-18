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
library(MASS)       #distribrutional analysis
library(markovchain)# for MCMC modeling
library(Evapotranspiration)
library(lubridate)  # dealing with dates
library(parallel)   # parallel computing
library(rslurm)     # parallel computing
library(tidyverse)  # for data wrangling

# 1c. Define data directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
backup_dir<-"/nfs/WHC-data/Regional_Analysis/backup/"
results_dir<-"/nfs/WHC-data/Figure Generation/"

####################################################################################
# Step 2:Input data for regional analysis-------------------------------------------
####################################################################################
#2.1 Load regional analysis function
source("R/regional_analysis.R")

#2.2 Delmarva
source("analysis/HUC12_Model_Input_Delmarva.R")
save.image(paste0(backup_dir,"Delmarva_Input.RData"))
remove(list=ls()[ls()!= 'backup_dir' & ls()!= 'results_dir' & ls()!= 'regional_analysis'])

#2.3 ppr
source("analysis/HUC12_Model_Input_PPR.R")
save.image(paste0(backup_dir,"PPR_Input.RData"))
remove(list=ls()[ls()!= 'backup_dir' & ls()!= 'results_dir' & ls()!= 'regional_analysis'])

#2.4 Florida
source("analysis/HUC12_Model_Input_Florida.R")
save.image(paste0(backup_dir,"Florida_Input.RData"))
remove(list=ls()[ls()!= 'backup_dir' & ls()!= 'results_dir' & ls()!= 'regional_analysis'])

####################################################################################
# Step 3: Regional simulations------------------------------------------------------
####################################################################################
#3.1 Define global simulation options-----------------------------------------------
cluster_name<-"sesync"
time_limit<-"12:00:00"
n.years<-10
n.nodes<-3
n.cpus<-8

#define functions from file 
source("R/WHC_2.R")

#3.2 Delmarva----------------------------------------------------------------------
#a load input data
load(paste0(backup_dir,"Delmarva_Input.RData"))

#b Create wrapper function 
dmv_fun<-function(ID){
   regional_analysis(WetID=ID,n.years, pet.VAR,precip.VAR,wetlands.shp,HUC12.shp, 
                     catchments.shp, flowlines.shp,fac.grd, soils.shp, dem.grd, 
                     nfw_centroid.shp, rootdepth.grd)
}

#c run using SLURM (this will take ~2.5 hrs)
sopts <- list(partition = cluster_name, time = time_limit)
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

#3.3 PPR----------------------------------------------------------------------
#a remove previous data
remove(list=ls()[ls()!=   'ppr' &
                 ls()!= 'delmarva'&
                 ls()!= 'cluster_name' &
                 ls()!= 'time_limit' &   
                 ls()!= 'n.years' &
                 ls()!= 'n.nodes' &
                 ls()!= 'n.cpus'  &
                 ls()!= 'regional_analysis' &
                 ls()!= 'wetland.hydrology' & 
                 ls()!= 'backup_dir' &
                 ls()!= 'results_dir'])

#b load ppr input data
load(paste0(backup_dir,"PPR_Input.RData"))

#c Create wrapper function 
ppr_fun<-function(ID){
  regional_analysis(WetID=ID,n.years, pet.VAR,precip.VAR,wetlands.shp,HUC12.shp, 
                    catchments.shp, flowlines.shp,fac.grd, soils.shp, dem.grd, 
                    nfw_centroid.shp, rootdepth.grd)}

#d run using SLURM (this will take ~2.5 hrs)
sopts  <- list(partition = cluster_name, time = time_limit)
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

#3.4 Florida----------------------------------------------------------------------
#a remove previous data
remove(list=ls()[  ls()!= 'ppr' &
                   ls()!= 'delmarva'&
                   ls()!= 'cluster_name' &
                   ls()!= 'time_limit' &   
                   ls()!= 'n.years' &
                   ls()!= 'n.nodes' &
                   ls()!= 'n.cpus'  &
                   ls()!= 'regional_analysis' &
                   ls()!= 'wetland.hydrology' & 
                   ls()!= 'backup_dir' &
                   ls()!= 'results_dir'])


#b load input data
load(paste0(backup_dir,"Florida_Input.RData"))

#c Create wrapper function 
florida_fun<-function(ID){
  regional_analysis(WetID=ID,n.years, pet.VAR,precip.VAR,wetlands.shp,HUC12.shp, 
                    catchments.shp, flowlines.shp,fac.grd, soils.shp, dem.grd, 
                    nfw_centroid.shp, rootdepth.grd)}

#d run using SLURM (this will take ~2.5 hrs)
sopts   <- list(partition = cluster_name, time = time_limit)
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

#3.5 Gather Data-------------------------------------------------------------------
#a check job status
print_job_status(delmarva)
print_job_status(ppr)
print_job_status(florida)

#a gather results
results_delmarva <- get_slurm_out(delmarva, outtype = "table")
results_ppr      <- get_slurm_out(ppr,      outtype = "table")
results_florida  <- get_slurm_out(florida,  outtype = "table")

#b write to output folder
write.csv(results_delmarva,  paste0(results_dir,"delmarva.csv"))
write.csv(results_ppr,       paste0(results_dir,"ppr.csv"))
write.csv(results_florida,   paste0(results_dir,"florida.csv"))

#clean up file
cleanup_files(delmarva)
cleanup_files(ppr)
cleanup_files(florida)






