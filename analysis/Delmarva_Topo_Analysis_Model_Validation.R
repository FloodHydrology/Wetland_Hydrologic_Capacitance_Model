###################################################################################
#Name: Baltimore Corner Topographic Analysis
#Coder: C. Nathan Jones
#Date: 1/22/2018
#Purpose: Topographic analysis for Baltimore Corner Wetlands
##################################################################################

#Clear Memory
rm(list=ls(all=TRUE))

#Define Master Working Directory
dir<-"//nfs/WHC-data/Validation_Modeling/WHC_BaltimoreCorner"

#Load Required Packages
library(raster)   #spatial analysis
library(rgdal)    #spatial analysis
library(rgeos)    #spatial analysis
library(maptools) #spatial analysis
library(plyr)     #data processing
library(dplyr)    #data processing

#Save to backup folder
load(paste0(dir, "/Backup/input_dat.RData"))


#Call Function DEM Processing F(x) into Memory
source("/R/DEM_Processing.R")

#Run Terrain Analysis Subroutine
DEM_Processing.fun(dem, pp, wd)

#Backup
save.image(paste0(dir,"/Backup/DEM_Processing.RData"))
