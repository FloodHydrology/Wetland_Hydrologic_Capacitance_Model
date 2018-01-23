###################################################################################
#Name: Baltimore Corner Topographic Analysis
#Coder: C. Nathan Jones
#Date: 1/22/2018
#Purpose: Topographic analysis for Baltimore Corner Wetlands
##################################################################################
#Clear Memory
rm(list=ls(all=TRUE))

#Define working and data directories
wd<-"C:\\Users/njones/Desktop/workspace/Wetland_Hydrologic_Capacitance_Model"
dir<-"X:/Validation_Modeling/WHC_BaltimoreCorner"

#Load Required Packages
library(RPyGeo)
library(raster)   #spatial analysis
library(rgdal)    #spatial analysis
library(rgeos)    #spatial analysis
library(maptools) #spatial analysis
library(plyr)     #data processing
library(dplyr)    #data processing

#Call Function DEM Processing F(x) into Memory
source(paste0(wd,"/R/DEM_Processing.R"))

#Save to backup folder
load(paste0(dir, "/Backup/input_data.RData"))

#Re-define working and data directories
wd<-"C:\\Users/njones/Desktop/workspace/Wetland_Hydrologic_Capacitance_Model"
dir<-"X:/Validation_Modeling/WHC_BaltimoreCorner"

#Run Terrain Analysis Subroutine
DEM_Processing.fun(dem, 
                   pp.shp,
                   wetland.shp,
                   python.path="C:\\Python27\\ArcGIS10.4\\", 
                   scratchspace="C:\\ScratchWorkspace")

#Plot for shits and giggles
plot(basin.shp)
plot(dem.grd, add=T)
plot(basin.shp, add=T)
plot(pp.shp, pch=19, col="red", cex=3, add=T)
plot(wetland.shp, pch=19, col="green", cex=3, add=T)

#Backup
save.image(paste0(dir,"/Backup/DEM_Processing.RData"))
