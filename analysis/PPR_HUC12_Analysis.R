###################################################################################
# Name: Sensitivity Analysis Dem
# Coder: C. Nathan Jones
# Date: 8 March 2018
# Purpose: Demo sensitivity analysis for individual wetland manuscript
##################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
# 1a. Clear Memory  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls(all=TRUE))

# 1b. Load Packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sp)          # for spatial analysis
library(raster)      # for spatial analysis
library(rgdal)       # for spatial analysis
library(rgeos)       # for spatial analysis
library(dplyr)       # for data processing
library(rslurm)      # parallel computing
library(geosphere)   # for spatial analysis

# 1c. Define Master Working Directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wd<-"//nfs/WHC-data/Regional_Analysis/PPR/"    # Define working directory for later reference
setwd(wd)                                      # Set working directory

# 1d. Obtain Model Inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download Model Inputs (and put everything in a uniform projection)
load("Model_inputs/climate.Rdata")                                      # Climate data from PPR_HUC12_Climate_Data                                             
wetlands.shp<-readOGR("Model_inputs/.","wetlands")                      # Import wetland shapefile
  wetlands.shp@data$WetID<-seq(1,length(wetlands.shp))                  # add (?) WetID column to shp file using sequence from 1:total num of wetlands
  p<-wetlands.shp@proj4string                                           # Define projection of interest
HUC12.shp<-readOGR("Model_inputs/.","HUC12")                            # Import HUC 12 shapefile
  HUC12.shp<-spTransform(HUC12.shp, wetlands.shp@proj4string)           # transform HUC 12 shapefile's projection into same as wetlands
catchments.shp<-readOGR("Model_inputs/.","catchments")                  # import NHD catchments
  catchments.shp<-spTransform(catchments.shp, p)                        # transform NHD catchment's projection into same as wetlands
flowlines.shp<-readOGR("Model_inputs/.","flow_net")                     # import NHD flowlines
  flowlines.shp<-spTransform(flowlines.shp, p)                          # transform NHD flowline's projection into same as wetlands
fac.grd<-raster("Model_inputs/fac2")                                    # import flow accumulation raster
soils.shp<-readOGR("Model_inputs/.","soils")                            # import soils shapefile
  soils.shp<-spTransform(soils.shp, p)                                  # transform soil file's projection into same as wetlands
  soils<-read.csv("Model_inputs/WHC_Soils_Input_PPR.csv")                   # import existing soil parameters csv
  soils.shp@data<-merge(soils.shp@data,soils, by.x='MUKEY', by.y="MUID")# append soils csv data into soils shapefile by matching MUKEY and MUID
  remove(soils)                                                         # delete soils df
dem.grd<-raster("Model_inputs/dem_cm")                                  # import DEM file
  mask<-spTransform(catchments.shp, dem.grd@crs)                        # transform file's projection
  dem.grd<-crop(dem.grd, mask)                                          # crop out portion of dem to keep relevant areas only
  remove(mask)                                                          # remove dummy variable mask
  dem.grd<-projectRaster(dem.grd, crs = p)                              # project raster
nfw_centroid.shp<-readOGR("Model_inputs/.","NFW_centroids")             # NFW from Lane and D'Amico 2016
  nfw_centroid.shp<-spTransform(nfw_centroid.shp, p)                    # project polygon
rootdepth.grd<-raster("Model_inputs/rootdepth")                         # Rooting depth estimate from Fan et al., 2017 [PNAS]
  rootdepth.grd<-projectRaster(rootdepth.grd, crs = p)                  # project raster
  rootdepth.grd<-raster::resample(rootdepth.grd, dem.grd)               # resample raster so it is same size

# 1e. Delineate "Isolated" Wetlands  (Use NFW from Lane and D'Amico 2016) ~~~~~~~~~~~~~
wetlands.shp$dist2stream<-apply(gDistance(flowlines.shp,wetlands.shp, byid=TRUE),1,min) # Calculate euclidean distance to stream network
wetlands.shp$area_m2<-gArea(wetlands.shp,byid=T)                                        # Caclculate wetland area
wetlands.shp<-wetlands.shp[nfw_centroid.shp,]

# 1f. Alternate distances ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wetlandCentroid <- centroid(wetlands.shp)                                   # get centroid of each wetland
distMat <- pointDistance(wetlandCentroid, lonlat = FALSE)                   # get the distances to each other point
diag(distMat) <- NA                                                         # replace the diagonal's zeros with NAs
wetlands.shp$dist2NearWet <- apply(distMat, 1, min, na.rm=TRUE)/2           # put the half distances into variable

# 1g. Plot for funzies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(HUC12.shp, lwd=2, border="grey10")                                     # plot overall HUC12 outline 
plot(catchments.shp, lwd=0.5, border="grey50", add=T)                       # plot interior catchments
plot(flowlines.shp, lwd=0.5, col="dark blue", add=T)                        # plot NHD flowlines
plot(wetlands.shp, col="dark blue", add=T)                                  # plot wetlands

# 1h. Save Image ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save.image("backup/Inputs.Rdata")                                           # save workspace for use later

####################################################################################
# Step 2: Create Function to run individual wetlands--------------------------------
####################################################################################
# 2.1 Set Up workspace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
remove(list=ls())                                           # clear environment
wd<-"//nfs/WHC-data/Regional_Analysis/PPR/"                 # Define working directory for later reference
setwd(wd)                                                   # Set working directory
load("backup/Inputs.Rdata")                                 # load inputs from previous processing
source("~/Wetland_Hydrologic_Capacitance_Model/R/WHC_2.R")  # compile WHC function 
source("~/Wetland_Hydrologic_Capacitance_Model/R/get_yc.R") # compile yc 

# 2.2 Calculate flowpath lengths for individual wetlands~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.2.1 Create function to cacluate distance from wetland edge to "groundwater inflection" point
divide_dist.fun<-function(cat){
  
  #a. Identify catchment of interest
  cat.shp<-catchments.shp[cat,]
  
  #b. Break catchment into 100 equally spaced points
  cat.grd<-raster(cat.shp)
  extent(cat.grd)<-extent(cat.shp)
  cellsize<-2*pi*((gArea(cat.shp)/pi)^.5)/100 #paremeter/1000
  res(cat.grd)<-c(cellsize,cellsize)
  cat.grd<-rasterize(cat.shp, cat.grd, 1)
  cat.grd<-boundaries(cat.grd, type="inner")
  cat.grd[cat.grd==0]<-NA
  cat.pnt<-rasterToPoints(cat.grd)
  cat.pnt<-SpatialPoints(cat.pnt)
  
  #c. Clip relevant wetlands
  wetlands.shp<-wetlands.shp[cat.shp,]
  
  #d. Create inner function to define median distance for each wetland
  fun<-function(WetID){
    
    #i. wetland of interest centroid
    ind_cent<-data.frame(centroid(wetlands.shp[wetlands.shp$WetID==WetID,]))
    colnames(ind_cent)<-c("x","y")
    ind_pnt<-ind_cent
    coordinates(ind_cent) <- ~x + y 
    ind_cent<-SpatialPoints(ind_cent)
    
    #ii. All centroids
    wet_cent<-data.frame(centroid(wetlands.shp))
    colnames(wet_cent)<-c("x","y")
    coordinates(wet_cent) <- ~x + y 
    wet_cent<-SpatialPoints(wet_cent)
    
    #iii. Calculate distances between wetlands
    wetlands.shp@data$dist2wet<-c(gDistance(ind_cent, wet_cent, byid=T))
    wetlands.shp<-wetlands.shp[rank(wetlands.shp$dist2wet)>1,]
    
    #iii. Calculate the dLe for each wetland
    wetlands.shp$dLe<-0
    for(i in 1:length(wetlands.shp$dist2wet)){
      #Create line between centroid
      flowpath <- data.frame(centroid(wetlands.shp[i,]))
      colnames(flowpath)<-c("x","y")
      flowpath<-rbind(flowpath, ind_pnt)
      coordinates(flowpath) <- ~x + y 
      flowpath.shp<-SpatialPoints(flowpath)
      flowpath.shp<-SpatialLines(list(Lines(Line(flowpath.shp), ID="a")))
      
      #Clip based on area wetlands
      flowpath.shp<-gDifference(flowpath.shp, wetlands.shp)
      
      #Calculate length
      if(length(flowpath.shp)==0){
        dLe<-0
      }else{
        dLe<-gLength(flowpath.shp)/2
      }
      
      #Define "delta' distance
      wetlands.shp$dLe[i]<-dLe
    }
    
    #iv. Deterimine length to polygon boundary
    watershed_boundary<-min(gDistance(ind_cent,cat.pnt, byid=T))/2
    
    #v. Export dLe
    c(WetID,median(c(wetlands.shp$dLe, watershed_boundary), na.rm=T))
  }
  
  #e. Run function (if there are enough wetlands)
  if(length(wetlands.shp)>1){
    
    #run inner function for each wetland
    output<-lapply(wetlands.shp$WetID, fun)
    output<-do.call(rbind, output)
    
    #Export output
    output
  }else{
    c(0,0)
  }
}

# 2.2.2 Run on the slurm server
sopts <- list(partition = "sesync", time= "1:00:00")
params<-data.frame(cat=seq(1,length(catchments.shp)))
flowpath_job<- slurm_apply(divide_dist.fun, params,
                           add_objects = c("catchments.shp","wetlands.shp"),
                           nodes = 4, cpus_per_node=8,
                           pkgs=c('sp','raster','rgdal','rgeos','geosphere'),
                           slurm_options = sopts)
results <- get_slurm_out(flowpath_job, outtype = "table")
cleanup_files(flowpath_job)

# 2.2.3 collect results and clean output
#Print results and clean
results<-data.frame(results)
colnames(results)<-c("WetID","dLe")
results<-results[results$WetID!=0,]
results<-results[!duplicated(results[,1]),]

#merge with wetlands.shp
wetlands.shp<-wetlands.shp[(wetlands.shp$WetID %in% results$WetID),]
wetlands.shp@data<-merge(wetlands.shp@data, results, by="WetID")

# 2.3 Create function to process data and run WHC for a wetland of interest~~~~~~~~~
fun<-function(WetID){ #
  
  #for testing
  #WetID<-wetlands.shp@data$WetID[1]
  
  # 2.3.1 Gather Required Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.3.1a Identify wetland of interest
  main_wetland.shp<-wetlands.shp[wetlands.shp$WetID==WetID,]
  main_wetland<-WetID
  
  # 2.3.1b Select spatial layers based on WetID  
  catchment_temp.shp<-catchments.shp[main_wetland.shp,]
  if(length(catchment_temp.shp)>1){catchment_temp.shp<-catchment_temp.shp[1,]}
  wetlands_temp.shp<-wetlands.shp[catchment_temp.shp,]
  soils_temp.shp<-soils.shp[catchment_temp.shp,]
  soils_temp.shp<-crop(soils_temp.shp,catchment_temp.shp)
  dem_temp.grd<-crop(dem.grd, catchment_temp.shp)
  dem_temp.grd<-mask(dem_temp.grd, catchment_temp.shp)
  dem_temp.grd<-dem_temp.grd/100 #Convert to m
  fac_temp.grd<-crop(fac.grd, catchment_temp.shp)
  fac_temp.grd<-mask(fac_temp.grd, catchment_temp.shp)
  root_temp.grd<-crop(rootdepth.grd, catchment_temp.shp)
  root_temp.grd<-mask(rootdepth.grd, catchment_temp.shp)
  
  # 2.3.1c For now, add catchment aggregate soils data to soils with missing data
  soils_temp.shp@data$y_cl[is.na(soils_temp.shp@data$y_cl)]<-mean(soils_temp.shp@data$y_cl, na.rm=T)
  soils_temp.shp@data$y_rd[is.na(soils_temp.shp@data$y_rd)]<-mean(soils_temp.shp@data$y_rd, na.rm=T)
  soils_temp.shp@data$s_fc[is.na(soils_temp.shp@data$s_fc)]<-mean(soils_temp.shp@data$s_fc, na.rm=T)
  soils_temp.shp@data$s_w[is.na(soils_temp.shp@data$s_w)]<-mean(soils_temp.shp@data$s_w, na.rm=T)
  soils_temp.shp@data$n[is.na(soils_temp.shp@data$n)]<-mean(soils_temp.shp@data$n, na.rm=T)
  soils_temp.shp@data$clay[is.na(soils_temp.shp@data$clay)]<-mean(soils_temp.shp@data$clay, na.rm=T)
  soils_temp.shp@data$ksat[is.na(soils_temp.shp@data$ksat)]<-mean(soils_temp.shp@data$ksat, na.rm=T)
  
  # 2.3.2 Estimate area and volume to stage relationships~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.3.2a Create vectors to store area and volume to stage relationships
  area<-matrix(0, ncol=length(wetlands_temp.shp), nrow=100)
  volume<-matrix(0, ncol=length(wetlands_temp.shp), nrow=100)
  
  # 2.3.2b Use loop to calculate based on previously published relationships
  for(i in 1:length(wetlands_temp.shp)){
    #Define TempID
    TempID<-wetlands_temp.shp$WetID[i]
    
    #Define Amax
    Amax<-wetlands_temp.shp$area_m2[wetlands_temp.shp$WetID==TempID]
    
    #Define rmax
    rmax<-(Amax/pi)^.5
    
    #Define Vmax (Using Wu and Lane 2016)
    Vmax<-(0.25*((Amax/10000)^1.4742))*10000 #m3
    
    #Estimate hmax using Hayashi and Kamp 2000
    hmax<-Vmax*1.5/Amax
    
    #create functions to calculate area and volume
    area.fun<-function(h){pi*((rmax*((h/hmax)^.25))^2)}
    volume.fun<-function(h){0.25*((area.fun(h)/10000)^1.4742)*10000}  
    
    #Apply functions 
    n.col<-which(wetlands_temp.shp$WetID==TempID)
    n.rows<-length(area.fun(seq(0,hmax,0.05)))
    area[1:n.rows,n.col]<-area.fun(seq(0,hmax,0.05))
    area[(n.rows+1):100,n.col]<-area.fun(hmax)
    volume[1:n.rows,n.col]<-volume.fun(seq(0,hmax,0.05))
    volume[(n.rows+1):100,n.col]<-volume.fun(hmax)
  }
  
  # 2.3.3 Populate GIW.INFO Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.3.3a Create input variables
  giw.INFO<-c("giw.ID","WetID","area_watershed","area_wetland","invert","vol_ratio", "dL", "dLe", "dz",  # geometric characteristics
              "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", "Sy",                        # soil characteristics
              "y_w_0", "s_t_0"                                                                           # initial conditions
  )
  
  # 2.3.3b Create giw.INFO matrix
  giw.INFO<-matrix(0, nrow=length(wetlands_temp.shp), ncol=length(giw.INFO), 
                   dimnames = list(seq(1,length(wetlands_temp.shp)), c(giw.INFO)))
  
  # 2.3.3c Define wetland variables (units == mm)
  for(i in 1:length(wetlands_temp.shp)){
    
    # i. Isolate soils data (if for osme reason we don't have ovlerlap, just use basin agregate)
    if(gIntersects(soils_temp.shp,wetlands_temp.shp[i,])==TRUE){
      temp_soils<-crop(soils_temp.shp,wetlands_temp.shp[i,])
    }else{
      temp_soils<-soils_temp.shp
    }
    
    # ii. Agregrate parameters based on space
    temp_soils$area<-gArea(temp_soils, byid=T)
    temp_soils<-temp_soils@data
    temp_soils<-colSums(temp_soils[,c("y_cl","y_rd","s_fc","s_w","n","clay","ksat")]*temp_soils[,"area"])/sum(temp_soils$area, na.rm=T)
    
    # iii. Isolate fac data and calculate ratio of upland drainage area to total drainage area
    temp_fac<-crop(fac_temp.grd, wetlands_temp.shp[i,])
    temp_fac<-mask(temp_fac, wetlands_temp.shp[i,])
    temp_fac<-ifelse(cellStats(temp_fac*0+1, sum)>0,
                     cellStats(temp_fac,max)/cellStats(fac_temp.grd,max),
                     0)
    
    # iv. Calculate relative wetland datum
    dz<-cellStats(mask(dem_temp.grd, wetlands_temp.shp[i,]), median, na.rm=T)- #median elevation of wetland elevation
      cellStats(mask(dem_temp.grd, wetlands_temp.shp[wetlands_temp.shp$WetID!=WetID,]), median, na.rm=T) #median elevation of all wetlands
    dz<-ifelse(is.na(dz)==TRUE, 0, dz)
    dz<-ifelse(dz==Inf, 0, dz)
    dz<-dz*1000
    
    # v. Calcluate wetland rooting depth
    temp_rd<-crop(root_temp.grd, wetlands_temp.shp[i,])
    temp_rd<-mask(temp_rd, wetlands_temp.shp[i,])
    temp_rd<-ifelse(cellStats(temp_rd*0+1, sum)>0,
                     cellStats(temp_rd,max)/cellStats(root_temp.grd,max),
                     0)

    # vi. Populate giw.INFO table
    giw.INFO[i,"giw.ID"]<-         i  #ID used in the model, note this is differnt than the WetID
    giw.INFO[i,"WetID"]<-          wetlands_temp.shp$WetID[i]
    giw.INFO[i,"area_watershed"]<- gArea(catchment_temp.shp)*(10^6)
    giw.INFO[i,"area_wetland"]<-   gArea(wetlands_temp.shp[i,], byid=T)*(10^6)
    giw.INFO[i,"invert"]<-         -1000*0.05*(length(area[,i][area[,i]!=max(area[,i], na.rm=T)]))   
    giw.INFO[i,"n"]<-              temp_soils["n"]
    giw.INFO[i,"s_fc"]<-           temp_soils["s_fc"]/100  
    giw.INFO[i,"psi"]<-            -16662*(temp_soils["n"]^7.8831) 
    giw.INFO[i,"y_cl"]<-           -1*temp_soils["y_cl"]   
    giw.INFO[i,"y_c"]<-            - temp_soils["y_rd"]/2                       #critical depth (mm) 
    giw.INFO[i,"s_wilt"]<-         temp_soils["s_w"]/100                       # soil moisture at permanent wilting point
    giw.INFO[i,"k_sat"]<-          -temp_soils["ksat"]*24                      # saturated condcuctivity (mm/day)
    giw.INFO[i,"RD"]<-             ifelse((temp_rd*1000)<temp_soils["y_cl"],   # Rooting Depth (mm)
                                      -1*temp_rd*1000,
                                      -1*temp_soils["y_cl"])
    giw.INFO[i,"b"]<-              12.524*(temp_soils["clay"]/100)+3.6907 
    giw.INFO[i,"Sy"]<-             giw.INFO[i,"n"]*(1-giw.INFO[i,"s_fc"])
    giw.INFO[i,"y_w_0"]<-          0
    giw.INFO[i,"s_t_0"]<-          giw.INFO[i,"s_fc"]
    giw.INFO[i,"vol_ratio"]<-      temp_fac                                   # ratio of upstream wetland volume
    giw.INFO[i, "dL"] <-           wetlands_temp.shp$dist2NearWet[i]*1000
    giw.INFO[i, "dLe"] <-          wetlands_temp.shp$dLe[i]*1000
    giw.INFO[i,"dz"]<-             dz
  }
  
  # 2.3.4 Populate land.INFO table~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.3.4a Create Land Info Table
  land.INFO<-c("area","invert",                                              # geometric characteristics
               "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", # soil characteristics
               "slope", "kb",                                                # larger watershed charactersitics
               "y_wt_0", "s_t_0","GW_bf_0",                                  # initial conditions
               "Sy",                                                         # calculated terms
               "wetland_invert","wetland_area","volume_max"                  # lumped wetland information
  )
  
  # 2.3.4b Create land.INFO matrix
  land.INFO<-matrix(0, nrow=1, ncol=length(land.INFO), dimnames = list(c(1), c(land.INFO)))
  
  # 2.3.4c Aggregate soils data by area
  temp_soils<-soils_temp.shp
  temp_soils$area<-gArea(temp_soils, byid=T)
  temp_soils<-temp_soils@data
  temp_soils<-colSums(temp_soils[,c("y_cl","y_rd","s_fc","s_w","n","clay","ksat")]*temp_soils[,"area"])/sum(temp_soils$area, na.rm=T)
  
  # 2.3.4d Rooting depth calculation
  RD<- cellStats(root_temp.grd, mean, na.rm=T)*1000
  RD<- ifelse(RD<temp_soils["y_cl"], -1*RD, -1*temp_soils["y_cl"])
  
  # 2.3.4e Populate land.INFO matrix (lenghth units in mm)
  land.INFO[,"area"]<-            gArea(catchment_temp.shp)*(10^6)            # area in mm^2
  land.INFO[,"n"]<-               temp_soils["n"]                             # porisity
  land.INFO[,"s_fc"]<-            temp_soils["s_fc"]/100                      # soil moisture at field capacity 
  land.INFO[,"psi"]<-             -16662*(temp_soils["n"]^7.8831)             # air entry pressure head (mm) --Relationship developed from Clapp and Hornberger, 1978
  land.INFO[,"y_cl"]<-            -1*temp_soils["y_cl"]                       # confining layer depth (mm) from SSURGO
  land.INFO[,"y_c"]<-             -temp_soils["y_rd"]/2                       # critical depth (mm) 
  land.INFO[,"s_wilt"]<-          temp_soils["s_w"]/100                       # soil moisture at permanent wilting point
  land.INFO[,"k_sat"]<-           -temp_soils["ksat"]*24                      # saturated condcuctivity (mm/day)
  land.INFO[,"RD"]<-              RD                                          # Rooting Depth (mm)
  land.INFO[,"b"]<-               12.524*(temp_soils["clay"]/100)+3.6907      # Presssure Head Power-Law Coefficient (b) --Relationship developed from Clapp and Hornberger, 1978
  land.INFO[,"y_wt_0"]<-          0
  land.INFO[,"s_t_0"]<-           land.INFO[,"s_fc"]
  land.INFO[,"GW_bf_0"]<-         0
  land.INFO[,"Sy"]<-              land.INFO[,"n"]*(1-land.INFO[,"s_fc"])
  land.INFO[,"wetland_invert"]<-  -1000*0.05*(length(rowSums(area)[rowSums(area)!=max(rowSums(area), na.rm=T)])+1)
  land.INFO[,"wetland_area"]<-    max(rowSums(area))*(1000^2) 
  land.INFO[,"volume_max"]<-      max(rowSums(volume))*(1000^3)
  land.INFO[,"kb"]<-              0.046                                         #Defined fromo literature. (Going to do with Gauge analysis eventually)
  
  # 2.3.5 Populate lumped.INFO table~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.3.5a #Create lumped.INFO table
  lumped.INFO<-c("r_w","dL", "dLe") #geometric characteristics
  
  # 2.3.5b Create lumped.INFO matrix
  lumped.INFO<-matrix(0, nrow=length(wetlands_temp.shp$WetID), ncol=3, dimnames = list(c(1:length(wetlands_temp.shp$WetID)), c(lumped.INFO)))
  
  # 2.3.5c Populate lumped.INFO matrix (length in mm); data for the wetlands in the lumped upland
  lumped.INFO[, "dL"] <- wetlands_temp.shp$dist2NearWet                  # 
  lumped.INFO[,"r_w"] <- (((wetlands_temp.shp$SHAPE_Area) / pi) ^ 0.5 )  # derive radius from area assuming circular geometry, convert
  lumped.INFO[,"dLe"] <- wetlands_temp.shp$dLe                           #
  lumped.INFO <- lumped.INFO *1000                                  # convert entire matrix from m to mm
  
  # 2.3.6 Run the model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.3.6a Define the wetland
  giw.ID<-giw.INFO[,"giw.ID"][giw.INFO[,"WetID"]==main_wetland]
  
  # 2.3.6b Create function to execute WHC and organize output terms
  n.years<-1000
  execute<-function(n.years){
    #i. Run WHC Model w/ tryCatch
    output<-tryCatch(wetland.hydrology(giw.INFO,land.INFO, lumped.INFO, precip.VAR, pet.VAR, n.years, area, volume, giw.ID),
                     error = function(e) data.frame(matrix(0,ncol=390,nrow=1)))
    
    #ii. Organize output 
    if(is.list(output)==T){
      #Attach list elements to enviornment
      attach(output)
      
      #Isolate SW-GW fluxes to and from wetland
      SW_GW<-GW_local.VAR[,3]
      GWin<-ifelse(SW_GW>0, SW_GW, 0)
      GWout<-ifelse(SW_GW<0, abs(SW_GW), 0)
      
      #Isolate giw.INFO for wetland of interest
      giw.INFO<-matrix(giw.INFO[giw.ID,],nrow=1,  dimnames = list(c(1), colnames(giw.INFO)))
      
      #Calcutlate waterbalance components 
      water_balance<-data.frame(precip=sum(precip.VAR)/n.years,
                                PET=sum(pet.VAR)/n.years,
                                ET=(sum(ET_lm.VAR[,3])+sum(ET_wt.VAR[,3]))/n.years,
                                runoff_in=sum(runoff_vol.VAR[,1]/giw.INFO[,"area_wetland"]*giw.INFO[,"vol_ratio"])/n.years,
                                SW_out=sum(spill_vol.VAR[runoff_vol.VAR[,3]==0,3])/n.years/giw.INFO[,"area_wetland"],
                                GW_out=sum(SW_GW[SW_GW<0])/giw.INFO[,"area_wetland"]/n.years,
                                GW_in=sum(SW_GW[SW_GW>0])/giw.INFO[,"area_wetland"]/n.years)
      
      #Calculate mean water level for each calander day
      hydrograph<-data.frame(day=rep(seq(1,365),n.years), y_w=y_w.VAR[1:(n.years*365),3])
      hydrograph$y_w<- hydrograph$y_w + abs(giw.INFO[,"invert"])
      hydrograph$y_w<-ifelse(hydrograph$y_w>0, 
                             hydrograph$y_w/abs(giw.INFO[,"invert"]), 
                             hydrograph$y_w/abs(giw.INFO[,"y_cl"]))
      hydrograph<- hydrograph %>% group_by(day) %>% summarise(y_w = mean(y_w))
      y_w<-data.frame(t(hydrograph$y_w))
      colnames(y_w)<-hydrograph$day
      
      #Estimate duration and magnitudes
      precip_vol<-sum(precip.VAR[1:(n.years*365)])*giw.INFO[,"area_wetland"]
      shift<-precip.VAR[1:(n.years*365)]+c(0,precip.VAR[1:(n.years*365-1)])
      precip.VAR<-c(precip.VAR,0)
      shift<-c(shift,0)
      duration<-data.frame(#Quick Flow
                           QFin_duration   = length(spill_vol.VAR[shift!=0,2])/n.years/2,
                           QFin_magnitude  = (sum(spill_vol.VAR[shift!=0,2])*giw.INFO[,"vol_ratio"] +
                                              sum(runoff_vol.VAR[shift!=0,1]*giw.INFO[,"vol_ratio"]) +
                                              sum(GWin[shift!=0]))/precip_vol,
                           QFout_duration  = length(spill_vol.VAR[shift!=0 & spill_vol.VAR[,3]!=0,3])/n.years, 
                           QFout_magnitude = sum(spill_vol.VAR[shift!=0,3])/precip_vol,
                           
                           #Surface water fluxes
                           SWin_duration   = length(spill_vol.VAR[shift==0 & spill_vol.VAR[,2]>0])/n.years, 
                           SWin_magnitude  = (sum(spill_vol.VAR[shift==0,2])+sum(runoff_vol.VAR[shift==0,1]))*giw.INFO[,"vol_ratio"]/precip_vol,
                           SWout_duration  = length(spill_vol.VAR[shift==0 & spill_vol.VAR[,3]>0])/n.years,
                           SWout_magnitude = sum(spill_vol.VAR[shift==0,3])/precip_vol,
                           
                           #Groundwater Fluxes
                           GWin_duration   = length(GWin[GWin>0 & shift==0])/n.years,
                           GWin_magnitude  = sum(GWin[GWin>0 & shift==0])/precip_vol,
                           GWout_duration  = length(GWout[GWout>0])/n.years,
                           GWout_magnitude = sum(GWout[GWout>0])/precip_vol)
      #Combine data
      output<-cbind(giw.INFO, water_balance, duration, y_w)
    }
  }
  
  # 2.3.6c Run function
  WHC<-execute(n.years)
  
  # 2.3.6d Name columns
  colnames(WHC)<-c(
    #GIW Info
    "giw.ID","WetID","area_watershed","area_wetland","invert",
    "vol_ratio", "dL", "dLe", "dz","n","s_fc","psi","y_cl" ,"y_c","s_wilt","k_sat","RD", "b","Sy", "y_w_0" , "s_t_0",
    #Water Balance
    "precip","PET","ET","runoff_in","SW_out","GW_out","GW_in", 
    #duration
    "QFin_duration","QFin_magnitude","QFout_duration","QFout_magnitude","SWin_duration","SWin_magnitude","SWout_duration","SWout_magnitude","GWin_duration","GWin_magnitude","GWout_duration","GWout_magnitude",
    #Normalized Flow
    seq(1,365))
  
  # 2.3.6e Add WetID info
  WHC$WetID=main_wetland
  
  # 2.3.6f Export WHC Results
  WHC
}

# 2.4 Save output
save.image("backup/Model_Setup.Rdata")  

####################################################################################
# Step 3: Run function -------------------------------------------------------------
####################################################################################
# 3.1 Set Up workspace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
remove(list=ls())                                # clear environment
wd<-"//nfs/WHC-data/Regional_Analysis/Delmarva"  # Define working directory for later reference
setwd(wd)                                        # Set working directory
load("backup/Model_Setup.Rdata")                 # load inputs from previous processing

# 3.2 Run the Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t0<-Sys.time()
sopts <- list(partition = "sesync", time = "12:00:00" )
params<-data.frame(WetID=wetlands.shp$WetID)
delmarva<- slurm_apply(fun, params,
                       add_objects = c(
                         #Functions
                         "wetland.hydrology",
                         #Spatial data
                         "fac.grd","catchments.shp","flowlines.shp",
                         "soils.shp","wetlands.shp","dem.grd",
                         #Climate data
                         "precip.VAR", "pet.VAR"),
                       nodes = 4, cpus_per_node=8,
                       pkgs=c('sp','raster','rgdal','rgeos','dplyr'),
                       slurm_options = sopts)

# 3.3 Check job status and collect results
print_job_status(delmarva)

# 3.4 Retrieve results
results <- get_slurm_out(delmarva, outtype = "table")
cleanup_files(delmarva)
tf<-Sys.time()
tf-t0

#Export results
write.csv(results, "results.csv")
save.image("backup/results.RData")

####################################################################################
# Step 4: Plot----------------------------------------------------------------------
####################################################################################
#Setup workspace
remove(list=ls())
load("backup/results.Rdata")

#Plot time series of wetland water level~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Transform data
level<-data.frame(t(results[,41:405]))

#Filter Level
level<-level[,level[1,]<0.975]

#Setup Plotting space
par(mar=c(3.1,3.5,0.35,0.35))
par(mgp=c(1.8,0.5,0))
par(ps=12)
par(cex.axis=10/12)
par(cex.lab=14/12)

#start plot
plot(level[,1], type="n", 
     ylim=c(-0.5,1), ylab = "Normalized Water Level [m/m]", 
     xlab="Julian Day"
)

#Plot water level
abline(h=0, lwd=2, lty=2)

#Plot lines from individual wetlands
for(i in 1:length(level[1,])){
  points(level[,i], type="l",lwd=0.1, col="grey40")
}
points(rowMeans(level), type="l", col="black", lwd=4)

#Plot relevant fluxes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df<-data.frame(SW_out = results$SW_out-results$runoff_in,
               GW_out = results$GW_out,
               GW_in =  results$GW_in)
df<-df/results$precip

#Correct Data
df<-df[df$GW_out>-1000,]
df$GW_out<-abs(df$GW_out)


#Setup Plotting space
par(mar=c(3.1,3.5,0.35,0.35))
par(mgp=c(1.8,0.5,0))
par(ps=12)
par(cex.axis=10/12)
par(cex.lab=14/12)

#Start Boxplot
boxplot(df, col="grey90", pch=19, outcol="#7f7f7f7D", cex=0.5, 
        ylim=c(0,0.5),ylab="Normalized Annual Flowpath Flux [year]", 
        names = c("Net SW Outflow","GW Outflow", "GW Inflow")
)

#Plot histogram of water balance [closure]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
water_balance<-results$precip+results$SW_in+results$GW_in-results$ET-results$SW_out-results$GW_out
hist(water_balance)

#Plot byplots of duration and magnitude~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(3,2))
plot(results$QFin_duration,  results$QFin_magnitude,  log="y", xlab="QF Inflow Duration [days]", ylab = "Magnitude [Flux/Precip]", pch=19, cex=0.5, col="grey30")
plot(results$QFout_duration, results$QFout_magnitude, log="y", ylim=c(1e-6,50), xlim=c(1,365), xlab="QF outflow Duration [days]", ylab = "Magnitude [Flux/Precip]", pch=19)
plot(results$SWin_duration,  results$SWin_magnitude,  log="y", ylim=c(1e-6,50), xlim=c(1,365), xlab="SW Inflow Duration [days]", ylab = "Magnitude [Flux/Precip]", pch=19)
plot(results$SWout_duration, results$SWout_magnitude, log="y", ylim=c(1e-6,50), xlim=c(1,365), xlab="SW Outflow Duration [days]", ylab = "Magnitude [Flux/Precip]", pch=19)
plot(results$GWout_duration, results$GWout_magnitude, log="y", ylim=c(1e-6,50), xlim=c(1,365), xlab="GW Inflow Duration [days]", ylab = "Magnitude [Flux/Precip]", pch=19)
plot(results$GWin_duration,  results$GWin_magnitude,  log="y", ylim=c(1e-6,50), xlim=c(1,365), xlab="GW Outflow Duration [days]", ylab = "Magnitude [Flux/Precip]", pch=19)

