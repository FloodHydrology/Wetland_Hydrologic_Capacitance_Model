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
library(sp)                                                                 # for spatial analysis
library(raster)                                                             # for spatial analysis
library(rgdal)                                                              # for spatial analysis
library(rgeos)                                                              # for spatial analysis
library(dplyr)                                                              # for data processing
library(rslurm)                                                             # parallel computing
library(geosphere)                                                          # for spatial analysis

# 1c. Define Master Working Directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wd<-"//nfs/WHC-data/Regional_Analysis/Delmarva"                             # Define working directory for later reference
setwd(wd)                                                                   # Set working directory

# 1d. Obtain Model Inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download Model Inputs (and put everything in a uniform projection)
load("inputs/climate.Rdata")                                                # Climate data from Delmarva_HUC12_Climate_Data
wetlands.shp<-readOGR("inputs/.","NWI")                                     # Import wetland shapefile
    wetlands.shp@data$WetID<-seq(1,length(wetlands.shp))                    # add (?) WetID column to shp file using sequence from 1:total num of wetlands
HUC12.shp<-readOGR("inputs/.","HUC12")                                      # Import HUC 12 shapefile
    HUC12.shp<-spTransform(HUC12.shp, wetlands.shp@proj4string)             # transform HUC 12 shapefile's projection into same as wetlands
catchments.shp<-readOGR("inputs/.","NHD_catchments")                        # import NHD catchments
    catchments.shp<-spTransform(catchments.shp, wetlands.shp@proj4string)   # transform NHD catchment's projection into same as wetlands
flowlines.shp<-readOGR("inputs/.","NHDplus")                                # import NHD flowlines
    flowlines.shp<-spTransform(flowlines.shp, wetlands.shp@proj4string)     # transform NHD flowline's projection into same as wetlands
fac.grd<-raster("inputs/fac")                                               # import flow accumulation raster
soils.shp<-readOGR("inputs/.","soils")                                      # import soils shapefile
    soils.shp<-spTransform(soils.shp, wetlands.shp@proj4string)             # transform soil file's projection into same as wetlands
    soils<-read.csv("inputs/WHC_Soils_Input.csv")                           # import existing soil parameters csv
    soils.shp@data<-merge(soils.shp@data,soils, by.x='MUKEY', by.y="MUID")  # append soils csv data into soils shapefile by matching MUKEY and MUID
    remove(soils)                                                           # delete soils df
dem.grd<-raster("inputs/NHDPlus02/Elev_Unit_a/elev_cm")                     # import DEM file
    mask<-spTransform(catchments.shp, dem.grd@crs)                          # transform file's projection
    dem.grd<-crop(dem.grd, mask)                                            # crop out portion of dem to keep relevant areas only
    remove(mask)                                                            # remove dummy variable mask
    dem.grd<-projectRaster(dem.grd, crs=catchments.shp@proj4string)         # project raster

# 1e. Delineate "Isolated" Wetlands  (e.g. atleast 1000 ft from the stream) ~~~~~~~~~~~
wetlands.shp$dist2stream<-apply(gDistance(flowlines.shp,wetlands.shp, byid=TRUE),1,min) # Calculate euclidean distance to stream network
wetlands.shp<-wetlands.shp[wetlands.shp$dist2stream>304.8,]                             # Remove wetlands <1000 ft from stream
wetlands.shp$area_m2<-gArea(wetlands.shp,byid=T)                                        # Caclculate wetland area
    
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
save.image("Inputs.Rdata")                                                  # save workspace for use later

####################################################################################
# Step 2: Create Function to run individual wetlands--------------------------------
####################################################################################
# 2.1 Set Up workspace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
remove(list=ls())                                                           # clear environment
load("Inputs.Rdata")                                                        # load inputs from previous processing
source("~/Wetland_Hydrologic_Capacitance_Model/R/WHC_2.R")                  # compile WHC function 
source("~/Wetland_Hydrologic_Capacitance_Model/R/get_yc.R")                 # compile yc 

# 2.2 Create function to process data and run WHC for a wetland of interest~~~~~~~~~
fun<-function(WetID){ #
  
  # 2.2.1 Gather Required Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.2.1a Identify wetland of interest
  main_wetland.shp<-wetlands.shp[wetlands.shp$WetID==WetID,]
  main_wetland<-WetID
  
  # 2.2.1b Select spatial layers based on WetID  
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
  
  # 2.2.1c For now, add catchment aggregate soils data to soils with missing data
  soils_temp.shp@data$y_cl[is.na(soils_temp.shp@data$y_cl)]<-mean(soils_temp.shp@data$y_cl, na.rm=T)
  soils_temp.shp@data$y_rd[is.na(soils_temp.shp@data$y_rd)]<-mean(soils_temp.shp@data$y_rd, na.rm=T)
  soils_temp.shp@data$s_fc[is.na(soils_temp.shp@data$s_fc)]<-mean(soils_temp.shp@data$s_fc, na.rm=T)
  soils_temp.shp@data$s_w[is.na(soils_temp.shp@data$s_w)]<-mean(soils_temp.shp@data$s_w, na.rm=T)
  soils_temp.shp@data$n[is.na(soils_temp.shp@data$n)]<-mean(soils_temp.shp@data$n, na.rm=T)
  soils_temp.shp@data$clay[is.na(soils_temp.shp@data$clay)]<-mean(soils_temp.shp@data$clay, na.rm=T)
  soils_temp.shp@data$ksat[is.na(soils_temp.shp@data$ksat)]<-mean(soils_temp.shp@data$ksat, na.rm=T)
  
  # 2.2.2 Estimate area and volume to stage relationships~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.2.2a Create vectors to store area and volume to stage relationships
  area<-matrix(0, ncol=length(wetlands_temp.shp), nrow=100)
  volume<-matrix(0, ncol=length(wetlands_temp.shp), nrow=100)
  
  # 2.2.2b Use loop to calculate based on previously published relationships
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
  
  # 2.2.3 Populate GIW.INFO Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.2.3a Create input variables
  giw.INFO<-c("giw.ID","WetID","area_watershed","area_wetland","invert","vol_ratio", "dL", "dLe",  # geometric characteristics
              "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", "Sy",                  # soil characteristics
              "y_w_0", "s_t_0"                                                                     # initial conditions
  )
  
  # 2.2.3b Create giw.INFO matrix
  giw.INFO<-matrix(0, nrow=length(wetlands_temp.shp), ncol=length(giw.INFO), 
                   dimnames = list(seq(1,length(wetlands_temp.shp)), c(giw.INFO)))
  
  # 2.2.3c Define wetland variables (units == mm)
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
    
    # iv. Populate giw.INFO table
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
    giw.INFO[i,"RD"]<-             -temp_soils["y_rd"]                         # Rooting Depth (mm)
    giw.INFO[i,"b"]<-              12.524*(temp_soils["clay"]/100)+3.6907 
    giw.INFO[i,"Sy"]<-             giw.INFO[i,"n"]*(1-giw.INFO[i,"s_fc"])
    giw.INFO[i,"y_w_0"]<-          0
    giw.INFO[i,"s_t_0"]<-          giw.INFO[i,"s_fc"]
    giw.INFO[i,"vol_ratio"]<-      temp_fac                                   # ratio of upstream wetland volume
    giw.INFO[i, "dL"] <-           wetlands_temp.shp$dist2NearWet[i]*1000
    giw.INFO[i, "dLe"] <-         200*1000
  }
  
  # 2.2.4 Populate land.INFO table~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.2.4a Create Land Info Table
  land.INFO<-c("area","invert",                                              # geometric characteristics
               "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", # soil characteristics
               "slope", "kb",                                                # larger watershed charactersitics
               "y_wt_0", "s_t_0","GW_bf_0",                                  # initial conditions
               "Sy",                                                         # calculated terms
               "wetland_invert","wetland_area","volume_max"                  # lumped wetland information
  )
  
  # 2.2.4b Create land.INFO matrix
  land.INFO<-matrix(0, nrow=1, ncol=length(land.INFO), dimnames = list(c(1), c(land.INFO)))
  
  # 2.2.4c Aggregate soils data by area
  temp_soils<-soils_temp.shp
  temp_soils$area<-gArea(temp_soils, byid=T)
  temp_soils<-temp_soils@data
  temp_soils<-colSums(temp_soils[,c("y_cl","y_rd","s_fc","s_w","n","clay","ksat")]*temp_soils[,"area"])/sum(temp_soils$area, na.rm=T)
  
  # 2.2.4d Populate land.INFO matrix (lenghth units in mm)
  land.INFO[,"area"]<-            gArea(catchment_temp.shp)*(10^6)            # area in mm^2
  land.INFO[,"n"]<-               temp_soils["n"]                             # porisity
  land.INFO[,"s_fc"]<-            temp_soils["s_fc"]/100                      # soil moisture at field capacity 
  land.INFO[,"psi"]<-             -16662*(temp_soils["n"]^7.8831)             # air entry pressure head (mm) --Relationship developed from Clapp and Hornberger, 1978
  land.INFO[,"y_cl"]<-            -1*temp_soils["y_cl"]                       # confining layer depth (mm) from SSURGO
  land.INFO[,"y_c"]<-             -temp_soils["y_rd"]/2                       # critical depth (mm) 
  land.INFO[,"s_wilt"]<-          temp_soils["s_w"]/100                       # soil moisture at permanent wilting point
  land.INFO[,"k_sat"]<-           -temp_soils["ksat"]*24                      # saturated condcuctivity (mm/day)
  land.INFO[,"RD"]<-              -temp_soils["y_rd"]                         # Rooting Depth (mm)
  land.INFO[,"b"]<-               12.524*(temp_soils["clay"]/100)+3.6907      # Presssure Head Power-Law Coefficient (b) --Relationship developed from Clapp and Hornberger, 1978
  land.INFO[,"y_wt_0"]<-          0
  land.INFO[,"s_t_0"]<-           land.INFO[,"s_fc"]
  land.INFO[,"GW_bf_0"]<-         300
  land.INFO[,"Sy"]<-              land.INFO[,"n"]*(1-land.INFO[,"s_fc"])
  land.INFO[,"wetland_invert"]<-  -1000*0.05*(length(rowSums(area)[rowSums(area)!=max(rowSums(area), na.rm=T)])+1)
  land.INFO[,"wetland_area"]<-    max(rowSums(area))*(1000^2) 
  land.INFO[,"volume_max"]<-      max(rowSums(volume))*(1000^3)
  land.INFO[,"kb"]<-              0.046                                         #Defined fromo literature. (Going to do with Gauge analysis eventually)

  
  # 2.2.5 Populate lumped.INFO table~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.2.5a #Create lumped.INFO table
  lumped.INFO<-c("r_w","dL", "dLe") #geometric characteristics
  
  # 2.2.5b Create lumped.INFO matrix
  lumped.INFO<-matrix(0, nrow=length(wetlands.shp$WetID), ncol=3, dimnames = list(c(1:length(wetlands.shp$WetID)), c(lumped.INFO)))
  
  # 2.2.5c Populate lumped.INFO matrix (length in mm); data for the wetlands in the lumped upland
  lumped.INFO[, "dL"] <- wetlands.shp$dist2NearWet                  # 
  lumped.INFO[,"r_w"] <- (((wetlands.shp$SHAPE_Area) / pi) ^ 0.5 )  # derive radius from area assuming circular geometry, convert
  lumped.INFO[,"dLe"] <- 200                 # please change me!
  lumped.INFO <- lumped.INFO *1000           # convert entire matrix from m to mm
  
  # 2.2.6 Run the model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.2.6a Define the wetland
  giw.ID<-giw.INFO[,"giw.ID"][giw.INFO[,"WetID"]==main_wetland]
  
  # 2.2.6b Create function to execute WHC and organize output terms
  n.years<-1000
  execute<-function(n.years){
    #i. Run WHC Model w/ tryCatch
    output<-tryCatch(wetland.hydrology(giw.INFO,land.INFO, lumped.INFO, precip.VAR, pet.VAR, n.years, area, volume, giw.ID),
                     error = function(e) error = function(e) data.frame(matrix(0,ncol=390,nrow=1)))
    
    #ii. Organize output 
    if(is.list(output)==T){
      #Attach list elements to enviornment
      attach(output)
      
      #Isolate SW-GW fluxes to and from wetland
      SW_GW<-GW_local.VAR[,3]
      
      #Isolate giw.INFO for wetland of interest
      giw.INFO<-matrix(giw.INFO[giw.ID,],nrow=1,  dimnames = list(c(1), colnames(giw.INFO)))
      
      #Calcutlate waterbalance components 
      water_balance<-data.frame(precip=sum(precip.VAR)/1000,
                                PET=sum(pet.VAR)/1000,
                                ET=(sum(ET_lm.VAR[,3])+sum(ET_wt.VAR[,3]))/1000,
                                SW_out=sum(spill_vol.VAR[,3])/1000/giw.INFO[,"area_wetland"],
                                SW_in=sum(runoff_vol.VAR[,1]/giw.INFO[,"area_wetland"]*giw.INFO[,"vol_ratio"])/1000,
                                GW_out=sum(SW_GW[SW_GW<0])/giw.INFO[,"area_wetland"]/1000,
                                GW_in=sum(SW_GW[SW_GW>0])/giw.INFO[,"area_wetland"]/1000)
      
      #Calculate mean water level for each calander day
      hydrograph<-data.frame(day=rep(seq(1,365),1000), y_w=y_w.VAR[1:365000,3])
      hydrograph$y_w<- hydrograph$y_w + abs(giw.INFO[,"invert"])
      hydrograph$y_w<-ifelse(hydrograph$y_w>0, 
                             hydrograph$y_w/abs(giw.INFO[,"invert"]), 
                             hydrograph$y_w/abs(giw.INFO[,"y_cl"]))
      hydrograph<- hydrograph %>% group_by(day) %>% summarise(y_w = mean(y_w))
      y_w<-data.frame(t(hydrograph$y_w))
      colnames(y_w)<-hydrograph$day
      
      #Combine data
      output<-cbind(giw.INFO, water_balance, y_w)
    }
  }
  
  # 2.2.6c Run function
  WHC<-execute(1000)
  
  # 2.2.6d Name columns
  colnames(WHC)<-c(
    #GIW Info
    "giw.ID","WetID","area_watershed","area_wetland","invert",
    "vol_ratio", "dL", "dLe", "n","s_fc","psi","y_cl" ,"y_c","s_wilt","k_sat","RD", "b","Sy", "y_w_0" , "s_t_0",
    #Water Balance
    "precip","PET","ET","SW_out","SW_in","GW_out","GW_in", 
    #Normalized Flow
    seq(1,365))
  
  # 2.2.6e Add WetID info
  WHC$WetID=main_wetland
  
  # 2.2.6f Export WHC Results
  WHC
}

####################################################################################
# Step 3: Run function -------------------------------------------------------------
####################################################################################
# for testing
plot(catchments.shp)
plot(catchments.shp[53,1], col = 'red', add = T)

catchments.shp <- catchments.shp[53,]
wetlands.shp <- wetlands.shp[catchments.shp,]

plot(wetlands.shp, add = T)

#----------------------------------------------------------------------------------
#Run function using SLURM server!
#Try running on server
sopts <- list(partition = "sesync")
params<-data.frame(WetID=wetlands.shp$WetID)
delmarva3<- slurm_apply(fun, params,
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

# Check job status and collect results
print_job_status(delmarva3)
save.image("results_intermediate.RData")

#------------------------------------------------------------------------
#Load data
remove(list=ls())
load("results_intermediate.RData")

#Get job
#delmarva<-slurm_job(jobname = "3a9025fe1fbb",nodes=4,jobid=291519)
results <- get_slurm_out(delmarva3, outtype = "table")
write.csv(results, "results.csv")
save.image("results.RData")
cleanup_files(delmarva2)

####################################################################################
# Step 4: Plot-------- -------------------------------------------------------------
####################################################################################
#Setup workspace
remove(list=ls())
load("results.Rdata")

#Plot time series of wetland water level~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Transform data
level<-data.frame(t(results[,26:390]))

#Setup Plotting space
par(mar=c(3.1,3.5,0.35,0.35))
par(mgp=c(1.8,0.5,0))
par(ps=12)
par(cex.axis=10/12)
par(cex.lab=14/12)

#start plot
plot(level[,1], type="n", 
     ylim=c(-0.5,0.5), ylab = "Normalized Water Level [m/m]", 
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
#Correct Data
results<-results[results$GW_out>-1000,]
results$GW_out<-abs(results$GW_out)


#Setup Plotting space
par(mar=c(3.1,3.5,0.35,0.35))
par(mgp=c(1.8,0.5,0))
par(ps=12)
par(cex.axis=10/12)
par(cex.lab=14/12)

#Start Boxplot
boxplot(results[,24:27]/results$precip, col="grey90", pch=19, outcol="#7f7f7f7D", cex=0.5,
        ylim=c(0,6), ylab="Normalized Annual Flowpath Flux [mm/year]", 
        names = c("SW Outflow","SW Inflow", "GW Outflow", "GW Inflow")
)



water_balance<-results$precip+results$SW_in+results$GW_in-results$ET-results$SW_out-results$GW_out

