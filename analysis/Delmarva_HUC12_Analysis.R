###################################################################################
#Name: Sensitivity Analysis Dem
#Coder: C. Nathan Jones
#Date: 5/18/2017
#Purpose: Demo sensitivity analysis for individual wetland manuscript
##################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
#Clear Memory
rm(list=ls(all=TRUE))

#Define Master Working Directory
wd<-"//nfs/WHC-data/Regional_Analysis/Delmarva"
setwd(wd)

#Load Required Packages
library(sp)       #spatial analysis
library(raster)   #spatial analysis
library(rgdal)    #spatial analysis
library(rgeos)    #spatial analysis
library(dplyr)    #data processing
library(rslurm)   #paralel computing

#Download Model Inputs (and put everything in a uniform projection)
load("inputs/climate.Rdata") #Climate data
wetlands.shp<-readOGR("inputs/.","NWI")
  wetlands.shp@data$WetID<-seq(1,length(wetlands.shp))
HUC12.shp<-readOGR("inputs/.","HUC12")
  HUC12.shp<-spTransform(HUC12.shp, wetlands.shp@proj4string)
catchments.shp<-readOGR("inputs/.","NHD_catchments")
  catchments.shp<-spTransform(catchments.shp, wetlands.shp@proj4string)
flowlines.shp<-readOGR("inputs/.","NHDplus")
  flowlines.shp<-spTransform(flowlines.shp, wetlands.shp@proj4string)
fac.grd<-raster("inputs/fac")
soils.shp<-readOGR("inputs/.","soils")
  soils.shp<-spTransform(soils.shp, wetlands.shp@proj4string)
  soils<-read.csv("inputs/WHC_Soils_Input.csv")
  soils.shp@data<-merge(soils.shp@data,soils, by.x='MUKEY', by.y="MUID")
  remove(soils)
dem.grd<-raster("inputs/NHDPlus02/Elev_Unit_a/elev_cm")
  mask<-spTransform(catchments.shp, dem.grd@crs)
  dem.grd<-crop(dem.grd, mask)
  remove(mask)
  dem.grd<-projectRaster(dem.grd, crs=catchments.shp@proj4string)

#Delineate "Isolated" Wetlands  (e.g. atleast 1000 ft from the stream)
wetlands.shp$dist2stream<-apply(gDistance(flowlines.shp,wetlands.shp, byid=TRUE),1,min) #Calculate euclidean distance to stream network
wetlands.shp<-wetlands.shp[wetlands.shp$dist2stream>304.8,] #Remove wetlands <1000 ft from stream
wetlands.shp$area_m2<-gArea(wetlands.shp,byid=T)  #Caclculate wetland area

#Plot for funzies 
plot(HUC12.shp, lwd=2, border="grey10")
plot(catchments.shp, lwd=0.5, border="grey50", add=T)
plot(flowlines.shp, lwd=0.5, col="dark blue", add=T)
plot(wetlands.shp, col="dark blue", add=T)

#Save Image
save.image("Inputs.Rdata")

####################################################################################
# Step 2: Create Function to run individual wetlands--------------------------------
####################################################################################
#Setup workspace
remove(list=ls())
load("Inputs.Rdata")
source("~/Wetland_Hydrologic_Capacitance_Model/R/WHC_2.R")

#Initiate Function 
fun<-function(WetID){
  
  #Gather Required data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Identify wetland of interest
  main_wetland.shp<-wetlands.shp[wetlands.shp$WetID==WetID,]
  main_wetland<-WetID
  
  #Select spatial layers based on WetID  
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
  
  #For now, add catchment aggregate soils data to soils with missing data
  soils_temp.shp@data$y_cl[is.na(soils_temp.shp@data$y_cl)]<-mean(soils_temp.shp@data$y_cl, na.rm=T)
  soils_temp.shp@data$y_rd[is.na(soils_temp.shp@data$y_rd)]<-mean(soils_temp.shp@data$y_rd, na.rm=T)
  soils_temp.shp@data$s_fc[is.na(soils_temp.shp@data$s_fc)]<-mean(soils_temp.shp@data$s_fc, na.rm=T)
  soils_temp.shp@data$s_w[is.na(soils_temp.shp@data$s_w)]<-mean(soils_temp.shp@data$s_w, na.rm=T)
  soils_temp.shp@data$n[is.na(soils_temp.shp@data$n)]<-mean(soils_temp.shp@data$n, na.rm=T)
  soils_temp.shp@data$clay[is.na(soils_temp.shp@data$clay)]<-mean(soils_temp.shp@data$clay, na.rm=T)
  soils_temp.shp@data$ksat[is.na(soils_temp.shp@data$ksat)]<-mean(soils_temp.shp@data$ksat, na.rm=T)
  
  #Estimate area and volume to stage relationships~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create vectors to store area and volume to stage relationships
  area<-matrix(0, ncol=length(wetlands_temp.shp), nrow=100)
  volume<-matrix(0, ncol=length(wetlands_temp.shp), nrow=100)
  
  #Use loop to calculate based on previously published relationships
  for(i in 1:length(wetlands_temp.shp)){
    #Define WetID
    WetID<-wetlands_temp.shp$WetID[i]
    
    #Define Amax
    Amax<-wetlands_temp.shp$area_m2[wetlands_temp.shp$WetID==WetID]
    
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
    n.col<-which(wetlands_temp.shp$WetID==WetID)
    n.rows<-length(area.fun(seq(0,hmax,0.05)))
    area[1:n.rows,n.col]<-area.fun(seq(0,hmax,0.05))
    area[(n.rows+1):100,n.col]<-area.fun(hmax)
    volume[1:n.rows,n.col]<-volume.fun(seq(0,hmax,0.05))
    volume[(n.rows+1):100,n.col]<-volume.fun(hmax)
  }
  
  
  #Populate GIW.INFO Table~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #create input variables
  giw.INFO<-c("giw.ID","WetID","area_watershed","area_wetland","invert","vol_ratio", #geometric characteristics
              "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", "Sy", #soil characteristics
              "y_w_0", "s_t_0" #initial conditions
  )
  
  #Create giw.INFO matrix
  giw.INFO<-matrix(0, nrow=length(wetlands_temp.shp), ncol=length(giw.INFO), 
                   dimnames = list(seq(1,length(wetlands_temp.shp)), c(giw.INFO)))
  
  #Define wetland variables (units == mm)
  for(i in 1:length(wetlands_temp.shp)){
    
    #Isolate soils data (if for osme reason we don't have ovlerlap, just use basin agregate)
    if(gIntersects(soils_temp.shp,wetlands_temp.shp[i,])==TRUE){
      temp_soils<-crop(soils_temp.shp,wetlands_temp.shp[i,])
    }else{
      temp_soils<-soils_temp.shp
    }
    
    #Agregrate parameters based on space
    temp_soils$area<-gArea(temp_soils, byid=T)
    temp_soils<-temp_soils@data
    temp_soils<-colSums(temp_soils[,c("y_cl","y_rd","s_fc","s_w","n","clay","ksat")]*temp_soils[,"area"])/sum(temp_soils$area, na.rm=T)
    
    #Isolate fac data and calculate ratio of upland drainage area to total drainage area
    temp_fac<-crop(fac_temp.grd, wetlands_temp.shp[i,])
    temp_fac<-mask(temp_fac, wetlands_temp.shp[i,])
    temp_fac<-ifelse(cellStats(temp_fac*0+1, sum)>0,
                     cellStats(temp_fac,max)/cellStats(fac_temp.grd,max),
                     0)
    
    #Populate giw.INFO table
    giw.INFO[i,"giw.ID"]<-         i  #ID used in the model, note this is differnt than the EetID
    giw.INFO[i,"WetID"]<-          wetlands_temp.shp$WetID[i]
    giw.INFO[i,"area_watershed"]<- gArea(catchment_temp.shp)*(10^6)
    giw.INFO[i,"area_wetland"]<-   gArea(wetlands_temp.shp[i,], byid=T)*(10^6)
    giw.INFO[i,"invert"]<-         -1000*0.05*(length(area[,i][area[,i]!=max(area[,i], na.rm=T)]))   
    giw.INFO[i,"n"]<-              temp_soils["n"]
    giw.INFO[i,"s_fc"]<-           temp_soils["s_fc"]/100  
    giw.INFO[i,"psi"]<-            -16662*(temp_soils["n"]^7.8831) 
    giw.INFO[i,"y_cl"]<-           -1*temp_soils["y_cl"]   
    giw.INFO[i,"y_c"]<-           -temp_soils["y_rd"]/2                       #critical depth (mm) 
    giw.INFO[i,"s_wilt"]<-        temp_soils["s_w"]/100                       #soil moisture at permanent wilting point
    giw.INFO[i,"k_sat"]<-         -temp_soils["ksat"]*24                      #saturated condcuctivity (mm/day)
    giw.INFO[i,"RD"]<-            -temp_soils["y_rd"]                         #Rooting Depth (mm)
    giw.INFO[i,"b"]<-              12.524*(temp_soils["clay"]/100)+3.6907 
    giw.INFO[i,"Sy"]<-             giw.INFO[i,"n"]*(1-giw.INFO[i,"s_fc"])
    giw.INFO[i,"y_w_0"]<-          0
    giw.INFO[i,"s_t_0"]<-          giw.INFO[i,"s_fc"]
    giw.INFO[i,"vol_ratio"]<-      temp_fac #ratio of upstream wetland volume
  }
  
  #Populate land.INFO table~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  land.INFO<-c("area","invert", #geometric characteristics
               "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", #soil characteristics
               "slope", "kb",#larger watershed charactersitics
               "y_wt_0", "s_t_0","GW_bf_0", #initial conditions
               "Sy", #calculated terms
               "wetland_invert","wetland_area","volume_max" #lumped wetland information
  )
  
  #Create land.INFO matrix
  land.INFO<-matrix(0, nrow=1, ncol=length(land.INFO), dimnames = list(c(1), c(land.INFO)))
  
  #Aggregate soils data by area
  temp_soils<-soils_temp.shp
  temp_soils$area<-gArea(temp_soils, byid=T)
  temp_soils<-temp_soils@data
  temp_soils<-colSums(temp_soils[,c("y_cl","y_rd","s_fc","s_w","n","clay","ksat")]*temp_soils[,"area"])/sum(temp_soils$area, na.rm=T)
  
  #Populate land.INFO matrix (lenghth units in mm)
  land.INFO[,"area"]<-            gArea(catchment_temp.shp)*(10^6)         #area in mm^2
  land.INFO[,"n"]<-               temp_soils["n"]                             #porisity
  land.INFO[,"s_fc"]<-            temp_soils["s_fc"]/100                      #soil moisture at field capacity 
  land.INFO[,"psi"]<-             -16662*(temp_soils["n"]^7.8831)             #air entry pressure head (mm) --Relationship developed from Clapp and Hornberger, 1978
  land.INFO[,"y_cl"]<-            -1*temp_soils["y_cl"]                       #confining layer depth (mm) from SSURGO
  land.INFO[,"y_c"]<-             -temp_soils["y_rd"]/2                       #critical depth (mm) 
  land.INFO[,"s_wilt"]<-          temp_soils["s_w"]/100                       #soil moisture at permanent wilting point
  land.INFO[,"k_sat"]<-           -temp_soils["ksat"]*24                      #saturated condcuctivity (mm/day)
  land.INFO[,"RD"]<-              -temp_soils["y_rd"]                         #Rooting Depth (mm)
  land.INFO[,"b"]<-               12.524*(temp_soils["clay"]/100)+3.6907      #Presssure Head Power-Law Coefficient (b) --Relationship developed from Clapp and Hornberger, 1978
  land.INFO[,"y_wt_0"]<-          0
  land.INFO[,"s_t_0"]<-           land.INFO[,"s_fc"]
  land.INFO[,"GW_bf_0"]<-         300
  land.INFO[,"Sy"]<-              land.INFO[,"n"]*(1-land.INFO[,"s_fc"])
  land.INFO[,"wetland_invert"]<-  -1000*0.05*(length(rowSums(area)[rowSums(area)!=max(rowSums(area), na.rm=T)])+1)
  land.INFO[,"wetland_area"]<-    max(rowSums(area))*(1000^2) 
  land.INFO[,"volume_max"]<-      max(rowSums(volume))*(1000^3)
  land.INFO[,"kb"]<-              0.046                                         #Defined fromo literatture. (Going to odo with Gauge analysis eventually)
  
  #Run the model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Define the wetland
  giw.ID<-giw.INFO[,"giw.ID"][giw.INFO[,"WetID"]==main_wetland]
  
  #run model  
  n.years<-1000
  execute<-function(n.years){
    tryCatch(wetland.hydrology(giw.INFO,land.INFO, precip.VAR, pet.VAR, n.years, area, volume, giw.ID),
             error = function(e) data.frame(matrix(0,ncol=390,nrow=1)))
  }
  WHC<-execute(1000)
  
  #Name columns
  colnames(WHC)<-c(#GIW Info
    "giw.ID","WetID","area_watershed","area_wetland","invert",
    "vol_ratio","n","s_fc","psi","y_cl" ,"y_c","s_wilt","k_sat","RD", "b","Sy", "y_w_0" , "s_t_0",
    #Water Balance
    "precip","PET","ET","SW_out","SW_in","GW_out","GW_in", 
    #Normalized Flow
    seq(1,365))
  
  #Add WetID info
  WHC$WetID=main_wetland
  
  #Export WHC Results
  WHC
}

####################################################################################
# Step 3: Run function -------------------------------------------------------------
####################################################################################
#----------------------------------------------------------------------------------
#Run function using SLURM server!
#Try running on server
sopts <- list(partition = "sesync")
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

# Check job status and collect results
print_job_status(delmarva)
save.image("results_intermediate.RData")

#------------------------------------------------------------------------
#Load data
remove(list=ls())
load("results_intermediate.RData")

#Get job
#delmarva<-slurm_job(jobname = "3a9025fe1fbb",nodes=4,jobid=291519)
results <- get_slurm_out(delmarva, outtype = "table")
write.csv(results, "results.csv")
save.image("results.RData")
cleanup_files(delmarva)

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
boxplot(results[,22:25]/results$precip, col="grey90", pch=19, outcol="#7f7f7f7D", cex=0.5,
        ylim=c(0,6), ylab="Normalized Annual Flowpath Flux [mm/year]", 
        names = c("SW Outflow","SW Inflow", "GW Outflow", "GW Inflow")
)



water_balance<-results$precip+results$SW_in+results$GW_in-results$ET-results$SW_out-results$GW_out

