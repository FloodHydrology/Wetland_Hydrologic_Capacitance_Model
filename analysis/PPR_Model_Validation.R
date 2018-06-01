###################################################################################
#Name: Delmarva Wetland Simulations (Model Validation)
#Coder: C. Nathan Jones
#Date: 1/22/2018
#Purpose: Validate the WHC model at Baltimore Corner TNC site
##################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
#Clear Memory
rm(list=ls(all=TRUE))

#Define working and data directories
wd<-"Wetland_Hydrologic_Capacitance_Model"
dir<-"//nfs/WHC-data/Validation_Modeling/WHC_BaltimoreCorner"

#Load Required Packages
library(raster)   #spatial analysis
library(rgdal)    #spatial analysis
library(rgeos)    #spatial analysis
library(maptools) #spatial analysis
library(plyr)     #data processing
library(dplyr)    #data processing
library(Evapotranspiration)

#Define projection for the project
proj<-CRS("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") 

#Obtain Model Inputs
dem<-raster(readGDAL(paste0(dir,"/Model Inputs/DEM/dem")))      #DEM (from USDA-Amir Sharifi)
soils.shp<-readOGR(dsn = paste0(dir,'/Model Inputs/Soils'), layer='SSURGO')  #SSURGO data
climate<-read.csv(paste0(dir,"/Model Inputs/Climate/ncdc.csv"))

#Convert all data sources into correct format
dem<-projectRaster(dem, crs=proj)
soils.shp<-spTransform(soils.shp, proj)

#Add pour point for watershed  
pp<-data.frame(x=947025.897981, y=4335592.83542)
coordinates(pp)<-~x+y
projection(pp)<-proj
pp.shp<-SpatialPoints(pp)
pp.shp<-SpatialPointsDataFrame(pp, data.frame(x=947025.897981, y=4335592.83542))

#Identify wetland of interest 
wetland<-data.frame(x=947057.170661, y=4335552.88131)
coordinates(wetland)<-~x+y
projection(wetland)<-proj
wetland.shp<-SpatialPoints(wetland)
wetland.shp<-SpatialPointsDataFrame(wetland, data.frame(x=947025.897981, y=4335592.83542))

#Plot to make sure there is overlap
plot(dem)
plot(soils.shp, border="grey60", cex=0.25, add=T)
plot(pp.shp, pch=19, col="red", cex=1, add=T)
plot(wetland.shp, pch=19, col="green", cex=1, add=T)

#Save to backup folder
save.image(paste0(dir, "/Backup/input_data.RData"))
  
####################################################################################
# Step 2: Topographic Analysis-------------------------------------------------------
####################################################################################
#To complete the topographic analysis, we use the RPyGeo package. It is a wrapper for 
#ArcGIS and Python geoprocessing capabilities.  Unfortunately, ArcGIS is not available
#on the SESYNC server, so this step has to be run on your local machine. Use the steps below. 
  #1) Identify python path (see rpygeo.build.env function help)
  #2) Map network drive to WHC-data folder
  #3) Clone GitHub repository to your machine [eg "C://GIT_Workspace"]
  #4) Create a "scratch workspace" folder to complete geoprocessing [e.g., "C://ScratchWorkspace"]
  #5) Run "Topo_Analysis" function
  #6) Be cautious of having simultaneous "git" folders. 

####################################################################################
# Step 3: Soils Analysis------------------------------------------------------------
####################################################################################
#Clear Memory
rm(list=ls(all=TRUE))

#Define working and data directories
wd3<-"~/Wetland_Hydrologic_Capacitance_Model"
dir3<-"//nfs/WHC-data/Validation_Modeling/WHC_BaltimoreCorner"

#Load soil database
soil.data<-read.csv(paste0(dir3,"/Model Inputs/Soils/WHC_Soils_Input.csv"))

#Load data from step 2
load(paste0(dir3, "/Backup/DEM_Processing.RData"))

#manipulate soils data
soils.shp<-spTransform(soils.shp, dem@crs) #transform the layer
soils.shp<-raster::intersect(soils.shp,watershed.shp) #clip to watershed
soils.shp@data$area_m2<-gArea(soils.shp, byid=T) #calculate area

#merge soils data with SSURGO database
soils<-soils.shp@data 
soils<-soils[,c("MUKEY","area_m2")]
soil.data$MUKEY<-paste(soil.data$MUID)
soils<-left_join(soils, soil.data)
remove(soil.data)

#Aggregate Soils Data
fun<-function(x){sum(soils[,x]*soils[,"area_m2"])/sum(soils[,"area_m2"])}
soils<-na.omit(soils)
soils<-c(fun("y_cl"), fun("y_rd"), fun("s_fc"), fun("s_w"), fun("n"), fun("clay"), fun("ksat"))
soils<-data.frame(soils)
rownames(soils)<-c("y_cl","y_rd","s_fc","s_w", "n","clay","ksat")
soils<-data.frame(t(soils))

#specific yeild (no depth-Sy relationship for now...maybe later?)
soils$Sy_soil<-soils$n*(1-soils$s_fc/100)

####################################################################################
# Step 4: Climate Data--------------------------------------------------------------
####################################################################################
#format climate data
climate$year<-as.numeric(substr(climate$DATE,1,4))
climate<-climate[climate$year>2007,]

#precip data
precip.VAR<-climate$PRCP
precip.VAR[precip.VAR<0]<-0

#process ET data
#create climatedata (first input file)
climatedata<-data.frame(substr(climate$DATE,1,4))
colnames(climatedata)<-"Year"
climatedata$Month<-as.numeric(substr(climate$DATE,5,6))
climatedata$Day<-as.numeric(substr(climate$DATE,7,8))
climatedata$Tmax<-(climate$TMAX)
climatedata$Tmin<-(climate$TMIN)
climatedata$RHmax<-0
climatedata$RHmin<-0
climatedata$Rs<-0
#create timeseries input file (second input file)
input<-ReadInputs(varnames = colnames(climatedata),
                  climatedata=climatedata,
                  stopmissing=c(10,10,3), 
                  timestep = "daily", 
                  interp_missing_days = F, 
                  interp_missing_entries = F,
                  interp_abnormal = F,
                  missing_method="DoY average",
                  abnormal_method="DoY average")
#create constants input file (third and final file)
data("constants")
constants$Elev<-10
constants$lat_rad<-38.8846*pi/180
#calculate ET
df<-ET.HargreavesSamani(data=input, 
                        constants=constants, 
                        ts="daily")
pet.VAR<-df$ET.Daily
pet.VAR[pet.VAR<0]<-0
pet.VAR[is.na(pet.VAR)==T]<-0

####################################################################################
# Step 5: Populate land.INFO and giw.INFO tables------------------------------------
####################################################################################
#create upland input variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
land.INFO<-c("area","invert", #geometric characteristics
             "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", #soil characteristics
             "slope", "kb",#larger watershed charactersitics
             "y_wt_0", "s_t_0","GW_bf_0", #initial conditions
             "Sy", #calculated terms
             "wetland_invert","wetland_area","volume_max" #lumped wetland information
)

#Create giw.INFO matrix
land.INFO<-matrix(0, nrow=1, ncol=length(land.INFO), dimnames = list(c(1), c(land.INFO)))

#Populate land.INFO matrix (lenghth units in mm)
land.INFO[,"area"]<-            gArea(watershed.shp)*(10^6)         #area in mm^2
land.INFO[,"n"]<-               soils$n                             #porisity
land.INFO[,"s_fc"]<-            soils$s_fc/100                      #soil moisture at field capacity 
land.INFO[,"psi"]<-             -16662*(soils$n^7.8831)             #air entry pressure head (mm) --Relationship developed from Clapp and Hornberger, 1978
land.INFO[,"y_cl"]<-            -1*soils$y_cl                       #confining layer depth (mm) from SSURGO
land.INFO[,"y_c"]<-             -soils$y_rd/2                       #critical depth (mm) 
land.INFO[,"s_wilt"]<-          soils$s_w/100                       #soil moisture at permanent wilting point
land.INFO[,"k_sat"]<-           -soils$ksat*24                      #saturated condcuctivity (mm/day)
land.INFO[,"RD"]<-              -soils$y_rd                         #Rooting Depth (mm)
land.INFO[,"b"]<-               12.524*(soils$clay/100)+3.6907      #Presssure Head Power-Law Coefficient (b) --Relationship developed from Clapp and Hornberger, 1978
land.INFO[,"y_wt_0"]<-          0
land.INFO[,"s_t_0"]<-           land.INFO[,"s_fc"]
land.INFO[,"GW_bf_0"]<-         300
land.INFO[,"Sy"]<-              soils$Sy
land.INFO[,"wetland_invert"]<-  -1000*0.05*(length(rowSums(area)[rowSums(area)!=max(rowSums(area), na.rm=T)])+1)
land.INFO[,"wetland_area"]<-    max(rowSums(area))*(1000^2) 
land.INFO[,"volume_max"]<-      max(rowSums(volume[,-giw.ID]))*(1000^3)
land.INFO[,"kb"]<-              0.046                               #Defined fromo literatture. (Going to odo with Gauge analysis eventually)

#create wetland input variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define number of wetlands
n.wetlands<-1

#create input variables
giw.INFO<-c("giw.ID","area_watershed","area_wetland","invert","vol_ratio", #geometric characteristics
            "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", "Sy", #soil characteristics
            "y_w_0", "s_t_0" #initial conditions
)

#Create giw.INFO matrix
giw.INFO<-matrix(0, nrow=n.wetlands, ncol=length(giw.INFO), dimnames = list(seq(1,n.wetlands,1), c(giw.INFO)))

#Populate giw.INFO matrix (length units in mm)
giw.INFO[,"giw.ID"]<-          giw.ID
giw.INFO[,"area_watershed"]<-  (gArea(basin.shp, byid = T)[giw.ID])*(1000^2)
giw.INFO[,"area_wetland"]<-    max(area[,giw.ID])*(1000^2)
giw.INFO[,"invert"]<-         -1000*0.05*(length(area[,giw.ID][area[,giw.ID]!=max(area[,giw.ID], na.rm=T)]))   
giw.INFO[,"n"]<-              soils$n
giw.INFO[,"s_fc"]<-           soils$s_fc/100  
giw.INFO[,"psi"]<-            -16662*(soils$n^7.8831) 
giw.INFO[,"y_cl"]<-           -1*soils$y_cl   
giw.INFO[,"y_c"]<-           -soils$y_rd/2                       #critical depth (mm) 
giw.INFO[,"s_wilt"]<-        soils$s_w/100                       #soil moisture at permanent wilting point
giw.INFO[,"k_sat"]<-         -soils$ksat*24                      #saturated condcuctivity (mm/day)
giw.INFO[,"RD"]<-            -soils$y_rd                         #Rooting Depth (mm)
giw.INFO[,"b"]<-              12.524*(soils$clay/100)+3.6907 
giw.INFO[,"Sy"]<-             giw.INFO[,"n"]*(1-giw.INFO[,"s_fc"])
giw.INFO[,"y_w_0"]<-          0
giw.INFO[,"s_t_0"]<-          giw.INFO[,"s_fc"]
giw.INFO[,"vol_ratio"]<-      0 #ratio of upstream wetland volume

#save image
save.image(paste0(dir3,"/Backup/model_parameters.RData"))

####################################################################################
# Step 6: Run WHC-------------------------------------------------------------------
####################################################################################
#Setup workspace~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
rm(list=ls(all=TRUE))

#Define working and data directories
wd5<-"~/Wetland_Hydrologic_Capacitance_Model"
dir5<-"//nfs/WHC-data/Validation_Modeling/WHC_BaltimoreCorner"

#Load WHC functiono
source(paste0(wd5,"/R/WHC_2.R"))

#Load previous workspace
load((paste0(dir5,"/Backup/model_parameters.RData")))

#Define time period for simulations
n.years<-length(pet.VAR)/365

#Calibrations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Remove coastal precipitation events
date<-strptime(climate$DATE,"%Y%m%d")
precip.VAR[as.yearmon(date)=="Sep 2008"]<-0
precip.VAR[as.yearmon(date)=="Oct 2008"]<-0
precip.VAR[date==paste(as.Date("2010-08-18"))]<-0
precip.VAR[as.yearmon(date)=="Oct 2010"]<-0

#Go Time!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Run Model
WHC<-wetland.hydrology(giw.INFO, 
                       land.INFO, 
                       precip.VAR, 
                       pet.VAR, 
                       n.years, 
                       area, 
                       volume,
                       giw.ID)

#Initial Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plotting Parameters
par(mar=c(3,3,0,0)+0.25)
par(mgp=c(2,0.6,0)) 
par(ps=12)
par(cex.lab=14/12)
par(cex.axis=10/12)

#plot blank plot
plot(watershed$y_wt, 
     #Plot Type
     type="n", 
     #Axes labels
     xlab="Simulation Day", ylab="Relative Elevation [mm]")

#Plot Upland Water Level
points(watershed$y_wt, type="l", lty=2, lwd=2, col="darkorange")

#Plot lumped wetland water level
points(watershed$y_w, type="l", lty=2, lwd=2, col="navyblue")

#plot legend
legend("bottomleft",
       lty=2,, lwd=2, col=c("navyblue","darkorange"),
       c("Upland Wetland Stage", "Upland Water Table"),
       bty="n")


