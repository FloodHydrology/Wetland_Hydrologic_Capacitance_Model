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

#Define Master Working Directory
dir<-"//nfs/WHC-data/Validation_Modeling/WHC_BaltimoreCorner"

#Load Required Packages
library(raster)   #spatial analysis
library(rgdal)    #spatial analysis
library(rgeos)    #spatial analysis
library(maptools) #spatial analysis
library(plyr)     #data processing
library(dplyr)    #data processing
library(Evapotranspiration)

#Obtain Model Inputs
dem<-raster(readGDAL(paste0(dir,"/Model Inputs/DEM/dem")))      #DEM (from USDA-Amir Sharifi)
soils.shp<-readOGR(dsn = paste0(dir,'/Model Inputs/Soils'), layer='SSURGO')  #SSURGO data
climate<-read.csv(paste0(dir,"/Model Inputs/Climate/ncdc.csv"))
giw.ID<-63 #Need to update with coordinates of the wetland in question.

#Add pp for watershed  (For now define this manually)
pp<-data.frame(x=500674.396458, y=155144.28711)
coordinates(pp)<-~x+y
projection(pp)<-dem@crs
pp.shp<-SpatialPoints(pp)
pp.shp<-SpatialPointsDataFrame(pp, data.frame(x=500674.396458, y=155144.28711))

#Convert all data sources into correct format
proj<-CRS("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") 
dem<-projectRaster(dem, crs=proj)
soils.shp<-spTransform(soils.shp, proj)
pp.shp<-spTransform(pp.shp, proj)

#Plot to make sure there is overlap
plot(dem)
plot(soils.shp, border="grey60", cex=0.25, add=T)
plot(pp.shp, pch=19, col="red", cex=3, add=T)

#Save to backup folder
save.image(paste0(dir, "/Backup/input_dat.RData"))
  
####################################################################################
# Step 2: Topographic Analysis-------------------------------------------------------
####################################################################################
#The topographic analysis requires connection to both ArcGIS and SAGA GIS.  This software
#are not currently available on the SESYNC server, so this step is run using the 
#Topo_Analysis_Model_Validation.R file. A few inmportant points:
   #1) You will need to save Topo_Analysis_Model_Validation.R on YOUR CPU
   #2) Note the "dir" variable will need to change from "/nfs" to the defined dir on your computer
   #3) I would suggest cloaning the github repository to a workspace folder (eg "C:\\Workspace").  

#Call Function DEM Processing F(x) into Memory
source(paste0(wd, "/R_Functions/DEM_Processing.R"))



#Run Terrain Analysis Subroutine
DEM_Processing.fun(dem, pp, wd)

#Backup
save.image(paste0(wd,"/Backup/DEM_Processing.RData"))

####################################################################################
# Step 3: Soils Analysis-------------------------------------------------------------
####################################################################################
#Rest wd
setwd(wd)

#manipulate soils data
soils.shp<-spTransform(soils.shp, dem@crs) #transform the layer
soils.shp<-raster::intersect(soils.shp,watershed.shp) #clip to watershed
soils.shp@data$area_m2<-gArea(soils.shp, byid=T) #calculate area

#merge soils data with SSURGO database
soils<-soils.shp@data 
soils<-soils[,c("MUKEY","area_m2")]
soil.data<-read.csv(paste0(wd,"/R_Functions/WHC_Soils_Input.csv"))
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
climatedata$Tmax.daily<-(climate$TMAX)
climatedata$Tmin.daily<-(climate$TMIN)
climatedata$RHmax.daily<-0
climatedata$RHmin.daily<-0
climatedata$Rs.daily<-0
#create timeseries input file (second input file)
input<-ReadInputs(climatedata,
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
df<-ET.HargreavesSamani(input, constants, ts="daily")
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

#define wetland id (visually for now)
plot(basin.shp)
invisible(text(getSpPPolygonsLabptSlots(basin.shp), 
               labels=as.character(basin.shp@data$ID), 
               cex=0.9))
giw.ID<-63

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
giw.INFO[,"area_wetland"]<-    max(area[,63])*(1000^2)
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
save.image(paste0(wd,"/Backup/inputs.RData"))

####################################################################################
# Step 6: Run WHC-------------------------------------------------------------------
####################################################################################
#Load WHC functiono
source(paste0(wd,"/R_Functions/WHC_2.R"))

#Define time period for simulations
ramp<-3
n.years<-length(pet.VAR)/365*ramp

#Calibrations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Surface water volume to stage relationship
load(paste0(wd,"/Backup/BC_stage_volume.RData"))
var<-14
area[,giw.ID]<-df$area[var]
area[1:var,giw.ID]<-df$area[1:var]
volume[,giw.ID]<-df$volume[var]
volume[1:var,giw.ID]<-df$volume[1:var]
giw.INFO[,"invert"]<- -1000*0.05*(length(area[,giw.ID][area[,giw.ID]!=max(area[,giw.ID], na.rm=T)]))

#Remove coastal precipitation events
date<-strptime(climate$DATE,"%Y%m%d")
#Fall 2008
precip.VAR[as.yearmon(date)=="Sep 2008"]<-0
precip.VAR[as.yearmon(date)=="Oct 2008"]<-0
#Fall 2010
precip.VAR[date==paste(as.Date("2010-08-18"))]<-0
#precip.VAR[as.yearmon(date)=="Sep 2010"]<-0
precip.VAR[as.yearmon(date)=="Oct 2010"]<-0

#Go Time!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Run Model
WHC<-wetland.hydrology(giw.INFO, 
                       land.INFO, 
                       rep(precip.VAR,ramp), 
                       rep(pet.VAR,ramp), 
                       n.years, 
                       area, 
                       volume,
                       giw.ID)

####################################################################################
# Step 8: Analyze outputs-----------------------------------------------------------
####################################################################################
#download water level data
data<-read.csv("Validation Data/baltimore_corner.csv")
data<-data[data$y_w>-1.5,] #remove data days
data<-aggregate(data$y_w, list(paste(data$date)), FUN='mean')
colnames(data)<-c("date", "y_w")
data$date<-strptime(data$date, format="%m/%d/%Y")
data<-data[order(data$date),]

#creatae date for simulated data
date<-strptime(climate$DATE,"%Y%m%d")
giw<-giw[(length(date)+1):(length(date)*2),]
giw$date<-date
watershed<-watershed[length(date)+1:length(date)*2,]
watershed$date<-date

#combine datasets
df<-data.frame(paste(giw$date), as.numeric(paste(giw$y_w)))
colnames(df)<-c('date','mod')
df<-merge(df, data, by='date')
colnames(df)<-c('date','mod','obs')
df$obs<-df$obs*1000-700

#Initial Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jpeg("M:\\McLaughlin_Lab/Jones/workspace/GIW Modeling/USDA Water Balance/WHC_BaltimoreCorner/Model Output/validation.jpg", res=300, width=4, height=2.5, units="in")
par(mar=c(3.25,3.25,0.35,0.35))
par(mgp=c(2,0.5,0))
plot(date,giw$y_w-giw.INFO[,"invert"],type="n", 
     ylim=c(-500,750), 
     ps=12, cex.axis=10/12, cex.lab=14/12,
     ylab="Water Level [mm]",xlab="[Year]",
     panel.first=grid(nx=8, ny=6, lty=2, lwd=0.9, col="grey80")
)
points(data$date, (data$y_w*1000), type="l",col="grey40", lty=2, lwd=1.1)
points(giw$date, giw$y_w-giw.INFO[,"invert"], type="l",lwd=1.25, col="grey10")
legend("bottom", c("Observed","Modeled"), lty=c(2,1), lwd=2, col=c("grey60","grey30"), cex=10/12, box.col="white")#,bty="n", box.col="white")
box()
dev.off()

#Calculate RSQ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Combine Daily
daily<-summary(lm(df[,2]~df[,3]))$r.squared

#monthly
df<-data.frame(aggregate(df$mod, list(format(as.Date(df$date), format="%m-%Y")), mean),
               aggregate(df$obs, list(format(as.Date(df$date), format="%m-%Y")), mean)[,2])
colnames(df)<-c("date","mod","obs")
monthly<-summary(lm(df$mod~df$obs))$r.squared    

#annual
df<-data.frame(aggregate(df$mod, list(substr(df$date,4,8)), mean),
               aggregate(df$obs, list(substr(df$date,4,8)), mean)[,2])
colnames(df)<-c("date","mod","obs")
annual<-summary(lm(df$mod~df$obs))$r.squared  

#print
c(daily, monthly, annual)

####################################################################################
# Step 9: Export Observed and Modeled
####################################################################################
write.csv(df,"C:\\Users/cnjones/Google Drive/Research Projects/Delmarva_GIW/Presentations/201612_AGU/ValidationPlot/Delmarva.csv", row.names=F)
