#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Name: PPR Sandbox
# Coder: C. Nathan Jones & Fred Cheng
# Date: 11/12/2019  
# Purpose: Explore Individual Fluxes from PPR Wetlands for Troubleshooting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup Workspace------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1a. Clear Memory  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls(all=TRUE))

# 1b. Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Spatial Packages
library(sf)         # for spatial analysis
library(raster)     # for spatial analysis
library(rgdal)      # for spatial analysis
library(rgeos)      # for spatial analysis
library(geosphere)  # for spatial analysis
#Computational Packages
library(MASS)       #distribrutional analysis
library(markovchain)# for MCMC modeling
library(Evapotranspiration)
#Data munging packages
library(lubridate)  # dealing with dates
library(dplyr)      # for data wrangling
library(tidyr)
library(tibble)
#Parralell Processing Packagaes
library(parallel)   # parallel computing
library(rslurm)     # parallel computing

#define functions from file 
source("R/WHC_4.R")
source("R/regional_analysis.R")

# 1c. Define data directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
backup_dir<-"/nfs/WHC-data/Regional_Analysis/backup/"
results_dir<-"/nfs/WHC-data/Figure Generation/PPR Sandbox/"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2:Create Test Wetland----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Gather spatial inputs for HUC----------------------------------------------
# source("analysis/HUC12_Model_Input_PPR.R")
# save.image(paste0(backup_dir,"PPR_Input.RData"))

#2.2 Condence to one wetland----------------------------------------------------
#Load spatial data
load(paste0(backup_dir,"PPR_Input.RData"))

#Determine median wetland and catchment charactersitics
target_size<-wetlands.shp %>% st_as_sf() %>% st_drop_geometry() %>% select(area_m2) %>% summarise(median(area_m2,na.rm=T)) %>% pull()
WetID <-wetlands.shp$WetID[wetlands.shp$area_m2==target_size]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3: Function to run WHC Model---------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create function to process data and run WHC for a wetland of interest~~~~~~~~~
model_function <-  function(WetID,
                            n.years,
                            pet.VAR,
                            precip.VAR,
                            snowmelt.VAR,
                            wetlands.shp,
                            HUC12.shp,
                            catchments.shp,
                            flowlines.shp,
                            fac.grd,
                            soils.shp,
                            dem.grd,
                            nfw_centroid.shp,
                            rootdepth.grd){ #
  
  # 1 tidyr::gather Required Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1.1 Identify wetland of interest
  main_wetland.shp    <- wetlands.shp[wetlands.shp$WetID==WetID,]
  main_wetland        <- WetID
  
  # 1.2 Select spatial layers based on WetID  
  catchment_temp.shp  <- catchments.shp[main_wetland.shp,]
  if(length(catchment_temp.shp)>1){catchment_temp.shp<-catchment_temp.shp[1,]}
  wetlands_temp.shp   <- wetlands.shp[catchment_temp.shp,]
  soils_temp.shp      <- soils.shp[catchment_temp.shp,]
  soils_temp.shp      <- raster::crop(soils_temp.shp,catchment_temp.shp)
  dem_temp.grd        <- raster::crop(dem.grd, catchment_temp.shp)
  dem_temp.grd        <- mask(dem_temp.grd, catchment_temp.shp)
  dem_temp.grd        <- dem_temp.grd/100 #Convert to m
  fac_temp.grd        <- raster::crop(fac.grd, catchment_temp.shp)
  fac_temp.grd        <- mask(fac_temp.grd, catchment_temp.shp)
  root_temp.grd       <- raster::crop(rootdepth.grd, catchment_temp.shp)
  root_temp.grd       <- mask(rootdepth.grd, catchment_temp.shp)
  
  # 1.3 For now, add catchment aggregate soils data to soils with missing data
  soils_temp.shp$y_c[is.na(soils_temp.shp$y_c)]   <- mean(soils_temp.shp$y_c, na.rm=T)
  soils_temp.shp$y_cl[is.na(soils_temp.shp$y_cl)] <- mean(soils_temp.shp$y_cl, na.rm=T)
  soils_temp.shp$s_fc[is.na(soils_temp.shp$s_fc)] <- mean(soils_temp.shp$s_fc, na.rm=T)
  soils_temp.shp$s_w[is.na(soils_temp.shp$s_w)]   <- mean(soils_temp.shp$s_w, na.rm=T)
  soils_temp.shp$n[is.na(soils_temp.shp$n)]       <- mean(as.numeric(paste(soils_temp.shp$n)), na.rm=T)
  soils_temp.shp$clay[is.na(soils_temp.shp$clay)] <- mean(soils_temp.shp$clay, na.rm=T)
  soils_temp.shp$ksat[is.na(soils_temp.shp$ksat)] <- mean(soils_temp.shp$ksat, na.rm=T)
  
  # 2 Estimate area and volume to stage relationships~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.1 Create vectors to store area and volume to stage relationships
  area     <- matrix(0, ncol=length(wetlands_temp.shp), nrow=100)
  volume   <- matrix(0, ncol=length(wetlands_temp.shp), nrow=100)
  
  # 2.2 Use loop to calculate based on previously published relationships
  for(i in 1:length(wetlands_temp.shp)){
    #Define shape factor (for now hard code, maybe make it a var later)
    p<-2
    
    #Define TempID
    TempID<-wetlands_temp.shp$WetID[i]
    
    #Define Amax
    Amax<-wetlands_temp.shp$area_m2[wetlands_temp.shp$WetID==TempID]
    
    #Define Vmax (Using Wu and Lane 2016)
    Vmax<-(0.01725*(Amax^1.30086)) #m3
    
    #Estimate hmax using Hayashi and Kamp 2000  
    hmax<-Vmax*(1+(2/p))/Amax
    
    #create functions to calculate area and volume
    area.fun<-function(h){Amax*(h/hmax)^(2/p)}
    volume.fun<-function(h){Amax*(h^(1+2/p))/(1+2/p)/(hmax^(2/p))}
    
    #Apply functions
    n.col   <- which(wetlands_temp.shp$WetID==TempID)
    n.rows  <- length(area.fun(seq(0,hmax,0.05)))
    area[1:n.rows,n.col]         <- area.fun(seq(0,hmax,0.05))
    area[(n.rows+1):100,n.col]   <- area.fun(hmax)
    volume[1:n.rows,n.col]       <- volume.fun(seq(0,hmax,0.05))
    volume[(n.rows+1):100,n.col] <- volume.fun(hmax)
  }
  
  # 3 Populate GIW.INFO Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 3.1 Create input variables
  giw.INFO<-c("giw.ID","WetID","area_watershed","area_wetland","invert","vol_ratio", "dL", "dLe", "dz",  # geometric characteristics
              "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", "Sy",                        # soil characteristics
              "y_w_0", "s_t_0"                                                                           # initial conditions
  )
  
  # 3.2 Create giw.INFO matrix
  giw.INFO<-matrix(0, nrow=length(wetlands_temp.shp), ncol=length(giw.INFO),
                   dimnames = list(seq(1,length(wetlands_temp.shp)), c(giw.INFO)))
  
  # 3.3 Define wetland variables (units == mm)
  for(i in 1:length(wetlands_temp.shp)){
    
    # a. Isolate soils data (if for osme reason we don't have ovlerlap, just use basin agregate)
    if(gIntersects(soils_temp.shp,wetlands_temp.shp[i,])==TRUE){
      temp_soils<-raster::crop(soils_temp.shp,wetlands_temp.shp[i,])
    }else{
      temp_soils<-soils_temp.shp
    }
    
    # b. Agregrate parameters based on space
    temp_soils$area<-gArea(temp_soils, byid=T)
    temp_soils<-temp_soils@data
    temp_soils<-colSums(temp_soils[,c("y_c","y_cl","s_fc","s_w","n","clay","ksat")]*temp_soils[,"area"])/sum(temp_soils$area, na.rm=T)
    
    # c. Isolate fac data and calculate ratio of upland drainage area to total drainage area
    temp_fac<-raster::crop(fac_temp.grd, wetlands_temp.shp[i,])
    temp_fac<-mask(temp_fac, wetlands_temp.shp[i,])
    temp_fac<-ifelse(cellStats(temp_fac*0+1, sum)>0,
                     cellStats(temp_fac,max)/cellStats(fac_temp.grd,max),
                     0)
    
    # d. Calculate relative wetland datum
    # dz<-cellStats(mask(dem_temp.grd, wetlands_temp.shp[i,]), median, na.rm=T)- #median elevation of wetland elevation
    #   cellStats(mask(dem_temp.grd, wetlands_temp.shp[wetlands_temp.shp$WetID!=WetID,]), median, na.rm=T) #median elevation of all wetlands
    # dz<-ifelse(is.na(dz)==TRUE, 0, dz)
    # dz<-ifelse(dz==Inf, 0, dz)
    # dz<-dz*1000
    
    # e. Calcluate wetland rooting depth
    temp_rd<-raster::crop(root_temp.grd, wetlands_temp.shp[i,])
    temp_rd<-mask(temp_rd, wetlands_temp.shp[i,])
    temp_rd<-ifelse(cellStats(temp_rd*0+1, sum)>0,
                    cellStats(temp_rd,max)/cellStats(root_temp.grd,max),
                    0.1)
    
    # f. Populate giw.INFO table
    giw.INFO[i,"giw.ID"]<-         i  #ID used in the model, note this is differnt than the WetID
    giw.INFO[i,"WetID"]<-          wetlands_temp.shp$WetID[i]
    giw.INFO[i,"area_watershed"]<- gArea(catchment_temp.shp)*(10^6)
    giw.INFO[i,"area_wetland"]<-   gArea(wetlands_temp.shp[i,], byid=T)*(10^6)
    giw.INFO[i,"invert"]<-         -1000*0.05*(length(area[,i][area[,i]!=max(area[,i], na.rm=T)]))
    giw.INFO[i,"n"]<-              temp_soils["n"]
    giw.INFO[i,"s_fc"]<-           temp_soils["s_fc"]/100
    giw.INFO[i,"psi"]<-            -16662*(temp_soils["n"]^7.8831)
    giw.INFO[i,"y_cl"]<-           -1*temp_soils["y_cl"]
    giw.INFO[i,"y_c"]<-            - temp_soils["y_c"]                       #critical depth (mm)
    giw.INFO[i,"s_wilt"]<-         temp_soils["s_w"]/100                       # soil moisture at permanent wilting point
    giw.INFO[i,"k_sat"]<-          temp_soils["ksat"]                      # saturated condcuctivity (mm/day)
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
    giw.INFO[i,"dz"]<-             wetlands_temp.shp$dz[i]
  }
  
  # 4 Populate land.INFO table~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 4.1 Create Land Info Table
  land.INFO<-c("area","invert",                                              # geometric characteristics
               "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", # soil characteristics
               "slope", "kb",                                                # larger watershed charactersitics
               "y_wt_0", "s_t_0","GW_bf_0",                                  # initial conditions
               "Sy",                                                         # calculated terms
               "wetland_invert","wetland_area","volume_max"                  # lumped wetland information
  )
  
  # 4.2 Create land.INFO matrix
  land.INFO<-matrix(0, nrow=1, ncol=length(land.INFO), dimnames = list(c(1), c(land.INFO)))
  
  # 4.3 Aggregate soils data by area
  temp_soils<-soils_temp.shp
  temp_soils$area<-gArea(temp_soils, byid=T)
  temp_soils<-temp_soils@data
  temp_soils<-colSums(temp_soils[,c("y_c","y_cl","s_fc","s_w","n","clay","ksat")]*temp_soils[,"area"])/sum(temp_soils$area, na.rm=T)
  
  # 4.4 Rooting depth calculation
  RD<- cellStats(root_temp.grd, mean, na.rm=T)
  RD<- ifelse(RD<temp_soils["y_cl"], -1*RD, -1*temp_soils["y_cl"])
  
  # 4.5 Populate land.INFO matrix (lenghth units in mm)
  land.INFO[,"area"]<-            gArea(catchment_temp.shp)*(10^6)            # area in mm^2
  land.INFO[,"n"]<-               temp_soils["n"]                             # porisity
  land.INFO[,"s_fc"]<-            temp_soils["s_fc"]/100                      # soil moisture at field capacity
  land.INFO[,"psi"]<-             -16662*(temp_soils["n"]^7.8831)             # air entry pressure head (mm) --Relationship developed from Clapp and Hornberger, 1978
  land.INFO[,"y_cl"]<-            -1*temp_soils["y_cl"]                       # confining layer depth (mm) from SSURGO
  land.INFO[,"y_c"]<-             -temp_soils["y_c"]                          # critical depth (mm)
  land.INFO[,"s_wilt"]<-          temp_soils["s_w"]/100                       # soil moisture at permanent wilting point
  land.INFO[,"k_sat"]<-           temp_soils["ksat"]                          # saturated condcuctivity (mm/day)
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
  
  #4.5 Define clay layer for both land.INFO and giw.INFO
  giw.INFO[,"y_cl"]<-land.INFO[,"y_cl"]
  
  # 5 Populate lumped.INFO table~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 5.1 #Create lumped.INFO table
  lumped.INFO<-c("r_w","dL", "dLe") #geometric characteristics
  
  # 5.2 Create lumped.INFO matrix
  lumped.INFO<-matrix(0, nrow=length(wetlands_temp.shp$WetID), ncol=3, dimnames = list(c(1:length(wetlands_temp.shp$WetID)), c(lumped.INFO)))
  
  # 5.3 Populate lumped.INFO matrix (length in mm); data for the wetlands in the lumped upland
  lumped.INFO[, "dL"] <- wetlands_temp.shp$dist2NearWet                  #
  lumped.INFO[,"r_w"] <- (((wetlands_temp.shp$SHAPE_Area) / pi) ^ 0.5 )  # derive radius from area assuming circular geometry, convert
  lumped.INFO[,"dLe"] <- wetlands_temp.shp$dLe                           #
  lumped.INFO <- lumped.INFO *1000                                  # convert entire matrix from m to mm
  
  # 6 Run the model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 6.1 Define the wetland
  giw.ID<-giw.INFO[,"giw.ID"][giw.INFO[,"WetID"]==main_wetland]
  
  # 6.2 Create function to execute WHC and organize output terms
  execute<-function(n.years){
    
    # a. Run WHC Model w/ tryCatch
    HUC12.shp$HUC_12<-HUC12.shp$HUC12
    output<-tryCatch(wetland.hydrology(giw.INFO,
                                       land.INFO,
                                       lumped.INFO,
                                       snowmelt.VAR,
                                       precip.VAR,
                                       pet.VAR,
                                       n.years,
                                       area,
                                       volume,
                                       giw.ID),
                     error = function(e) tibble::tibble(var = "error",
                                                        value = -9999,
                                                        day = 0,
                                                        HUC12=HUC12.shp$HUC_12,
                                                        WetID = WetID,
                                                        scale = "wetland"))
    
    #Gather info (brute force, cause I just gotta get shit done...)
    output<-tibble(
      giw_ET_lm = output$ET_lm.VAR[,3],
      giw_ET_wt = output$ET_wt.VAR[,3],
      giw_exfil = output$exfil.VAR[,3],
      giw_GW_local = output$GW_local.VAR[,3],
      giw_loss_lm = output$loss_lm.VAR[,3],
      giw_s_ex = output$s_ex.VAR[,3], 
      giw_s_lim = output$s_lim.VAR[,3],
      giw_runoff= output$runoff_vol.VAR[,3],
      giw_spill = output$spill_vol.VAR[,3],
      giw_y_hm = output$y_hm.VAR[,3],
      giw_y_sat = output$y_sat.VAR[,3],
      giw_y_wt = output$y_w.VAR[,3],
      giw_y_w = output$y_w.VAR[,3],
      upland_ET_lm = output$ET_lm.VAR[,1],
      upland_ET_wt = output$ET_wt.VAR[,1],
      upland_exfil = output$exfil.VAR[,1],
      upland_GW_local = output$GW_local.VAR[,1],
      upland_GW_baseflow = output$GW_bf.VAR[,1],
      upland_loss_lm = output$loss_lm.VAR[,1],
      upland_s_ex = output$s_ex.VAR[,1], 
      upland_s_lim = output$s_lim.VAR[,1],
      upland_runoff= output$runoff_vol.VAR[,1],
      upland_spill = output$spill_vol.VAR[,2],
      upland_y_hm = output$y_hm.VAR[,1],
      upland_y_sat = output$y_sat.VAR[,1],
      upland_y_wt = output$y_wt.VAR[,1],
      upland_y_w = output$y_w.VAR[,2]
    )  
    
    #Export
    output
  }
  
  # 6.3 Run function
  WHC<-execute(n.years)
  
  # 6.6 Export WHC Results
  WHC
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 4: Run Function----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define remaining input vars
n.years<-3

#Run function
df<-model_function(WetID,
                   n.years,
                   pet.VAR,
                   precip.VAR,
                   snowmelt.VAR,
                   wetlands.shp,
                   HUC12.shp,
                   catchments.shp,
                   flowlines.shp,
                   fac.grd,
                   soils.shp,
                   dem.grd,
                   nfw_centroid.shp,
                   rootdepth.grd)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 5: Create plots----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#5.1 Wetland--------------------------------------------------------------------
#a. Wetland Fluxes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Identify wetland area
wet_area<-wetlands.shp %>% st_as_sf() %>% st_drop_geometry %>% rename(ID=WetID) %>% filter(ID==WetID) %>% select(area_m2) %>% pull

#Plot Wetland Fluxes
df %>% 
  #Select giw fluxes
  select(contains('giw')) %>% 
  #Remove level and storage data
  select(-contains('_y_')) %>% 
  select(-contains('_s_')) %>% 
  #Divide by wetland area where needed
  mutate(giw_GW_local=giw_GW_local/wet_area,
        giw_runoff = giw_runoff/wet_area) %>% 
  #Turn into long format (for plotting) %>% 
  mutate(day = rownames(.), 
         day = as.numeric(day)) %>% 
  pivot_longer(-day) %>% 
  #Normalize by wetland area
  mutate(value = value) %>% 
  #Start ggplot
  ggplot() +
    #add line  
    geom_line(aes(x=day, y=value)) +
    #Add Facets
    facet_grid(name~., scale='free') +
    #Add verticle lines
    geom_vline(aes(xintercept = c(365)), lty=2, col="grey30") +
    geom_vline(aes(xintercept = c(730)), lty=2, col="grey30") +
    #add black/white theme
    theme_bw()+
      ggtitle('Wetland Fluxes')+
      ylab("Flow Rate [mm/day]")

#Print plot
ggsave(paste0(results_dir,"wetland_fluxes.png"), width = 7, height = 9)

#b. Wetland Storages~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plot Wetland Storages
df %>% 
  #Select giw fluxes
  select(contains('giw')) %>% 
  #Remove level and storage data
  select(matches('_y_|_s_')) %>% 
  #Turn into long format (for plotting) %>% 
  mutate(day = rownames(.), 
         day = as.numeric(day)) %>% 
  pivot_longer(-day) %>% 
  #Normalize by wetland area
  mutate(value = value) %>% 
  #Remove final day
  filter(day<1095) %>% 
  #Start ggplot
  ggplot() +
  #add line  
  geom_line(aes(x=day, y=value)) +
  #Add verticle lines
  geom_vline(aes(xintercept = c(365)), lty=2, col="grey30") +
  geom_vline(aes(xintercept = c(730)), lty=2, col="grey30") +
  #Add Facets
  facet_grid(name~., scale='free') +
  #add black/white theme
  theme_bw()+
  ggtitle('Wetland Storages')+
  ylab("Storage [mm -OR - % Sat]")

#Print plot
ggsave(paste0(results_dir,"wetland_storages.png"), width = 7, height = 9)

#5.2 Upland---------------------------------------------------------------------
#a. Upland Fluxes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Identify watershed area
wet_area<-catchments.shp[wetlands.shp[wetlands.shp$WetID==WetID,],] %>% st_as_sf %>% st_drop_geometry() %>% mutate(area_m2 = 1000000*AREASQKM) %>% select(area_m2) %>% pull

#Plot Watershed Fluxes
df %>% 
  #Select giw fluxes
  select(contains('upland')) %>% 
  #Remove level and storage data
  select(-contains('_y_')) %>% 
  select(-contains('_s_')) %>% 
  #Convert to mm where needed
  mutate(upland_GW_local=upland_GW_local/wet_area,
         upland_runoff = upland_runoff/wet_area) %>% 
  #Turn into long format (for plotting) %>% 
  mutate(day = rownames(.), 
         day = as.numeric(day)) %>% 
  pivot_longer(-day) %>% 
  #Normalize by wetland area
  mutate(value = value) %>% 
  #Start ggplot
  ggplot() +
  #add line  
  geom_line(aes(x=day, y=value)) +
  #Add verticle lines
  geom_vline(aes(xintercept = c(365)), lty=2, col="grey30") +
  geom_vline(aes(xintercept = c(730)), lty=2, col="grey30") +
  #Add Facets
  facet_grid(name~., scale='free') +
  #add black/white theme
  theme_bw()+
  ggtitle('Watershed Fluxes')+
  ylab("Flow Rate [mm/day]")

#Print plot
ggsave(paste0(results_dir,"watershed_fluxes.png"), width = 7, height = 9)

#b. Watershed Storages~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plot Watershed Storages
df %>% 
  #Select giw fluxes
  select(contains('upland')) %>% 
  #Remove level and storage data
  select(matches('_y_|_s_')) %>% 
  #Turn into long format (for plotting) %>% 
  mutate(day = rownames(.), 
         day = as.numeric(day)) %>% 
  pivot_longer(-day) %>% 
  #Normalize by wetland area
  mutate(value = value) %>% 
  #Remove final day
  filter(day<1095) %>% 
  #Start ggplot
  ggplot() +
  #add line  
  geom_line(aes(x=day, y=value)) +
  #Add verticle lines
  geom_vline(aes(xintercept = c(365)), lty=2, col="grey30") +
  geom_vline(aes(xintercept = c(730)), lty=2, col="grey30") +
  #Add Facets
  facet_grid(name~., scale='free') +
  #add black/white theme
  theme_bw()+
  ggtitle('Watershed Storages')+
  ylab("Storage [mm -OR - % Sat]")

#Print plot
ggsave(paste0(results_dir,"watershed_storages.png"), width = 7, height = 9)
