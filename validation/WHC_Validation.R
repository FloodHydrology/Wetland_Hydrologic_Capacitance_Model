# ==================================================================================
# Name: HUC12 Analysis
# Coder: C. Nathan Jones & Fred Cheng
# Date: 10 Jan 2019
# Purpose: Main function to prepare inputs and run Model for Individual Wetland
# ==================================================================================
rm(list=ls(all=TRUE))
# 1. Set up workspace --------------------------------------------------------------
library(MASS)       #distribrutional analysis
library(markovchain)# for MCMC modeling
library(Evapotranspiration)
library(lubridate)  # dealing with dates
library(tidyverse)  # for data wrangling

setwd("~/Wetland_Hydrologic_Capacitance_Model")
# setwd('D:/OneDrive - University of Waterloo/Research/WHC Model/Wetland_Hydrologic_Capacitance_Model/')
source('R/WHC_Individual.R')

# 2. Bring in external data --------------------------------------------------------
data <- read.csv('validation/precip.csv') 
precip.VAR <- as.matrix(data$precip_mm)

data <- read.csv('validation/PET.csv') 
pet.VAR <- as.matrix(data$PET)

remove(data)

# 3. User defined model variables ---------------------------------------------------
n.years <- 5
giw.ID  <- 1

# 4. User defined site characteristics =============================================

giw.INFO<-c("giw.ID","WetID","area_watershed","area_wetland","invert","vol_ratio", "dL", "dLe", "dz",  # geometric characteristics
            "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", "Sy",                        # soil characteristics
            "y_w_0", "s_t_0"                                                                           # initial conditions
            )

# 4.1 Create giw.INFO matrix -----------------------------------------------------
giw.INFO<-matrix(0, nrow=1, ncol=length(giw.INFO),
                 dimnames = list(1, c(giw.INFO)))

giw.INFO[1,"giw.ID"]<-         1  #ID used in the model, note this is differnt than the WetID
giw.INFO[1,"WetID"]<-          1
giw.INFO[1,"area_watershed"]<- (5800*80) * 1000^2         # in mm^2
giw.INFO[1,"area_wetland"]<-   (5800) * 1000^2            # in mm^2
giw.INFO[1,"invert"]<-         -0.2 * 1000                # in mm^2
giw.INFO[1,"n"]<-              0.52
giw.INFO[1,"s_fc"]<-           0.36
giw.INFO[1,"psi"]<-            -16662*(giw.INFO[1,"n"]^7.8831)
giw.INFO[1,"y_cl"]<-           -0.33 * 1000
giw.INFO[1,"y_c"]<-            -0.75 * 1000                   # critical depth (mm)
giw.INFO[1,"s_wilt"]<-         0.17                           # soil moisture at permanent wilting point
giw.INFO[1,"k_sat"]<-          -0.27*1000                      # saturated condcuctivity (mm/day)
giw.INFO[1,"RD"]<-             0.10*1000                      # for wetland??

clay <-                        15                             # guesstimate
giw.INFO[1,"b"]<-              12.524*(clay/100)+3.6907
giw.INFO[1,"Sy"]<-             giw.INFO[1,"n"]*(1-giw.INFO[1,"s_fc"])
giw.INFO[1,"y_w_0"]<-          0
giw.INFO[1,"s_t_0"]<-          giw.INFO[1,"s_fc"]
giw.INFO[1,"vol_ratio"]<-      1                              # ratio of upstream wetland volume FYC:????
giw.INFO[1, "dL"] <-           300 *1000
giw.INFO[1, "dLe"] <-          300 *1000
giw.INFO[1,"dz"]<-             1                              # FYC:?????

# 4.2 Create Land Info Table ---------------------------------------------------
land.INFO<-c("area","invert",                                              # geometric characteristics
             "n","s_fc","psi","y_cl", "y_c", "s_wilt", "k_sat", "RD", "b", # soil characteristics
             "slope", "kb",                                                # larger watershed charactersitics
             "y_wt_0", "s_t_0","GW_bf_0",                                  # initial conditions
             "Sy",                                                         # calculated terms
             "wetland_invert","wetland_area","volume_max"                  # lumped wetland information
)

land.INFO<-matrix(0, nrow=1, ncol=length(land.INFO), dimnames = list(c(1), c(land.INFO)))

# 4.5 Populate land.INFO matrix (lenghth units in mm)
land.INFO[,"area"]<-            (5800*80) * 1000^2                          # area in mm^2
land.INFO[,"n"]<-               0.52                                        # porisity
land.INFO[,"s_fc"]<-            0.36                                        # soil moisture at field capacity
land.INFO[,"psi"]<-             -16662*(land.INFO[,"n"]^7.8831)             # air entry pressure head (mm) --Relationship developed from Clapp and Hornberger, 1978
land.INFO[,"y_cl"]<-            -0.33 * 1000                                # confining layer depth (mm) from SSURGO
land.INFO[,"y_c"]<-             -0.75 * 1000                                # critical depth (mm)
land.INFO[,"s_wilt"]<-          0.17                                        # soil moisture at permanent wilting point
land.INFO[,"k_sat"]<-           -0.27*1000                                  # saturated condcuctivity (mm/day)
land.INFO[,"RD"]<-              0.4*1000                                    # Rooting Depth (mm)
land.INFO[,"b"]<-               12.524*(clay/100)+3.6907      # Presssure Head Power-Law Coefficient (b) --Relationship developed from Clapp and Hornberger, 1978
land.INFO[,"y_wt_0"]<-          0
land.INFO[,"s_t_0"]<-           land.INFO[,"s_fc"]
land.INFO[,"GW_bf_0"]<-         0
land.INFO[,"Sy"]<-              land.INFO[,"n"]*(1-land.INFO[,"s_fc"])
land.INFO[,"wetland_invert"]<-  -1000*0.05
land.INFO[,"wetland_area"]<-    0.1
land.INFO[,"volume_max"]<-      0.1
land.INFO[,"kb"]<-              0.046       

# 4.3 Populate lumped.INFO table -------------------------------------------------------------------
# 5.1 #Create lumped.INFO table
lumped.INFO<-c("r_w","dL", "dLe") #geometric characteristics

# 5.2 Create lumped.INFO matrix
lumped.INFO<-matrix(0, nrow=1, ncol=3, dimnames = list(c(1:1), c(lumped.INFO)))

# 5.3 Populate lumped.INFO matrix (length in mm); data for the wetlands in the lumped upland
lumped.INFO[, "dL"] <- 1
lumped.INFO[,"r_w"] <- 1
lumped.INFO[,"dLe"] <- 1
lumped.INFO <- lumped.INFO *1000                                  # convert entire matrix from m to mm




# 5. Placeholder area, volume relationships -----------------------------------------------------
# 2.1 Create vectors to store area and volume to stage relationships
area<-matrix(0, ncol=(2), nrow=100)
volume<-matrix(0, ncol=(2), nrow=100)

# 2.2 Use loop to calculate based on previously published relationships
for(i in 1:2){
  #Define TempID
  TempID<-1
  
  #Define Amax
  Amax<-     giw.INFO[1,"area_wetland"]           # in mm^2
  
  #Define rmax
  rmax<-(Amax/pi)^.5
  
  #Define Vmax (Using Wu and Lane 2016)
  Vmax<-(0.25*((Amax/10000)^1.4742))*10000 #m3
  
  #Estimate hmax using Hayashi and Kamp 2000
  hmax<- 5
  
  #create functions to calculate area and volume
  area.fun<-function(h){pi*((rmax*((h/hmax)^.25))^2)}
  volume.fun<-function(h){0.25*((area.fun(h)/10000)^1.4742)*10000}
  
  #Apply functions
  n.col<-2
  n.rows<-length(seq(0,hmax-0.05,0.05))
  area[1:n.rows,]<- Amax
  
  volume[1:n.rows,]<-(seq(0,hmax-0.05,0.05)) * Amax
  
}




# 6. Run WHC -------------------------------------------------------------------------------------
results <- wetland.hydrology(giw.INFO,
                  land.INFO,
                  lumped.INFO,
                  precip.VAR,
                  pet.VAR,
                  n.years,
                  area,
                  volume,
                  giw.ID)








