# ==================================================================================
# Name:    HUC12 Analysis
# Coder:   C. Nathan Jones & Fred Cheng
# Date:    10 Jan 2019
# Purpose: Main function to prepare inputs and run Model for Individual Wetland
# NOTE:    THIS IS FOR A GENERIC SANDBOX MODEL TO TEST PARAMETERS AND STUFF; 
#          NOT FOR EXISTING PPR SITE THAT WE HAVE DATA FOR
# ==================================================================================
rm(list=ls(all=TRUE))

# 1. Set up workspace --------------------------------------------------------------
library(MASS)       # distribrutional analysis
library(markovchain)# for MCMC modeling
library(Evapotranspiration)
library(lubridate)  # dealing with dates
library(tidyverse)  # for data wrangling

WetID <- 26
n.years <- 5
giw.ID <- WetID

# 2. Load Necessary Data ----------------------------------------------------------
  # Run all of HUC12_Model_Input_PPR.R except dLe is assumed to be dL for now
  # wetlands.shp$dLe <- wetlands.shp$dist2NearWet

# 3. Generate .INFO Files ---------------------------------------------------------
  # Run all of regional_analysis.R using WetID <- 26
  # Run until lumped.INFO is run
  # setwd("/nfs/WHC-data/Sandbox Model")
  # save.image(file = 'PPR_ProcessedInputs_ID26_all.RData')
  # save(HUC12.shp, area, giw.INFO, land.INFO, lumped.INFO, pet.VAR, precip.VAR, volume,
  #      file = 'PPR_ProcessedInputs_ID26_selected.RData')

load(file = 'PPR_ProcessedInputs_ID26_selected.RData')

# 4. Source and Run WHC Model ----------------------------------------------------
source("~/Wetland_Hydrologic_Capacitance_Model/R/WHC_Sandbox.R")

output <- wetland.hydrology(giw.INFO,
                            land.INFO,
                            lumped.INFO,
                            precip.VAR,
                            pet.VAR,
                            n.years,
                            area,
                            volume,
                            giw.ID)

if(is.list(output)==T){
  giw.INFO<-giw.INFO[giw.INFO[,"giw.ID"]==giw.ID,]
  
  #i Wetland Scale Estimates
  #Calculate wetland scale annual water balance
  wetland_balance<-
    tibble::tibble(precip     = sum(precip.VAR)/(length(precip.VAR)/365),
                   pet        = sum(pet.VAR)/(length(pet.VAR)/365),
                   et         = (sum(output$ET_lm.VAR[,3])+sum(output$ET_wt.VAR[,3]))/n.years,
                   qf_in      = sum(output$runoff_vol.VAR[,3])/giw.INFO["area_wetland"]/n.years,
                   qf_out     = sum(output$spill_vol.VAR[output$runoff_vol.VAR[,3]!=0,3])/n.years/giw.INFO["area_wetland"],
                   sw_out     = sum(output$spill_vol.VAR[output$runoff_vol.VAR[,3]==0,3])/n.years/giw.INFO["area_wetland"],
                   gw_out     = sum(output$GW_local.VAR[output$GW_local.VAR<0])/giw.INFO["area_wetland"]/n.years,
                   gw_in      = sum(output$GW_local.VAR[,3][output$GW_local.VAR>0])/giw.INFO["area_wetland"]/n.years) %>%
    tidyr::gather(key="var") %>%
    dplyr::mutate(day=0)
  
  #Calculate mean annual duration of fluxes
  wetland_duration<-
    tibble::tibble(rain_day   = length(precip.VAR[precip.VAR!=0])/(length(precip.VAR)/365),
                   gw_in_day  = length(output$GW_local.VAR[output$GW_local.VAR[,3]>0,3])/n.years,
                   gw_out_day = length(output$GW_local.VAR[output$GW_local.VAR[,3]<0,3])/n.years,
                   qf_in_day  = length(output$runoff_vol.VAR[output$runoff_vol.VAR[,3]>0,3])/n.years,
                   qf_out_day = length(output$spill_vol.VAR[output$runoff_vol.VAR[,3]!=0 & 
                                                              output$spill_vol.VAR[,3]!=0,3])/n.years,
                   sw_out_day = length(output$spill_vol.VAR[output$runoff_vol.VAR[,3]==0 & 
                                                              output$spill_vol.VAR[,3]!=0,3])/n.years) %>%
    tidyr::gather(key="var") %>%
    dplyr::mutate(day=0)
  
  #Calculate mean daily fluxes at wetland scale
  wetland_fluxes<-
    tibble::tibble(
      day       = c(rep(seq(1,365),n.years),1), 
      y_w       = output$y_w.VAR[,3], 
      gw_local  = output$GW_local.VAR[,3]/giw.INFO["area_wetland"], 
      spill_out = output$spill_vol.VAR[,3]/giw.INFO["area_wetland"], 
      runoff_in = output$runoff_vol.VAR[,3]/giw.INFO["area_wetland"]) %>%
    dplyr::mutate(y_w = y_w + abs(giw.INFO["invert"]), 
                  y_w = dplyr::if_else(y_w>0,
                                       y_w/abs(giw.INFO["invert"]),
                                       y_w/abs(giw.INFO["y_cl"]))) %>%
    dplyr::group_by(day) %>%
    dplyr::summarise_all(., mean) %>%
    tidyr::gather(var,value,-day)
  
  #Create wetland output
  output_wetland<-dplyr::bind_rows(wetland_balance, wetland_duration, wetland_fluxes) %>%
    dplyr::mutate(HUC12 = HUC12.shp$HUC12,
                  WetID = WetID, 
                  scale = "weltand")
  
  #ii. Catchement scale estimates
  #Calcluate catchment scale annual water balance
  catchment_balance<-
    tibble::tibble(precip     = sum(precip.VAR)/(length(precip.VAR)/365),
                   pet        = sum(pet.VAR)/(length(pet.VAR)/365),
                   et         = (sum(output$ET_lm.VAR[,1])+sum(output$ET_wt.VAR[,1]))/n.years,
                   sw_out     = sum(output$spill_vol.VAR[,2])/land.INFO[,"area"]/n.years,
                   gw_out     = sum(-1*output$GW_bf.VAR[,1])/n.years) %>%
    tidyr::gather(key="var") %>%
    dplyr::mutate(day=0)
  
  #Calculate mean annual duration of catchment fluxes
  catchment_duration<-
    tibble::tibble(rain_day   = length(precip.VAR[precip.VAR!=0])/(length(precip.VAR)/365),
                   sw_out_day = length(output$spill_vol.VAR[output$spill_vol.VAR[,2]!=0,1])/n.years,
                   gw_out_day = length(output$GW_bf.VAR[output$GW_bf.VAR[,1]<0,1])/n.years) %>%
    tidyr::gather(key="var") %>%
    dplyr::mutate(day=0)
  
  #Calculate mean daily fluxes at catchment scale
  catchment_fluxes<-
    tibble::tibble(
      day       = c(rep(seq(1,365),n.years),1), 
      y_w       = output$y_w.VAR[,2], 
      y_wt      = output$y_wt.VAR[,1],
      bf_out    = output$GW_bf.VAR[,1]/land.INFO[,"area"], 
      spill_out = output$spill_vol.VAR[,2]/land.INFO[,"area"]) %>%
    dplyr::mutate(y_wt = y_wt/abs(land.INFO[,"y_cl"]), 
                  y_w =  y_w/abs(land.INFO[,"wetland_invert"])) %>%
    dplyr::group_by(day) %>%
    dplyr::summarise_all(., mean) %>%
    tidyr::gather(var,value,-day)
  
  #Create catchment output
  output_catchment<-dplyr::bind_rows(catchment_balance, catchment_duration, catchment_fluxes) %>%
    dplyr::mutate(HUC12 = HUC12.shp$HUC12,
                  WetID = WetID, 
                  scale = "catchment")
  
  #Combine output in long form
  output2<-dplyr::bind_rows(output_wetland,output_catchment)
}

plot(output$y_w.VAR[,3], cex = 0.3, ylim=c(-500, 500))
lines(output$y_wt.VAR[,1], col = 'blue')
legend(0, 0, legend=c("Wetland Stage (y_w)", "Upland Water Table (y_wt)"),
       col=c("black", "blue"), lty=1:1, cex=0.8)

plot(output$GW_local.VAR[,3])
