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

WetID   <- 1
n.years <- 100
giw.ID  <- WetID

# 2. Load Necessary Data ----------------------------------------------------------
# Run all of HUC12_Model_Input_PPR.R
# setwd("~/Wetland_Hydrologic_Capacitance_Model")
# source("~/Wetland_Hydrologic_Capacitance_Model/analysis/HUC12_Model_Input_PPR.R")

# 3. Generate .INFO Files ---------------------------------------------------------
  # Run all of regional_analysis.R using WetID <- 26
  # Run until lumped.INFO is run

WetID <- 60

setwd("/nfs/WHC-data/Sandbox Model")
# save.image(file = '20190401_PPR_ProcessedInputs_ID26_all.RData')
# save(HUC12.shp, area, giw.INFO, land.INFO, lumped.INFO, pet.VAR, precip.VAR, volume,  snowmelt.VAR,
#       file = '20190401_PPR_ProcessedInputs_ID26_selected.RData')
# load(file = '20190401_PPR_ProcessedInputs_ID26_selected.RData')

# save.image(file = '20190417_PPR_ProcessedInputs_ID60_all.RData')
# save(HUC12.shp, area, giw.INFO, land.INFO, lumped.INFO, pet.VAR, precip.VAR, volume,  snowmelt.VAR,
#      file = '20190417_PPR_ProcessedInputs_ID60_selected.RData')
load(file = '20190417_PPR_ProcessedInputs_ID60_selected.RData')


# land.INFO[,'Sy'] <- 1
land.INFO[,'RD'] <- land.INFO[,'RD'] /4
# precip.VAR <- precip.VAR
# pet.VAR <- pet.VAR /2
land.INFO[,"kb"]<-              0.046 /1000
# land.INFO[,"y_c"]<-         -100
# giw.INFO[,"y_c"]<-         -100
land.INFO[,"y_cl"]<-         -5000
# giw.INFO[,"y_cl"]<-         -500
# land.INFO[,"GW_bf_0"] <- 0.01

# 4. Source and Run WHC Model ----------------------------------------------------
# source("~/Wetland_Hydrologic_Capacitance_Model/R/WHC_Sandbox.R")
source("~/Wetland_Hydrologic_Capacitance_Model/R/WHC_3g.R")

output <- wetland.hydrology(giw.INFO, 
                            land.INFO, 
                            lumped.INFO, 
                            snowmelt.VAR, 
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
                   R_wt       = sum(output$R.VAR[,1])/n.years,
                   R_lmz      = sum(output$R_lm.VAR[,1])/n.years,
                   pet        = sum(pet.VAR)/(length(pet.VAR)/365),
                   et         = (sum(output$ET_lm.VAR[,1])+sum(output$ET_wt.VAR[,1]))/n.years,
                   et_lm      = sum(output$ET_lm.VAR[,1])/n.years,
                   et_hm      = sum(output$exfil.VAR[,1])/n.years,
                   et_wt      = sum(output$ET_wt.VAR[,1])/n.years,
                   sw_out     = sum(output$spill_vol.VAR[,2])/land.INFO[,"area"]/n.years,
                   gw_out     = sum(-1*output$GW_bf.VAR[,1])/n.years,
                   closure    = sum(precip.VAR)/(length(precip.VAR)/365) -
                                (sum(output$ET_lm.VAR[,1])+sum(output$ET_wt.VAR[,1]))/n.years -
                                sum(output$spill_vol.VAR[,2])/land.INFO[,"area"]/n.years -
                                sum(-1*output$GW_bf.VAR[,1])/n.years) %>%
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



# plot(output$GW_local.VAR[,3])
# plot(output$GW_bf.VAR[,1])

start_year <- 15
end_year <- 20
start_ind <- start_year*365
end_ind <- end_year*365

# Climate Plots -----------------------------------------------------------------------
par(mfrow=c(3,1))
plot(precip.VAR[start_ind:end_ind], type = 'line', main = 'Precip')
rect((90 + 365*seq(from=0,to=end_year-1)), 0, 330+365*seq(from=0,to=end_year-1), 60, 
     col= rgb(0.1,0,1,alpha=0.2), border = NA)
plot(snowmelt.VAR[start_ind:end_ind], type = 'line', main = 'Snowmelt')
rect((90 + 365*seq(from=0,to=end_year-1)), 0, 330+365*seq(from=0,to=end_year-1), 60, 
     col= rgb(0.1,0,1,alpha=0.2), border = NA)
plot(pet.VAR[start_ind:end_ind], type = 'line', main = 'PET')
rect((90 + 365*seq(from=0,to=end_year-1)), 0, 330+365*seq(from=0,to=end_year-1), 60, 
     col= rgb(0.1,0,1,alpha=0.2), border = NA)

# Flux Plots -----------------------------------------------------------------------
par(mfrow=c(3,2), mar = c(2,3,1,1)) 
# plot(output$R.VAR[,"land"]+output$loss_lm.VAR[,"land"]*(1-exp(-land.INFO[,"kb"])))
# plot(output$y_wt.VAR[,"land"] - output$land.INFO[, "y_cl"])
plot(output$R.VAR[start_ind:end_ind,'land'], type = "line",
     main = 'R', ylab="R")
plot(output$R_lm.VAR[start_ind:end_ind,'land'], type = "line",
    main ='R_lm', ylab='R_l')
plot(output$loss_lm.VAR[start_ind:end_ind,'land'], type = "line",
    main ='loss_lm', ylab='loss_lm')
# plot(output$GW_bf.VAR[start_ind:end_ind,'land'], type = "line",
#      main ='baseflow', ylab='baseflow')
plot(output$GW_local.VAR[start_ind:end_ind,'land'], type = "line",
     main ='GW local', ylab='GW Local')
plot(output$ET_wt.VAR[start_ind:end_ind,'land'], type = "line", 
     main ='ET Components')
  lines(output$exfil.VAR[start_ind:end_ind,'land'], type = "line", col = 'blue')
  lines(output$ET_lm.VAR[start_ind:end_ind,'land'], type = "line", col = 'red')
  legend('topright', legend = c('ET_WT', 'ET_HMZ', 'ET_LMZ'), lty = 1, col = c('black', 'blue', 'red'))
plot(output$y_sat.VAR[start_ind:end_ind,'land'], col = 'blue', type='l')
  title('upland WT and HMZ/LMZ boundary')
  lines(output$y_hmz.VAR[start_ind:end_ind,'land'], col = 'red', type='l', lty = 3)
  legend('topright', legend = c('WT', 'HMZ/LMZ Boundary'), 
          lty = 1, col = c('blue', 'red'))
  
# par(mfrow=c(1,1)) 
# plot(output$s_ex.VAR[start_ind:end_ind,'land'], type = "line")
# title('s_ex')

# Water Table ------------------------------------------------------------------------
par(mfrow=c(3,1))
plot(output$y_sat.VAR[start_ind:end_ind,'land'], col = 'blue', type='l')
  title('upland WT and HMZ/LMZ boundary')
  lines(output$y_hm.VAR[start_ind:end_ind,'land'], col = 'red', type='l', lty = 3)
  legend('topright', legend = c('WT', 'HMZ/LMZ Boundary'), 
       lty = 1, col = c('blue', 'red'))
plot(output$y_w.VAR[start_ind:end_ind,'wetland'], col = 'blue', type='l')
  title('lumped wetland stage')
plot(output$y_w.VAR[start_ind:end_ind,3], col = 'blue', type='l')
  lines(output$y_hm.VAR[start_ind:end_ind,3], col = 'red', type='l', lty = 3)
  title('wetland stage')
  abline(h = giw.INFO['invert'])
  rect((90 + 365*seq(from=0,to=end_year-1)), giw.INFO['y_cl'], 330+365*seq(from=0,to=end_year-1), 0, 
       col= rgb(0.1,0,1,alpha=0.2), border = NA)

# plot(output$y_w.VAR[start_ind:end_ind,3], col = 'red', type='l')
# # lines(output$y_w.VAR[(3*365):(5*365),2], col = 'red')
# title('upland WT')
# legend(0, 0, legend=c("Wetland of Interest Stage (y_w)", "Upland Water Table (y_wt)", "Lumped Wetland"),
#        col=c("black", "blue", "red"), lty=1:1, cex=0.8, xpd = T)
# abline(h = land.INFO[,'y_c'])
# abline(h = land.INFO[,'y_cl'])


  