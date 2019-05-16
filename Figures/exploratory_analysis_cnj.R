#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Exploratory Analysis
#Coder: C. Nathan Jones (njones@sesync.org)
#Date: 5/15/2019
#Purpose: Explore initial outputs from WHC modeling effort
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1. Setup workspace---------------------------------------------------------------------------------
#Clear environmnet
remove(list=ls())

#Call required libraries
library(tidyverse)
library(gridExtra)

#define dictories of interest
results_dir<-"/nfs/WHC-data/Figure Generation/"

#2. Organize Data-----------------------------------------------------------------------------------
#download data
dmv<-read_csv(paste0(results_dir,"delmarva.csv"))
ppr<-read_csv(paste0(results_dir,"ppr.csv"))
florida<-read_csv(paste0(results_dir,"florida.csv"))

#Add landscape info
dmv$landscape<-'dmv'
ppr$landscape<-'ppr'
florida$landscape<-'florida'

#bind rows
dmv$HUC12<-as.numeric(paste(dmv$HUC12))
ppr$HUC12<-as.numeric(paste(ppr$HUC12))
florida$HUC12<-as.numeric(paste(florida$HUC12))
df<-bind_rows(dmv,florida,ppr) %>% filter(!stringr::str_detect(var,'Error'))
remove(dmv,ppr,florida)

#3.Wetland Scale Plot-------------------------------------------------------------------------------
#Wetland water level~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
waterLevel<-df %>% 
  filter(var=='y_w', scale=='wetland') %>%
  group_by(day, scale, landscape) %>%
  summarise(med=median(value, na.rm=T), 
            upr = quantile(value, 0.75, na.rm=T), 
            lwr = quantile(value, 0.25, na.rm=T)) %>%
  ggplot(aes(x=day, y=med)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), bg="grey70") +
  geom_line(lty=2) + 
  facet_grid(rows=vars(landscape), 
             scales='fixed') +
  theme_bw() +
  ylab("Median Water Level [m/m]") +
  xlab(NULL)

#Wetland water balance~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boxplots<-df %>% 
  #Select Relevant Data
  dplyr::filter(var %in% c("precip","et","sw_in","sw_out","gw_out","gw_in"), 
         scale=='wetland')  %>%
  dplyr::select(WetID, var, value, landscape) %>%
  #Normalize to precip
  tidyr::spread(var,value) %>%
  dplyr::mutate(
         et     = et/precip,
         sw_in  = sw_in/precip,
         sw_out = sw_out/precip,
         gw_in  = gw_in/precip,
         gw_out = gw_out/precip) %>%
  dplyr::select(-precip) %>%
  tidyr::gather(var,val, -c(WetID,landscape)) %>%
  #plot 
  ggplot(aes(var,val)) +
    facet_grid(rows=vars(landscape)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() + 
      ylim(0,2) + ylab("Proportion of Precip [mm/mm]") + xlab(NULL)
  

#Plot and save~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output<-grid.arrange(waterLevel,boxplots,nrow=1, top = "Wetland Scale")
ggsave(paste0(results_dir, "Output/wetland_scale.png"), 
       device = "png",width = 7, height=8, units="in")


#4.Catchment Scale Plot-------------------------------------------------------------------------------
#Catchment water level~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
waterLevel<-df %>% 
  filter(var=='y_wt', scale=='catchment') %>%
  group_by(day, scale, landscape) %>%
  summarise(med=median(value, na.rm=T), 
            upr = quantile(value, 0.75, na.rm=T), 
            lwr = quantile(value, 0.25, na.rm=T)) %>%
  ggplot(aes(x=day, y=med)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), bg="grey70") +
  geom_line(lty=2) + 
  facet_grid(rows=vars(landscape), 
             scales='fixed') +
  theme_bw() +
  ylab("Median Water Tabel Elevation [m/m]") +
  xlab(NULL)

#Wetland water balance~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boxplots<-df %>% 
  #Select Relevant Data
  dplyr::filter(var %in% c("precip","et","sw_out","gw_out"), 
                scale=='catchment')  %>%
  dplyr::select(WetID, var, value, landscape) %>%
  #Normalize to precip
  tidyr::spread(var,value) %>%
  dplyr::mutate(
    et     = et/precip,
    sw_out = sw_out/precip,
    gw_out = gw_out/precip) %>%
  dplyr::select(-precip) %>%
  tidyr::gather(var,val, -c(WetID,landscape)) %>%
  #plot 
  ggplot(aes(var,val)) +
    facet_grid(rows=vars(landscape)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() + 
      ylim(0,0.5) + ylab("Proportion of Precip [mm/mm]") + xlab(NULL)


#Plot and save~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid.arrange(waterLevel,boxplots,nrow=1, top="Catchment Scale")
ggsave(paste0(results_dir, "Output/catchment_scale.png"), 
       device = "png",width = 7, height=8, units="in")
