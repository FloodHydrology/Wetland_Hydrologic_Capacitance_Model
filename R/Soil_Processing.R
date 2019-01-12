###################################################################################
#Name: SSURGO Soils Conversion
#Coder: C. Nathan Jones
#Date: 5/19/2016
#Purpose: Convert SSURGO Soils databse into WHC model inputs
##################################################################################

#to add new soils, use spreadsheet from this website: https://hydrologicou.wordpress.com/software/swat-soil/

####################################################################################
#Step 1: Setup Workspace
####################################################################################
#Clear Memory
rm(list=ls(all=TRUE))

#add appropriate libarary
library(parallel)
library(foreign)
library(data.table)
library(dplyr)
library(magrittr)

#Download Fred's get_yc function!
source("R/get_yc.R")

####################################################################################
#Step 2: Create function to calculate values---------------------------------------
####################################################################################
#Create function 
fun<-function(n){

  #download appropriate tables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #download component files
  component<-fread(paste0(files[n],"/tabular/comp.txt"), sep="|")
  colnames(component)<-c('comppct_l',	'comppct_r',	'comppct_h',	'compname',	'compkind',	'majcompflag',	'otherph',	'localphase',	'slope_l',	'slope_r',	'slope_h',	'slopelenusle_l',	'slopelenusle_r',	'slopelenusle_h',	'runoff',	'tfact',	'wei',	'weg',	'erocl',	'earthcovkind1',	'earthcovkind2',	'hydricon',	'hydricrating',	'drainagecl',	'elev_l',	'elev_r',	'elev_h',	'aspectccwise',	'aspectrep',	'aspectcwise',	'geomdesc',	'albedodry_l',	'albedodry_r',	'albedodry_h',	'airtempa_l',	'airtempa_r',	'airtempa_h',	'map_l',	'map_r',	'map_h',	'reannualprecip_l',	'reannualprecip_r',	'reannualprecip_h',	'ffd_l',	'ffd_r',	'ffd_h',	'nirrcapcl',	'nirrcapscl',	'nirrcapunit',	'irrcapcl',	'irrcapscl',	'irrcapunit',	'cropprodindex',	'constreeshrubgrp',	'wndbrksuitgrp',	'rsprod_l',	'rsprod_r',	'rsprod_h',	'foragesuitgrpid',	'wlgrain',	'wlgrass',	'wlherbaceous',	'wlshrub',	'wlconiferous',	'wlhardwood',	'wlwetplant',	'wlshallowwat',	'wlrangeland',	'wlopenland',	'wlwoodland',	'wlwetland',	'soilslippot',	'frostact',	'initsub_l',	'initsub_r',	'initsub_h',	'totalsub_l',	'totalsub_r',	'totalsub_h',	'hydgrp',	'corcon',	'corsteel',	'taxclname',	'taxorder',	'taxsuborder',	'taxgrtgroup',	'taxsubgrp',	'taxpartsize',	'taxpartsizemod',	'taxceactcl',	'taxreaction',	'taxtempcl',	'taxmoistscl',	'taxtempregime',	'soiltaxedition',	'castorieindex',	'flecolcomnum',	'flhe',	'flphe',	'flsoilleachpot',	'flsoirunoffpot',	'fltemik2use',	'fltriumph2use',	'indraingrp',	'innitrateleachi',	'misoimgmtgrp',	'vasoimgtgrp',	'mukey',	'cokey')
  
  #download chorizon files
  chorizon<-fread(paste0(files[n],"/tabular/chorizon.txt"), sep="|")
  colnames(chorizon)<-c('hzname',	'desgndisc',	'desgnmaster',	'desgnmasterprime',	'desgnvert',	'hzdept_l',	'hzdept_r',	'hzdept_h',	'hzdepb_l',	'hzdepb_r',	'hzdepb_h',	'hzthk_l',	'hzthk_r',	'hzthk_h',	'fraggt10_l',	'fraggt10_r',	'fraggt10_h',	'frag3to10_l',	'frag3to10_r',	'frag3to10_h',	'sieveno4_l',	'sieveno4_r',	'sieveno4_h',	'sieveno10_l',	'sieveno10_r',	'sieveno10_h',	'sieveno40_l',	'sieveno40_r',	'sieveno40_h',	'sieveno200_l',	'sieveno200_r',	'sieveno200_h',	'sandtotal_l',	'sandtotal_r',	'sandtotal_h',	'sandvc_l',	'sandvc_r',	'sandvc_h',	'sandco_l',	'sandco_r',	'sandco_h',	'sandmed_l',	'sandmed_r',	'sandmed_h',	'sandfine_l',	'sandfine_r',	'sandfine_h',	'sandvf_l',	'sandvf_r',	'sandvf_h',	'silttotal_l',	'silttotal_r',	'silttotal_h',	'siltco_l',	'siltco_r',	'siltco_h',	'siltfine_l',	'siltfine_r',	'siltfine_h',	'claytotal_l',	'claytotal_r',	'claytotal_h',	'claysizedcarb_l',	'claysizedcarb_r',	'claysizedcarb_h',	'om_l',	'om_r',	'om_h',	'dbtenthbar_l',	'dbtenthbar_r',	'dbtenthbar_h',	'dbthirdbar_l',	'dbthirdbar_r',	'dbthirdbar_h',	'dbfifteenbar_l',	'dbfifteenbar_r',	'dbfifteenbar_h',	'dbovendry_l',	'dbovendry_r',	'dbovendry_h',	'partdensity',	'ksat_l',	'ksat_r',	'ksat_h',	'awc_l',	'awc_r',	'awc_h',	'wtenthbar_l',	'wtenthbar_r',	'wtenthbar_h',	'wthirdbar_l',	'wthirdbar_r',	'wthirdbar_h',	'wfifteenbar_l',	'wfifteenbar_r',	'wfifteenbar_h',	'wsatiated_l',	'wsatiated_r',	'wsatiated_h',	'lep_l',	'lep_r',	'lep_h',	'll_l',	'll_r',	'll_h',	'pi_l',	'pi_r',	'pi_h',	'aashind_l',	'aashind_r',	'aashind_h',	'kwfact',	'kffact',	'caco3_l',	'caco3_r',	'caco3_h',	'gypsum_l',	'gypsum_r',	'gypsum_h',	'sar_l',	'sar_r',	'sar_h',	'ec_l',	'ec_r',	'ec_h',	'cec7_l',	'cec7_r',	'cec7_h',	'ecec_l',	'ecec_r',	'ecec_h',	'sumbases_l',	'sumbases_r',	'sumbases_h',	'ph1to1h2o_l',	'ph1to1h2o_r',	'ph1to1h2o_h',	'ph01mcacl2_l',	'ph01mcacl2_r',	'ph01mcacl2_h',	'freeiron_l',	'freeiron_r',	'freeiron_h',	'feoxalate_l',	'feoxalate_r',	'feoxalate_h',	'extracid_l',	'extracid_r',	'extracid_h',	'extral_l',	'extral_r',	'extral_h',	'aloxalate_l',	'aloxalate_r',	'aloxalate_h',	'pbray1_l',	'pbray1_r',	'pbray1_h',	'poxalate_l',	'poxalate_r',	'poxalate_h',	'ph2osoluble_l',	'ph2osoluble_r',	'ph2osoluble_h',	'ptotal_l',	'ptotal_r',	'ptotal_h',	'excavdifcl',	'excavdifms',	'cokey',	'chkey')
  
  #Aggregate values~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Calculate horizon thickness
  chorizon$thickness<-chorizon$hzdepb_r-chorizon$hzdept_r
  
  #Clean up particle desnsity
  chorizon$partdensity[is.na(chorizon$partdensity)]<-2.65
  
  #Summorize depth weighted values
  depth_average_horizons<- chorizon        %>% 
                           group_by(cokey) %>%
                           summarise(y_cl = ifelse(min(ksat_r, na.rm=T)<1.25,min(hzdepb_r[ksat_r<1.25], na.rm=T), max(hzdepb_r, na.rm=T)), #clay layer depth [cm]
                                     clay = sum(claytotal_r*thickness, na.rm=T)/sum(thickness, na.rm=T),   # % clay
                                     sand = sum(sandtotal_r*thickness, na.rm=T)/sum(thickness, na.rm=T),
                                     ksat = sum(ksat_r*thickness, na.rm=T)/sum(thickness, na.rm=T),        # ksat (um/s)
                                     s_fc = sum(wthirdbar_r*thickness, na.rm=T)/sum(thickness, na.rm=T),   # field capacity (%)
                                     s_w  = sum(wfifteenbar_r*thickness, na.rm=T)/sum(thickness, na.rm=T), # wilting point  (%)
                                     n    = sum((1 - (dbthirdbar_r/partdensity))*thickness, na.rm=T)/sum(thickness, na.rm=T) #porosity (%)
                                     )
  component<-merge(component, depth_average_horizons, by="cokey")
  
  #Summarize area weighted values
  area_average_mu<- component       %>%
                    group_by(mukey) %>%
                    summarise(y_cl = sum(y_cl*comppct_r, na.rm=T)/100, #clay layer depth [cm]
                              clay = sum(clay*comppct_r, na.rm=T)/100, # % clay
                              sand = sum(sand*comppct_r, na.rm=T)/100, # % sand
                              ksat = sum(ksat*comppct_r, na.rm=T)/100, # ksat (um/s)
                              s_fc = sum(s_fc*comppct_r, na.rm=T)/100, # field capacity (%)
                              s_w  = sum(s_w*comppct_r, na.rm=T)/100,  # wilting point  (%)
                              n    = sum(n*comppct_r, na.rm=T)/100     #porosity (%)
                    )
  
  #Estimate yc values using get_yc from Fred
  get_yc<-function(mukey){
    clay<-area_average_mu$clay[area_average_mu$mukey==mukey]
    sand<-area_average_mu$sand[area_average_mu$mukey==mukey]
    y_c<-get.yc(clay, sand)
    c(mukey, y_c)
  }
  y_c<-lapply(area_average_mu$mukey, get_yc)
  y_c<-data.frame(do.call(rbind,y_c))
  colnames(y_c)<-c("mukey", "y_c")
  area_average_mu<-left_join(area_average_mu, y_c)
  
  #Export results output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  area_average_mu
}

####################################################################################
#Step 3: Execute Function-----------------------------------------------------------
####################################################################################
#Temporarly change wd
temp<-"/nfs/WHC-data/Data/SSURGO"

#create list of files containing soils data
files<-list.files(temp, full.names = T)

#Create function to catch errors
execute<-function(n){tryCatch(fun(n), error=function(e) c(files[n], rep(0,6)))}

#Run 
t0<-Sys.time()
n.cores<-detectCores()
x<-mclapply(seq(1,length(files)), execute, mc.cores=n.cores)
tf<-Sys.time()
tf-t0

#Gather data
output<-do.call(rbind,x)

#Export Ouput
write.csv(output, "/nfs/WHC-data/WHC_Soils_Input.csv")
write.csv(output, "data/soils_lookup.csv")

