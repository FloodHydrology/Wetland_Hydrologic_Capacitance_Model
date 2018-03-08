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

#Set Working Directory
wd<-"/nfs/WHC-data/Soil Database (SSURGO Conversion)"
setwd(paste0(wd))

#add appropriate libarary
library(parallel)
library(foreign)
library(data.table)
library(dplyr)
library(magrittr)

#Download SWAT SSURGO database
soils<-fread("SWAT_SSURGO_Database_PPR.csv")
#colnames(soils)<-c("OBJECTID","MUID","SEQN","SNAM","S5ID","CMPPCT","NLAYERS","HYDGRP","SOL_ZMX","ANION_EXCL","SOL_CRK","TEXTURE","SOL_Z1","SOL_BD1","SOL_AWC1","SOL_K1","SOL_CBN1","CLAY1","SILT1","SAND1","ROCK1","SOL_ALB1","USLE_K1","SOL_EC1","SOL_Z2","SOL_BD2","SOL_AWC2","SOL_K2","SOL_CBN2","CLAY2","SILT2","SAND2","ROCK2","SOL_ALB2","USLE_K2","SOL_EC2","SOL_Z3","SOL_BD3","SOL_AWC3","SOL_K3","SOL_CBN3","CLAY3","SILT3","SAND3","ROCK3","SOL_ALB3","USLE_K3","SOL_EC3","SOL_Z4","SOL_BD4","SOL_AWC4","SOL_K4","SOL_CBN4","CLAY4","SILT4","SAND4","ROCK4","SOL_ALB4","USLE_K4","SOL_EC4","SOL_Z5","SOL_BD5","SOL_AWC5","SOL_K5","SOL_CBN5","CLAY5","SILT5","SAND5","ROCK5","SOL_ALB5","USLE_K5","SOL_EC5","SOL_Z6","SOL_BD6","SOL_AWC6","SOL_K6","SOL_CBN6","CLAY6","SILT6","SAND6","ROCK6","SOL_ALB6","USLE_K6","SOL_EC6","SOL_Z7","SOL_BD7","SOL_AWC7","SOL_K7","SOL_CBN7","CLAY7","SILT7","SAND7","ROCK7","SOL_ALB7","USLE_K7","SOL_EC7","SOL_Z8","SOL_BD8","SOL_AWC8","SOL_K8","SOL_CBN8","CLAY8","SILT8","SAND8","ROCK8","SOL_ALB8","USLE_K8","SOL_EC8","SOL_Z9","SOL_BD9","SOL_AWC9","SOL_K9","SOL_CBN9","CLAY9","SILT9","SAND9","ROCK9","SOL_ALB9","USLE_K9","SOL_EC9","SOL_Z10","SOL_BD10","SOL_AWC10","SOL_K10","SOL_CBN10","CLAY10","SILT10","SAND10","ROCK10","SOL_ALB10","USLE_K10","SOL_EC10","SOL_CAL1","SOL_CAL2","SOL_CAL3","SOL_CAL4","SOL_CAL5","SOL_CAL6","SOL_CAL7","SOL_CAL8","SOL_CAL9","SOL_CAL10","SOL_PH1","SOL_PH2","SOL_PH3","SOL_PH4","SOL_PH5","SOL_PH6","SOL_PH7","SOL_PH8","SOL_PH9","SOL_PH10","HYDGRP_ORIG")

####################################################################################
#Step 2: Add missing variables
####################################################################################
#Download soils data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Temporarly change wd
temp<-"/nfs/njones-data/Research Projects/GIW_Modeling/Data/SSURGO/"
setwd(paste0(temp))

#Collect chorizon file paths
files<-"ND093" #Compile list of unzipped files

#download comp files using data.table package
comp<-fread(paste0(files,"/tabular/comp.txt"), sep="|")
colnames(comp)<-c('comppct_l',	'comppct_r',	'comppct_h',	'compname',	'compkind',	'majcompflag',	'otherph',	'localphase',	'slope_l',	'slope_r',	'slope_h',	'slopelenusle_l',	'slopelenusle_r',	'slopelenusle_h',	'runoff',	'tfact',	'wei',	'weg',	'erocl',	'earthcovkind1',	'earthcovkind2',	'hydricon',	'hydricrating',	'drainagecl',	'elev_l',	'elev_r',	'elev_h',	'aspectccwise',	'aspectrep',	'aspectcwise',	'geomdesc',	'albedodry_l',	'albedodry_r',	'albedodry_h',	'airtempa_l',	'airtempa_r',	'airtempa_h',	'map_l',	'map_r',	'map_h',	'reannualprecip_l',	'reannualprecip_r',	'reannualprecip_h',	'ffd_l',	'ffd_r',	'ffd_h',	'nirrcapcl',	'nirrcapscl',	'nirrcapunit',	'irrcapcl',	'irrcapscl',	'irrcapunit',	'cropprodindex',	'constreeshrubgrp',	'wndbrksuitgrp',	'rsprod_l',	'rsprod_r',	'rsprod_h',	'foragesuitgrpid',	'wlgrain',	'wlgrass',	'wlherbaceous',	'wlshrub',	'wlconiferous',	'wlhardwood',	'wlwetplant',	'wlshallowwat',	'wlrangeland',	'wlopenland',	'wlwoodland',	'wlwetland',	'soilslippot',	'frostact',	'initsub_l',	'initsub_r',	'initsub_h',	'totalsub_l',	'totalsub_r',	'totalsub_h',	'hydgrp',	'corcon',	'corsteel',	'taxclname',	'taxorder',	'taxsuborder',	'taxgrtgroup',	'taxsubgrp',	'taxpartsize',	'taxpartsizemod',	'taxceactcl',	'taxreaction',	'taxtempcl',	'taxmoistscl',	'taxtempregime',	'soiltaxedition',	'castorieindex',	'flecolcomnum',	'flhe',	'flphe',	'flsoilleachpot',	'flsoirunoffpot',	'fltemik2use',	'fltriumph2use',	'indraingrp',	'innitrateleachi',	'misoimgmtgrp',	'vasoimgtgrp',	'mukey',	'cokey')

#download chorizon files
chorizon<-fread(paste0(files,"/tabular/chorizon.txt"), sep="|")
colnames(chorizon)<-c('hzname',	'desgndisc',	'desgnmaster',	'desgnmasterprime',	'desgnvert',	'hzdept_l',	'hzdept_r',	'hzdept_h',	'hzdepb_l',	'hzdepb_r',	'hzdepb_h',	'hzthk_l',	'hzthk_r',	'hzthk_h',	'fraggt10_l',	'fraggt10_r',	'fraggt10_h',	'frag3to10_l',	'frag3to10_r',	'frag3to10_h',	'sieveno4_l',	'sieveno4_r',	'sieveno4_h',	'sieveno10_l',	'sieveno10_r',	'sieveno10_h',	'sieveno40_l',	'sieveno40_r',	'sieveno40_h',	'sieveno200_l',	'sieveno200_r',	'sieveno200_h',	'sandtotal_l',	'sandtotal_r',	'sandtotal_h',	'sandvc_l',	'sandvc_r',	'sandvc_h',	'sandco_l',	'sandco_r',	'sandco_h',	'sandmed_l',	'sandmed_r',	'sandmed_h',	'sandfine_l',	'sandfine_r',	'sandfine_h',	'sandvf_l',	'sandvf_r',	'sandvf_h',	'silttotal_l',	'silttotal_r',	'silttotal_h',	'siltco_l',	'siltco_r',	'siltco_h',	'siltfine_l',	'siltfine_r',	'siltfine_h',	'claytotal_l',	'claytotal_r',	'claytotal_h',	'claysizedcarb_l',	'claysizedcarb_r',	'claysizedcarb_h',	'om_l',	'om_r',	'om_h',	'dbtenthbar_l',	'dbtenthbar_r',	'dbtenthbar_h',	'dbthirdbar_l',	'dbthirdbar_r',	'dbthirdbar_h',	'dbfifteenbar_l',	'dbfifteenbar_r',	'dbfifteenbar_h',	'dbovendry_l',	'dbovendry_r',	'dbovendry_h',	'partdensity',	'ksat_l',	'ksat_r',	'ksat_h',	'awc_l',	'awc_r',	'awc_h',	'wtenthbar_l',	'wtenthbar_r',	'wtenthbar_h',	'wthirdbar_l',	'wthirdbar_r',	'wthirdbar_h',	'wfifteenbar_l',	'wfifteenbar_r',	'wfifteenbar_h',	'wsatiated_l',	'wsatiated_r',	'wsatiated_h',	'lep_l',	'lep_r',	'lep_h',	'll_l',	'll_r',	'll_h',	'pi_l',	'pi_r',	'pi_h',	'aashind_l',	'aashind_r',	'aashind_h',	'kwfact',	'kffact',	'caco3_l',	'caco3_r',	'caco3_h',	'gypsum_l',	'gypsum_r',	'gypsum_h',	'sar_l',	'sar_r',	'sar_h',	'ec_l',	'ec_r',	'ec_h',	'cec7_l',	'cec7_r',	'cec7_h',	'ecec_l',	'ecec_r',	'ecec_h',	'sumbases_l',	'sumbases_r',	'sumbases_h',	'ph1to1h2o_l',	'ph1to1h2o_r',	'ph1to1h2o_h',	'ph01mcacl2_l',	'ph01mcacl2_r',	'ph01mcacl2_h',	'freeiron_l',	'freeiron_r',	'freeiron_h',	'feoxalate_l',	'feoxalate_r',	'feoxalate_h',	'extracid_l',	'extracid_r',	'extracid_h',	'extral_l',	'extral_r',	'extral_h',	'aloxalate_l',	'aloxalate_r',	'aloxalate_h',	'pbray1_l',	'pbray1_r',	'pbray1_h',	'poxalate_l',	'poxalate_r',	'poxalate_h',	'ph2osoluble_l',	'ph2osoluble_r',	'ph2osoluble_h',	'ptotal_l',	'ptotal_r',	'ptotal_h',	'excavdifcl',	'excavdifms',	'cokey',	'chkey')

#Manipulate soils data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Use dplyr package to select major component of each map unit
comp<-data.frame(comp)
comp<-comp[,c("mukey","cokey","comppct_r")]
comp<-comp %>%
  group_by(mukey) %>%
  filter(comppct_r==max(comppct_r)) %>%
  filter(mukey==unique(mukey)) %>%
  filter(row_number()==1) %>%
  arrange(mukey,cokey,comppct_r)

#Define wthrdbar_r variable for each cokey
df1<-chorizon %>% group_by(cokey) %>% summarise(wthirdbar_r1=wthirdbar_r[order(wthirdbar_r)==1]) 
df2<-chorizon %>% group_by(cokey) %>% summarise(wthirdbar_r2=wthirdbar_r[order(wthirdbar_r)==2])
df_wthirdbar_r<-left_join(df1,df2)

#Define wthrdbar_r variable for each cokey
df1<-chorizon %>% group_by(cokey) %>% summarise(wfifteenbar_r1=wfifteenbar_r[order(wfifteenbar_r)==1]) 
df2<-chorizon %>% group_by(cokey) %>% summarise(wfifteenbar_r2=wfifteenbar_r[order(wfifteenbar_r)==2])
df_wfifteenbar_r<-left_join(df1,df2)

#Merge soil information
chorizon<-left_join(df_wthirdbar_r,df_wfifteenbar_r)
comp<-left_join(comp,chorizon)
soils$mukey<-soils$MUID
soils<-data.frame(soils)
soils<-merge(soils,comp, by.x="Mukey",by.y="mukey")
remove(list=ls()[ls()!="soils"])

####################################################################################
#Step 3: Vertically integrate relevant terms
####################################################################################
#prep soils data 
soils$ZMX<-soils$SOL_ZMX
soils$SOL_ZMX<-NULL
soils$SOL_X1<-soils$SOL_Z1
soils$SOL_X2<-soils$SOL_Z2-soils$SOL_Z1
soils$SOL_X3<-soils$SOL_Z3-soils$SOL_Z1
soils$SOL_X4<-soils$SOL_Z4-soils$SOL_Z1
soils$SOL_X5<-soils$SOL_Z5-soils$SOL_Z1
soils$SOL_X6<-soils$SOL_Z6-soils$SOL_Z1
soils$SOL_X7<-soils$SOL_Z7-soils$SOL_Z1
soils$SOL_X8<-soils$SOL_Z8-soils$SOL_Z1
soils$SOL_X9<-soils$SOL_Z9-soils$SOL_Z1
soils$SOL_X10<-soils$SOL_Z10-soils$SOL_Z1

#define soil properties
soils$p_bd<-rowSums(soils[,grep("SOL_BD",colnames(soils))]*soils[,grep("SOL_X", colnames(soils))],na.rm=T)/rowSums(soils[,grep("SOL_X", colnames(soils))],na.rm=T)
soils$n<-soils$p_bd/2.65
soils$s_fc<-rowSums(soils[,grep("wthirdbar",colnames(soils))]*soils[,168:169], na.rm = T)/rowSums(soils[,168:169],na.rm=T)
soils$s_w<-rowSums(soils[,grep("wfifteenbar",colnames(soils))]*soils[,168:169], na.rm = T)/rowSums(soils[,168:169],na.rm=T)
soils$clay<-rowSums(soils[,grep("CLAY",colnames(soils))]*soils[,grep("SOL_X", colnames(soils))], na.rm = T)/rowSums(soils[,grep("SOL_X", colnames(soils))],na.rm=T)
soils$sand<-rowSums(soils[,grep("SAND",colnames(soils))]*soils[,grep("SOL_X", colnames(soils))], na.rm = T)/rowSums(soils[,grep("SOL_X", colnames(soils))],na.rm=T)
soils$ksat<-rowSums(soils[,grep("SOL_K",colnames(soils))]*soils[,grep("SOL_X", colnames(soils))], na.rm = T)/rowSums(soils[,grep("SOL_X", colnames(soils))],na.rm=T)
soils$y_rd<-soils[,"SOL_Z1"]

#Calculate confining layer depth
#high ksat soils
fun<-function(x){length(x[is.na(x)==F])}
soils[soils==0]<-NA
soils<-soils %>% rowwise() %>% mutate(index=fun(c(SOL_K1,SOL_K2,SOL_K3, SOL_K4,SOL_K5, SOL_K6, SOL_K7, SOL_K8,SOL_K9, SOL_K10)))
#low ksat soils
soils<-soils %>% rowwise() %>% mutate(ksat_min=min(c(SOL_K1,SOL_K2,SOL_K3, SOL_K4,SOL_K5, SOL_K6, SOL_K7, SOL_K8,SOL_K9, SOL_K10),na.rm=T))
fun<-function(x){which(x<=4.5)[1]}
soils[soils$ksat_min<4.5,]<-soils[soils$ksat_min<4.5,] %>% rowwise() %>% mutate(index=fun(c(SOL_K1,SOL_K2,SOL_K3, SOL_K4,SOL_K5, SOL_K6, SOL_K7, SOL_K8,SOL_K9, SOL_K10)))
fun<-function(x){t(soils[x,grep("SOL_Z", colnames(soils))])[soils$index[x]]}
soils$y_cl<-sapply(seq(1,nrow(soils)), fun)
soils$y_cl<-as.numeric(soils$y_cl)

####################################################################################
#Step 4: Critical depth calculation [FRED]
####################################################################################


####################################################################################
#Step 5: Export Soils Database!
####################################################################################
wd<-"/nfs/njones-data/Research Projects/GIW_Modeling/Model Development/Soil Database (SSURGO Conversion)/"
setwd(paste0(wd))
soils$MUID<-soils$Mukey
soils$n<-abs(soils$n)
soils$clay[soils$clay<0]<-1
soils$clay[soils$clay>100]<-100
soils$sand[soils$sand<0]<-1
soils$sand[soils$sand>100]<-100
soils$ksat<-abs(soils$ksat)
soils<-soils[,c("MUID","y_cl","y_rd","s_fc","s_w","n","clay","ksat")]
write.csv(soils, "WHC_Soils_Input_PPR.csv")
