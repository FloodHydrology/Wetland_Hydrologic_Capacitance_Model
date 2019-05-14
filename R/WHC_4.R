###################################################################################
#Name: Wetland Hydrologic Capacitance Model 
#Coder: C. Nathan Jones
#Date: 5/13/2016
#Purpose: Apply the model to a REAL wetland!!!  (Baltimore Corner)
##################################################################################

# MODEL 4 -> Revised Laio Vadose Zone for Lumped upland

####################################################################################
# Step 1: Setup Workspace-----------------------------------------------------------
####################################################################################
#Start Function
wetland.hydrology<-function(giw.INFO, land.INFO, lumped.INFO, snowmelt.VAR, precip.VAR, pet.VAR, n.years, area, volume, giw.ID){
  
  # Need to touch up later!!!
  s_star <- (land.INFO[, 's_fc'] - land.INFO[, 's_wilt'] ) * 0.673 + land.INFO[, 's_wilt']  # very rough estimation; interpolated relationship from data in paper
  psi_fc <- land.INFO[,"psi"]*land.INFO[,"s_fc"]^(-4.05) #see eq. 26 in Laio 2011
  
  ####################################################################################
  # Step 2: Create stage-storage and stage-area relationships------------------------
  ####################################################################################
  #Locate correct GIW.INFO info
  giw.INFO<-matrix(giw.INFO[giw.ID,],nrow=1,  dimnames = list(c(1), colnames(giw.INFO)))
  
  #Create relationships for GIW module~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #isolate area and volume matrices
  temp<-(length(area[,giw.ID][area[,giw.ID]!=max(area[,giw.ID], na.rm=T)])+1)
  temp<-data.frame(seq(giw.INFO[,"invert"],0,50),
                   area[1:temp,giw.ID]*1000^2,
                   volume[1:temp,giw.ID]*1000^3)
  colnames(temp)<-c("stage", "area", "volume")
  
  #create volume matrix
  vol_giw.INFO<-matrix(0, ncol=2, nrow=(10+length(temp[,1])))
  colnames(vol_giw.INFO)<-c("stage","vol_mm^3")
  vol_giw.INFO[,"stage"]<-c(seq(giw.INFO[,"y_cl"],giw.INFO[,"invert"], length.out=10), temp[,1])
  
  #Calculate volume for each step
  for(i in 2:length(vol_giw.INFO[,1])){
    #Define wetland ID (for now just 1)
    wet.INFO<-1
    
    #Define Stage
    y<-vol_giw.INFO[i,"stage"]
    
    #Estimate Subsurface Storage~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(y<=giw.INFO[,"invert"]){
      #calculate depth of saturated zone
      y_sat<-y-giw.INFO[,"psi"] #saturated zone
      
      #estimate depth of HMZ & LMZ (Laio et al 2011)
      if(y_sat>giw.INFO[wet.INFO,"y_c"]){
        #high moisture zone is above the ground
        y_hm<-giw.INFO[wet.INFO,"invert"]
      }else{
        psi_fc<-giw.INFO[wet.INFO,"psi"]*giw.INFO[wet.INFO,"s_fc"]^(-4.05) #see eq. 26 in Laio 2011
        A<-(psi_fc-giw.INFO[wet.INFO,"psi"]-giw.INFO[wet.INFO,"y_c"])/(psi_fc-giw.INFO[wet.INFO,"psi"]-giw.INFO[wet.INFO,"y_c"]-(5*abs(giw.INFO[wet.INFO,"RD"])))
        if(y_sat>(-5*abs(giw.INFO[wet.INFO,"RD"])+psi_fc-giw.INFO[wet.INFO,"psi"])){
          A<-abs(A) #Change for now.
          y_hm<-(1-A^(3/4))*(y_sat-giw.INFO[wet.INFO,"y_c"])-(A^2*(1-A^(-1/4)))/(-giw.INFO[wet.INFO,"y_c"]+psi_fc-giw.INFO[wet.INFO,"psi"])*(y_sat-giw.INFO[wet.INFO,"y_c"])^2
        }else{
          y_hm<-y_sat-psi_fc+giw.INFO[wet.INFO,"psi"]
        }
      }
      
      #Calculate Storage
      vol_giw.INFO[i,"vol_mm^3"]<-(y-giw.INFO[,"y_cl"])*giw.INFO[,"Sy"]*giw.INFO[,"area_wetland"]+
        (y_hm-y)*giw.INFO[,"s_fc"]*giw.INFO[,"area_wetland"]
    }
    
    #Estimate Surface Storage~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(y>giw.INFO[,"invert"]){
      subsurface_vol<-vol_giw.INFO[vol_giw.INFO[,"stage"]==giw.INFO[,"invert"],2][1]
      vol_giw.INFO[i,"vol_mm^3"]<-temp$volume[which(abs(temp$stage-y)==min(abs(temp$stage-y)))]+subsurface_vol
    }
  }
  
  #Create interpolation functions
  stage2vol_giw.fun<-approxfun(vol_giw.INFO[,1], vol_giw.INFO[,2])
  vol2stage_giw.fun<-approxfun(vol_giw.INFO[,2], vol_giw.INFO[,1])
  stage2area_giw.fun<-approxfun(c(c(giw.INFO[,"y_cl"]),seq(to=0,from=giw.INFO[,"invert"], by=50)),
                                c(0,area[,giw.ID][area[,giw.ID]!=max(area[,giw.ID])], max(area[,giw.ID]))
  )
  
  
  #Create relationship for upland module~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Note, this is only for the surface water component
  #Remove individual GIW from volume and area df
  volume<-volume[,-giw.ID]
  area<-area[,-giw.ID]
  
  #Calculate Depth of wetland
  if(length(volume)!=100){
    depth<-0.05*(length(rowSums(area)[rowSums(area)!=max(rowSums(area), na.rm=T)])+1)
    
    #create dataframe of depth and surface area (in m)
    area<-data.frame(c(0,rowSums(area[1:(depth/0.05),])), seq(-depth,0, 0.05))
    colnames(area)<-c("area_m2","y_w")
    
    #create dataframe of depth and volume (in m)
    volume<-data.frame(c(0,rowSums(volume[1:((depth/0.05)),])), seq(-depth,0, 0.05))
    colnames(volume)<-c("volume_m3","y_w")
    #volume$volume_m3<-volume$volume_m3+volume$volume_m3[3]
    
    #create interpolation functions
    stage2area.fun<-approxfun(area[,2]*1000, area[,1]*(1000^2))
    stage2volume.fun<-approxfun(volume[,2]*1000, volume[,1]*(1000^3))
    volume2stage.fun<-approxfun(volume[,1]*(1000^3),volume[,2]*1000, yleft = min(volume[,2]*1000), yright=max(volume[,2]*1000))
  }else{
    #Create Dataframe for depth and area
    depth<-0.05*length(area[area!=max(area)]+1)
    area<-data.frame(area_m2 = c(0,area[1:(depth/0.05)]),
                     y_w     = seq(-depth, 0, 0.05))
    volume<-data.frame(volume_m3 = c(0,volume[1:((depth/0.05))]),
                       y_w       = seq(-depth, 0, 0.05))
    
    #create interpolation functions
    stage2area.fun<-approxfun(area[,2]*1000, area[,1]*(1000^2))
    stage2volume.fun<-approxfun(volume[,2]*1000, volume[,1]*(1000^3))
    volume2stage.fun<-approxfun(volume[,1]*(1000^3),volume[,2]*1000, yleft = min(volume[,2]*1000), yright=max(volume[,2]*1000))
    
  }
  
  ####################################################################################
  # Step 3: Create dynamic variables (.VAR)-------------------------------------------
  ####################################################################################
  #For this iteration of the model, set n.wetlands =1
  n.wetlands<-1
  
  if(giw.INFO==0){giw.INFO<-matrix(0, ncol=2, dimnames = list(c(1),c("giw.ID", "area_wetland")))}
  
  #create vector of dynamic variable names
  dynamic_variables<-c(
    "y_wt",       #Water table elevation (mm)
    "y_sat",      #elevation of saturated zone (mm)
    "y_hmz",      #elevation of HMZ/LMZ boundary (mm)
    "y_w",        #Surface water elevation(mm)
    "V_w",        #Volume of wetland water (mm^3)
    "dy_wt",      #change in water table elevation (mm)
    "dV_w",       #change in surface water volume (mm^3)
    "r_a",        #Radius of the effective aquifer
    "R",          #Recharge (mm)
    "exfil",      #Exfiltration from WT, ie ET_HMZ (mm)
    "R_lm",       #infiltration to low moisture zone (mm)
    "GW_local",   #Local ground water exchange (mm^3/day)
    "ET_wt",      #ET from the water table (mm/day)
    "GW_bf",      #GW lost to "Baseflow" (mm/day)
    "ET_lm",      #ET in the low moisture zone (mm/day)
    "loss_lm",    #Loss from over saturation
    "y_hm",       #elevation of hm/lm boundary (mm)
    "ds",         #change in volumetric water content
    "s_ex",       #volumetric water content excess of s_fc
    "s_lim",      #volumetric water content with limit of s_fc
    "ET_lm",      #ET in the low moisture zone
    "dh",         #change in head between water table and wetland surface water
    "Ax",         #surface area of GW-Wetland boundary (mm^2)
    "As",         #surface area of lumped wetlands (mm^2)
    "runoff_vol", #saturation excess runoff in mm^3
    "spill_vol"   #spill from lumped wetland module in mm^3
  )
  
  #create matrix for each variable (GO LOOPS)
  n.row<-as.integer(n.years*365)+1
  for(i in 1:length(dynamic_variables)){
    assign(paste0(dynamic_variables[i], ".VAR"),
           matrix(0, ncol=(n.wetlands+2),
                  nrow=n.row,
                  dimnames=list(seq(1,n.row, 1), c("land","wetland",giw.ID))))
  }
  
  ####################################################################################
  # Step 4: Create functions to calculate SW-GW interactions (at daily timestep)------
  ####################################################################################
  #Wetland SW-GW dynamics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  giw.FUN<-function(day, wetland){
    #Define wetland location in .VAR and .INFO matricies matricies~~~~~~~~~~~~
    wet.VAR  <- wetland+2
    wet.INFO <- wetland
    
    #Characterize the Vadose Zone~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Water Table Elevation
    y_wt.VAR[day,wet.VAR] <<- y_w.VAR[day,wet.VAR]
    if(y_wt.VAR[day,wet.VAR]>giw.INFO[wet.INFO,"invert"]){y_wt.VAR[day,wet.VAR]<<-giw.INFO[wet.INFO,"invert"]}
    y_sat.VAR[day,wet.VAR] <<- y_wt.VAR[day,wet.VAR]+abs(giw.INFO[wet.INFO,"psi"])
    
    #High Moisture Zone
    if(y_sat.VAR[day,wet.VAR]>giw.INFO[wet.INFO,"y_c"]){
      #high moisture zone is above the ground
      y_hm.VAR[day,wet.VAR]<<-giw.INFO[wet.INFO,"invert"]
    }else{
      psi_fc <- giw.INFO[wet.INFO,"psi"]*giw.INFO[wet.INFO,"s_fc"]^(-4.05) #see eq. 26 in Laio 2011
      A      <- (psi_fc-giw.INFO[wet.INFO,"psi"]-giw.INFO[wet.INFO,"y_c"])/
                (psi_fc-giw.INFO[wet.INFO,"psi"]-giw.INFO[wet.INFO,"y_c"]-(5*(abs(giw.INFO[wet.INFO,"RD"]))))
      if(y_sat.VAR[day,wet.VAR]>(-5*abs(giw.INFO[wet.INFO,"RD"])+psi_fc-giw.INFO[wet.INFO,"psi"])){
        A    <- abs(A)
        y_hm.VAR[day,wet.VAR] <<- (1-A^(3/4))*(y_sat.VAR[day,wet.VAR]-giw.INFO[wet.INFO,"y_c"]) - (A^2*(1-A^(-1/4))) / 
                                  (-giw.INFO[wet.INFO,"y_c"]+psi_fc-giw.INFO[wet.INFO,"psi"])*(y_sat.VAR[day,wet.VAR]-giw.INFO[wet.INFO,"y_c"])^2
      }else{
        y_hm.VAR[day,wet.VAR] <<- y_sat.VAR[day,wet.VAR]-psi_fc+giw.INFO[wet.INFO,"psi"]
      }
    }
    
    #Estimate Evapotransporation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #ET from the water table
    ET_wt.VAR[day, wet.VAR] <<- ifelse(y_sat.VAR[day, wet.VAR]>giw.INFO[wet.INFO,"y_c"],
                                       pet.VAR[day],
                                       pet.VAR[day]*exp(y_hm.VAR[day,wet.VAR]/abs(giw.INFO[wet.INFO, "RD"]))) #Equation 5 McLaughlin et al 2014
    
    #PET in LMZ
    PET_lm<<-pet.VAR[day]-ET_wt.VAR[day,wet.VAR]
    
    #ET from the low moisture zone
    if(s_lim.VAR[day,wet.VAR]>s_star){
      ET_lm.VAR[day,wet.VAR]<<-PET_lm
    }else{
      if(s_lim.VAR[day,wet.VAR]>giw.INFO[wet.INFO,"s_wilt"]){
        ET_lm.VAR[day,wet.VAR]<<-PET_lm*(s_lim.VAR[day, wet.VAR]-giw.INFO[wet.INFO,"s_wilt"])/(giw.INFO[wet.INFO,"s_fc"]-giw.INFO[wet.INFO,"s_wilt"])
      }else{
        ET_lm.VAR[day,wet.VAR]<<-0
      }
    }
    
    #Calculate Change in Soil Moisture~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Recharge
    R_lm.VAR[day, wet.VAR]<<-ifelse(precip.VAR[day]>giw.INFO[wet.INFO,"n"]*abs(y_hm.VAR[day,wet.VAR])*(1-giw.INFO[wet.INFO,"s_fc"]),
                                    giw.INFO[wet.INFO,"n"]*abs(y_hm.VAR[day,wet.VAR])*(1-giw.INFO[wet.INFO,"s_fc"]),
                                    precip.VAR[day])
    
    #loss
    loss_lm.VAR[day,wet.VAR]<<-giw.INFO[wet.INFO,"n"]*(abs(y_hm.VAR[day,wet.VAR]))*(s_ex.VAR[day,wet.VAR]-s_lim.VAR[day, wet.VAR])
    
    #change in soil moisture
    ds.VAR[day,wet.VAR]<<-ifelse(abs(y_hm.VAR[day,wet.VAR])>0,
                                 (R_lm.VAR[day,wet.VAR]-ET_lm.VAR[day,wet.VAR]-loss_lm.VAR[day,wet.VAR])/(abs(y_hm.VAR[day,wet.VAR])*giw.INFO[wet.INFO, "n"]),
                                 0)
    s_ex.VAR[day+1,wet.VAR]<<-s_ex.VAR[day,wet.VAR]+ds.VAR[day,wet.VAR]
    s_lim.VAR[day+1,wet.VAR]<<-ifelse(s_ex.VAR[day+1,wet.VAR]>giw.INFO[wet.INFO, "s_fc"], giw.INFO[wet.INFO, "s_fc"], s_ex.VAR[day+1,wet.VAR])
    
    
    #Characterize the change in Volume~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Recharge
    vadose_storage<-ifelse(y_hm.VAR[day, wet.VAR]<giw.INFO[wet.INFO,"y_c"],
                           giw.INFO[wet.INFO, "n"]*abs(giw.INFO[wet.INFO,"invert"]-y_hm.VAR[day, wet.VAR])*(giw.INFO[wet.INFO,"s_fc"]-s_lim.VAR[day,wet.VAR]),
                           0)
    R.VAR[day,wet.VAR]<<-ifelse(precip.VAR[day]>vadose_storage,
                                precip.VAR[day]-vadose_storage,
                                0)
    
    #Exchange with the lumped module
    if((y_w.VAR[day, wet.VAR]-y_wt.VAR[day,1])==0 | day==1){
      GW_local.VAR[day,wet.VAR]<<- 0
    }else{
      r_w   <- (stage2area_giw.fun(vol2stage_giw.fun(V_w.VAR[day, wet.VAR]))/pi)^0.5
      r_w   <- ifelse(r_w>0, r_w, (giw.INFO[,"area_wetland"]/pi)^0.5)
      r_ws  <- giw.INFO[wet.INFO, "dLe"] + r_w
      y_w   <- y_w.VAR[day, wet.VAR] + giw.INFO[wet.INFO, "dz"]
      GW_local.VAR[day,wet.VAR]<<- pi*giw.INFO[wet.INFO,"k_sat"] * 
        ((y_wt.VAR[day, 1]- land.INFO[,"y_cl"])^2-(y_w- land.INFO[,"y_cl"])^2) / 
        log(r_ws/r_w) 
    }
    
    #Adjust for differences in contribruting watershed area
    runoff_vol.VAR[day,wet.VAR]<<- runoff_vol.VAR[day,"land"]*giw.INFO[,"vol_ratio"]
    
    #Calculate the change in volume of wetland water table
    dV_w.VAR[day, wet.VAR]<<- R.VAR[day,wet.VAR]*giw.INFO[wet.INFO,"area_wetland"]-
                              ET_wt.VAR[day, wet.VAR]*giw.INFO[wet.INFO,"area_wetland"]+
                              GW_local.VAR[day,wet.VAR]+
                              runoff_vol.VAR[day,wet.VAR]+
                              spill_vol.VAR[day,"wetland"]*giw.INFO[wet.INFO,"vol_ratio"]
    
    #Characterize the change in water level
    #Calculate the next days volume
    V_w.VAR[day+1, wet.VAR]<<-ifelse((V_w.VAR[day, wet.VAR]+dV_w.VAR[day, wet.VAR])<0,
                                     0,
                                     V_w.VAR[day, wet.VAR]+dV_w.VAR[day, wet.VAR])
    #Calculate Spill Volume
    if(V_w.VAR[day+1, wet.VAR]>max(vol_giw.INFO[,2],na.rm=T)){
      spill_vol.VAR[day+1, wet.VAR]<<-V_w.VAR[day+1, wet.VAR]-max(vol_giw.INFO[,2],na.rm=T)
      V_w.VAR[day+1, wet.VAR]<<-max(vol_giw.INFO[,2], na.rm=T)
    }
    
    #Convert to water level
    y_w.VAR[day+1, wet.VAR]<<-vol2stage_giw.fun(V_w.VAR[day+1, wet.VAR])
    
    #The end for now
  }
  
  #Landscape GW dynamics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  upland.FUN<-function(day){
    # =========================================================================
    # 5. Calculate change in wetland water depth
    # =========================================================================
    #Estimate wetland surface area
    As.VAR[day,"wetland"]<<-stage2area.fun(y_w.VAR[day,"wetland"])
    
    #Calculate GW_local (mm^3, assume water flowing out of the wetland is +)
    GW_local_mat <- pi*land.INFO[,"k_sat"] * 
      ((y_sat.VAR[day, "land"]- land.INFO[,"y_cl"])^2-(y_w.VAR[day, "wetland"]- land.INFO[,"y_cl"])^2) /
      log((lumped.INFO[,'dLe'] + lumped.INFO[, 'r_w'])/lumped.INFO[,"r_w"])
    GW_local.VAR[day, "wetland"] <<- sum(GW_local_mat) 
    
     #change in wetland storage (mm^3)
    dV_w.VAR[day,"wetland"]<<-precip.VAR[day]  * land.INFO[,"wetland_area"]+
                              snowmelt.VAR[day]* land.INFO[,"area"]-
                              pet.VAR[day]     * As.VAR[day,"wetland"]+
                              GW_local.VAR[day,"wetland"]+
                              ((1-giw.INFO[,"vol_ratio"])*runoff_vol.VAR[day,"land"])
    V_w.VAR[day+1,"wetland"]<<-V_w.VAR[day,"wetland"]+dV_w.VAR[day,"wetland"]
    
    #If Wetland is dry
    if(V_w.VAR[day+1,"wetland"]<=0){
      V_w.VAR[day+1,"wetland"]<<- 0
      GW_local.VAR[day, "wetland"] <<- 0
    }
    
    #Calculate Surface water spillage (if there is any)
    if(V_w.VAR[day+1,"wetland"]>land.INFO[,"volume_max"]){
      spill_vol.VAR[day+1,"wetland"]<<-V_w.VAR[day+1,"wetland"]-land.INFO[,"volume_max"]
      V_w.VAR[day+1,"wetland"]<<-land.INFO[,"volume_max"]
    }
    
    #Convert Volume to stage
    y_w.VAR[day+1,"wetland"]<<-volume2stage.fun(V_w.VAR[day+1,"wetland"])
    
    # =========================================================================
    # 6. Calculate change in upland water table elevation
    # =========================================================================
    #Characterize the Vadose Zone ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A_up <-(psi_fc - land.INFO[,"psi"] - land.INFO[,"y_c"])/
           (psi_fc - land.INFO[,"psi"] - land.INFO[,"y_c"]-(5*abs(land.INFO[,"RD"])))


    if(y_sat.VAR[day,"land"] < land.INFO[,"y_c"]){ # Case 1: Deep Water Table (LMZ Forms)
      
      # Determine location of LMZ/HMZ boundary ----------------------------------
      if(y_sat.VAR[day,"land"] < (-5*abs(land.INFO[,"RD"]) + psi_fc - land.INFO[,"psi"])){ # Check if very deep water table (cond 3 in eqn 26 of Laio 2009)
        y_hmz.VAR[day,"land"] <<- y_sat.VAR[day,"land"] - psi_fc + land.INFO[,"psi"]
      }else{ # Otherwise Cond 2 in eqn 26 of Laio et al 2009 ---------------------
        y_hmz.VAR[day,"land"] <<- (1-A_up^0.75) * (y_sat.VAR[day,"land"] - land.INFO[,"y_c"]) - A_up^2 * (1-A_up^-0.25) /
                                  (-land.INFO[,"y_c"] + psi_fc - giw.INFO[wet.INFO,"psi"]) *
                                  (y_sat.VAR[day,"land"] - land.INFO[,"y_c"])^2
      } # end of if statement for LMZ/HMZ Boundary
      
      
      # 6a. ET and Recharge dynamics-----------------------------------------------
      # 6a.1 high moisture zone dynamics-------------------------------------------
      ET_wt.VAR[day,"land"]<<- pet.VAR[day]* exp(y_sat.VAR[day,"land"]/abs(land.INFO[,"RD"]))
      exfil.VAR[day,"land"]<<- pet.VAR[day]*(exp(y_hmz.VAR[day,"land"]/abs(land.INFO[,"RD"])) -
                                             exp(y_sat.VAR[day,"land"]/abs(land.INFO[,"RD"])))
      R.VAR[day,"land"]    <<- 0
      
    
      # 6a.2 low moisture zone dynamics -------------------------------------------
      vadose_storage <- land.INFO[,"n"]*abs(y_hmz.VAR[day,"land"])*(1-s_ex.VAR[day,"land"])
      PET_lm         <- max(0, pet.VAR[day] - ET_wt.VAR[day,"land"] - exfil.VAR[day,"land"])  
      
      if(s_ex.VAR[day,"land"] > s_star){  # If LMZ greater than _star 
        ET_lm.VAR[day,"land"] <<- PET_lm
      }else if (s_ex.VAR[day,"land"] > land.INFO[,"s_wilt"]){
        ET_lm.VAR[day,"land"] <<- PET_lm*(s_ex.VAR[day,"land"] - land.INFO[,"s_wilt"])/
                                         (s_star - land.INFO[,"s_wilt"])
      }else{
        ET_lm.VAR[day,"land"] <<- 0
      }
      
      R_lm.VAR[day,'land'] <<- precip.VAR[day]
      
      beta <- 2*land.INFO[,"b"]+4
      
      if (s_ex.VAR[day,"land"] >= land.INFO[,"s_fc"]){
        loss_lm.VAR[day,"land"] <<- land.INFO[,"n"]*(s_ex.VAR[day,"land"]-s_lim.VAR[day,"land"])*abs(y_hm.VAR[day,"land"])
      } else {
        loss_lm.VAR[day,"land"] <<- 0
      }
  
      
      # Calculate soil moisture balance ----------------------------------------------
      ds.VAR[day,"land"]<<-ifelse( y_hmz.VAR[day,"land"] < 0,
                                  (R_lm.VAR[day,"land"] - ET_lm.VAR[day,"land"] - loss_lm.VAR[day,"land"]) /
                                   abs(y_hmz.VAR[day,"land"]) / land.INFO[,"n"],
                                  0)
      s_ex.VAR[day+1, "land"]<<-max(0, s_ex.VAR[day,"land"]+ds.VAR[day,"land"])
      s_lim.VAR[day+1,"land"]<<-min(s_ex.VAR[day+1, "land"] , land.INFO[,"s_fc"])
      
    } else {# Case 2: Shallow Water Table (No LMZ) -------------------------------------
      y_hmz.VAR[day,"land"]   <<- 0
      R_lm.VAR[day,"land"]    <<- 0
      loss_lm.VAR[day,"land"] <<- 0
      
      R.VAR[day,"land"]    <<- precip.VAR[day]
      
      ET_wt.VAR[day,"land"]<<- pet.VAR[day]* exp(y_sat.VAR[day,"land"]/abs(land.INFO[,"RD"]))
      
      exfil.VAR[day,"land"]<<- ifelse(y_sat.VAR[day,"land"] == 0, 0,pet.VAR[day]*(1 - exp(y_sat.VAR[day,"land"]/abs(land.INFO[,"RD"])))  )
      
      
    } # End if statement for ET and Recharge calculations
    
    # OLD OLD OLD 

    #Estimate Change in Water Table Depth~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Local GW Flux (Assume water flowing out of the wetlands is positive)
    GW_local.VAR[day,"land"]<<--1*GW_local.VAR[day,"wetland"]+GW_local.VAR[day,3]
    GW_local<-GW_local.VAR[day,"land"]/(land.INFO[,"area"]-sum(giw.INFO[,"area_wetland"]))
    
    # print(GW_local.VAR[day,"wetland"])
    # print(GW_local.VAR[day,])
    
    #Flux out of the watershed (e.g. baseflow from SWAT manual)
    GW_bf.VAR[day,"land"]<<-ifelse(day==1,
                                   land.INFO[,"GW_bf_0"],
                                   ifelse(GW_bf.VAR[day-1,"land"]*exp(-land.INFO[,"kb"])+(R.VAR[day,"land"]+loss_lm.VAR[day,"land"]-ET_wt.VAR[day,"land"])*(1-exp(-land.INFO[,"kb"]))<0,
                                          (GW_bf.VAR[day-1,"land"]*exp(-land.INFO[,"kb"])+(R.VAR[day,"land"]+loss_lm.VAR[day,"land"]-ET_wt.VAR[day,"land"])*(1-exp(-land.INFO[,"kb"]))),
                                          0))
    
    #Water Table Depth
    y_sat.VAR[day+1,"land"]<<- y_sat.VAR[day,"land"] + (1/land.INFO[,"Sy"])*(R.VAR[day,"land"]+
                                                                              loss_lm.VAR[day,"land"]-
                                                                              exfil.VAR[day,"land"]-
                                                                              ET_wt.VAR[day,"land"]+ 
                                                                              GW_bf.VAR[day,"land"]-
                                                                              GW_local) #change in water table elevation
    
    y_wt.VAR[day+1,"land"] <<- y_sat.VAR[day+1,"land"] + land.INFO[,"psi"]
    
    if(y_sat.VAR[day+1,"land"] > 0){ # runoff calculations
      runoff_vol.VAR[day+1,"land"] <<- y_sat.VAR[day+1,"land"] *(land.INFO[,"area"]-land.INFO[,"wetland_area"]) * land.INFO[,"Sy"]
      y_sat.VAR[day+1,"land"]      <<- 0
      
      if(y_wt.VAR[day+1,"land"] > 0){
        y_wt.VAR[day+1,"land"] <<- 0
      }
    }else{
      if((y_wt.VAR[day,"land"]+dy_wt.VAR[day,"land"])<land.INFO[,"y_cl"]){
        #Calculate deficit
        #divide deficit by the number of fluxes
        #substract modified deficit from all fluxes 
        
        
        
        y_wt.VAR[day+1,"land"]<<-land.INFO[,"y_cl"]
        
        
        
      }else{
        # y_wt.VAR[day+1,"land"]<<-y_wt.VAR[day,"land"]+dy_wt.VAR[day,"land"]
      }
    } #END OF IF STATEMENT FOR Y_WT
  }
  
  ####################################################################################
  # Step 5: Model wetland/upland hydrology--------------------------------------------
  ####################################################################################
  #Run model with wetlands~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #set random seed
  set.seed(1)
  
  #set initial conditions
  V_w.VAR[1,"wetland"]<-land.INFO[,"volume_max"]
  land.INFO[,"GW_bf_0"]<-10
  y_wt.VAR[1,]<-land.INFO[,"y_wt_0"]
  s_ex.VAR[1,]<-land.INFO[,"s_t_0"]
  s_lim.VAR[1,]<-land.INFO[,"s_t_0"]
  y_w.VAR[1,]<-land.INFO[,"y_wt_0"]
  
  #run this bad boy
  for(i in 1:((365*n.years)-1)){
    # print(i)
    giw.FUN(day=i,wetland = 1)
    upland.FUN(i)
  }
  
  ####################################################################################
  # Step 6: Summarize Results---------------------------------------------------------
  ####################################################################################
  #Create list to pass to outter function
  output<-c(ls(pattern=".VAR"), ls(pattern=".INFO"))
  output<-mget(output)
  
  #Close Function
  output
}