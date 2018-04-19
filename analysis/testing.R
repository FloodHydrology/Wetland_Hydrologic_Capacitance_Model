#Testin sandbox

#Section 1: Save workspace ------------------
remove(list=ls()[ls()!= 'giw.INFO'    &
                 ls()!= 'land.INFO'   &
                 ls()!= 'lumped.INFO' &
                 ls()!= 'precip.VAR'  &
                 ls()!= 'pet.VAR'     &
                 ls()!= 'n.years'     &
                 ls()!= 'area'        &
                 ls()!= 'volume'      &
                 ls()!= 'giw.ID'])
save.image("backup/test.RData")                  

#Run model -----------------------------------
#Setup workspace
remove(list=ls())
load("backup/test.RData")
source("~/Wetland_Hydrologic_Capacitance_Model/R/WHC_2.R") 

#Run model
output<-wetland.hydrology(giw.INFO,land.INFO, lumped.INFO, precip.VAR, pet.VAR, 10, area, volume, giw.ID)
attach(output)

#Plot Individual Wetland Module Fluxes------------
giw.INFO<-matrix(giw.INFO[giw.ID,],nrow=1,  dimnames = list(c(1), colnames(giw.INFO)))
giw<-data.frame(timestep  = seq(1,  3651, 1),
                y_w       = y_w.VAR[,3],
                recharge  = R.VAR[,3],
                ET        = ET_wt.VAR[,3], 
                GW_local  = GW_local.VAR[,3]/giw.INFO[1,"area_wetland"],
                runoff    = runoff_vol.VAR[,3]/giw.INFO[1,"area_wetland"],
                spill     = spill_vol.VAR[,3]*giw.INFO[1,"vol_ratio"]
)


#For GW Local
plot(giw$y_w+giw.INFO[1,"dz"], col="blue", type="l", ylim=c(-400,0))
#points(y_w.VAR[1000:1300,"wetland"], col="dark green", type="l", lty=2)
points(y_wt.VAR[,1], col="brown", type="l")
abline(h=giw.INFO[,"invert"])


plot((giw$y_w+giw.INFO[1,"dz"])-y_wt.VAR[,1], type="l")
abline(h=0)

plot(giw$GW_local, type="l")



plot(giw$timestep, giw$y_w, type= "l")
plot(giw$timestep, giw$recharge, type="l")
plot(giw$timestep, giw$GW_local, type="l")
plot(giw$timestep, giw$runoff)
