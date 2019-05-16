#Run these bad boys 
source("R/WHC_4.R")
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
attach(output)

#Plot y_wt water balance components-----------------------------------------------------------------
par(mfrow=c(7,1)) 
par(mar=c(2,5,1,1))
start<-2
stop<-2000
plot(y_wt.VAR[start:stop,"land"], type="l", ylab="Water Table Elevation")
plot(R.VAR[start:stop,"land"], type="l", ylab="Direct Recharge")
plot(loss_lm.VAR[start:stop,"land"], type="l", ylab="Recharge from LMZ")
plot(exfil.VAR[start:stop,"land"], type="l", ylab="Exfiltration")
plot(ET_wt.VAR[start:stop,"land"], type="l", ylab="ET from Water Table")
plot(GW_bf.VAR[start:stop,"land"], type="l", ylab="Baseflow")
plot(GW_local.VAR[start:stop,"land"]/(land.INFO[,"area"]-sum(giw.INFO[,"area_wetland"])), type="l", ylab="SW-GW Exchange")

#Plot s water balance components--------------------------------------------------------------------
dev.off()
start<-1000
stop<-2000
par(mfrow=c(2,1)) 
par(mar=c(2,5,1,1))
plot(y_hmz.VAR[start:stop,"land"], type="l")
plot(s_ex.VAR[start:stop,"land"], type="l")
abline(h=land.INFO[,"s_fc"])
abline(h=land.INFO[,"s_wilt"])

#Inviestigate loss term
beta <- 2*land.INFO[,"b"]+4
if(s_ex.VAR[day,"land"] >= land.INFO[,"s_fc"]){
  loss_lm.VAR[day,"land"] <- land.INFO[,"n"]*(s_ex.VAR[day,"land"]-s_lim.VAR[day,"land"])*abs(y_hm.VAR[day,"land"])
  s_ex.VAR[day,"land"] <- land.INFO[,"s_fc"]
} 