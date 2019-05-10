#Plot 
load("backup_WHC_3e.RData")
plot(output$y_wt.VAR[(365*10):(365*15),1], type="l", ylab="Water Tabel Elevation [mm]", xlab="Day")
for(i in 1:10){
  abline(v=365*i)
}

plot(output$y_w.VAR[(365*10):(365*11),2], type="l", ylab="Lumped Wetland Water Level [mm]", xlab="Day")
for(i in 1:10){
  abline(v=365*i)
}

#3f
load("backup_WHC_3f.RData")
plot(output$y_wt.VAR[(365*10):(365*11),1], type="l", ylab="Water Tabel Elevation [mm]", xlab="Day")
plot(output$y_w.VAR[(365*3):(365*5),2], type="l", ylab="Lumped Wetland Water Level [mm]", xlab="Day")
