dz.fun<-function(WetID){

  #Isolate wetland of interest
  wetland.shp<-wetlands.shp[wetlands.shp$WetID==WetID,]
  
  #Define wetland centroid
  wet_cent<-data.frame(centroid(wetland.shp))
  colnames(wet_cent)<-c("x","y")
  wet_pnt<-wet_cent
  coordinates(wet_cent) <- ~x + y 
  wet_cent<-SpatialPoints(wet_cent)
  
  #Define wetland and watershed radius 
  r_wetland<-(gArea(wetland.shp)/pi)^0.5
  r_watershed<-r_wetland + wetland.shp$dLe
  if(is.na(r_watershed)){
    r_watershed<-median(wetlands.shp$dLe, na.rm=T)
  }
  
  #define upland dem 
  upland_mask<-gBuffer(wet_cent, width=r_watershed)
  upland_dem<-crop(dem.grd, upland_mask)
  upland_dem<-mask(upland_dem, upland_mask)
  
  #Define wetland dem
  wetland_dem<-crop(dem.grd, wetland.shp)

  #Estimate difference (mm)
  dz<- cellStats(wetland_dem, median) - cellStats(upland_dem, median)  
     #If the median wetland elevation is 8, and the median upland elevation is 10, dz should be -2
  #Export 
  c(WetID, dz)
}
