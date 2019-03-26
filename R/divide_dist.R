divide_dist.fun<-function(cat){
  
  #a. Identify catchment of interest
  cat.shp<-catchments.shp[cat,]
  
  #b. Break catchment boundary into 100 equally spaced points
  cat.grd<-raster(cat.shp)
  extent(cat.grd)<-extent(cat.shp)
  cellsize<-2*pi*((gArea(cat.shp)/pi)^.5)/100 #paremeter/1000
  res(cat.grd)<-c(cellsize,cellsize)
  cat.grd<-rasterize(cat.shp, cat.grd, 1)
  cat.grd<-boundaries(cat.grd, type="inner")
  cat.grd[cat.grd==0]<-NA
  cat.pnt<-rasterToPoints(cat.grd)
  if(nrow(cat.pnt)>0){cat.pnt<-SpatialPoints(cat.pnt)}
  
  #c. Clip relevant wetlands
  wetlands.shp<-wetlands.shp[cat.shp,]
  
  #d. Create inner function to define median distance for each wetland
  fun<-function(WetID){
    
    #i. wetland of interest centroid
    ind_cent<-data.frame(centroid(wetlands.shp[wetlands.shp$WetID==WetID,]))
    colnames(ind_cent)<-c("x","y")
    ind_pnt<-ind_cent
    coordinates(ind_cent) <- ~x + y 
    ind_cent<-SpatialPoints(ind_cent)
    
    #ii. All centroids
    wet_cent<-data.frame(centroid(wetlands.shp))
    colnames(wet_cent)<-c("x","y")
    coordinates(wet_cent) <- ~x + y 
    wet_cent<-SpatialPoints(wet_cent)
    
    #iii. Calculate distances between wetlands
    wetlands.shp@data$dist2wet<-c(gDistance(ind_cent, wet_cent, byid=T))
    wetlands.shp<-wetlands.shp[rank(wetlands.shp$dist2wet)>1,]
    
    #iii. Calculate the dLe for each wetland
    wetlands.shp$dLe<-0
    for(i in 1:length(wetlands.shp$dist2wet)){
      #Create line between centroid
      flowpath <- data.frame(centroid(wetlands.shp[i,]))
      colnames(flowpath)<-c("x","y")
      flowpath<-rbind(flowpath, ind_pnt)
      coordinates(flowpath) <- ~x + y 
      flowpath.shp<-SpatialPoints(flowpath)
      flowpath.shp<-SpatialLines(list(Lines(Line(flowpath.shp), ID="a")))
      
      #Clip based on area wetlands
      flowpath.shp<-gDifference(flowpath.shp, wetlands.shp)
      
      #Calculate length
      if(length(flowpath.shp)==0){
        dLe<-0
      }else{
        dLe<-gLength(flowpath.shp)/2
      }
      
      #Define "delta' distance
      wetlands.shp$dLe[i]<-dLe
    }
    
    #iv. Deterimine length to polygon boundary
    watershed_boundary<-min(gDistance(ind_cent,cat.pnt, byid=T))/2
    
    #v. Export dLe
    c(WetID,median(c(wetlands.shp$dLe, watershed_boundary), na.rm=T))
  }
  
  #e. Run function (if there are enough wetlands)
  if(length(wetlands.shp)>1){
    
    #run inner function for each wetland
    output<-lapply(wetlands.shp$WetID, fun)
    output<-do.call(rbind, output)
    
    #Export output
    output
  }else{
    c(0,0)
  }
}