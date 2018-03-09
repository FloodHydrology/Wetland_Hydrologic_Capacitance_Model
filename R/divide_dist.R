#To do list
#-add to wetland.shp
#-add to lump.INFO
#-remove "for testings" WetID 

#Setup workspace
rm(list=ls(all=TRUE))
load("Inputs.Rdata")   

#Create function to cacluate distance from wetland edge to "groundwater inflection" point
divide_dist.fun<-function(WetID){

  #wetland of interest centoird
  ind_cent<-data.frame(centroid(wetlands.shp[wetlands.shp$WetID==WetID,]))
  colnames(ind_cent)<-c("x","y")
  coordinates(ind_cent) <- ~x + y 
  ind_cent<-SpatialPoints(ind_cent)
  
  #All centroids
  wet_cent<-data.frame(centroid(wetlands.shp))
  colnames(wet_cent)<-c("x","y")
  coordinates(wet_cent) <- ~x + y 
  wet_cent<-SpatialPoints(wet_cent)

  #Calculate distances between wetlands
  wetlands.shp@data$dist2wet<-c(gDistance(ind_cent, wet_cent, byid=T))
  wetlands.shp<-wetlands.shp[rank(wetlands.shp$dist2wet)<3,]

  #Create line between centroid
  flowpath <- data.frame(centroid(wetlands.shp))
  colnames(flowpath)<-c("x","y")
  coordinates(flowpath) <- ~x + y 
  flowpath.shp<-SpatialPoints(flowpath)
  flowpath.shp<-SpatialLines(list(Lines(Line(flowpath.shp), ID="a")))
  
  #Clip based on area wetlands
  centroid.shp<-gDifference(flowpath.shp, wetlands.shp)
  
  #Define "delta' distance
  gLength(flowpath.shp)/2
}
wetlands.shp$dLe<-0
for(i in 1:length(wetlands.shp$WetID)){
  print(i)
  wetlands.shp$dLe[i]<-divide_dist.fun(wetlands.shp$WetID[i])
}
