#To do list
#-add to wetland.shp
#-add to lump.INFO
#-remove "for testings" WetID 

#Setup workspace
rm(list=ls(all=TRUE))
load("Inputs.Rdata")   

#Create function to cacluate distance from wetland edge to "groundwater inflection" point
divide_dist.fun<-function(cat){

  #Identify catchment of interest
  cat.shp<-catchments.shp[cat,]
  
  #Clip relevant wetlands
  wetlands.shp<-wetlands.shp[cat.shp,]
  
  #Create inner function to define median distance for each wetland
  fun<-function(WetID){
  
    #wetland of interest centroid
    ind_cent<-data.frame(centroid(wetlands.shp[wetlands.shp$WetID==WetID,]))
    colnames(ind_cent)<-c("x","y")
    ind_pnt<-ind_cent
    coordinates(ind_cent) <- ~x + y 
    ind_cent<-SpatialPoints(ind_cent)
    
    #All centroids
    wet_cent<-data.frame(centroid(wetlands.shp))
    colnames(wet_cent)<-c("x","y")
    coordinates(wet_cent) <- ~x + y 
    wet_cent<-SpatialPoints(wet_cent)
  
    #Calculate distances between wetlands
    wetlands.shp@data$dist2wet<-c(gDistance(ind_cent, wet_cent, byid=T))
    wetlands.shp<-wetlands.shp[rank(wetlands.shp$dist2wet)>1,]
  
    #Calculate the dLe for each wetland
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
    
    #Export dLe
    median(wetlands.shp$dLe, na.rm=T)
 }
  
  #Run function (if there are enough wetlands)
  if(length(wetlands.shp)>2){
    
    #run inner function for each wetland
    output<-sapply(wetlands.shp$WetID, fun)

    #Export median of output
    median(output, na.rm=T)
  }else{
    0
  }
}

#run bp_analysis in parralel
library(parallel)
n.cores<-detectCores()
cl <- makePSOCKcluster(n.cores) #Make 4 processors
clusterEvalQ(cl, library(sp))  
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(rgeos))
clusterEvalQ(cl, library(rgdal))
clusterEvalQ(cl, library(geosphere))
clusterExport(cl, c('catchments.shp',"wetlands.shp"), env=.GlobalEnv)  #Send Clusters function with the execute function
t0<-Sys.time()
a<-seq(1,length(catchments.shp))
x<-parLapply(cl, a, divide_dist.fun) #Run execute Function
tf<-Sys.time()
stopCluster(cl)  #Turn clusters off
tf-t0

save.image("radius.RData")
