#To do list
#-add to wetland.shp
#-add to lump.INFO
#-remove "for testings" WetID 

#Setup workspace
rm(list=ls(all=TRUE))
load("backup/Inputs.Rdata")  




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

save.image("backup/radius.RData")
