###################################################################################
#Name: Topographic Analysis
#Coder: C. Nathan Jones
#Date: 8/5/2016
#Purpose: Provide Function for Topographic analysis 
##################################################################################

DEM_Processing.fun<-function(dem.grd,pp.shp, 
                             python.path="C:\\Python27\\ArcGIS10.4\\", 
                             scratchspace="C:\\ScratchWorkspace"){
  
  ####################################################################################
  #Step 1: Setup Workspace------------------------------------------------------------
  ####################################################################################
  #load appropriate libarary
  library(RPyGeo)
  library(raster)
  library(maptools)
  library(rgeos)
  library(rgdal)
  library(dplyr)
  library(magrittr)

  ####################################################################################
  #Step 2: Identify "sinks" and basins (ArcPyGeo)-------------------------------------
  ####################################################################################
  #Setup Rpygeo environment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create python environment
  py.env<-rpygeo.build.env(python.path = python.path, 
                           workspace = scratchspace,
                           overwriteoutput = 1) 
  
  #Set R directoy to python scratchworkspace
  setwd(scratchspace)
  
  #Clear Scratch Space
  file.remove(list.files())

  #Find Sinks and associated watersheds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Filter the DEM
  dem_filter.grd<- focal(dem.grd, w=focalWeight(dem.grd, 5, "Gauss"))
  dem_filter.grd<- focal(dem_filter.grd, w=focalWeight(dem_filter.grd, 5, "Gauss"))
  dem_filter.grd<- focal(dem_filter.grd, w=focalWeight(dem_filter.grd, 5, "Gauss"))
  dem_filter.grd@crs<-dem.grd@crs
  writeRaster(dem_filter.grd, file="dem_filter.asc", overwrite=T)
  
  #fill sinks 
  rpygeo.geoprocessor(fun="Fill_sa",
                      c("dem_filter.asc", "dem_fill"),
                      env=py.env
  )
  
  #Compute Flow Direction
  rpygeo.FlowDirection.sa("dem_fill","fdr_esri", env=py.env)
  
  #Compute Flow Accumulation
  rpygeo.FlowAccumulation.sa("fdr_esri","fac_esri", env=py.env)
  fac.grd<-raster("fac_esri")
  
  #Define pour point
  pp_buffer.shp<-gBuffer(pp.shp, width=30, byid=T)
  snap_pp.grd<-mask(fac.grd, pp_buffer.shp)
  snap_pp.grd[snap_pp.grd!=cellStats(snap_pp.grd, max)]<-NA
  snap_pp.grd<-snap_pp.grd*0+1
  writeRaster(snap_pp.grd, file="snap.asc", overwrite=T)
  
  #Delineate Watershed
  rpygeo.geoprocessor("Watershed_sa",
                      list("fdr_esri","snap.asc", "Value"),
                      env=py.env)
  rpygeo.geoprocessor(fun="RasterToPolygon_conversion",
                      c("value", "watershed.shp", "NO_SIMPLIFY","VALUE"),
                      env=py.env)
  watershed.shp<-readOGR(".","watershed")
  
  #Extract DEM by mask
  dem_mask.grd<-mask(dem.grd, gBuffer(watershed.shp,width=15))
  dem_save.grd<-dem_mask.grd
  dem_mask.grd<-crop(dem_mask.grd, gBuffer(watershed.shp,width=15))
  dem_mask.grd<-focal(dem_mask.grd, w=focalWeight(dem_mask.grd, 5, "Gauss"))
  writeRaster(dem_mask.grd,"dem_mask.asc", overwrite=T)
  
  
  #fill sinks 
  rpygeo.geoprocessor(fun="Fill_sa",
                      c("dem_mask.asc", "dem_mask",0.1),
                      env=py.env
  )
  
  #Redo flow direction
  rpygeo.FlowDirection.sa("dem_mask","fdr_esri", env=py.env)
  
  #Identify Sinks in watershed
  rpygeo.Sink.sa("fdr_esri", "sink", env=py.env)
  
  #run basin tool
  rpygeo.geoprocessor(fun="Basin_sa", 
                      c("fdr_esri","basin"),
                      env=py.env)
  
  #Convert sink and basin to polygons
  rpygeo.geoprocessor(fun="RasterToPoint_conversion",
                      c("sink", "sink.shp", "VALUE"),
                      env=py.env)
  rpygeo.geoprocessor(fun="RasterToPolygon_conversion",
                      c("basin", "basin.shp", "NO_SIMPLIFY","VALUE"),
                      env=py.env)
  
 
  #Bring everything into the R environment
  sink.shp<-readOGR(".","sink")
  basin.shp<-readOGR(".","basin")
  
  ####################################################################################
  #Step 3: Estimate Spill elevations (Raster)-----------------------------------------
  ####################################################################################
  #rename basin id's
  basin.shp@data$ID<-seq(1:length(basin.shp))
  
  #Create Inundate Function
  inundate.fun<-function(basin.id){
    #select basin
    temp.shp<-basin.shp[basin.id,]
    
    #mask dem
    temp.grd<-mask(dem.grd, temp.shp)
    
    #Create Minimum Raster
    temp_min.grd<-mask(dem.grd, temp.shp)
    temp_min.grd<-temp.grd*0+minValue(temp_min.grd)
    
    #Create function to return conditional raster 
    Con<-function(condition, trueValue, falseValue){
      return(condition * trueValue + (!condition)*falseValue)
    }
    
    #Create dataframe to house information
    df<-data.frame(matrix(0, ncol=4, nrow=100))
    colnames(df)<-c("ele", "area","volume", "outflow")
    df$ele<-seq(0.05,5,5/100)
    
    #inundate
    outflow<-0
    i<-1
    while(outflow==0){
      
      #Inundation Elevation (relative)
      z<-0.05*i
      
      #Inundation Area Calculation
      area<-Con(temp.grd>(temp_min.grd+z),0,1)
      df$area[i]<-cellStats(area, 'sum')*res(area)[1]*res(area)[2]
      volume<-(((z+temp.grd)-temp.grd)*area)*res(area)[1]*res(area)[2]
      df$volume[i]<-cellStats(volume, 'sum')
      
      #Determine if the wetland is "spilling" (e.g. outer cell is saturated)
      outflow<-cellStats(area*boundaries(temp_min.grd, type="inner"), 'sum')
      
      #Next step (terminate if more than 100 steps)
      i<-i+1
      if(i>99){outlflow<-9999}
    }
    
    #Remove zeros
    df$area[df$area==0]<-max(df$area, na.rm=T)
    df$volume[df$volume==0]<-max(df$volume, na.rm=T)
    
    #Add to raster
    extent(area)<-inundate.grd@extent
    area[is.na(area)]<-0
    inundate.grd<-inundate.grd+area
    inundate.grd[is.na(inundate.grd)]<-0
    assign('inundate.grd',inundate.grd,env=.GlobalEnv)
    
    #Export Volume
    c(df$area, df$volume)
  }
  
  #Create raster to populate inundation
  inundate.grd<-dem.grd[[1]]*0
  
  #apply inundate function to basin
  df<-sapply(seq(1:length(basin.shp)), inundate.fun)
  df<-data.frame(matrix(unlist(df), ncol=length(basin.shp), byrow=F))
  
  #Split output into area vs volume data.frames
  area<-df[1:100,]
  volume<-df[101:200,]
  
  ####################################################################################
  #Step 4: Export Data----------------------------------------------------------------
  ####################################################################################
  #Assign data to Global Env
  assign('watershed.shp',watershed.shp, env=.GlobalEnv)
  assign('sink.shp',sink.shp, env=.GlobalEnv)
  assign('basin.shp',basin.shp, env=.GlobalEnv)
  assign('area',area, env=.GlobalEnv)
  assign('volume',volume, env=.GlobalEnv)
  assign('dem.grd',dem_save.grd, env=.GlobalEnv)
}
