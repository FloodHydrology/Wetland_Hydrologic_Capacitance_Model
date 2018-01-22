###################################################################################
#Name: Topographic Analysis
#Coder: C. Nathan Jones
#Date: 8/5/2016
#Purpose: Provide Function for Topographic analysis 
##################################################################################

DEM_Processing.fun<-function(dem,pp,wd, 
                             saga.path="C:\\Program Files/saga-6.2.0_x64", 
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
  
  #Define Variables
  dem.grd<-dem
  pp.shp<-pp
  
  ####################################################################################
  #Step 2: Watershed Delineation (RSAGA)----------------------------------------------
  ####################################################################################
  #SetupSaga Environment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Filter dem
  focal<-focalWeight(x=dem.grd, type="Gauss")
  dem_filter.grd<- focal(dem.grd, w=focalWeight(dem.grd, 5, "Gauss"))
  dem_filter.grd<- focal(dem_filter.grd, w=focalWeight(dem_filter.grd, 5, "Gauss"))
  dem_filter.grd<- focal(dem_filter.grd, w=focalWeight(dem_filter.grd, 5, "Gauss"))
  dem_filter.grd@crs<-dem.grd@crs
  
  #Create Watershed Delineation function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  delineate.fun<-function(outlet.shp,dem_filter.grd,fdr.grd, saga.path, scratchspace){
    #define working directory
    setwd(scratchspace)
    
    #define saga working environment
    saga<-rsaga.env(path=saga.path)
    
    #remove na from dem
    dem.grd[is.na(dem.grd)]<-1
    
    #convert dem to sgrid
    writeRaster(dem.grd, file="dem.sdat", overwrite=T)
    
    #fill dem
    rsaga.fill.sinks('dem.sgrd','dem_fill.sgrd', minslope=0.0001,env=saga)
    
    #Create catchment area grid
    rsaga.topdown.processing('dem_fill.sgrd', 
                             out.carea = 'dem_ca.sgrd', 
                             out.flowpath = 'flowpath.sgrd',
                             method=3,
                             env=saga)
    rsaga.sgrd.to.esri('dem_ca.sgrd','dem_ca.asc', env=saga)
    dem_ca.grd<-raster('dem_ca.asc')
    rsaga.sgrd.to.esri('flowpath.sgrd','fdr.asc', env=saga)
    fdr.grd<-raster('fdr.asc')
    
    
    #Snap pour point
    buffer <-   raster::extract(dem_ca.grd, outlet.shp, buffer=5,cellnumbers = T)[[1]] %>%  as.data.frame
    snap_loc <- buffer$cell[which.max(buffer$value)]
    snap_loc <- xyFromCell(dem_ca.grd, snap_loc)
    
    #Delineate Watershed
    rsaga.geoprocessor(lib='ta_hydrology',4,  
                       param=list(TARGET_PT_X = snap_loc[1,1],
                                  TARGET_PT_Y = snap_loc[1,2], 
                                  ELEVATION = 'dem_fill.sgrd', 
                                  AREA= 'dem_ca.sgrd',
                                  METHOD=0),
                       env=saga)
    
    #Convert to polygon
    rsaga.geoprocessor(lib = 'shapes_grid', 6,
                       param = list(GRID = 'dem_ca.sgrd',
                                    POLYGONS = 'watershed.shp',
                                    CLASS_ALL = 0,
                                    CLASS_ID = 100,
                                    SPLIT = 0), 
                       env=saga)
    watershed.shp<-readShapeSpatial('watershed.shp')
    
    #Write watershed.shp to global environment
    watershed.shp
  }
  
  #Delineate
  watershed.shp<-delineate.fun(pp.shp, dem.grd,fdr.grd,saga.path, scratchspace)
  watershed.shp@proj4string<-dem.grd@crs
  
  ####################################################################################
  #Step 3: Identify "sinks" and basins (ArcPyGeo)-------------------------------------
  ####################################################################################
  #Just for reference: http://sdmg.forestry.oregonstate.edu/node/89
  
  #Setup Rpygeo environment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create python environment
  py.env<-rpygeo.build.env(python.path = "C:\\Python27\\ArcGIS10.2\\", 
                           workspace = scratchspace,
                           overwriteoutput = 1) 
  
  #Set R directoy to python scratchworkspace
  setwd("C:\\Python27\\ScratchWorkspace\\")
  
  #Clear Scratch Space
  file.remove(list.files())
  
  #Copy pertinent files to python scratchspace
  writeRaster(dem.grd, "dem.asc", overwrite=T)
  writeOGR(watershed.shp, ".", "watershed",driver="ESRI Shapefile")
  
  #Find Sinks and associated watersheds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Filter the DEM
  dem_filter.grd<- focal(dem.grd, w=focalWeight(dem.grd, 5, "Gauss"))
  dem_filter.grd<- focal(dem_filter.grd, w=focalWeight(dem_filter.grd, 5, "Gauss"))
  dem_filter.grd<- focal(dem_filter.grd, w=focalWeight(dem_filter.grd, 5, "Gauss"))
  dem_filter.grd@crs<-dem.grd@crs
  writeRaster(dem_filter.grd, file="dem_filter.asc", overwrite=T)
  
  #fill sinks (up to 0.1m)
  rpygeo.geoprocessor(fun="Fill_sa",
                      c("dem_filter.asc", "dem_fill",0.1),
                      env=py.env
  )
  
  #Compute Flow Direction
  rpygeo.FlowDirection.sa("dem_fill","fdr_esri", env=py.env)
  
  #Identify Sink
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
  
  #Clip everything to delineated watershed
  rpygeo.geoprocessor(fun="Clip_analysis ",
                      c("basin.shp","watershed.shp","basin_clip"),
                      env=py.env)
  rpygeo.geoprocessor(fun="Clip_analysis ",
                      c("sink.shp","watershed.shp","sink_clip"),
                      env=py.env)
  
  #Bring everything into the R environment
  sink.shp<-readOGR(".","sink_clip")
  basin.shp<-readOGR(".","basin_clip")
  
  ####################################################################################
  #Step 4: Estimate Spill elevations (Raster)-----------------------------------------
  ####################################################################################
  #rename basin id's
  basin.shp@data$ID<-seq(1:length(basin.shp))
  
  #set working directory
  setwd(paste0(wd,"/ScratchWorkspace"))
  
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
  #Step 5: Export Data----------------------------------------------------------------
  ####################################################################################
  #Assign data to Global Env
  assign('watershed.shp',watershed.shp, env=.GlobalEnv)
  assign('sink.shp',sink.shp, env=.GlobalEnv)
  assign('basin.shp',basin.shp, env=.GlobalEnv)
  assign('area',area, env=.GlobalEnv)
  assign('volume',volume, env=.GlobalEnv)
}
