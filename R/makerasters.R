#' Generate maps of Pleistocene island extents
#' 
#' This function generates a series of map outputs based on the interval bins contained in the interval 
#' file (see getintervals_time() and getintervals_sealvl()). This function takes a user-supplied input bathymetry raster and 
#' generates a series of island extents based on the mean sea levels specified in the interval file. This function generates maps in 
#' three formats: ESRI shapefile format, a flat raster format (with no topography), and a topographic raster format that preserves 
#' the original bathymetric elevations of each island pixel, reprojected to the user-specified projected coordinate system. 
#' 
#' @param inputraster An input bathymetry raster in ASCII (.asc) format. Although PleistoDist should theoretically be able to 
#' use any time of ASCII-formatted bathymetry grid as input, this tool has been tested specifically with data from the General
#' Bathymetric Chart of the Oceans (GEBCO: https://www.gebco.net). Locality-specific bathymetric maps can be downloaded from 
#' https://download.gebco.net/. 
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map 
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system. 
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's 
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result 
#' in distance matrices calculated in decimal degrees rather than in distance units. 
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl() 
#' function. By default, this function will use the "intervals.csv" file stored in the output folder, but users can also specify their 
#' own custom interval file (with nice round mean sea level values, for example) although the file will need to have a column named MeanDepth. 
#' 
#' @export
makemaps <- function(inputraster,epsg,intervalfile="output/intervals.csv") {
  intervalfile = read.csv(intervalfile)
  inputraster = terra::rast(inputraster)

  numintervals = max(intervalfile$Interval)
  
  for (x in 0:numintervals) {
    
    message("Generating maps for interval ",x,", sea level = ",intervalfile$MeanDepth[x+1])
    
    #reclassify values below interval sea level as NA
    convmat <- cbind(as.numeric(terra::global(inputraster,fun="min")),intervalfile$MeanDepth[x+1],NA)
    outraster <- terra::classify(inputraster,convmat,right=FALSE)
    
    #reproject outraster using bilinear method (since raw input is continuous)
    outraster_projected <- terra::project(outraster,paste("EPSG:",epsg,sep=""),method="bilinear")
    
    #check to see if raster folder exists, create folder if it doesn't already exist
    if (base::dir.exists("output/raster_topo") == FALSE) {
      base::dir.create("output/raster_topo")
    }
    
    #write projected raster to raster folder
    terra::writeRaster(outraster_projected,paste("output/raster_topo/interval",x,".tif",sep=""),filetype="GTiff",overwrite=TRUE)
    
    #reclassify all values as 1 to create flat raster
    convmat <- cbind(intervalfile$MeanDepth[x+1],as.numeric(global(inputraster,fun="max")),1)
    outraster_flat <- terra::classify(outraster_projected,convmat)
    
    #write flat raster to raster_flat folder
    if (base::dir.exists("output/raster_flat") == FALSE) {
      base::dir.create("output/raster_flat")
    }
    terra::writeRaster(outraster_flat,paste("output/raster_flat/interval",x,".tif",sep=""),filetype="GTiff",overwrite=TRUE)
    
    outvector <- terra::as.polygons(outraster_flat)
    outvector <- terra::disagg(outvector)
    
    #check to see if shapefile folder exists, create folder if it doesn't already exist
    if (base::dir.exists("output/shapefile") == FALSE) {
      base::dir.create("output/shapefile")
    }
    
    #write reprojected vector to shapefile folder
    terra::writeVector(outvector,paste("output/shapefile/interval",x,".shp",sep=""),filetype="ESRI Shapefile",overwrite=TRUE)
  }
  message("Done!")
}