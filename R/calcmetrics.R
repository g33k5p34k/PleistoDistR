#' Calculate Euclidean distance between points
#'
#' This function calculates the straight-line distance ("as the crow flies) between user-specified points.
#' Because these points are invariant across space and time and independent of sea level change, this calculation
#' is only performed once instead of repeatedly across different sea level or time intervals.
#'
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @return This function outputs a pairwise distance matrix of inter-point Euclidean distances in long format.
#' @examples
#' #load example data (AntPops.shp, adapted from Darwell et al (2020))
#' #note that no interval file input is needed (since Euclidean distances are time-invariant)
#' pleistodist_euclidean(points="inst/extdata/AntPops.shp",epsg=3141)
#' @export
pleistodist_euclidean <- function(points,epsg) {
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
    p2_names <- points_transformed$Name

  } else {
    p1_names <- points_transformed$FID
    p2_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing Euclidean distance matrix
  euclideandist <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names,include.equals=T)[,2],
    Island2 = expand.grid(p1_names,p2_names,include.equals=T)[,1]
  )

  euclidean <- sf::st_distance(points_transformed)
  blankmatrix <- expand.grid(1:nrow(points_transformed),1:nrow(points_transformed))
  blankmatrix$interval0 <- tidyr::gather(as.data.frame(euclidean))$value
  #blankmatrix <- blankmatrix %>%
    #ilter(Var1>=Var2)
  euclideandist$interval0 <- as.numeric(blankmatrix$interval0)
  write.csv(euclideandist,"output/island_euclideandist.csv")
}

#' Calculate inter-island centroid-to-centroid distance over time
#'
#' This function calculates the centroid-to-centroid distance between island landmasses over Pleistocene time. For each pairwise combination
#' of source points provided by the user, this function selects the landmasses that correspond to each point, and calculates the distance between
#' the centroids of the selected landmasses for each time/sea level interval. If both points lie on the same landmass, a distance value of 0 is returned. If one or both islands
#' is/are underwater, then a value of NA is returned. After calculating the pairwise island distances for each interval, this function calculates
#' the time-weighted average centroid-to-centroid distance between each island pair.
#'
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean distance calculations.
#' @return This function outputs a pairwise distance matrix of pairwise inter-island centroid-to-centroid distances in long format, with
#' one column per interval.
#' @examples
#' pleistodist_centroid(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv")
#' @export
pleistodist_centroid <- function(points,epsg,intervalfile="output/intervals.csv") {
  intervalfile = read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)
  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
    p2_names <- points_transformed$Name

  } else {
    p1_names <- points_transformed$FID
    p2_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing centroid-to-centroid distance matrix
  centroiddist <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names,include.equals=T)[,2],
    Island2 = expand.grid(p1_names,p2_names,include.equals=T)[,1]
  )

  #create dataframe for storing relative island width matrix
  relativeislandwidth <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names)[,2],
    Island2 = expand.grid(p1_names,p2_names)[,1]
  )

  #for each interval...
  for (i in 0:numintervals) {
    #load interval shapefile
    invector <- sf::st_read(paste("output/shapefile/interval",i,".shp",sep="")) #load interval shapefile
    intvlname <- paste("interval",i,sep="")

    #create empty lists for storing distance values at each interval
    centroidd <- c()
    iwidth <- c()

    message(paste("Analysing island distances for interval ",i,"...",sep=""))
    #start of first order for-loop to select island 1
    for (f in 1:nrow(points_transformed)) {

      #start second order for-loop to select island 2
      for (t in 1:nrow(points_transformed)) {

        #extract coordinates for island 1
        p1 <- points_transformed$geometry[f]
        #extract coordinates for island 2
        p2 <- points_transformed$geometry[t]

        message(paste("Calculating distance between points ",p1_names[f]," and ",p2_names[t],"... ",sep=""))

        if (f == t) { #write NA for self-comparisons
          centroidd <- c(centroidd,NA)
          iwidth <- c(iwidth,NA)
        } else if ((nrow(invector[p1, ,op=sf::st_contains]) || nrow(invector[p2, ,op=sf::st_contains])) == 0) { #check to see if island1 or island2 is submerged under water
          centroidd <- c(centroidd,NA)
          iwidth <- c(iwidth,NA)
          message("One or both islands are underwater during this interval, writing a value of NA")
        } else {
          combinedpoints <- sf::st_sfc(rbind(p1,p2),crs=epsg)
          islandpair <- invector[combinedpoints, ,op=sf::st_contains]
          if (nrow(islandpair) == 1) { #check to see if both points are on the same island, if TRUE write distance of 0 for island-based comparisons
            centroidd <- c(centroidd,0)
            iwidth <- c(iwidth,NA)
            message("Both points are on the same island, writing a value of 0 for centroid-to-centroid distance")
          } else {
            island1 <- invector[p1, ,op=sf::st_contains]
            island2 <- invector[p2, ,op=sf::st_contains]
            island1_centroid <- sf::st_centroid(island1)
            island2_centroid <- sf::st_centroid(island2)
            centroidline <- sf::st_linestring(rbind(sf::st_coordinates(island1_centroid),sf::st_coordinates(island2_centroid)))

            m = -1/(line(centroidline)$coefficients[2])
            c = sf::st_coordinates(island1_centroid)[,2]-(m*sf::st_coordinates(island1_centroid)[,1])

            #calculate 2 points on widthline at edges of bounding box
            x1 = sf::st_bbox(island1)$xmax
            y1 = m*x1+c

            x2 = sf::st_bbox(island1)$xmin
            y2 = m*x2+c

            coord1 <- sf::st_point(c(x1,y1))
            coord2 <- sf::st_point(c(x2,y2))
            widthline <- sf::st_sfc(sf::st_linestring(c(coord1,coord2)),crs=epsg)
            #because complex island shapes can slice widthlines into multiple fragments, convert widthline into coordinates,
            #generate a new line, and calculate length from this new line
            islandwidth<- sf::st_length(sf::st_linestring(sf::st_coordinates(sf::st_intersection(widthline,island1))))
            iwidth = c(iwidth,islandwidth)

            cdistance <- as.numeric(sf::st_length(centroidline))
            centroidd <- c(centroidd,cdistance)
            message(paste("Distance between island centroids is ",cdistance,"m",sep=""))
          }
        }
      }
    }
    centroiddist[intvlname] <- centroidd
    relativeislandwidth[intvlname] <- iwidth
  }
  centroiddist$mean <- apply(centroiddist[3:ncol(centroiddist)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  relativeislandwidth$mean <- apply(relativeislandwidth[3:ncol(relativeislandwidth)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(centroiddist,"output/island_centroiddist.csv")
  write.csv(relativeislandwidth,"output/island_relativewidth.csv")
}

#' Calculate least inter-island shore-to-shore distance over time
#'
#' This function calculates the least shore-to-shore distance between island landmasses over Pleistocene time. For each pairwise combination
#' of source points provided by the user, this function selects the landmasses that correspond to each point, and calculates the minimum
#' shore-to-shore distance between the selected landmasses for each time/sea level interval. If both points lie on the same landmass, a distance value of 0 is returned. If one or both islands
#' is/are underwater, then a value of NA is returned. After calculating the pairwise island distances for each interval, this function calculates
#' the time-weighted average least shore-to-shore distance between each island pair.
#'
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean distance calculations.
#' @return This function outputs a pairwise distance matrix of pairwise inter-island least shore-to-shore distances in long format, with
#' one column per interval.
#' @examples
#' pleistodist_leastshore(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv")
#' @export
pleistodist_leastshore <- function(points,epsg,intervalfile="output/intervals.csv") {
  intervalfile = read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)
  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
    p2_names <- points_transformed$Name

  } else {
    p1_names <- points_transformed$FID
    p2_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing shore-to-shore distance matrix
  shoredist <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names,include.equals=T)[,2],
    Island2 = expand.grid(p1_names,p2_names,include.equals=T)[,1]
  )

  #create dataframe for storing relative island width matrix
  relativeislandwidth <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names)[,2],
    Island2 = expand.grid(p1_names,p2_names)[,1]
  )

  #for each interval...
  for (i in 0:numintervals) {
    #load interval shapefile
    invector <- sf::st_read(paste("output/shapefile/interval",i,".shp",sep="")) #load interval shapefile
    intvlname <- paste("interval",i,sep="")
    #create empty lists for storing distance values at each interval
    shored <- c()
    iwidth <- c()
    message(paste("Analysing island distances for interval ",i,"...",sep=""))
    #start of first order for-loop to select island 1
    for (f in 1:nrow(points_transformed)) {
      #start second order for-loop to select island 2
      for (t in 1:nrow(points_transformed)) {
        #extract coordinates for island 1
        p1 <- points_transformed$geometry[f]
        #extract coordinates for island 2
        p2 <- points_transformed$geometry[t]
        message(paste("Calculating distance between points ",p1_names[f]," and ",p2_names[t],"... ",sep=""))
        if (f == t) { #write NA for self-comparisons
          shored <- c(shored,NA)
          iwidth <- c(iwidth,NA)
        } else if ((nrow(invector[p1, ,op=sf::st_contains]) || nrow(invector[p2, ,op=sf::st_contains])) == 0) { #check to see if island1 or island2 is submerged under water
          shored <- c(shored,NA)
          iwidth <- c(iwidth,NA)
          message("One or both islands are underwater during this interval, writing a value of NA")
        } else {
          combinedpoints <- sf::st_sfc(rbind(p1,p2),crs=epsg)
          islandpair <- invector[combinedpoints, ,op=sf::st_contains]
          if (nrow(islandpair) == 1) { #check to see if both points are on the same island, if TRUE write distance of 0 for island-based comparisons
            shored <- c(shored,0)
            iwidth <- c(iwidth,0)
            message("Both points are on the same island, writing a value of 0 for least shore-to-shore distance")
          } else {
            island1 <- invector[p1, ,op=sf::st_contains]
            island2 <- invector[p2, ,op=sf::st_contains]
            island1_centroid <- sf::st_centroid(island1)
            island2_centroid <- sf::st_centroid(island2)
            centroidline <- sf::st_linestring(rbind(sf::st_coordinates(island1_centroid),sf::st_coordinates(island2_centroid)))

            m = -1/(line(centroidline)$coefficients[2])
            c = sf::st_coordinates(island1_centroid)[,2]-(m*sf::st_coordinates(island1_centroid)[,1])

            #calculate 2 points on widthline at edges of bounding box
            x1 = sf::st_bbox(island1)$xmax
            y1 = m*x1+c

            x2 = sf::st_bbox(island1)$xmin
            y2 = m*x2+c

            coord1 <- sf::st_point(c(x1,y1))
            coord2 <- sf::st_point(c(x2,y2))
            widthline <- sf::st_sfc(sf::st_linestring(c(coord1,coord2)),crs=epsg)
            #because complex island shapes can slice widthlines into multiple fragments, convert widthline into coordinates,
            #generate a new line, and calculate length from this new line
            islandwidth<- sf::st_length(sf::st_linestring(sf::st_coordinates(sf::st_intersection(widthline,island1))))
            iwidth = c(iwidth,islandwidth)

            sdistance <- as.numeric(sf::st_distance(island1,island2))
            shored <- c(shored,sdistance)
            message(paste("Least shore-to-shore distance between islands is ",sdistance,"m",sep=""))
          }
        }
      }
    }
    shoredist[intvlname] <- shored
    relativeislandwidth[intvlname] <- iwidth
  }
  shoredist$mean <- apply(shoredist[3:ncol(shoredist)],1,weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  relativeislandwidth$mean <- apply(relativeislandwidth[3:ncol(relativeislandwidth)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(shoredist,"output/island_leastshoredist.csv")
  write.csv(relativeislandwidth,"output/island_relativewidth.csv")
}

#' Calculate least-cost distance between points over time
#'
#' This function calculates the least-cost overland distance between each pairwise combination of user-specified points. By assigning exposed
#' land a resistance value of 1 and submerged areas a resistance value of 999,999, this function calculates the shortest route between each pair of points
#' that minimises the number of overwater crossings. If both points lie on the same landmass, a distance value of 0 is returned. If one or both islands
#' is/are underwater, then a value of NA is returned. After calculating the pairwise least-cost distances for each interval, this function calculates
#' the time-weighted average least-cost distance between each pair of points.
#'
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean distance calculations.
#' @return This function outputs a pairwise distance matrix of pairwise inter-island least shore-to-shore distances in long format, with
#' one column per interval.
#' @examples
#' pleistodist_leastcost(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv")
#' @export
pleistodist_leastcost <- function(points,epsg,intervalfile="output/intervals.csv") {
  intervalfile = read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)
  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
    p2_names <- points_transformed$Name

  } else {
    p1_names <- points_transformed$FID
    p2_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing point-to-point least cost distance matrix
  leastcostdist <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names,include.equals=T)[,2],
    Island2 = expand.grid(p1_names,p2_names,include.equals=T)[,1]
  )

  #create dataframe for storing relative island width matrix
  relativeislandwidth <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names)[,2],
    Island2 = expand.grid(p1_names,p2_names)[,1]
  )

  #for each interval...
  for (i in 0:numintervals) {
    #load interval shapefile
    invector <- sf::st_read(paste("output/shapefile/interval",i,".shp",sep="")) #load interval shapefile
    inraster <- raster::raster(paste("output/raster_flat/interval",i,".tif",sep="")) #load flat raster for least cost path calculation
    inraster[is.na(inraster[])] <- 1/9999999 #convert raster NA values (water) to very low conductance values
    trraster <- gdistance::transition(inraster,mean,directions=8) #convert inraster into a transition matrix
    trraster_corr <- gdistance::geoCorrection(trraster,type="c",scl=FALSE) #apply geographical correction to the transition matrix (for least cost path)
    intvlname <- paste("interval",i,sep="")

    #create empty lists for storing distance values at each interval
    leastcostd <- c()
    iwidth <- c()
    message(paste("Analysing island distances for interval ",i,"...",sep=""))
    #start of first order for-loop to select island 1
    for (f in 1:nrow(points_transformed)) {
      #start second order for-loop to select island 2
      for (t in 1:nrow(points_transformed)) {
        #extract coordinates for island 1
        p1 <- points_transformed$geometry[f]
        #extract coordinates for island 2
        p2 <- points_transformed$geometry[t]
        message(paste("Calculating distance between points ",p1_names[f]," and ",p2_names[t],"... ",sep=""))
        if (f == t) { #write NA for self-comparisons
          leastcostd <- c(leastcostd,NA)
          iwidth <- c(iwidth,NA)
        } else if ((nrow(invector[p1, ,op=sf::st_contains]) || nrow(invector[p2, ,op=sf::st_contains])) == 0) { #check to see if island1 or island2 is submerged under water
          leastcostd <- c(leastcostd,NA)
          iwidth <- c(iwidth,NA)
          message("One or both islands are underwater during this interval, writing a value of NA")
        } else {
          combinedpoints <- sf::st_sfc(rbind(p1,p2),crs=epsg)
          islandpair <- invector[combinedpoints, ,op=sf::st_contains]
          if (nrow(islandpair) == 1) { #check to see if both points are on the same island, if TRUE write distance of 0 for island-based comparisons
            lcdistance <- as.numeric(sp::SpatialLinesLengths(gdistance::shortestPath(trraster_corr,sf::st_coordinates(p1),sf::st_coordinates(p2),output="SpatialLines")))
            leastcostd <- c(leastcostd,lcdistance) #for least cost distance, calculate least cost path
            iwidth <- c(iwidth,0)
            message(paste("Least cost distance between points is ",lcdistance,"m",sep=""))
          } else {
            island1 <- invector[p1, ,op=sf::st_contains]
            island2 <- invector[p2, ,op=sf::st_contains]
            island1_centroid <- sf::st_centroid(island1)
            island2_centroid <- sf::st_centroid(island2)
            centroidline <- sf::st_linestring(rbind(sf::st_coordinates(island1_centroid),sf::st_coordinates(island2_centroid)))

            m = -1/(line(centroidline)$coefficients[2])
            c = sf::st_coordinates(island1_centroid)[,2]-(m*sf::st_coordinates(island1_centroid)[,1])

            #calculate 2 points on widthline at edges of bounding box
            x1 = sf::st_bbox(island1)$xmax
            y1 = m*x1+c

            x2 = sf::st_bbox(island1)$xmin
            y2 = m*x2+c

            coord1 <- sf::st_point(c(x1,y1))
            coord2 <- sf::st_point(c(x2,y2))
            widthline <- sf::st_sfc(sf::st_linestring(c(coord1,coord2)),crs=epsg)
            #because complex island shapes can slice widthlines into multiple fragments, convert widthline into coordinates,
            #generate a new line, and calculate length from this new line
            islandwidth<- sf::st_length(sf::st_linestring(sf::st_coordinates(sf::st_intersection(widthline,island1))))
            iwidth = c(iwidth,islandwidth)
            lcdistance <- as.numeric(sp::SpatialLinesLengths(gdistance::shortestPath(trraster_corr,sf::st_coordinates(p1),sf::st_coordinates(p2),output="SpatialLines")))
            leastcostd <- c(leastcostd,lcdistance)
            message(paste("Least cost distance between points is ",lcdistance,"m",sep=""))
          }
        }
      }
    }
    leastcostdist[intvlname] <- leastcostd
    relativeislandwidth[intvlname] <- iwidth
  }
  leastcostdist$mean <- apply(leastcostdist[3:ncol(leastcostdist)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  relativeislandwidth$mean <- apply(relativeislandwidth[3:ncol(relativeislandwidth)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(leastcostdist,"output/island_leastcostdist.csv")
  write.csv(relativeislandwidth,"output/island_relativewidth.csv")
}

#' Calculate inter-island mean shore-to-shore distance over time
#'
#' This function calculates the mean shore-to-shore distance between island landmasses over Pleistocene time. This function is slightly different
#' from other PleistoDist functions because the resultant distance matrix is asymmetrical. This is due to the fact that the mean shore-to-shore distance from
#' a large island to a small one will be greater relative to the reverse, since larger islands have longer and more complex shorelines, so the mean dispersal distance
#' from any point along the larger island's shoreline facing the smaller island is likely to be greater than in the reverse case. To calculate the mean
#' shore-to-shore distance, PleistoDist generates a number of points at equal intervals (1 point every 100m by default) along the shoreline of the source island at each time/sea level interval, filters out points not directly
#' facing the receiving island, and calculates the average minimum distance between each point along the source island shoreline and the shoreline of the
#' receiving island. Because the number of shoreline points can run into the thousands for very large islands, thereby increasing computational time, users
#' can set a cap on the maximum number of shoreline points to sample (the default cap is 1000 points). As with other PleistoDist functions, if both points lie on the same landmass, a distance value of 0 is returned. If one or both islands
#' is/are underwater, then a value of NA is returned. After calculating the pairwise island distances for each interval, this function calculates
#' the time-weighted average mean shore-to-shore distance between each island pair.
#'
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean distance calculations.
#' @param maxsamp Maximum number of points to sample along the donor island shoreline, mainly to limit computational time. If users are
#' ok with long runtimes, then maxsamp can be set to NA.
#' @return This function outputs an asymmetric pairwise distance matrix of mean shore-to-shore distances in long format, with
#' one column per interval. Because this distance matrix is asymmetric, the directionality of the distance calculation matters.
#' @examples
#' #load example data, with shoreline sampling points capped at 1000
#' pleistodist_meanshore(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv",maxsamp=1000)
#' #load example data, with no cap on the number of sampling points along the shoreline (will lead to very long runtimes for large datasets)
#' pleistodist_meanshore(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv",maxsamp=NA)
#' @export
pleistodist_meanshore <- function(points,epsg,intervalfile="output/intervals.csv",maxsamp=1000) {
  intervalfile = read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)
  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
    p2_names <- points_transformed$Name

  } else {
    p1_names <- points_transformed$FID
    p2_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing mean shore distance matrix
  meanshoredist <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names)[,2],
    Island2 = expand.grid(p1_names,p2_names)[,1]
  )

  #create dataframe for storing relative island width matrix
  relativeislandwidth <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names)[,2],
    Island2 = expand.grid(p1_names,p2_names)[,1]
  )

  #for each interval...
  for (i in 0:numintervals) {
    #load interval shapefile
    invector <- sf::st_read(paste("output/shapefile/interval",i,".shp",sep="")) #load interval shapefile
    intvlname <- paste("interval",i,sep="")
    meanshored <- c()
    iwidth <- c()

    #start of first order for-loop to select island 1
    for (f in 1:nrow(points_transformed)) {

      #start second order for-loop to select island 2
      for (t in 1:nrow(points_transformed)) {

        #extract coordinates for island 1
        p1 <- points_transformed$geometry[f]
        #extract coordinates for island 2
        p2 <- points_transformed$geometry[t]

        message(paste("Calculating distance between points ",p1_names[f]," and ",p2_names[t],"... ",sep=""))

        if (f == t) { #write NA for self-comparisons
          meanshored <- c(meanshored,NA)
          iwidth <- c(iwidth,NA)
        } else if ((nrow(invector[p1, ,op=sf::st_contains]) || nrow(invector[p2, ,op=sf::st_contains])) == 0) { #check to see if island1 or island2 is submerged under water
          meanshored <- c(meanshored,NA)
          iwidth <- c(iwidth,NA)
          message("One or both islands are underwater during this interval, writing a value of NA")
        } else {
          combinedpoints <- sf::st_sfc(rbind(p1,p2),crs=epsg)
          islandpair <- invector[combinedpoints, ,op=sf::st_contains]
          if (nrow(islandpair) == 1) { #check to see if both points are on the same island, if TRUE write distance of 0 for island-based comparisons
            meanshored <- c(meanshored,0)
            iwidth <- c(iwidth,NA)
            message("Both points are on the same island, writing a value of 0")
          } else {
            island1 <- invector[p1, ,op=sf::st_contains]
            island2 <- invector[p2, ,op=sf::st_contains]
            island1_centroid <- sf::st_centroid(island1)
            island2_centroid <- sf::st_centroid(island2)
            centroidline <- sf::st_linestring(rbind(sf::st_coordinates(island1_centroid),sf::st_coordinates(island2_centroid)))

            m = -1/(line(centroidline)$coefficients[2])
            c = sf::st_coordinates(island1_centroid)[,2]-(m*sf::st_coordinates(island1_centroid)[,1])

            #calculate 2 points on widthline at edges of bounding box
            x1 = sf::st_bbox(island1)$xmax
            y1 = m*x1+c

            x2 = sf::st_bbox(island1)$xmin
            y2 = m*x2+c

            coord1 <- sf::st_point(c(x1,y1))
            coord2 <- sf::st_point(c(x2,y2))
            widthline <- sf::st_sfc(sf::st_linestring(c(coord1,coord2)),crs=epsg)
            #because complex island shapes can slice widthlines into multiple fragments, convert widthline into coordinates,
            #generate a new line, and calculate length from this new line
            islandwidth<- sf::st_length(sf::st_linestring(sf::st_coordinates(sf::st_intersection(widthline,island1))))
            iwidth = c(iwidth,islandwidth)

            island1_shoreline <- sf::st_cast(island1,"LINESTRING")
            #pick the line segment that is closest to the centroid of the receiving island (since there may be inland lakes)
            lineseg <- which.min(sf::st_distance(island1_shoreline,island2))
            island1_shoreline <- island1_shoreline$geometry[lineseg]

            #generate points at 100m intervals along shoreline
            island1_shorepoints <- sf::st_line_sample(island1_shoreline,density=1/100)
            island1_shorepoints <- sf::st_cast(island1_shorepoints,"POINT")

            #downsample the number of shorepoints if the island is very big (to speed up computations)
            if (length(island1_shorepoints) > maxsamp && !is.na(maxsamp)) {
              island1_shorepoints <- sf::st_line_sample(island1_shoreline,n=maxsamp)
              island1_shorepoints <- sf::st_cast(island1_shorepoints,"POINT")
            }

            #draw a line between each shorepoint and the receiving island, and exclude points with lines crossing the shore multiple times (indicating that it faces away from the receiving island)
            island1_shorepoints_filter1 <- sf::st_crosses(sf::st_nearest_points(island1_shorepoints,island2),island1_shoreline,sparse=FALSE)
            #subset points that face receiving island
            island1_shorepoints_filtered <- island1_shorepoints[which(island1_shorepoints_filter1==FALSE)]
            #apply a second filter to eliminate points that cross the island more than once
            island1_shorepoints_filter2 <- sf::st_overlaps(sf::st_nearest_points(island1_shorepoints_filtered,island2),island1_shoreline,sparse=FALSE)
            island1_shorepoints_final <- island1_shorepoints_filtered[which(island1_shorepoints_filter2==FALSE)]

            msd <- as.numeric(mean(sf::st_distance(island1_shorepoints_final,island2)))
            meanshored <- c(meanshored,msd)
            #meanshored_stdev <- c(meanshored_stdev,as.numeric(sd(sf::st_distance(island1_shorepoints,islandpair$geometry[2]))))
            message(paste("The mean shore-to-shore distance between islands is ",msd,"m",sep=""))
          }
        }
      }
    }
    meanshoredist[intvlname] <- meanshored
    relativeislandwidth[intvlname] <- iwidth
  }
  meanshoredist$mean <- apply(meanshoredist[3:ncol(meanshoredist)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  relativeislandwidth$mean <- apply(relativeislandwidth[3:ncol(relativeislandwidth)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(meanshoredist,"output/island_meanshoredist.csv")
  write.csv(relativeislandwidth,"output/island_relativewidth.csv")
}

#' Calculate width of islands relative to each other
#'
#' This function calculates the width of source islands relative to their recipient islands for every pairwise combination of islands as specified
#' by the user-defined source points. The island width is calculated by drawing a line between the centroids of the two islands, and plotting
#' a line perpendicular to the centroid-centroid line extending from the centroid of the source island to the island edges. This widthline changes
#' relative to the orientation of the receiving island, and is mostly needed for calculating the expected net inter-island migration, as defined
#' in chapter 6 of MacArthur and Wilson (1967).
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean distance calculations.
#' @return This function outputs an asymmetric pairwise distance matrix of relative island widths in long format, with
#' one column per interval. Because this distance matrix is asymmetric, the directionality of the distance calculation matters.
#' @examples
#' pleistodist_relativewidth(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv")
#' @export
pleistodist_relativewidth <- function (points, epsg,intervalfile="output/intervals.csv") {
  intervalfile = read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)
  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
    p2_names <- points_transformed$Name
  } else {
    p1_names <- points_transformed$FID
    p2_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing relative island width matrix
  relativeislandwidth <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names)[,2],
    Island2 = expand.grid(p1_names,p2_names)[,1]
  )

  #for each interval...
  for (i in 0:numintervals) {
    #load interval shapefile
    invector <- sf::st_read(paste("output/shapefile/interval",i,".shp",sep="")) #load interval shapefile
    intvlname <- paste("interval",i,sep="")

    #create empty lists for storing distance values at each interval
    iwidth <- c()

    message(paste("Analysing island width for interval ",i,"...",sep=""))
    #start of first order for-loop to select island 1
    for (f in 1:nrow(points_transformed)) {

      #start second order for-loop to select island 2
      for (t in 1:nrow(points_transformed)) {

        #extract coordinates for island 1
        p1 <- points_transformed$geometry[f]
        #extract coordinates for island 2
        p2 <- points_transformed$geometry[t]

        message(paste("Calculating distance between points ",p1_names[f]," and ",p2_names[t],"... ",sep=""))

        if (f == t) { #write NA for self-comparisons
          iwidth <- c(iwidth,NA)
        } else if ((nrow(invector[p1, ,op=sf::st_contains]) || nrow(invector[p2, ,op=sf::st_contains])) == 0) { #check to see if island1 or island2 is submerged under water
          iwidth <- c(iwidth,NA)
          message("One or both islands are underwater during this interval, writing a value of NA")
        } else {
          combinedpoints <- sf::st_sfc(rbind(p1,p2),crs=epsg)
          islandpair <- invector[combinedpoints, ,op=sf::st_contains]
          if (nrow(islandpair) == 1) { #check to see if both points are on the same island, if TRUE write distance of 0 for island-based comparisons
            iwidth <- c(iwidth,NA)
            message("Both points are on the same island, writing a value of 0 for centroid-to-centroid distance")
          } else {
            island1 <- invector[p1, ,op=sf::st_contains]
            island2 <- invector[p2, ,op=sf::st_contains]
            island1_centroid <- sf::st_centroid(island1)
            island2_centroid <- sf::st_centroid(island2)
            centroidline <- sf::st_linestring(rbind(sf::st_coordinates(island1_centroid),sf::st_coordinates(island2_centroid)))

            m = -1/(line(centroidline)$coefficients[2])
            c = sf::st_coordinates(island1_centroid)[,2]-(m*sf::st_coordinates(island1_centroid)[,1])

            #calculate 2 points on widthline at edges of bounding box
            x1 = sf::st_bbox(island1)$xmax
            y1 = m*x1+c

            x2 = sf::st_bbox(island1)$xmin
            y2 = m*x2+c

            coord1 <- sf::st_point(c(x1,y1))
            coord2 <- sf::st_point(c(x2,y2))
            widthline <- sf::st_sfc(sf::st_linestring(c(coord1,coord2)),crs=epsg)
            #because complex island shapes can slice widthlines into multiple fragments, convert widthline into coordinates,
            #generate a new line, and calculate length from this new line
            islandwidth<- sf::st_length(sf::st_linestring(sf::st_coordinates(sf::st_intersection(widthline,island1))))
            iwidth = c(iwidth,islandwidth)
          }
        }
      }
    }
    relativeislandwidth[intvlname] <- iwidth
  }
  relativeislandwidth$mean <- apply(relativeislandwidth[3:ncol(relativeislandwidth)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(relativeislandwidth,"output/island_relativewidth.csv")
}

#' Calculate area, perimeter and surface area of islands over time
#'
#' This function calculates the area, perimeter, and surface area of user-specified islands at each time/sea level interval, and
#' performs the same calculations as the pleistoshape_area(), pleistoshape_perimeter, and pleistoshape_surfacearea() functions combined.
#' This function is nonetheless provided since it allows users to calculate all three shape metrics simultaneously, which is faster than running
#' each constituent function individually. PleistoDist will return a value of NA if the island is submerged underwater for that particular interval. After calculating shape metrics for each interval, this function calculates
#' the time-weighted average area, perimeter, and surface area of each selected island.
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean calculations.
#' @return This function outputs three matrices of island shape estimates in long format, with one column per interval in each matrix.
#' @examples
#' pleistoshape_all(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv")
#' @export
pleistoshape_all <- function (points,epsg,intervalfile="output/intervals.csv") {
  intervalfile <- read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)

  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
  } else {
    p1_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing island area
  islandarea <- tibble::tibble(
    Island = p1_names
  )
  #create dataframe for storing island perimeter
  islandperimeter <- tibble::tibble(
    Island = p1_names
  )
  #create dataframe for storing island perimeter
  islandsurfacearea <- tibble::tibble(
    Island = p1_names
  )

  #for each interval...
  for (i in 0:numintervals) {
    #load interval shapefile
    invector <- sf::st_read(paste("output/shapefile/interval",i,".shp",sep=""))
    inraster <- terra:rast(paste("output/raster_topo/interval",i,".tif",sep=""))-intervalfile$MeanDepth[i+1]
    intvlname <- paste("interval",i,sep="")
    ia <- c() #create blank list for temporarily storing area values (before offloading into the main dataframe)
    ip <- c() #create blank list for temporarily storing perimeter values (before offloading into the main dataframe)
    isa <- c() #create blank list for temporarily storing surface area values (before offloading into the main dataframe)

    #start of first order for-loop to select island 1
    for (f in 1:nrow(points_transformed)) {

      #extract coordinates for island
      p1 <- points_transformed$geometry[f]

      message("Calculating island area and perimeter for point ",p1_names[f],sep="")

      if (nrow(invector[p1, ,op=sf::st_contains]) == 0) { #check to see if island is submerged under water
        ia <- c(ia,0)
        ip <- c(ip,0)
        isa <- c(isa,0)
        message("Island is underwater during this interval, writing a value of 0")
      } else {
        iarea <- as.numeric(sf::st_area(invector[p1, ,op=sf::st_contains]))
        iperimeter <- as.numeric(lwgeom::st_perimeter(invector[p1, ,op=sf::st_contains]))
        islandmask <- terra::vect(invector[p1, ,op=sf::st_contains]) #use shapefile to create a mask layer
        islandtopo <- terra::mask(inraster,islandmask) #mask topo raster to isolate single island
        islandtopo_cropped <- terra::crop(islandtopo,islandmask) #crop raster extent to single island
        islandtopo_cropped_sp <- sf::as(raster::raster(islandtopo_cropped),"SpatialPixelsDataFrame") #convert spatraster to sp format
        isurfacearea <- sp::surfaceArea(islandtopo_cropped_sp) #calculate surface area (incorporating elevation)
        ia <- c(ia,iarea)
        ip <- c(ip,iperimeter)
        isa <- c(isa,isurfacearea)
        message(paste("Island perimeter is ",iperimeter,"m",sep=""))
        message(paste("Island area is ",iarea,"m^2",sep=""))
        message(paste("Island surface area is ",isurfacearea,"m^2",sep=""))
      }
    }
    islandarea[intvlname] <- ia
    islandperimeter[intvlname] <- ip
    islandsurfacearea[intvlname] <- isa
  }
  islandarea$mean <- apply(islandarea[2:ncol(islandarea)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  islandperimeter$mean <- apply(islandperimeter[2:ncol(islandperimeter)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  islandsurfacearea$mean <- apply(islandsurfacearea[2:ncol(islandsurfacearea)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(islandarea,"output/island_area.csv")
  write.csv(islandperimeter,"output/island_perimeter.csv")
  write.csv(islandsurfacearea,"output/island_surfacearea.csv")
}

#' Calculate island perimeter over time
#'
#' This function calculates the perimeter of user-specified islands at each time/sea level interval, and calculates the weighted average perimter over the user-defined time period.
#' PleistoDist will return a value of NA if the island is submerged underwater for that particular interval.
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean calculations.
#' @return This function outputs a matrix of island perimeter estimates in long format, with one column per interval.
#' @examples
#' pleistoshape_perimeter(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv")
#' @export
pleistoshape_perimeter <- function (points,epsg,intervalfile="output/intervals.csv") {
  intervalfile <- read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)
  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
  } else {
    p1_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }
  #create dataframe for storing island perimeter
  islandperimeter <- tibble::tibble(
    Island = p1_names
  )
  #for each interval...
  for (i in 0:numintervals) {
    #load interval shapefile
    invector <- sf::st_read(paste("output/shapefile/interval",i,".shp",sep=""))
    ip <- c() #create blank list for temporarily storing perimeter values (before offloading into the main dataframe)

    #start of first order for-loop to select island 1
    for (f in 1:nrow(points_transformed)) {

      #extract coordinates for island
      p1 <- points_transformed$geometry[f]

      message("Calculating island area and perimeter for point ",p1_names[f],sep="")

      if (nrow(invector[p1, ,op=sf::st_contains]) == 0) { #check to see if island is submerged under water
        ip <- c(ip,0)
        message("Island is underwater during this interval, writing a value of 0")
      } else {
        iperimeter <- as.numeric(lwgeom::st_perimeter(invector[p1, ,op=sf::st_contains]))
        ip <- c(ip,iperimeter)
        message(paste("Island perimeter is ",iperimeter,"m",sep=""))
      }
    }
    islandperimeter[intvlname] <- ip
  }
  islandperimeter$mean <- apply(islandperimeter[2:ncol(islandperimeter)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(islandperimeter,"output/island_perimeter.csv")
}

#' Calculate island area over time
#'
#' This function calculates the 2-dimensional area of user-specified islands at each time/sea level interval, and calculates the weighted average area over the user-defined time period.
#' PleistoDist will return a value of NA if the island is submerged underwater for that particular interval.
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean calculations.
#' @return This function outputs a matrix of 2D island area estimates in long format, with one column per interval.
#' @examples
#' pleistoshape_area(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv")
#' @export
pleistoshape_area <- function (points,epsg,intervalfile="output/intervals.csv") {
  intervalfile <- read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)
  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
  } else {
    p1_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing island area
  islandarea <- tibble::tibble(
    Island = p1_names
  )

  #for each interval...
  for (i in 0:numintervals) {
    #load interval shapefile
    invector <- sf::st_read(paste("output/shapefile/interval",i,".shp",sep=""))
    intvlname <- paste("interval",i,sep="")
    ia <- c() #create blank list for temporarily storing area values (before offloading into the main dataframe)

    #start of first order for-loop to select island 1
    for (f in 1:nrow(points_transformed)) {
      #extract coordinates for island
      p1 <- points_transformed$geometry[f]
      message("Calculating island area ",p1_names[f],sep="")

      if (nrow(invector[p1, ,op=sf::st_contains]) == 0) { #check to see if island is submerged under water
        ia <- c(ia,0)
        message("Island is underwater during this interval, writing a value of 0")
      } else {
        iarea <- as.numeric(sf::st_area(invector[p1, ,op=sf::st_contains]))
        ia <- c(ia,iarea)
        message(paste("Island area is ",iarea,"m^2",sep=""))
      }
    }
    islandarea[intvlname] <- ia
  }
  islandarea$mean <- apply(islandarea[2:ncol(islandarea)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(islandarea,"output/island_area.csv")
}

#' Calculate island surface area over time
#'
#' This function calculates the surface area of user-specified islands at each time/sea level interval, and calculates the weighted average surface area over the user-defined time period.
#' PleistoDist will return a value of NA if the island is submerged underwater for that particular interval.
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean calculations.
#' @return This function outputs a matrix of island surface area estimates in long format, with one column per interval.
#' @examples
#' pleistoshape_surfacearea(points="inst/extdata/AntPops.shp",epsg=3141,intervalfile="output/intervals.csv")
#' @export
pleistoshape_surfacearea <- function (points,epsg,intervalfile="output/intervals.csv") {
  intervalfile <- read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)
  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
  } else {
    p1_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing island perimeter
  islandsurfacearea <- tibble::tibble(
    Island = p1_names
  )

  #for each interval...
  for (i in 0:numintervals) {
    #load interval shapefile
    invector <- sf::st_read(paste("output/shapefile/interval",i,".shp",sep=""))
    inraster <- terra::rast(paste("output/raster_topo/interval",i,".tif",sep=""))-intervalfile$MeanDepth[i+1]
    intvlname <- paste("interval",i,sep="")
    isa <- c() #create blank list for temporarily storing surface area values (before offloading into the main dataframe)

    #start of first order for-loop to select island 1
    for (f in 1:nrow(points_transformed)) {

      #extract coordinates for island
      p1 <- points_transformed$geometry[f]
      message("Calculating island surface area for point ",p1_names[f],sep="")

      if (nrow(invector[p1, ,op=sf::st_contains]) == 0) { #check to see if island is submerged under water
        isa <- c(isa,0)
        message("Island is underwater during this interval, writing a value of 0")
      } else {
        islandmask <- terra::vect(invector[p1, ,op=sf::st_contains]) #use shapefile to create a mask layer
        islandtopo <- terra::mask(inraster,islandmask) #mask topo raster to isolate single island
        islandtopo_cropped <- terra::crop(islandtopo,islandmask) #crop raster extent to single island
        islandtopo_cropped_sp <- sf::as(raster::raster(islandtopo_cropped),"SpatialPixelsDataFrame") #convert spatraster to sp format
        isurfacearea <- sp::surfaceArea(islandtopo_cropped_sp) #calculate surface area (incorporating elevation)
        isa <- c(isa,isurfacearea)
        message(paste("Island surface area is ",isurfacearea,"m^2",sep=""))
      }
    }
    islandsurfacearea[intvlname] <- isa
  }
  islandsurfacearea$mean <- apply(islandsurfacearea[2:ncol(islandsurfacearea)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(islandsurfacearea,"output/island_surfacearea.csv")
}

#' Estimate net inter-island migration over time
#'
#' This function makes use of the equation from chapter 6 of MacArthur & Wilson (1967) that estimates the number of propagules from a
#' source island arriving at a receiving island as a function of inter-island distance, area, and relative width. Although estimating
#' the unidirectional rate of migration from one island to another requires several other coefficients (alpha: the coefficient of population density,
#' and lambda: the coefficient of overwater dispersal), calculating the ratio of migration rate for bidirectional migration causes
#' these coefficients to cancel out, allowing for net migration rate to be estimated using purely geomorphological and distance-based
#' metrics. The net migration rate can therefore be calculated using the equation (atan(t2/(2*d))*A1)/(atan(t1/(2*d))*A2), where t is
#' the island width relative to the other island, d is the inter-island distance, and A is the island area, for islands 1 and 2. If the
#' net migration rate is greater than 1, then the net movement of propagules occurs from island 1 to island 2. If the net migration rate
#' is less than 1, then the net movement of propagules occurs from island 2 to island 1.
#' @param points A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
#' labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
#' output distance matrix will use the attribute FID values instead.
#' @param epsg The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
#' projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
#' Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
#' associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
#' in distance matrices calculated in decimal degrees rather than in distance units.
#' @param disttype The type of inter-island distance matrix to use for this calculation, either "centroid", "leastshore", or "meanshore".
#' @param intervalfile This is the master control file generated using either the getintervals_time() or getintervals_sealvl()
#' function that defines the number of intervals, the sea level at each interval, and the duration of each interval. By default, this
#' function will use the "intervals.csv" file stored in the output folder, but users can also specify their own custom interval
#' file (with nice round mean sea level values, for example), although users need to ensure that the same column names are preserved, and
#' be aware that custom interval files may lead to inaccurate weighted mean calculations.
#' @return This function outputs a matrix of island surface area estimates in long format, with one column per interval.
#' @examples
#' pleistodist_netmig(points="inst/extdata/AntPops.shp",epsg=3141,disttype="centroid",intervalfile="output/intervals.csv")
#' pleistodist_netmig(points="inst/extdata/AntPops.shp",epsg=3141,disttype="meanshore",intervalfile="output/intervals.csv")
#' @export
pleistodist_netmig <- function(points,epsg,disttype,intervalfile="output/intervals.csv") {
  #check to make sure a valid disttype is supplied
  if (disttype %in% c("centroid","leastshore","meanshore")==FALSE) {
    stop("Please enter a valid distance matrix type (centroid, leastshore, or meanshore)")
  }
  #parse disttype input, reading the specified distance matrix file from the output folder
  if (file.exists(base::paste0("output/island_",disttype,"dist.csv"))==TRUE) {
    distfile = read.csv(base::paste0("output/island_",disttype,"dist.csv"))
    message(base::paste0("Reading distance matrix: output/island_",disttype,"dist.csv"))
  } else {
    #if PleistoDist doesn't detect the correct distance matrix, then call the associated function to generate the distance matrix
    message("No distance matrix file found in output folder, generating distance matrix from scratch (this might take a while...)")
    do.call(base::paste0("pleistodist_",disttype),list(intervalfile,points,epsg))
    distfile <- read.csv(base::paste0("output/island_",disttype,"dist.csv"))
  }
  widthfile <- read.csv("output/island_relativewidth.csv") #read relative island widths
  if ((disttype != "euclidean") && ncol(widthfile)<=4) {
    pleistodist_relativewidth(intervalfile,points,epsg)
    widthfile <- read.csv(base::paste0("output/island_relativewidth.csv"))
    message("Regenerating relative island width file...")
  }
  if (file.exists("output/island_area.csv")==TRUE) {
    areafile <- read.csv("output/island_area.csv")
  } else {
    pleistoshape_area(intervalfile,points,epsg)
    areafile <- read.csv("output/island_area.csv")
    message("No island area file found in output folder, generating island areas from scratch")
  }

  intervalfile <- read.csv(intervalfile)
  points = sf::st_read(points,fid_column_name="FID")

  #check projection on input points, set to WGS84 (EPSG:4326) if undefined
  if (is.na(sf::st_crs(points)$input)) {
    sf::st_crs(points) = 4326
  }

  #reproject input points to target projection
  points_transformed <- sf::st_transform(points,epsg)
  numintervals = max(intervalfile$Interval)

  #check to see if the points contains a column of unique names, otherwise use FID numbers as identifiers
  if ("Name" %in% colnames(points_transformed) && anyDuplicated(points_transformed$Name) == 0) {
    p1_names <- points_transformed$Name
    p2_names <- points_transformed$Name
  } else {
    p1_names <- points_transformed$FID
    p2_names <- points_transformed$FID
    message("No 'Name' column with unique identifiers detected in points file, defaulting to FID values instead")
  }

  #create dataframe for storing relative migration distance matrix
  relmig <- tibble::tibble(
    Island1 = expand.grid(p1_names,p2_names,include.equals=T)[,2],
    Island2 = expand.grid(p1_names,p2_names,include.equals=T)[,1],
  )
  for (i in 0:numintervals) {
    intvlname <- paste("interval",i,sep="")
    mr <- c()
    message(base::paste0("Calculating relative migration for ",intvlname))
    for (f in 1:nrow(points_transformed)) {
      for (t in 1:nrow(points_transformed)) {
        d12 <- distfile %>% dplyr::filter(Island1==p1_names[f],Island2==p2_names[t]) %>% dplyr::select(all_of(intvlname)) %>% as.numeric()
        d21 <- distfile %>% dplyr::filter(Island1==p2_names[t],Island2==p1_names[f]) %>% dplyr::select(all_of(intvlname)) %>% as.numeric()
        d = (d12+d21)/2
        if (is.na(d12)==TRUE) {
          mr <- c(mr,NA)
        } else if (d12 == 0) {
          mr <- c(mr,1)
        } else {
          A1 <- areafile %>% dplyr::filter(Island==p1_names[f]) %>% dplyr::select(all_of(intvlname)) %>% as.numeric()
          A2 <- areafile %>% dplyr::filter(Island==p2_names[t]) %>% dplyr::select(all_of(intvlname)) %>% as.numeric()
          t1 <- widthfile %>% dplyr::filter(Island1==p1_names[f],Island2==p2_names[t]) %>% dplyr::select(all_of(intvlname)) %>% as.numeric()
          t2 <- widthfile %>% dplyr::filter(Island1==p2_names[t],Island2==p1_names[f]) %>% dplyr::select(all_of(intvlname)) %>% as.numeric()
          migratio <- (atan(t2/(2*d))*A1)/(atan(t1/(2*d))*A2)
          mr <- c(mr,migratio)
        }
      }
    }
    relmig[intvlname] <- mr
  }
  relmig$mean <- apply(relmig[3:ncol(relmig)],1,stats::weighted.mean,intervalfile$TimeInterval,na.rm=TRUE)
  write.csv(relmig,base::paste0("output/island_netmigration_",disttype,".csv"))
}
