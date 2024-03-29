% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcmetrics.R
\name{pleistodist_meanshore}
\alias{pleistodist_meanshore}
\title{Calculate inter-island mean shore-to-shore distance over time}
\usage{
pleistodist_meanshore(
  points,
  epsg,
  intervalfile,
  mapdir,
  outdir,
  maxsamp = 1000
)
}
\arguments{
\item{points}{A user-generated multi-point shapefile (.shp) containing at least two points. If the shapefile attribute table contains a column
labelled 'Name', the output distance matrix will use the identifiers in this column to label each pairwise comparison. Otherwise, the
output distance matrix will use the attribute FID values instead.}

\item{epsg}{The projected coordinate system in EPSG code format. Because of the curvature of the Earth's surface, we need to apply a map
projection to accurately calculate straight-line distances between points instead of using the default WGS84 geographical coordinate system.
Users should specify a projected coordinate system appropriate to the geographic region being analysed using the projection's
associated EPSG code (https://epsg.org/home). Geographic coordinate system projections are not recommended as those will result
in distance matrices calculated in decimal degrees rather than in distance units.}

\item{intervalfile}{This is the master control file generated using either the getintervals_time() or getintervals_sealvl() functions, or custom generated by the user,
that defines the number of intervals, the sea level at each interval, and the duration of each interval.}

\item{mapdir}{The directory containing map outputs from the makemaps() function
(i.e. the directory containing the raster_flat, raster_topo, and shapefile folders).}

\item{outdir}{Output directory for the inter-island mean shore-to-shore distance matrix file.
If the specified output directory doesn't already exist, PleistoDist will create the output directory.}

\item{maxsamp}{Maximum number of points to sample along the donor island shoreline, mainly to limit computational time. If users are
ok with long runtimes, then maxsamp can be set to NA. Default number of sample points is 1000.}
}
\value{
This function outputs an asymmetric pairwise distance matrix of mean shore-to-shore distances in long format, with
one column per interval. Because this distance matrix is asymmetric, the directionality of the distance calculation matters.
}
\description{
This function calculates the mean shore-to-shore distance between island landmasses over Pleistocene time. This function is slightly different
from other PleistoDist functions because the resultant distance matrix is asymmetrical. This is due to the fact that the mean shore-to-shore distance from
a large island to a small one will be greater relative to the reverse, since larger islands have longer and more complex shorelines, so the mean dispersal distance
from any point along the larger island's shoreline facing the smaller island is likely to be greater than in the reverse case. To calculate the mean
shore-to-shore distance, PleistoDist generates a number of points at equal intervals (1 point every 100m by default) along the shoreline of the source island at each time/sea level interval, filters out points not directly
facing the receiving island, and calculates the average minimum distance between each point along the source island shoreline and the shoreline of the
receiving island. Because the number of shoreline points can run into the thousands for very large islands, thereby increasing computational time, users
can set a cap on the maximum number of shoreline points to sample (the default cap is 1000 points). As with other PleistoDist functions, if both points lie on the same landmass, a distance value of 0 is returned. If one or both islands
is/are underwater, then a value of NA is returned. After calculating the pairwise island distances for each interval, this function calculates
the time-weighted average mean shore-to-shore distance between each island pair.
}
\examples{
\donttest{
#create temp directory
path <- file.path(tempdir())
#create points shapefile
points <- sf::st_multipoint(rbind(c(178.68408,-16.82573), c(179.13585,-16.59487)))
#convert points to feature geometry
points <- sf::st_sfc(points)
points <- sf::st_cast(points, "POINT")
#set default projection (WGS84)
sf::st_crs(points) <- 4326
#save shapefile
sf::write_sf(points,paste0(path,"/points.shp"))
#load bathymetry file
fiji <- system.file("extdata","FJ.asc",package="PleistoDist")
#generate interval file for 1 interval and 20 kya cutoff time, binning by time
getintervals_time(time=20,intervals=1,outdir=path)
#generate maps based on interval file, projecting map using EPSG:3141
makemaps(inputraster=fiji,epsg=3141,intervalfile=paste0(path,"/intervals.csv"),outdir=path)
#calculate inter-island mean shore-to-shore distances, projecting points using EPSG:3141,
#with 500 sampling points.
pleistodist_meanshore(points=paste0(path,"/points.shp"),epsg=3141,
    intervalfile=paste0(path,"/intervals.csv"),mapdir=path,outdir=path,maxsamp=500)
}
}
