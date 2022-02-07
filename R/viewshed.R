#' Calculates if point is visible from another fixed point, given local topography
#' Code written by Barry Rowlingson (spacedman) with modifications by David Tan, posted on
#' https://stackoverflow.com/questions/21841387/r-code-that-evaluates-line-of-sight-los-between-two-lat-lon-points
#' @keywords internal
cansee <- function(r, xy1, xy2, h1=0, h2=0){
  ### can xy1 see xy2 on DEM r?
  ### r is a DEM in same x,y, z units
  ### xy1 and xy2 are 2-length vectors of x,y coords
  ### h1 and h2 are extra height offsets
  ###  (eg top of mast, observer on a ladder etc)
  xyz = rasterprofile(r, xy1, xy2)
  np = nrow(xyz)-1
  h1 = xyz$z[1] + h1
  h2 = xyz$z[np] + h2
  hpath = h1 + (0:np)*(h2-h1)/np
  return(!any(hpath < xyz$z,na.rm=TRUE))
}

#' Calculates if matrix of points is visible from another fixed point, given local topography
#' Code written by Barry Rowlingson (spacedman), posted on
#' https://stackoverflow.com/questions/21841387/r-code-that-evaluates-line-of-sight-los-between-two-lat-lon-points
#' @keywords internal
viewTo <- function(r, xy, xy2, h1=0, h2=0, progress="none"){
  ## xy2 is a matrix of x,y coords (not a data frame)
  plyr::aaply(xy2, 1, function(d){cansee(r,xy,d,h1,h2)}, .progress=progress)
}

#' Code written by Barry Rowlingson (spacedman) with modifications by David Tan, posted on
#' https://stackoverflow.com/questions/21841387/r-code-that-evaluates-line-of-sight-los-between-two-lat-lon-points
#' @keywords internal
rasterprofile <- function(r, xy1, xy2){
  ### sample a raster along a straight line between two points
  ### try to match the sampling size to the raster resolution
  dx = sqrt( (xy1[1]-xy2[1])^2 + (xy1[2]-xy2[2])^2 )
  nsteps = 1 + base::round(dx/ min(terra::res(r)))
  xc = xy1[1] + (0:nsteps) * (xy2[1]-xy1[1])/nsteps
  yc = xy1[2] + (0:nsteps) * (xy2[2]-xy1[2])/nsteps
  zc <- r[cellFromXY(r,cbind(xc,yc))]
  colnames(zc) <- NULL
  data.frame(x=xc, y=yc, z=zc)
}
