##################################################################################
##       Test script for Adaptive mesh on the globe using INLA                  ##
## We try to generate the mesh on a low resolution map of global coastlines.    ##
##                                                                              ##
## ***DATA***                                                                   ##
## * We use the GSHHG downloaded from:                                          ##
##   http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/                 ##
## * Current version of the data: Version 2.3.5 April 12, 2016.                 ##
## * For this script we used the crude resolution (C).                          ##
## * All shape files are stored under "shapefiles/GSHHS/"                       ##
## * CRS: geodetic longitude, latitude, locations on the WGS-84 ellipsoid       ##
##                                                                              ##
## ***Meshing procedure outline***                                              ##
## *1. Load shapefiles                                                          ##
## *2. Group polygons                                                           ##
## *3. Mesh within the polygons                                                 ##
## *4. Get the locations from the mesh                                          ##
## *5. Convert the longlat location to Cartesian                                ##
## *6  Use the location coordinates to mesh on the globe                        ##
##################################################################################

library(INLA)
library(rgdal)
library(maptools)
library(GEOmap)
library(rgl)

#### 1. Loadshape files
## Continents and ocean islands
COIlines <- readOGR(dsn = "shapefiles/GSHHS", layer = "GSHHS_c_L1")
plot(COIlines, col = "grey")
summary(COIlines)

## Remove the small ocean islands, say with area smaller than 20
COIareas <- sapply(slot(COIlines, "polygons"), function(x) slot(x, "area"))
Clands <- COIareas >  20
COIlines2 <- COIlines[Clands,]
plot(COIlines2)

## Antarctica grounding line
AnGlines <- readOGR(dsn= "Shapefiles/GSHHS", layer = "GSHHS_c_L6")
plot(AnGlines)
summary(AnGlines)

## check and remove small islands
AnGareas <- sapply(slot(AnGlines, "polygons"), function(x) slot(x, "area"))
AnGland <- AnGareas > 20
AnGlines2 <- AnGlines[AnGland,]
plot(AnGlines2)

## Plot Both maps
plot(COIlines2)
plot(AnGlines2, add =T)

## Conbine the two SpatialPolygonDataFrame objs
row.names(COIlines2) <- paste("C", row.names(COIlines2), sep="_")
row.names(AnGlines2) <- paste("A", row.names(AnGlines2), sep="_")
BLand <- spRbind(COIlines2, AnGlines2)
plot(BLand)

## Get the sampling points along the coastlines using the inla segment function
## Try a coarse mesh on the area defined by coastline
B.bdry <- inla.sp2segment(BLand) # get the boundary from the polygon object
B.points <- B.bdry$loc # get the coordinates of the coastline points (long lat)
B.xyz <- do.call(cbind, Lll2xyz(B.points[,2], B.points[,1]))
plot3d(B.xyz)

MeshB <- inla.mesh.2d(loc = B.xyz, cutoff = 0.005, max.edge = 0.2)
plot(MeshB, rgl = TRUE)
plot3d(B.xyz, type = "l", add = TRUE, lwd = 10, col = "red")

## Try to plot the lines of the polygons 
## Function for extract coordinates of the polygon and converted it to 3d xyz
poly2xyz <- function(polygs){
  coordsll <- slot(polygs, "Polygons")[[1]]@coords
  coordsxyz <- do.call(cbind, Lll2xyz(coordsll[,2], coordsll[,1]))
}
## Function for ploting the 3d coastlines from a SpatialPolygonDataFrame
plot.polys <- function(Polys){
  poly.list <- slot(Polys, "polygons")
  poly.coords <- lapply(poly.list, poly2xyz)
  lapply(poly.coords, plot3d, type="l", col = "red", lwd =4, add = TRUE)
}

plot(MeshB, rgl = TRUE)
plot.polys(BLand)


summary(MeshB)
MeshB$n # number of triangle cells
