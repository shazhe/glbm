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
## *2. Simplify: remove small islands                                           ##
## *3. Generate sampling points                                                 ##
## *4. Generate mesh on the globe                                               ##
## *5. Plot the results and save                                                ##
##################################################################################

#### 0. Load packages and set working directory
setwd("O:/globalmass_code/Mesh")
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

#### 2. Remove small islands
## To simplify the coastlines, remove all the islands with area smaller than 20 
AnGareas <- sapply(slot(AnGlines, "polygons"), function(x) slot(x, "area"))
AnGland <- AnGareas > 20
AnGlines2 <- AnGlines[AnGland,]
plot(AnGlines2)

## Conbine the two SpatialPolygonDataFrame objs
row.names(COIlines2) <- paste("C", row.names(COIlines2), sep="_")
row.names(AnGlines2) <- paste("A", row.names(AnGlines2), sep="_")
BLand <- spRbind(COIlines2, AnGlines2)
plot(BLand) #plot the combined map

#### 3. Get the sampling points as starting nodes for inla mesh
## use the points defined the coastlines as the sampling points
## Generate a boundary object from the polygons and extract the points defining the boundary
B.bdry <- inla.sp2segment(BLand) # get the boundary from the polygon object
B.points <- B.bdry$loc # get the coordinates of the coastline points (long lat)
B.xyz <- do.call(cbind, Lll2xyz(B.points[,2], B.points[,1])) # project longlat back to cartesian xyz
plot3d(B.xyz) # visualise these points in 3d

#### 4. Generate the mesh on the globe by using the sampling points
MeshB <- inla.mesh.2d(loc = B.xyz, cutoff = 0.005, max.edge = 0.2)
summary(MeshB)
MeshB$n # number of triangle cells

#### 5. Plot the result
plot(MeshB,rgl=TRUE) # plot the mesh grid in 3d
plot3d(B.xyz, add = TRUE, col = "red", cex = 2) # add the coastlines points
filename <- writeWebGL(dir = file.path(getwd(), "GlobeMesh"),  # save the plot as an html
                       width = 1000, reuse = TRUE)

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
  lapply(poly.coords, points3d, col = "red", cex=5, add = TRUE)
}

