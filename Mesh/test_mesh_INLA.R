## Test script for generating global mesh using package INLA

## If the packages is not installed, use
## install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
library(INLA)

## A simplest way to generate coarse and uniform mesh on the globe
mesh_g <- inla.mesh.create(globe = 10) # specify the number of sub-segments to use, when subdividing an icosahedron
plot(mesh_g)

## For a 3d interactive plot use rgl, if rgl not installed, run
## install.packages("rgl")
plot(mesh.g, rgl = TRUE)

## Now try to build an mesh using the Antarctic Data
## Load the map of the data
require(maptools)
AntCoast <- readShapePoly(fn = "ArcGIS_shapefiles/moa_coast_poly.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
AntGround <- readShapePoly(fn = "ArcGIS_shapefiles/moa_grounding_poly.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
AntIslands <- readShapePoly(fn = "ArcGIS_shapefiles/islands_poly.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
## Plot the map
plot(AntCoast, border = "grey") 
plot(AntGround, add = T)
plot(AntIslands, add = T)

## Try a coarse mesh on the area defined by coastline
AC.bdry <- inla.sp2segment(AntCoast) #get the boundary from the polygon object
mesh_A <- inla.mesh.2d(boundary = AC.bdry, cutoff = 1e5, max.edge = c(1e5, 5e5))
summary(mesh_A)
mesh_A$n # number of triangle cells
plot(mesh_A)

## Try iteratively generate adaptive mesh -- dense at the coast and sparse otherwise
## Got the idea and help from https://groups.google.com/forum/#!topic/r-inla-discussion-group/b3uinntrh9E
## An example from the tutorial p16 of https://drive.google.com/drive/folders/0B55S_W3X6C38WS13T1NtTVlvblU
## load the map of norway
library(maps); library(maptools); library(mapdata)
data("worldHiresMapEnv")
nor_coast_poly <- map("worldHires", "norway", fill=TRUE, col="transparent",
                      plot=FALSE, ylim=c(58,72))
IDs <- sapply(strsplit(nor_coast_poly$names, ":"), function(x) x[1])
nor_coast_poly_sp <- map2SpatialPolygons(nor_coast_poly, IDs=IDs,
                                         proj4string=CRS("+proj=longlat +datum=WGS84"))
plot(nor_coast_poly_sp)

## select a small area within norway along the coast
loc3 = matrix(c(4,60.6, 8,60.6, 4,61.5, 8,61.5), 4, 2, byrow = T)
max.edge.length = 0.04
square = bakka.square.polygon(xlim=range(loc3[ ,1]), ylim=range(loc3[ ,2]),
                              ret.SP = T)
square@proj4string=nor_coast_poly_sp@proj4string

## The coastline and the land
poly3 = gIntersection(nor_coast_poly_sp, square)
plot(poly3, xlim=range(loc3[ ,1]), ylim=range(loc3[ ,2]), axes=T) 

## In this example, we are interested in generating an denser grid for water, so we get the complimentary of poly3
poly3.water = gDifference(square, poly3)
plot(poly3.water)

## Step 1 build the dense mesh for water
poly3.water.seg <- inla.sp2segment(poly3.water)
mesh3 = inla.mesh.2d(boundary = poly3.water.seg, max.e = max.edge.length)
plot(mesh3)

## Step 2 build convex mesh along wrapping around the water -- so that
## we have small triagles at the coastlines and there is no "corners" (that can affect the finite element method)
mesh3 = inla.mesh.2d(loc=mesh3$loc, max.edge = max.edge.length*2,
                     offset=max.edge.length*2)
plot(mesh3, main="")
mesh3$n

## Step 3 build the outer mesh extension to correct the boundary effect of SPDE
mesh3 = inla.mesh.2d(loc=mesh3$loc, max.edge = max.edge.length*4, offset=0.5)
plot(mesh3, main="")
mesh3$n

## Next plot the grid and differentiating between water and land
mesh3 = DT.mesh.addon.posTri(mesh3) # Add central points to the generated triangles 
points = SpatialPoints(mesh3$posTri) # Get the central points
points@proj4string = poly3@proj4string # Set the projection of these points to be the same as the land polygon
barrier3 = over(nor_coast_poly_sp, points, returnList=T) # Find points on land
barrier3 = unlist(barrier3) 
Omega3 = DT.Omega(list(barrier3, 1:mesh3$t), mesh3) # Create a union domain from lists of subdomains
Omega.SP3 = DT.polygon.omega(mesh3, Omega3) # Construct spatial polygon for the subdomains

plot(mesh3, main="")
plot(Omega.SP3[[1]], add=T, col='grey')
plot(Omega.SP3[[2]], add=T, col='lightblue')
plot(mesh3, add=T)
points(loc3, col='black')

## We can try a similar approach to the Antarctica map
## Step 1 build a sparse mesh within the coastline
## Step 1 build the dense mesh for water
mesh_A <- inla.mesh.2d(boundary = AC.bdry, cutoff = 1e4, max.edge = 5e5)
plot(mesh_A, main="")


## Step 2 build convex mesh along wrapping around the water -- so that
## we have small triagles at the coastlines and there is no "corners" (that can affect the finite element method)
mesh_A = inla.mesh.2d(loc=mesh_A$loc, max.edge = 4e5,
                     offset=4e5)
plot(mesh_A, main="")
mesh3$n

## Step 3 build the outer mesh extension to correct the boundary effect of SPDE
mesh_A = inla.mesh.2d(loc=mesh_A$loc, max.edge = 8e5, offset=1e6)
plot(mesh_A, main="")
mesh3$n


