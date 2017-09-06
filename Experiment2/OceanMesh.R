#### Generate mesh over the global ocean with the land being boudary and exterior 
#### Based on example from barrier models.
#### 0. Load packages and set working directory
library(INLA)
library(rgdal)
library(maptools)
library(GEOmap)
library(rgl)

#### 1. Loadshape files
## Continents and ocean islands
COIlines <- readOGR(dsn = "GSHHS", layer = "GSHHS_c_L1")
plot(COIlines, col = "grey")
summary(COIlines)

## Remove the small ocean islands, say with area smaller than 20
COIareas <- sapply(slot(COIlines, "polygons"), function(x) slot(x, "area"))
Clands <- COIareas >  20
COIlines2 <- COIlines[Clands,]
plot(COIlines2)

## Antarctica grounding line
AnGlines <- readOGR(dsn= "GSHHS", layer = "GSHHS_c_L6")
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
meshR <- inla.mesh.create(globe = 20)
MeshB <- inla.mesh.2d(loc = meshR$loc, interior = B.bdry, cutoff = 0.05, max.edge = 0.2)
summary(MeshB)
MeshB$n # number of triangle cells
save(BLand, B.xyz, MeshB, file = "../Results/Mesh/MeshGlobe.RData")
