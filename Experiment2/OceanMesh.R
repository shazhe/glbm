#### Generate mesh over the global ocean with the land being boudary and exterior 
#### Based on example from barrier models.
#### 0. Load packages and set working directory
library(INLA)
library(rgdal)
library(maptools)
library(GEOmap)
library(rgl)

lands <- readOGR(dsn = "C:/ZSwork/glbm/Experiment2/ne_110m_land", layer = "ne_110m_land")
#landpoly <- lands@polygons[[1]]@Polygons
#landpolyL <- lapply(1:length(landpoly), function(x) Polygons(landpoly[x], x)) 
#spLand <- SpatialPolygons(landpolyL)
spLand <- SpatialPolygons(lands@polygons)
plot(spLand, col = "grey")
summary(spLand)

landsareas <- sapply(spLand@polygons, function(x) slot(x, "area"))
Llands <- landsareas > 20
spLand2 <- spLand[Llands]
plot(spLand2, col ="grey")

#### 3. Get the sampling points as starting nodes for inla mesh
## use the points defined the coastlines as the sampling points
## Generate a boundary object from the polygons and extract the points defining the boundary
B.bdry <- inla.sp2segment(spLand2) # get the boundary from the polygon object
B.points <- B.bdry$loc # get the coordinates of the coastline points (long lat)
B.xyz <- do.call(cbind, Lll2xyz(B.points[,2], B.points[,1])) # project longlat back to cartesian xyz
B.bdry$loc <- B.xyz
plot3d(B.xyz) # visualise these points in 3d

#### 4. Generate the mesh on the globe by using the sampling points
meshR <- inla.mesh.create(globe = 60)
MeshB <- inla.mesh.2d(loc = meshR$loc, interior = B.bdry, cutoff = 0.02, max.edge = 0.02)
MeshB2 <- inla.mesh.2d(loc = rbind(meshR$loc, B.xyz), cutoff = 0.05, max.edge = 0.05)
plot(MeshB, rgl = T)
summary(MeshB)



##### Smooth the coastline by moving average


### add 


