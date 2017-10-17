#### Generate mesh one the sphere except in one polygon
#### Use the world coase line
library(INLA)
library(rgdal)
library(maptools)
library(GEOmap)
library(rgl)

### 1 INLA mesh -- boundary and interior
loc.bnd = matrix(c(0,0, 1,0, 1,1, 0,1), 4, 2, byrow=TRUE)
loc.int = matrix(c(0.3,0.3, 0.8,0.4, 0.7,0.9, 0.2, 0.8, 0.3,0.3), 5, 2, byrow=TRUE)
segm.bnd = inla.mesh.segment(loc.bnd)
segm.int = inla.mesh.segment(loc.int, is.bnd=FALSE)
int.Poly <- SpatialPolygons(list(Polygons(list(Polygon(loc.int)), ID = "1")))


## Points to be meshed
loc = matrix(runif(10*2),10,2)*0.9+0.05
mesh = inla.mesh.2d(loc = loc, boundary = segm.bnd, interior = segm.int, max.edge = 0.05)
# Note that when using interior here offset and max.edge can only be scalars.
# Bondary defines the boundary of the domain of interest
# interior defines a set of segments inside the domain of interest that are desired to the triagle edges.
plot(mesh)
## From this we can build functions to separate polygons and the holes in it.

### 2 meshes for separate domains
## Use the above example we try to build continuous mesh with different shapes within the two domains
## Learn from the different terrain/ barrier model examples
## Idea: Find the triangles inside the polygons and mask them from the mesh when building up the spde model

## Find which triangles are inside the polygon -- centroid incide the polygon
# Find all the centroid of all triangles from the mesh
mesh_center <- function(mesh){ # this is exactly the same as Bakka's function "DT.mesh.addon.posTri()"
  mesh$nt <- nrow(mesh$graph$tv)
  mesh$t_center <- matrix(0, nrow = mesh$nt, ncol = 3)
  for (i in 1:mesh$nt){
    t_coord <- mesh$loc[mesh$graph$tv[i,],]
    mesh$t_center[i,] <- colMeans(t_coord)
  }
  return(mesh)
}

mesh2 <- mesh_center(mesh)
insidePoly <- over(int.Poly, SpatialPoints(mesh2$t_center[,1:2]), returnList=T)
plot(mesh)
points(mesh2$t_center[unlist(insidePoly),1:2], col = "blue")








# - checking which mesh triangles are inside the barrier area
barrier = unlist(barrier)
Omega = DT.Omega(list(barrier, 1:mesh$t), mesh)
Omega.SP = DT.polygon.omega(mesh, Omega)
# - creates polygons for the different areas: Barrier area and Normal area


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
AusPoly <- spLand2[4]



# for(i in 1:19){
#   spLand2@polygons[[i]]@ID <- "1"
# }

#### 3. Get the sampling points as starting nodes for inla mesh
## use the points defined the coastlines as the sampling points
## Generate a boundary object from the polygons and extract the points defining the boundary
B.bdry <- inla.sp2segment(AusPoly) # get the boundary from the polygon object
B.points <- B.bdry$loc # get the coordinates of the coastline points (long lat)
B.xyz <- do.call(cbind, Lll2xyz(B.points[,2], B.points[,1])) # project longlat back to cartesian xyz
B.bdry$loc <- B.xyz
plot3d(B.xyz) # visualise these points in 3d

#### 4. Generate the mesh on the globe by using the sampling points
meshR <- inla.mesh.create(globe = 40)
MeshB <- inla.mesh.2d(loc = meshR$loc, boundary = B.bdry, cutoff = 0.02, max.edge = 0.02)
MeshB2 <- inla.mesh.2d(loc = rbind(meshR$loc, B.xyz), cutoff = 0.05, max.edge = 0.05)
plot(MeshB, rgl = T)
summary(MeshB)



##### Smooth the coastline by moving average


### add 


