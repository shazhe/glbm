### GIA paper



#### Illstration of subset model

library(INLA); library(rgdal); library(maptools); library(GEOmap); library(rgl)
source("~/glbm/BHM_sphere/functions-barriers-dt-models-march2017.R")
source("~/glbm/BHM_sphere/functions.R")
set.seed(18)

# boundaries and interiors
loc.bnd = matrix(c(0,0, 1,0, 1,1, 0,1), 4, 2, byrow=TRUE)
loc.int = matrix(c(0.3,0.3, 0.8,0.3, 0.8,0.5, 0.3,0.5, 0.3,0.3), 5, 2, byrow=TRUE)
segm.bnd = inla.mesh.segment(loc.bnd)
segm.int = inla.mesh.segment(loc.int, is.bnd=FALSE)
int.Poly <- SpatialPolygons(list(Polygons(list(Polygon(loc.int)), ID = "in")))
out.Poly <- SpatialPolygons(list(Polygons(list(Polygon(loc.bnd)), ID = "out")))

# first create uniformaly dense points
loc0 <- spsample(out.Poly, n = 1000, type = "hexagonal", offset = c(0, 0))
## mesh for the entire region
mesh <- inla.mesh.2d(loc = loc0,  boundary = segm.bnd, cutoff = 0.02, max.edge = 0.5)
## separate the mesh
mesh = dt.mesh.addon.posTri(mesh)
TinP = unlist(over(int.Poly, SpatialPoints(mesh$posTri), returnList=T))
Omega = dt.Omega(list(TinP, 1:mesh$t), mesh)
Omega.SP = dt.polygon.omega(mesh, Omega)
plot(mesh, main ="Mesh and Omega", asp = 1)
plot(Omega.SP[[1]], add=T, col='grey')
plot(Omega.SP[[2]], add=T, col='lightblue')
plot(mesh, add=T)

## get the subset mesh
mesh_sub <- mesh.sub(mesh, Omega, 2, sphere = FALSE)
plot(mesh_sub)

spde1 <- inla.spde2.matern(mesh_sub, B.tau = matrix(c(-5, -1, 1),1,3), B.kappa = matrix(c(2, 0, -1), 1,3),
                  theta.prior.mean = c(0,0), theta.prior.prec = c(1,1))
Q <- inla.spde.precision(spde1, c(0,0))

data_loc <- matrix(c(0.5, 0.55,  0.5, 0.25, 0.8, 0.55), byrow = TRUE, ncol = 2)
Amat <- inla.spde.make.A(mesh = mesh_sub, loc = data_loc)
Sigma_data <- solve(Q)
