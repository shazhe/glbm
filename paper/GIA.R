### GIA paper

#### Illstration of subset model

library(INLA); library(rgdal); library(maptools); library(GEOmap); library(rgl)
source("C:/ZSwork/glbm/BHM_sphere/functions-barriers-dt-models-march2017.R")
source("C:/ZSwork/glbm/BHM_sphere/functions.R")
set.seed(18)

# boundaries and interiors
loc.bnd = matrix(c(0,0, 5,0, 5,5, 0,5), 4, 2, byrow=TRUE)
loc.int = matrix(c(1.5,2.3, 3.5,2.3, 3.5, 2.7, 1.5,2.7, 1.5,2.3), 5, 2, byrow=TRUE)
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

## Transform the parameters for the SPDE_GMRF approximation
lsigma0 <- log(2)
lrho0 <- log(1)
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

spde1 <- inla.spde2.matern(mesh_sub, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                  theta.prior.mean = c(0,0), theta.prior.prec = c(1,1))

spde2 <- inla.spde2.matern(mesh, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                           theta.prior.mean = c(0,0), theta.prior.prec = c(1,1))

Q1 <- inla.spde.precision(spde1, c(0,0))
Q2 <- inla.spde.precision(spde2, c(0,0))

data_loc <- matrix(c(2,3,  2, 2,  3,3), byrow = TRUE, ncol = 2)
plot(mesh_sub, asp = 1)
points(data_loc, col = 2)
Amat1 <- inla.spde.make.A(mesh = mesh_sub, loc = data_loc)
Amat2 <- inla.spde.make.A(mesh = mesh, loc = data_loc)
Sigma_data1 <- cov2cor(Amat1%*%solve(Q1, t(Amat1)))
Sigma_data2 <- cov2cor(Amat2%*%solve(Q2, t(Amat2)))
