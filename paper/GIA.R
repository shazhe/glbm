### GIA paper

#### Illstration of subset model
library(INLA); library(rgdal); library(maptools); library(GEOmap); library(rgl)
source("C:/ZSwork/glbm/BHM_sphere/functions-barriers-dt-models-march2017.R")
source("C:/ZSwork/glbm/BHM_sphere/functions.R")
set.seed(18)

# boundaries and interiors

loc.bnd = matrix(c(0,0, 5,0, 5,5, 0,5), 4, 2, byrow=TRUE)
segm.bnd = inla.mesh.segment(loc.bnd)
out.Poly <- SpatialPolygons(list(Polygons(list(Polygon(loc.bnd)), ID = "out")))
loc.int <- matrix( c(0.87,0.27, 0.13,0.13, 0.26,0.24, 0.07,0.97, 0.18,0.92, 0.31,0.32), byrow=T, 6, 2)*3 + 1
int.bound = inla.nonconvex.hull(loc.int, convex=0.6)
int.Poly <- SpatialPolygons(list(Polygons(list(Polygon(int.bound$loc)), ID = "in")))
plot(out.Poly); plot(int.Poly, add = T)



# first create uniformaly dense points
loc0 <- spsample(out.Poly, n = 1000, type = "hexagonal", offset = c(0, 0))
## mesh for the entire region
mesh <- inla.mesh.2d(loc = loc0,  boundary = segm.bnd, offset = c(0.5, 0.5), cutoff = 0.1, max.edge = 0.5)
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
plot(mesh_sub, asp = 1)

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

data_loc <- matrix(c(2.2, 2.8, 1.3,2.8, 3.1,2.8), byrow = TRUE, ncol = 2)
plot(mesh_sub, asp = 1, main = "")
points(data_loc, col = 2, pch = 19)
text(data_loc, pos = 3, c("A", "B", "C"), col = 2)
Amat1 <- inla.spde.make.A(mesh = mesh_sub, loc = data_loc)
Amat2 <- inla.spde.make.A(mesh = mesh, loc = data_loc)
Sigma_data1 <- cov2cor(Amat1%*%solve(Q1, t(Amat1)))
Sigma_data2 <- cov2cor(Amat2%*%solve(Q2, t(Amat2)))

## Space partition model example of correlation between A and B
# generate points inside the polygon
loc1 <- spsample(int.Poly, n = 1000, type = "hexagonal", offset = c(0, 0))
## generate mesh
mesh2 <- inla.mesh.2d(loc = rbind(loc0,loc1),  boundary = segm.bnd, offset = c(0.5, 0.5), cutoff = 0.03, max.edge = 0.5)
mesh2 <- inla.mesh.2d(loc = mesh2$loc,  boundary = segm.bnd, offset = c(0.5, 0.5), cutoff = 0.03, max.edge = 0.5)
## separate the mesh
mesh2 = dt.mesh.addon.posTri(mesh2)
TinP = unlist(over(int.Poly, SpatialPoints(mesh2$posTri), returnList=T))
Omega2 = dt.Omega(list(TinP, 1:mesh2$t), mesh)
Omega.SP2 = dt.polygon.omega(mesh2, Omega2)
plot(mesh2, main ="Mesh and Omega", asp = 1)
plot(Omega.SP2[[1]], add=T, col='grey')
plot(Omega.SP2[[2]], add=T, col='lightblue')
plot(mesh2, add=T)



## build the spdes for the partition model
## Create the precsion matrix function
Q.function = dt.create.Q(mesh2, Omega2)
ranges = c(log(0.5), log(1))
Q3 = Q.function(theta = c(log(5), ranges))

data_loc <- matrix(c(2.2,1.9, 2.2,2.8, 3.1, 2.8, 3.1, 1.9 ), byrow = TRUE, ncol = 2)
plot(mesh2, main ="Mesh and Omega", asp = 1)
plot(Omega.SP2[[1]], add=T, col='grey')
#plot(Omega.SP2[[2]], add=T, col='lightblue')
plot(mesh2, add=T)
points(data_loc, col = 2, pch = 19)
text(data_loc, pos = 3, c("A", "B", "C", "D"), col = 2)
Amat3 <- inla.spde.make.A(mesh = mesh2, loc = data_loc)
Sigma_data3 <- cov2cor(Amat3%*%solve(Q3, t(Amat3)))
