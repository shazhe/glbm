##### This script shows the effects of mesh on building the spde precision #####
library(rgl)
library(INLA)
library(fields)
wkdir <- getwd()
source("glbm/Experiment1a/Rscript/MVSTplus.R")

## Priors mean and variance for rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

## Setup for spde parameters
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]

lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

#### load GPS data
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
GPSX <- ifelse(GPS_obs$x_center > 180, GPS_obs$x_center-360, GPS_obs$x_center)
GPSY <- GPS_obs$y_center

GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obs$y_center, lon = GPS_obs$x_center))
#### load GIA prior data
GIA_ice6g <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GIA_truth.txt", header = T)
polycoords <- GIA_ice6g[,c(6:13, 6,7)]
plist <- lapply(GIA_ice6g$ID, function(x) Polygons(list(Polygon(cbind(lon = as.numeric(polycoords[x, c(1,3,5,7,9)]), 
                                                                      lat = as.numeric(polycoords[x, c(2,4,6,8,10)])))), ID = x))

Plist <- SpatialPolygons(plist, proj4string = CRS("+proj=longlat"))
GIA_ice6g_sp <- SpatialPolygonsDataFrame(Plist, data = GIA_ice6g[,c("trend", "std")])

## Treat the gridded data as point data sampled at the centre of the grid
## converte the LonLat to xyz coordinates
GIA_loc6 <- do.call(cbind, Lll2xyz(lat = GIA_ice6g$y_center, lon = GIA_ice6g$x_center))
GIA_mu6 <- GIA_ice6g$trend
GIA_sd6 <- GIA_ice6g$std

#### 1 Genreate a mesh by from the given GIA grid
Mesh_GIAs <- list()
Mesh_GIAs[[1]] <- inla.mesh.2d(loc = GIA_loc6, cutoff = 0.1, max.edge = 0.1) ## small
Mesh_GIAs[[2]] <- inla.mesh.2d(loc = GIA_loc6, cutoff = 0.045, max.edge = 0.1) ## Medium
Mesh_GIAs[[3]] <- inla.mesh.2d(loc = GIA_loc6, cutoff = 0.023, max.edge = 0.1) ## Large

pdf(file = paste0(getwd(), "/experimentBHM/prior_GIAmesh.pdf"), width = 8, height = 11)
mains <- paste("Prior marginal standard error -- ",  c("small mesh", "medium mesh", "large mesh"))
par(mfrow = c(3,1))
for (i in 1:3){
GIA_spde <- inla.spde2.matern(Mesh_GIAs[[i]], B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

Q00 <- inla.spde2.precision(GIA_spde, theta = c(0,0))
GIA_sprior <- sqrt(1/diag(Q00))
proj<- inla.mesh.projector(Mesh_GIAs[[i]], projection = "longlat", dims = c(361,181))
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_sprior)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = mains[i])
points(GPSX, GPSY, pch = 20)
}
dev.off()


#### 2 Genreate a mesh from a globe with equidistant points from equidistant lattice
mesh_regs <- list()
mesh_regs[[1]] <- inla.mesh.create(globe = 10)
mesh_regs[[2]] <- inla.mesh.create(globe = 30)
mesh_regs[[3]] <- inla.mesh.create(globe = 60)


pdf(file = paste0(getwd(), "/experimentBHM/prior_Globemesh.pdf"), width = 8, height = 11)
mains <- paste("Prior marginal standard error -- ",  c("small mesh", "medium mesh", "large mesh"))
par(mfrow = c(3,1))
for (i in 1:3){
  GIA_spde <- inla.spde2.matern(mesh_regs[[i]], B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                                theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))
  
  Q00 <- inla.spde2.precision(GIA_spde, theta = c(0,0))
  GIA_sprior <- sqrt(1/diag(Q00))
  proj<- inla.mesh.projector(mesh_regs[[i]], projection = "longlat", dims = c(361,181))
  image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_sprior)), col = topo.colors(40),
             xlab = "Longitude", ylab = "Latitude", main = mains[i])
  points(GPSX, GPSY, pch = 20)
}
dev.off()


#### 3 Genreate a mesh from GPS location
mesh_GPS <- list()
mesh_GPS[[1]] <- inla.mesh.2d(loc = GPS_loc, cutoff = 0.001,  max.edge = 0.1)
mesh_GPS[[2]] <- inla.mesh.2d(loc = GPS_loc, cutoff = 0.001, max.edge = 0.08)
mesh_GPS[[3]] <- inla.mesh.2d(loc = GPS_loc, cutoff = 0.001, max.edge = 0.04)
GPS_mu <- list()
pdf(file = paste0(getwd(), "/experimentBHM/prior_GPSmesh.pdf"), width = 8, height = 11)
mains <- paste("Prior marginal standard error -- ",  c("small mesh", "medium mesh", "large mesh"))
par(mfrow = c(3,1))
for (i in 1:3){
  GIA_spde <- inla.spde2.matern(mesh_GPS[[i]], B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                                theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))
  
  Q00 <- inla.spde2.precision(GIA_spde, theta = c(0,0))
  GIA_sprior <- sqrt(1/diag(Q00))
  proj<- inla.mesh.projector(mesh_GPS[[i]], projection = "longlat", dims = c(361,181))
  image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_sprior)), col = topo.colors(40),
             xlab = "Longitude", ylab = "Latitude", main = mains[i])
  points(GPSX, GPSY, pch = 20)
}
dev.off()

# # ## Plot the results 3D
# vals <- GIA_sprior
# open3d()
# par3d(windowRect = c(100, 100, 900, 900))
# layout3d(matrix(1:2, 2,1), heights = c(4,2), sharedMouse = TRUE)
# par3d(zoom = 0.8)
# t_lim <- range(vals)
# t_Clens <- round((t_lim[2] - t_lim[1])*100) + 1
# t_Cpal <- topo.colors(t_Clens, alpha=0)
# t_Cols<- t_Cpal[round((vals - t_lim[1])*100) + 1]
# plot3d(GPS_loc, pch = "20", size = 5, xlab = "", ylab = "", zlab = "", axe=FALSE)
# plot(mesh_GPS[[3]], rgl = TRUE, col= t_Cols, edge.color = rgb(0, 0.5, 0.6, alpha =0), add = TRUE)
# 
# ## plot the color bar
# next3d(reuse = FALSE)
# bgplot3d({z=matrix(1:t_Clens, nrow = t_Clens)
# y=1
# x=seq(t_lim[1],t_lim[2],len = t_Clens)
# par(cex = 1.5, fin = c(8, 1), mai = c(0,0, 0.5, 0), oma = c(1, 0, 0, 0))
# image(x,y,z,col = topo.colors(t_Clens),axes=FALSE,xlab="",ylab="")
# title("INLA Posterior ICE6G GIA, mu = 0")
# axis(1)})
