###################################################################
##                    GIA Data Screening                         ##
## This script read in and do screening of the GIA data produced ##
## from the ICE-6G forward model.                                ##
## We will use this data to:                                     ##
## 1. Generate a mesh with similar resolution                    ##
## 2. Use the data as prior mean of the GP(GMRF)                 ##
## 3. Learn initial values for lenght scale and variance         ##
###################################################################
library(rgl)
library(GEOmap)
library(INLA)
library(MVST)
#### 0 Read in GIA data
###################################################################
GIA_ice6g <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GIA_truth.txt", header = T)

## Treat the gridded data as point data sampled at the centre of the grid
## converte the LonLat to xyz coordinates
GIA_loc <- do.call(cbind, Lll2xyz(lat = GIA_ice6g$y_center, lon = GIA_ice6g$x_center))

## Plot the GIA prior
open3d()
par3d(windowRect = c(100, 100, 900, 900))
layout3d(matrix(1:2, 2,1), heights = c(4,2), sharedMouse = TRUE)
par3d(zoom = 0.8)
t_lim <- range(GIA_ice6g$trend)
t_Clens <- round((t_lim[2] - t_lim[1])*100) + 1
t_Cpal <- topo.colors(t_Clens, alpha=0)
t_Cols<- t_Cpal[round((GIA_ice6g$trend - t_lim[1])*100) + 1]
plot3d(GIA_loc, col= t_Cols)

## plot the color bar
next3d(reuse = FALSE)
bgplot3d({z=matrix(1:t_Clens, nrow = t_Clens)
y=1
x=seq(t_lim[1],t_lim[2],len = t_Clens)
par(cex = 1.5, fin = c(8, 1), mai = c(0,0, 0.5, 0), oma = c(1, 0, 0, 0))
image(x,y,z,col = topo.colors(t_Clens),axes=FALSE,xlab="",ylab="")
title("The ICE6G GIA")
axis(1)})

#### 1 Genreate a mesh by from the given points
###################################################################
## To get 1 by 1 degree resolution, we try to generate about 64800 triangles 
Mesh_GIA <- inla.mesh.2d(loc = GIA_loc, cutoff = 0.0167, max.edge = 0.038)
summary(Mesh_GIA)
plot(Mesh_GIA, rgl = T)

#### 2 Set up GIA priors
###################################################################
GIA_mu <- GIA_ice6g$trend
## Parameters needed for setting up the Matern covariance function
## nu -- mean square differentialbilty of the process -- poorly identified -- fix at 3/2
## var -- marginal variance
## kappa -- length scale -- range rho = sqrt(8nu/kappa) for correlation 0.1
## These can be determined by experts or data or both.

GIA_cov<- MVST::Matern(as.matrix(dist()), nu = 3/2, var = 4, kappa = 0.1) # create the Matern covmat

save(Mesh_GIA, GIA_mu, GIA_cov, file = "o:/glbm/Experiment1a/Mesh_GIA.RData")
