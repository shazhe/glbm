#########################################################
## Compare INLA and MVST functions for GMRF precsion   ##
#########################################################

#### 0. Set working directory, load data and packages
setwd("O:/globalmass_code")
load("Mesh/GlobeMesh/MeshGlobe.RData")
library(INLA)
library(rgdal)
library(maptools)
library(GEOmap)
library(rgl)

library(dplyr)
library(ggplot2)
library(MVST)

#### 1. Simulate process with Matern kernel
# alpha = nu + d/2
### a. MVST
mglb_tv <- MeshB$loc # extract the triaulation vertices
mglb_S1 <- Matern(as.matrix(dist(mglb_tv)), nu = 1/2, var=1, kappa = 0.1) # create the Matern covmat
mglb_Q1 <- GMRF(Q=as(chol2inv(chol(mglb_S1)),"dgCMatrix")) # the precmat
mglb_x1 <- sample_GMRF(mglb_Q1) # simulate the true processs


### b.INLA
spde1 <- inla.spde2.matern(MeshB, alpha=2)
theta <- c(-log(4*pi*1*0.1^2)/2, log(0.1)) # c(-log(4*pi*s2x*kappa^2), log(kappa))
spde_Q <- inla.spde2.precision(spde1, theta)
spde_samp <- inla.qsample(Q=spde_Q)



#### 1. list functions in INLA for building the precision for the GMRF using MeshB
## 1.1 Construct the element for the Gaussian weight matrices
W_inla <- inla.mesh.fem(MeshB)

## 1.2
mglb_tv <- MeshB$loc # extract the triaulation vertices
mglb_sc <- mesh.centroid(MeshB) # extract the centroids the cells

## Simulate the process based on the veritices
mglb_S1 <- Matern(as.matrix(dist(mglb_tv)), nu = 3/2, var=4, kappa = 0.1) # create the Matern covmat
mglb_Q1 <- GMRF(Q=as(chol2inv(chol(mglb_S1)),"dgCMatrix")) # the precmat
mglb_x1 <- sample_GMRF(mglb_Q1) # simulate the true processs
