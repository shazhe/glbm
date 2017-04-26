###################################################################
##          A Bayesian Gaussian Update                           ##
## This script performs a full Bayesian Gaussian Update          ##
## of the GIA processes                                          ##
## Given                                                         ##
## 1. prior assumptions of the parameters                        ##
## 2. subset of GPS obs                                          ##
## 3. GIA forward model priors                                   ##
###################################################################
#### INLA Approximation
library(rgl)
library(GEOmap)
library(INLA)
#INLA:::inla.dynload.workaround() 
library(MVST)
library(fields)
#load("experimentBHM/Mesh_GIA.RData")
#GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/Mesh_GIA.RData")
GPS_obs <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obsU$y_center, lon = GPS_obsU$x_center))

## Find the mapping between observations and processes basis
CMat <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
GIA_spde <- inla.spde2.matern(Mesh_GIA)
y <- as.vector(GPS_obsU$trend)
