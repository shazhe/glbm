###################################################################
##                    Gaussian Update                            ##
## This script performs a Gaussian Update of the GIA processes   ##
## Given                                                         ##
## 1. fixed parameter values                                     ##
## 2. subset of GPS obs                                          ##
## 3. GIA forward model prior                                    ##
###################################################################
library(rgl)
library(GEOmap)
library(INLA)
library(MVST)
load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/Mesh_GIA.RData")
GPS_obs <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GPS_synthetic.txt", header = T)
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obs$y_center, lon = GPS_obs$x_center))
CMat <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)

GIA_GMRF<- GMRF(mu = GIA_mu, Q = Q_GIA) 
GPS_obsDF <- Obs(df=data.frame(x = GPS_loc[,1], y = GPS_loc[,2], w = GPS_loc[,3], z = GPS_obs$trend, 
                          std = GPS_obs$std, t = rep(0, nrow(GPS_loc))))
L <- link(GIA_GMRF, GPS_obsDF, Cmat = CMat)
e <- new("link_list",list(L))
v <- new("block_list",list(G1 = GIA_GMRF, O = GPS_obsDF))
G <- new("Graph",e = e,v = v)
G_reduced <- compress(G)
Results <- Infer(G_reduced)

GIA_post <- Results$Post_GMRF@rep
GIA_mpost<- GIA_post$x_mean
## Get the posterior mean and variance
GIA_spost <- sqrt(GIA_post$x_margvar)
GIA_diff <- GIA_mpost - GIA_mu
## Plot the GIA prior
vals <- GIA_mpost
open3d()
par3d(windowRect = c(100, 100, 900, 900))
layout3d(matrix(1:2, 2,1), heights = c(4,2), sharedMouse = TRUE)
par3d(zoom = 0.8)
t_lim <- range(vals)
t_Clens <- round((t_lim[2] - t_lim[1])*100) + 1
t_Cpal <- topo.colors(t_Clens, alpha=0)
t_Cols<- t_Cpal[round((vals - t_lim[1])*100) + 1]
plot3d(GIA_Mloc, col= t_Cols)

## plot the color bar
next3d(reuse = FALSE)
bgplot3d({z=matrix(1:t_Clens, nrow = t_Clens)
y=1
x=seq(t_lim[1],t_lim[2],len = t_Clens)
par(cex = 1.5, fin = c(8, 1), mai = c(0,0, 0.5, 0), oma = c(1, 0, 0, 0))
image(x,y,z,col = topo.colors(t_Clens),axes=FALSE,xlab="",ylab="")
title("The ICE6G GIA")
axis(1)})

writeWebGL(dir = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a", 
           filename= "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/GIA_globe.html",  # save the plot as an html
           width = 1000, reuse = TRUE)

