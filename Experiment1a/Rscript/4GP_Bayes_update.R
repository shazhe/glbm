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
library(MVST)
load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/Mesh_GIA.RData")
GPS_obs <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obs$y_center, lon = GPS_obs$x_center))

## Find the mapping between observations and processes basis
CMat <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
## Setting pars in inla.matern
# d = 2 for S2
# alpha = nu + d/2, fractional operator, 0< alpha < 2, choose nu = 1, alpha = 1
# range -- initial set to be minimum of the mesh grid size
# log (tau) = theta1
# log (kappa) = theta2
#size0 <- min(c(diff(range(Mesh_GIA$loc[, 1])), diff(range(Mesh_GIA$loc[, 2])), diff(range(Mesh_GIA$loc[, 3]))))
range0 <- 0.0538
sigma0 <- sqrt(4.72)
kappa0 <- sqrt(8)/range0
tau0 <- 1/(4*pi*kappa0^2*sigma0^2)
GIA_spde <- inla.spde2.matern(Mesh_GIA)
GIA_Q0 <- inla.spde.precision(GIA_spde, theta=c(log(sqrt(tau0)), log(kappa0)))

x <- as.vector(GIA_mu)
y <- as.vector(GPS_obs$trend)
intc <- rep(NA, length(y))
covs <- rep(NA, length(y))
## Y = Ax + error 
mesh.index <- inla.spde.make.index(name = "field", n.spde = GIA_spde$n.spde, n.repl = 1)
st.est <- inla.stack(data = list(y=y), A = list(CMat,1), 
                     effects = list(c(list(intercept = GIA_mu), mesh.index),  list(cov = covs)), tag = "est")

st.pred <- inla.stack(data = list(y=NA), A = list(CMat), 
                     effects = list(c(list(intercept = GIA_mu), mesh.index)), tag = "pred")

stGIA <- inla.stack(st.est, st.pred)

formula = y ~ -1 + f(field, model = GIA_spde)
resultI <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde),
                family = "normal", control.predictor = list(A = inla.stack.A(stGIA), compute = TRUE))

res1 <- inla.spde2.result(resultI, "field", spde=GIA_spde)
par(mfrow = c(2,2))
plot(res1[["marginals.range.nominal"]][[1]], type = "l", main = "range")
plot(res1[["marginals.variance.nominal"]][[1]], type = "l", main = "variance")
plot(res1[["marginals.kappa"]][[1]], type = "l", main = "kappa")
plot(res1[["marginals.tau"]][[1]], type = "l", main = "tau")

GIA_mpostC <- resultI$summary.random$field$mean
GIA_mpost <- resultI$summary.random$field$mean + GIA_mu
GIA_spost <- resultI$summary.random$field$mean


## Plot the results
vals <- GIA_mpostC
open3d()
par3d(windowRect = c(100, 100, 900, 900))
layout3d(matrix(1:2, 2,1), heights = c(4,2), sharedMouse = TRUE)
par3d(zoom = 0.8)
t_lim <- range(vals)
t_Clens <- round((t_lim[2] - t_lim[1])*100) + 1
t_Cpal <- topo.colors(t_Clens, alpha=0)
t_Cols<- t_Cpal[round((vals - t_lim[1])*100) + 1]
plot3d(GPS_loc, cex = 2, xlab = "", ylab = "", zlab = "", axe=FALSE)
plot(Mesh_GIA, rgl = TRUE, col= t_Cols, edge.color = rgb(0, 0.5, 0.6, alpha =0), add = TRUE)

## plot the color bar
next3d(reuse = FALSE)
bgplot3d({z=matrix(1:t_Clens, nrow = t_Clens)
y=1
x=seq(t_lim[1],t_lim[2],len = t_Clens)
par(cex = 1.5, fin = c(8, 1), mai = c(0,0, 0.5, 0), oma = c(1, 0, 0, 0))
image(x,y,z,col = topo.colors(t_Clens),axes=FALSE,xlab="",ylab="")
title("Posterior ICE6G GIA variance -- mu = 0")
axis(1)})
