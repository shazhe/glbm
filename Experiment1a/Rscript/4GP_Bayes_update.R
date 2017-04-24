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
load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/Mesh_GIA.RData")
GPS_obs <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obsU$y_center, lon = GPS_obsU$x_center))

## Find the mapping between observations and processes basis
CMat <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
GIA_spde <- inla.spde2.matern(Mesh_GIA)
y <- as.vector(GPS_obsU$trend)


#### Test1: no GIA_ice6g prior mean info
st.est <- inla.stack(data = list(y=y), A = list(CMat,1), 
                     effects = list(mi = 1:GIA_spde$n.spde, intercept = rep(1, length(y))), tag = "est")
formula = y ~ 0 + intercept + f(mi, model = GIA_spde)
result.est <- inla(formula, data = inla.stack.data(st.est, spde = GIA_spde),
                   control.predictor=list(A=inla.stack.A(st.est)))

#### Test2: use GIA_ice6g as prior mean for offset 
st.est <- inla.stack(data = list(y=y), A = list(CMat,CMat), 
                     effects = list(mi = 1:GIA_spde$n.spde, offset = GIA_mu), tag = "est")
formula = y ~ 0 + offset + f(mi, model = GIA_spde)
result.est <- inla(formula, data = inla.stack.data(st.est, spde = GIA_spde),
                   control.predictor=list(A=inla.stack.A(st.est)))


#### Test3 defaut fixed
st.est <- inla.stack(data = list(y=y), A = list(CMat), 
                     effects = list(mi = 1:GIA_spde$n.spde), tag = "est")
formula = y ~  f(mi, model = GIA_spde)
result.est <- inla(formula, data = inla.stack.data(st.est, spde = GIA_spde),
                   control.predictor=list(A=inla.stack.A(st.est)))

#### Test4 no fixed no intercept
st.est2 <- inla.stack(data = list(y=y), A = list(CMat), 
                     effects = list(mi = 1:GIA_spde$n.spde), tag = "est")
formula2 = y ~  -1 + f(mi, model = GIA_spde)
result.est2 <- inla(formula2, data = inla.stack.data(st.est2, spde = GIA_spde),
                   control.predictor=list(A=inla.stack.A(st.est2)))







st.est <- inla.stack(data = list(y=y), A = list(CMat, CMat), 
                     effects = list(c(mi = 1:GIA_spde$n.spde, pmu = GIA_mu), tag = "est"))
formula = y ~ -1 + pmu + f(mi, model = GIA_spde)
result.est <- inla(formula, data = inla.stack.data(st.est, spde = GIA_spde))


st.pred <- inla.stack(data = list(y=NA), A = list(diag(1, 32546),diag(1,32546)), 
                     effects = list(mi=1:GIA_spde$n.spde, pmu = GIA_mu), tag = "pred")

stGIA <- inla.stack(st.est, st.pred)
formula = y ~ 0 + pmu + f(mi, model = GIA_spde)
result.pred <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde),
               control.predictor = list(A = inla.stack.A(stGIA), compute = TRUE))


par(mfrow = c(3,2))
res1 <- inla.spde2.result(result.est, "mi", GIA_spde)
plot(res1[["marginals.range.nominal"]][[1]], type = "l", main = "range")
plot(res1[["marginals.variance.nominal"]][[1]], type = "l", main = "variance")
plot(res1[["marginals.kappa"]][[1]], type = "l", main = "kappa")
plot(res1[["marginals.tau"]][[1]], type = "l", main = "tau")


plot(res1[["marginals.log.kappa"]][[1]], type = "l", main = "log.kappa")
plot(res1[["marginals.log.tau"]][[1]], type = "l", main = "log.tau")

GIA_mpost <- result.pred$summary.random$mi$mean
GIA_spost <- result.est$summary.random$mi$sd
GIA_mmpost <- result.est$marginals.random$mi$index.1

GIA_mpost2 <- result.pred$summary.linear.predictor$mean

## Plot the results
vals <- GIA_mpost
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
title("INLA Posterior ICE6G GIA, mu = 0")
axis(1)})
