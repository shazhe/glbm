###################################################################
##          A Bayesian Gaussian Update (INLA)                    ##
## This script performs a full Bayesian Gaussian Update          ##
## of the GIA processes                                          ##
## Use a lot of memories -- run on nuuk at least                 ##
## Given                                                         ##
## 1. prior assumptions of the parameters                        ##
## 2. subset of GPS obs                                          ##
## 3. GIA forward model mean adjustment                          ##
###################################################################
#### INLA Approximation
## require
library(GEOmap)
library(INLA)
if(runserver){
INLA:::inla.dynload.workaround()
  }#Use this on server with old gcc library
library(fields)
library(maps)
library(maptools)
source("glbm/Experiment1a/Rscript/MVSTplus.R")

## Setup for spde parameters
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]

lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

## Get the GPS locations
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obs$y_center, lon = GPS_obs$x_center))
GPSX <- ifelse(GPS_obs$x_center > 180, GPS_obs$x_center-360, GPS_obs$x_center)
GPSY <- GPS_obs$y_center

## Find the mapping between observations and processes basis
Ay <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)

## create new synthetic data, plot and compare
y0 <- as.vector(Ay %*% Mesh_GIA_sp@data$GIA_m)
ydata <- GPS_obs$trend - y0

#### 1: GIA_ice6g prior mean info
st.est <- inla.stack(data = list(y=ydata), A = list(Ay), 
                     effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
## Predict at the GPS location and mesh grid
Ip <- Matrix(0, GIA_spde$n.spde, GIA_spde$n.spde)
diag(Ip) <- 1
st.pred <- inla.stack(data = list(y=NA), A = list(rbind(Ay, Ip)), 
                      effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

hyper <- list(prec = list(fixed = TRUE, initial = 0))
formula = y ~ -1 + f(GIA, model = GIA_spde)

res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian", 
                 scale = c(rep((1/40),2*length(ydata)), rep(1, GIA_spde$n.spde)),
                 control.family = list(hyper = hyper),
                   control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))
summary(res_inla)

## Plot hyperparameters
pars_GIA <- inla.spde2.result(res_inla, "GIA", GIA_spde, do.transf=TRUE)
theta_mean <- pars_GIA$summary.theta$mean
theta_sd <- pars_GIA$summary.theta$sd

## plot the range parameter rho
lrho_mode <- pars_GIA$summary.log.range.nominal$mode
lrho_mean <- pars_GIA$summary.log.range.nominal$mean
lrho_sd <- pars_GIA$summary.log.range.nominal$sd
rho_post <- lognormT(logmu = lrho_mean, logv = (lrho_sd)^2)
rho_mean <- exp(pars_GIA$summary.log.range.nominal$mean)
rho_mode <- exp(lrho_mean - lrho_sd^2)


pdf(file = paste0(wkdir, exname, "hyperpar.pdf"), width = 6, height = 8)
par(mfrow = c(3,2))
## plot log(rho)
plot(pars_GIA$marginals.log.range.nominal[[1]], type = "l", 
     main = bquote(bold(log(rho)("mode") == .(round(lrho_mode, 4))))) # The posterior from inla output
xx <- pars_GIA$marginals.log.range.nominal[[1]][,1]
yy0 <- dnorm(xx, mean = lrho0, sd = theta2_s) # The prior
lines(x = xx, y = yy0, col = 2)

## plot rho
plot(pars_GIA$marginals.range.nominal[[1]], type = "l", 
     main = bquote(bold(rho("mode") == .(round(rho_mode, 4))))) # The posterior from inla output
xx <- pars_GIA$marginals.range.nominal[[1]][,1]
yy0 <- dlnorm(xx, meanlog = lrho0, sdlog = sqrt(theta2_s)) # The prior
lines(x = xx, y = yy0, col = 2)


## Plot the field variance parameter sigma
lsigma_mode <- pars_GIA$summary.log.variance.nominal$mode
lsigma_mean <- pars_GIA$summary.log.variance.nominal$mean
lsigma_sd <- pars_GIA$summary.log.variance.nominal$sd
sigma_post <- lognormT(logmu = lrho_mean, logv = (lrho_sd)^2)
sigma_mode <- exp(lsigma_mean - lsigma_sd^2)

## plot log(sigma)
plot(pars_GIA$marginals.log.variance.nominal[[1]], type = "l",
     main = bquote(bold(log(sigma)("mode") == .(round(lsigma_mode, 4))))) # The posterior from inla output
xx <- pars_GIA$marginals.log.variance.nominal[[1]][,1]
yy0 <- dnorm(xx, mean = lsigma0, sd = theta1_s) # The prior
lines(x = xx, y = yy0, col = 2)

## plot sigma
plot(pars_GIA$marginals.variance.nominal[[1]], type = "l", xlim = c(0, 20), 
     main = bquote(bold(sigma("mode") == .(round(sigma_mode, 4))))) # The posterior from inla output
xx <- seq(from = 0.2, to=20, length.out = 1000)
yy0 <- dlnorm(xx, meanlog = lsigma0, sdlog = sqrt(theta1_s)) # The prior
lines(x = xx, y = yy0, col = 2)

dev.off()




## Plot the predicted GIA field mean and variance
proj <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(361,181))

pidx <- inla.stack.index(stGIA, tag = "pred")
GIA_mpost <- res_inla$summary.linear.predictor$mean[pidx$data[-c(1:length(ydata))]] + Mesh_GIA_sp@data$GIA_m
GIA_spost <- res_inla$summary.linear.predictor$sd[pidx$data[-c(1:length(ydata))]]

GPS_mpost <- res_inla$summary.linear.predictor$mean[pidx$data[1:length(ydata)]]
GPS_spost <- res_inla$summary.linear.predictor$sd[pidx$data[1:length(ydata)]]

GPScol <- ifelse(ydata > 0, 2, 1)


pdf(file = paste0(wkdir, exname, "GIAfield.pdf"), width = 8, height = 10)
par(mfrow = c(2,1))
## The mean field
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_mpost)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterior marginals -- mean")
points(GPSX, GPSY, col = GPScol, pch = 20, cex = 0.8)

## The standard error
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_spost)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterior marginal standard error")
points(GPSX, GPSY, cex = (GPS_spost/mean(GPS_spost))^2, col = GPScol, pch = 1)

dev.off()



save(res_inla, GIA_spde, ydata, mu_r, v_r, mu_s, v_s, Mesh_GIA, file = paste0(wkdir, exname, "inla.RData"))










# ## Plot the results 3D
# vals <- GIA_spost
# open3d()
# par3d(windowRect = c(100, 100, 900, 900))
# layout3d(matrix(1:2, 2,1), heights = c(4,2), sharedMouse = TRUE)
# par3d(zoom = 0.8)
# t_lim <- range(vals)
# t_Clens <- round((t_lim[2] - t_lim[1])*100) + 1
# t_Cpal <- topo.colors(t_Clens, alpha=0)
# t_Cols<- t_Cpal[round((vals - t_lim[1])*100) + 1]
# plot3d(GPS_loc, pch = "20", size = 5, col = GPScol, xlab = "", ylab = "", zlab = "", axe=FALSE)
# plot(Mesh_GIA, rgl = TRUE, col= t_Cols, edge.color = rgb(0, 0.5, 0.6, alpha =0), add = TRUE)
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




