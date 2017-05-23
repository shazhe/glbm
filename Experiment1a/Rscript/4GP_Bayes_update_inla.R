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
INLA:::inla.dynload.workaround() #Use this on server with old gcc library
library(fields)
library(maps)
library(maptools)
source("experimentBHM/Rscript/MVSTplus.R")

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
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obsU$y_center, lon = GPS_obsU$x_center))
GPSX <- ifelse(GPS_obsU$x_center > 180, GPS_obsU$x_center-360, GPS_obsU$x_center)
GPSY <- GPS_obsU$y_center

## Find the mapping between observations and processes basis
Ay <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)

## create new synthetic data, plot and compare
y0 <- as.vector(Ay %*% Mesh_GIA_sp@data$GIA_m)
ydata <- y0 + err1
GPS_sd <- SpatialPointsDataFrame(coords = cbind(GPSX, GPSY), data = data.frame(GPS=y0))

## plot the "prior mean" and synthetic data
pdf(file = paste0(wkdir, exname, "data.pdf"), width = 6, height = 8)
par(mfrow = c(2,1))
map(interior = FALSE)
plot(GPS_sd, pch = 19, cex = sqrt(abs(y0)), col = sign(y0)+6, add = TRUE)

map(interior = FALSE)
plot(GPS_sd, pch = 19, cex = sqrt(abs(ydata)), col = sign(err1)+6, add = TRUE)
dev.off()


#### 1: GIA_ice6g prior mean info
st.est <- inla.stack(data = list(y=ydata), A = list(Ay,Ay), 
                     effects = list(GIA = 1:GIA_spde$n.spde, offset = Mesh_GIA_sp@data$GIA_m), tag = "est")
Ip <- Matrix(0, GIA_spde$n.spde, GIA_spde$n.spde)
diag(Ip) <- 1
st.pred <- inla.stack(data = list(y=NA), A = list(Ip, Ip), 
                      effects = list(mi=1:GIA_spde$n.spde, offset = Mesh_GIA_sp@data$GIA_m), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

formula = y ~ -1 + offset + f(GIA, model = GIA_spde)
res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde),
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


#pdf(file = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/inla1hyper.pdf", width = 8.5, height = 11)
par(mfrow = c(2,2))
## plot log(rho)
plot(pars_GIA$marginals.log.range.nominal[[1]], type = "l", 
     main = bquote(bold(log(rho)[mode] == .(round(lrho_mode, 4))))) # The posterior from inla output
xx <- pars_GIA$marginals.log.range.nominal[[1]][,1]
yy0 <- dnorm(xx, mean = lrho0, sd = theta2_s) # The prior
yy <- dnorm(xx, mean = theta_mean[2] +lrho0, sd = theta_sd[2]) # Posterior derived from theta2
lines(x = xx, y = yy0, col = 2)
points(x = xx, y = yy, col = 3)

## plot rho
plot(pars_GIA$marginals.range.nominal[[1]], type = "l", 
     main = bquote(bold(rho[mode] == .(round(rho_mode, 4))))) # The posterior from inla output
xx <- pars_GIA$marginals.range.nominal[[1]][,1]
yy0 <- dlnorm(xx, meanlog = lrho0, sdlog = sqrt(theta2_s)) # The prior
yy <- dlnorm(xx, mean = lrho0 + theta_mean[2], sd = theta_sd[2]) # Posterior derived from theta2
lines(x = xx, y = yy0, col = 2)
points(x = xx, y = yy, col = 3)


## Plot the field variance parameter sigma
lsigma_mode <- pars_GIA$summary.log.variance.nominal$mode
lsigma_mean <- pars_GIA$summary.log.variance.nominal$mean
lsigma_sd <- pars_GIA$summary.log.variance.nominal$sd
sigma_post <- lognormT(logmu = lrho_mean, logv = (lrho_sd)^2)
sigma_mode <- exp(lsigma_mean - lsigma_sd^2)

## plot log(sigma)
plot(pars_GIA$marginals.log.variance.nominal[[1]], type = "l",
     main = bquote(bold(log(sigma)[mode] == .(round(lsigma_mode, 4))))) # The posterior from inla output
xx <- pars_GIA$marginals.log.variance.nominal[[1]][,1]
yy0 <- dnorm(xx, mean = lsigma0, sd = theta1_s) # The prior
yy <- dnorm(xx, mean = theta_mean[2] - theta_mean[1] +lsigma0, sd = sqrt(theta_sd[1]^2 + theta_sd[2]^2)) # Posterior derived from theta2
lines(x = xx, y = yy0, col = 2)
points(x = xx, y = yy, col = 3)

## plot rho
plot(pars_GIA$marginals.variance.nominal[[1]], type = "l", xlim = c(0, 10), 
     main = bquote(bold(sigma[mode] == .(round(sigma_mode, 4))))) # The posterior from inla output
xx <- seq(from = 0, to=10, length.out = 1000)
yy0 <- dlnorm(xx, meanlog = lsigma0, sdlog = sqrt(theta1_s)) # The prior
yy <- dlnorm(xx, mean = lsigma0 + theta_mean[2], sd = theta_sd[2]) # Posterior derived from theta2
lines(x = xx, y = yy0, col = 2)


save(res_inla1, GIA_spde, file = "res_inlaC.RData")




##plot results
GIA1_mpost <- res_inla$summary.random$GIA$mean
GIA1_spost <- res_inla$summary.random$GIA$sd

pidx <- inla.stack.index(stGIA, tag = "pred")
GIA1_mpred <- res_inla$summary.fitted.values$mean[pidx$data]
GIA1_spred <- res_inla$summary.fitted.values$sd[pidx$data]
  
diff1 <- GIA1_mpred - Mesh_GIA_sp@data$GIA_m




proj1 <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(361,181))

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_mpost)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- mean")
points(GPSX, GPSY, pch = 20)
points(GPSX[idx0], GPSY[idx0], cex = 2, pch = 1, col = 2)


image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_mpred)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- mean")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(Mesh_GIA_sp@data$GIA_m)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "ICE6G GIA")
points(GPSX, GPSY, pch = 20)
dev.off()

#pdf(file = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/inla1Var.pdf", width = 8.5, height =11)
#par(mfrow = c(2,1))
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_spost)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginal standard error")
points(GPSX, GPSY, pch = 20, cex = abs(err1)/10, col = sign(err1)+3)
#points(GPSX[idx0], GPSY[idx0], cex = 2, pch = 1, col = 2)


image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_spred)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginal standard error")
points(GPSX, GPSY, pch = 20,cex = abs(err1), col = sign(err1)+3)


#dev.off()



































#### 2: no GIA_ice6g prior mean info
st.est1 <- inla.stack(data = list(y=y), A = list(Ay), 
                      effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
st.pred1 <- inla.stack(data = list(y=NA), A = Ip, 
                       effects = list(mi=1:GIA_spde$n.spde), tag = "pred")
stGIA1 <- inla.stack(st.est1, st.pred1)
formula1 = y ~  -1 + f(GIA, model = GIA_spde)
res_inla1 <- inla(formula1, data = inla.stack.data(stGIA1, spde = GIA_spde),
                  control.predictor=list(A=inla.stack.A(stGIA1), compute = TRUE))
save(res_inla1, file = "res_inla1.RData")
#load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/res_inla1.RData")




## Plot hyperparameters
res_GIA2 <- inla.spde2.result(res_inla2, "GIA", GIA_spde, do.transf=TRUE)
pdf(file = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/inla2hyper.pdf", width = 8.5, height = 11)
par(mfrow = c(3,2))
plot(res_GIA2$marginals.variance.nominal[[1]], type = "l", main = "variance")
plot(res_GIA2$marginals.range.nominal[[1]], type = "l", main = "range")
plot(res_GIA2$marginals.kappa[[1]], type = "l", main = "kappa")
plot(res_GIA2$marginals.tau[[1]], type = "l", main = "tau")
plot(res_inla2$marginals.hyperpar$`Precision for the Gaussian observations`, type = "l", main = "precision of error")
dev.off()

##plot results
index2 <- inla.stack.index(stGIA2, tag = "pred")
GIA2_mpost <- res_inla2$summary.random$GIA$mean
GIA2_spost <- res_inla2$summary.random$GIA$sd

GIA2_mpred <- res_inla2$summary.fitted.values$mean[index2$data]
GIA2_spred <- res_inla2$summary.fitted.values$sd[index2$data]
diff2 <- GIA2_mpred - GIA_mu

pdf(file = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/inla2Mean.pdf", width = 8.5, height = 11)
par(mfrow = c(3,1))
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA2_mpost)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginal mean")
points(GPSX, GPSY, pch = 20)
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA2_mpred)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Fitted response")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA_mu)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "ICE6G GIA")
points(GPSX, GPSY, pch = 20)
dev.off()

pdf(file = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/inla2Var.pdf", width = 8.5, height =11)
par(mfrow = c(3,1))
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA2_spost)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginal standard error")
points(GPSX, GPSY, pch = 20)
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA2_spred)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Fitted response standard error")
points(GPSX, GPSY, pch = 20)
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(diff2)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Difference between posterier marginal mean and ICE6G")
points(GPSX, GPSY, pch = 20)
dev.off()


## Plot the results 3D
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

#writeWebGL(filename= "GIA_globe.html",  # save the plot as an html
#           width = 1000, reuse = TRUE)




