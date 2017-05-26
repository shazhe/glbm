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
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obsU$y_center, lon = GPS_obsU$x_center))
GPSX <- ifelse(GPS_obsU$x_center > 180, GPS_obsU$x_center-360, GPS_obsU$x_center)
GPSY <- GPS_obsU$y_center

## Find the mapping between observations and processes basis
Ay <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)

## create new synthetic data, plot and compare
y0 <- as.vector(Ay %*% Mesh_GIA_sp@data$GIA_m)
ydata1 <- y0 + err1
ydata2 <- y0 + err2
GPS_sd <- SpatialPointsDataFrame(coords = cbind(GPSX, GPSY), data = data.frame(GPS=y0))



#### 1: GIA_ice6g prior mean info
st.est1 <- inla.stack(data = list(y=ydata1), A = list(Ay,Ay), 
                     effects = list(GIA = 1:GIA_spde$n.spde, offset = Mesh_GIA_sp@data$GIA_m), tag = "est")

st.est2 <- inla.stack(data = list(y=ydata2), A = list(Ay,Ay), 
                      effects = list(GIA = 1:GIA_spde$n.spde, offset = Mesh_GIA_sp@data$GIA_m), tag = "est")

Ip <- Matrix(0, GIA_spde$n.spde, GIA_spde$n.spde)
diag(Ip) <- 1

st.pred <- inla.stack(data = list(y=NA), A = list(Ip, Ip), 
                      effects = list(mi=1:GIA_spde$n.spde, offset = Mesh_GIA_sp@data$GIA_m), tag = "pred")
stGIA1 <- inla.stack(st.est1, st.pred)
stGIA2 <- inla.stack(st.est2, st.pred)

formula = y ~ -1 + offset + f(GIA, model = GIA_spde)

res_inla1 <- inla(formula, data = inla.stack.data(stGIA1, spde = GIA_spde),
                 control.predictor=list(A=inla.stack.A(stGIA1), compute =TRUE))

res_inla2 <- inla(formula, data = inla.stack.data(stGIA2, spde = GIA_spde),
                  control.predictor=list(A=inla.stack.A(stGIA2), compute =TRUE))


## Plot hyperparameters
pars_GIA1 <- inla.spde2.result(res_inla1, "GIA", GIA_spde, do.transf=TRUE)
pars_GIA2 <- inla.spde2.result(res_inla2, "GIA", GIA_spde, do.transf=TRUE)

pdf(file = paste0(wkdir, exname, "hyperpar.pdf"), width = 6, height = 8)
par(mfrow = c(3,2))
## plot log(rho)
plot(pars_GIA1$marginals.log.range.nominal[[1]], type = "l", 
     main =expression(bold(log(rho)))) # The posterior from inla output
lines(pars_GIA2$marginals.log.range.nominal[[1]], type = "l", col = 2) # The posterior from inla output

## plot rho
plot(pars_GIA1$marginals.range.nominal[[1]], type = "l", 
     main = expression(bold(log(rho)))) # The posterior from inla output
lines(pars_GIA2$marginals.range.nominal[[1]], type = "l", col = 2) 

## plot log(sigma)
plot(pars_GIA1$marginals.log.variance.nominal[[1]], type = "l",
     main = expression(bold(log(sigma)))) # The posterior from inla output
lines(pars_GIA2$marginals.log.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output

## plot sigma
plot(pars_GIA1$marginals.variance.nominal[[1]], type = "l", xlim = c(0, 20), 
     main = expression(bold(sigma))) # The posterior from inla output
lines(pars_GIA2$marginals.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output

## plot the measurement error sigma_e
sigma_e_marginals1 <- inla.tmarginal(function(x) 1/x, res_inla1$marginals.hyperpar[[1]])
sigma_e_marginals2 <- inla.tmarginal(function(x) 1/x, res_inla2$marginals.hyperpar[[1]])
plot(sigma_e_marginals1, 
     main = expression(bold({sigma[e]^2})), type = "l")
lines(sigma_e_marginals2, type = "l", col = 2) # The posterior from inla output

dev.off()


## Plot the predicted GIA field mean and variance
GIA_spost1 <- res_inla1$summary.random$GIA$sd
GIA_spost2 <- res_inla2$summary.random$GIA$sd
proj1 <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(361,181))

GPScol1 <- ifelse(err1 > 0, 2, 1)
GPScol2 <- ifelse(err2 > 0, 2, 1)

pdf(file = paste0(wkdir, exname, "GIAfield.pdf"), width = 8, height = 10)
par(mfrow = c(2,1))
## The standard error
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA_spost1)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterior marginal standard error -- Small error")
points(GPSX, GPSY, cex = sqrt(abs(err1)), col = GPScol1, pch = 1)
rect(xleft =-15 , ybottom =-85 , xright = 10, ytop = -60, lwd = 1.5, border = 2)


image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA_spost2)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterior marginal standard error -- Large error")
points(GPSX, GPSY, cex = sqrt(abs(err2)), col = GPScol2, pch = 1)
rect(xleft =-15 , ybottom =-85 , xright = 10, ytop = -60, lwd = 1.5, border = 2)

dev.off()



save(res_inla1, res_inla2, GIA_spde, ydata1, ydata2, mu_r, v_r, mu_s, v_s, Mesh_GIA, file = paste0(wkdir, exname, "inla.RData"))


