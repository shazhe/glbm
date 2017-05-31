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

library(GEOmap)
library(INLA)
if(runserver){
  INLA:::inla.dynload.workaround()
}#Use this on server with old gcc library
library(fields)
library(maps)
library(maptools)
source("glbm/Experiment1a/Rscript/MVSTplus.R")

## Get the GPS locations
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obsU$y_center, lon = GPS_obsU$x_center))
GPSX <- ifelse(GPS_obsU$x_center > 180, GPS_obsU$x_center-360, GPS_obsU$x_center)
GPSY <- GPS_obsU$y_center

#### load the results from different mesh sizes
res_inlas <- list()
projs <- list()
pars_GIAs <- list()
meshSize <- c("s", "m", "l")
errSize <- c("s", "L")
for (i in 1:3){
    exnames <- paste0("/experimentBHM/",meshSize[i], "Mesh_sErr_")
    load(paste0(wkdir, exnames, "inla.RData"))
    res_inlas[[i]] <- res_inla
    pars_GIAs[[i]] <- inla.spde2.result(res_inla, "GIA", GIA_spde, do.transf=TRUE)
    projs[[i]] <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(361,181))
}

for (i in 4:6){
  exnames <- paste0("/experimentBHM/",meshSize[i-3], "Mesh_LErr_")
  load(paste0(wkdir, exnames, "inla.RData"))
  res_inlas[[i]] <- res_inla
  pars_GIAs[[i]] <- inla.spde2.result(res_inla, "GIA", GIA_spde, do.transf=TRUE)
}

for (ee in 1:2){
  L <- ifelse(ee == 1, 0, 3)
## Plot hyperparameters -- small errors
pdf(file = paste0(wkdir, exname, ee, "hyperpar.pdf"), width = 6, height = 8)
par(mfrow = c(3,2))
## plot log(rho)
maxd <- max(c(pars_GIAs[[1+L]]$marginals.log.range.nominal[[1]][,2],
              pars_GIAs[[2+L]]$marginals.log.range.nominal[[1]][,2],
              pars_GIAs[[3+L]]$marginals.log.range.nominal[[1]][,2]))
plot(pars_GIAs[[1+L]]$marginals.log.range.nominal[[1]], type = "l", ylim = c(0,signif(maxd,2)),
     main =expression(bold(log(rho)))) # The posterior from inla output
lines(pars_GIAs[[2+L]]$marginals.log.range.nominal[[1]], type = "l", col = 2) # The posterior from inla output
lines(pars_GIAs[[3+L]]$marginals.log.range.nominal[[1]], type = "l", col = 4)

## plot rho
maxd <- max(c(pars_GIAs[[1+L]]$marginals.range.nominal[[1]][,2],
              pars_GIAs[[2+L]]$marginals.range.nominal[[1]][,2],
              pars_GIAs[[3+L]]$marginals.range.nominal[[1]][,2]))
plot(pars_GIAs[[1+L]]$marginals.range.nominal[[1]], type = "l", ylim = c(0,signif(maxd,2)),
     main = expression(bold(rho))) # The posterior from inla output
lines(pars_GIAs[[2+L]]$marginals.range.nominal[[1]], type = "l", col = 2) 
lines(pars_GIAs[[3+L]]$marginals.range.nominal[[1]], type = "l", col = 4)


## plot log(sigma)
maxd <- max(c(pars_GIAs[[1+L]]$marginals.log.variance.nominal[[1]][,2],
              pars_GIAs[[2+L]]$marginals.log.variance.nominal[[1]][,2],
              pars_GIAs[[3+L]]$marginals.log.variance.nominal[[1]][,2]))

plot(pars_GIAs[[1+L]]$marginals.log.variance.nominal[[1]], type = "l", ylim = c(0,signif(maxd,2)),
     main = expression(bold(log(sigma)))) # The posterior from inla output
lines(pars_GIAs[[2+L]]$marginals.log.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output
lines(pars_GIAs[[3+L]]$marginals.log.variance.nominal[[1]], type = "l", col = 4) # The posterior from inla output

## plot sigma
maxd <- max(c(pars_GIAs[[1+L]]$marginals.variance.nominal[[1]][,2],
              pars_GIAs[[2+L]]$marginals.variance.nominal[[1]][,2],
              pars_GIAs[[3+L]]$marginals.variance.nominal[[1]][,2]))
maxv <- max(c(pars_GIAs[[1+L]]$marginals.variance.nominal[[1]][,1],
              pars_GIAs[[2+L]]$marginals.variance.nominal[[1]][,1],
              pars_GIAs[[3+L]]$marginals.variance.nominal[[1]][,1]))

plot(pars_GIAs[[1+L]]$marginals.variance.nominal[[1]], type = "l",  xlim = c(0, maxv/20), ylim = c(0, signif(maxd, 2)),
     main = expression(bold(sigma))) # The posterior from inla output
lines(pars_GIAs[[2+L]]$marginals.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output
lines(pars_GIAs[[3+L]]$marginals.variance.nominal[[1]], type = "l", col = 4) # The posterior from inla output

## plot the measurement error sigma_e
sigma_e_marginals1 <- inla.tmarginal(function(x) sqrt(1/x), res_inlas[[1+L]]$marginals.hyperpar[[1]])
sigma_e_marginals2 <- inla.tmarginal(function(x) sqrt(1/x), res_inlas[[2+L]]$marginals.hyperpar[[1]])
sigma_e_marginals3 <- inla.tmarginal(function(x) sqrt(1/x), res_inlas[[3+L]]$marginals.hyperpar[[1]])
xxd <- range(c(sigma_e_marginals1[,1], sigma_e_marginals2[,1], sigma_e_marginals3[,1]))
plot(sigma_e_marginals1, xlim = round(xxd, digits = 4),
     main = expression(bold({sigma[e]})), type = "l")
lines(sigma_e_marginals2, type = "l", col = 2) # The posterior from inla output
lines(sigma_e_marginals3, type = "l", col = 4) 

dev.off()


## Plot the predicted GIA field mean and variance
GIA_spost1 <- res_inlas[[1+L]]$summary.random$GIA$sd
GIA_spost2 <- res_inlas[[2+L]]$summary.random$GIA$sd
GIA_spost3 <- res_inlas[[3+L]]$summary.random$GIA$sd

pdf(file = paste0(wkdir, exname, ee, "GIAfield.pdf"), width = 8, height = 10)
par(mfrow = c(3,1))
## The standard error
image.plot(projs[[1]]$x, projs[[1]]$y, inla.mesh.project(projs[[1]], as.vector(GIA_spost1)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterior marginal standard error -- Small mesh grid")
points(GPSX, GPSY)

image.plot(projs[[2]]$x, projs[[2]]$y, inla.mesh.project(projs[[2]], as.vector(GIA_spost2)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterior marginal standard error -- Medium mesh grid")
points(GPSX, GPSY)

image.plot(projs[[3]]$x, projs[[3]]$y, inla.mesh.project(projs[[3]], as.vector(GIA_spost3)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterior marginal standard error -- Large mesh grid")
points(GPSX, GPSY)
dev.off()
}
