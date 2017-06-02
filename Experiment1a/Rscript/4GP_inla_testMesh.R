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
res_inlas <- list()
projs <- list()
pars_GIAs <- list()
meshSize <- c("s", "m", "l")

for (i in 1:3){
  exname <- paste0("/experimentBHM/",meshSize[i],"Mesh")
  load(paste0("experimentBHM/Mesh_GIA", meshSize[i], ".RData"))
  source("glbm/Experiment1a/Rscript/4GP_Bayes_update_inla.R")
  res_inlas[[i]] <- res_inla
  pars_GIAs[[i]] <- pars_GIA
  projs[[i]] <- proj
}


## Plot hyperparameters -- small errors
pdf(file = paste0(wkdir, exname, "_hyperpar.pdf"), width = 6, height = 8)
par(mfrow = c(3,2))
## plot log(rho)
maxd <- max(c(pars_GIAs[[1]]$marginals.log.range.nominal[[1]][,2],
              pars_GIAs[[2]]$marginals.log.range.nominal[[1]][,2],
              pars_GIAs[[3]]$marginals.log.range.nominal[[1]][,2]))
plot(pars_GIAs[[1]]$marginals.log.range.nominal[[1]], type = "l", ylim = c(0,signif(maxd,2)),
     main =expression(bold(log(rho)))) # The posterior from inla output
lines(pars_GIAs[[2]]$marginals.log.range.nominal[[1]], type = "l", col = 2) # The posterior from inla output
lines(pars_GIAs[[3]]$marginals.log.range.nominal[[1]], type = "l", col = 4)

## plot rho
maxd <- max(c(pars_GIAs[[1]]$marginals.range.nominal[[1]][,2],
              pars_GIAs[[2]]$marginals.range.nominal[[1]][,2],
              pars_GIAs[[3]]$marginals.range.nominal[[1]][,2]))
plot(pars_GIAs[[1]]$marginals.range.nominal[[1]], type = "l", ylim = c(0,signif(maxd,2)),
     main = expression(bold(rho))) # The posterior from inla output
lines(pars_GIAs[[2]]$marginals.range.nominal[[1]], type = "l", col = 2) 
lines(pars_GIAs[[3]]$marginals.range.nominal[[1]], type = "l", col = 4)


## plot log(sigma)
maxd <- max(c(pars_GIAs[[1]]$marginals.log.variance.nominal[[1]][,2],
              pars_GIAs[[2]]$marginals.log.variance.nominal[[1]][,2],
              pars_GIAs[[3]]$marginals.log.variance.nominal[[1]][,2]))

plot(pars_GIAs[[1]]$marginals.log.variance.nominal[[1]], type = "l", ylim = c(0,signif(maxd,2)),
     main = expression(bold(log(sigma)))) # The posterior from inla output
lines(pars_GIAs[[2]]$marginals.log.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output
lines(pars_GIAs[[3]]$marginals.log.variance.nominal[[1]], type = "l", col = 4) # The posterior from inla output

## plot sigma
maxd <- max(c(pars_GIAs[[1]]$marginals.variance.nominal[[1]][,2],
              pars_GIAs[[2]]$marginals.variance.nominal[[1]][,2],
              pars_GIAs[[3]]$marginals.variance.nominal[[1]][,2]))
maxv <- max(c(pars_GIAs[[1]]$marginals.variance.nominal[[1]][,1],
              pars_GIAs[[2]]$marginals.variance.nominal[[1]][,1],
              pars_GIAs[[3]]$marginals.variance.nominal[[1]][,1]))

plot(pars_GIAs[[1]]$marginals.variance.nominal[[1]], type = "l",  xlim = c(0, maxv/20), ylim = c(0, signif(maxd, 2)),
     main = expression(bold(sigma))) # The posterior from inla output
lines(pars_GIAs[[2]]$marginals.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output
lines(pars_GIAs[[3]]$marginals.variance.nominal[[1]], type = "l", col = 4) # The posterior from inla output

dev.off()


## Plot the predicted GIA field mean and variance
GIA_spost1 <- res_inlas[[1]]$summary.random$GIA$sd
GIA_spost2 <- res_inlas[[2]]$summary.random$GIA$sd
GIA_spost3 <- res_inlas[[3]]$summary.random$GIA$sd

pdf(file = paste0(wkdir, exname, "_GIAfield.pdf"), width = 8, height = 10)
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
