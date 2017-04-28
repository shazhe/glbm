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
library(rgl)
library(GEOmap)
library(INLA)
INLA:::inla.dynload.workaround() #Use this on server with old gcc library
library(MVST)
library(fields)
load("experimentBHM/Mesh_GIA.RData")
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
#load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/Mesh_GIA.RData")
#GPS_obs <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obsU$y_center, lon = GPS_obsU$x_center))

## Find the mapping between observations and processes basis
Ay <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
y <- as.vector(GPS_obsU$trend)


#### 1: GIA_ice6g prior mean info
st.est1 <- inla.stack(data = list(y=y), A = list(Ay,Ay), 
                     effects = list(GIA = 1:GIA_spde$n.spde, offset = Mesh_GIA_sp@data$GIA_m), tag = "est")
Ip <- Matrix(0, GIA_spde$n.spde, GIA_spde$n.spde)
diag(Ip) <- 1
st.pred1 <- inla.stack(data = list(y=NA), A = list(Ip, Ip), 
                      effects = list(mi=1:GIA_spde$n.spde, offset = Mesh_GIA_sp@data$GIA_m), tag = "pred")
stGIA1 <- inla.stack(st.est1, st.pred1)

formula1 = y ~ -1 + offset + f(GIA, model = GIA_spde)
res_inla1 <- inla(formula1, data = inla.stack.data(stGIA1, spde = GIA_spde),
                   control.predictor=list(A=inla.stack.A(stGIA1), compute =TRUE))

save(res_inla1, file = "res_inla1a.RData")
#load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/res_inla2.RData")

## Plot hyperparameters
res_GIA1 <- inla.spde2.result(res_inla1, "GIA", GIA_spde, do.transf=TRUE)

pdf(file = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/inla1hyper.pdf", width = 8.5, height = 11)
par(mfrow = c(3,2))
plot(res_GIA1$marginals.variance.nominal[[1]], type = "l", main = "variance")
plot(res_GIA1$marginals.range.nominal[[1]], type = "l", main = "range")
plot(res_GIA1$marginals.kappa[[1]], type = "l", main = "kappa")
plot(res_GIA1$marginals.tau[[1]], type = "l", main = "tau")
plot(res_inla1$marginals.hyperpar$`Precision for the Gaussian observations`, type = "l", main = "precision of error")
dev.off()

##plot results
GIA1_mpost <- res_inla1$summary.random$GIA$mean
GIA1_spost <- res_inla1$summary.random$GIA$sd

pidx <- inla.stack.index(stGIA1, tag = "pred")
GIA1_mpred <- res_inla1$summary.fitted.values$mean[pidx$data]
GIA1_spred <- res_inla1$summary.fitted.values$sd[pidx$data]
  
diff1 <- GIA1_mpred - Mesh_GIA_sp@data$GIA_m

proj1 <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(361,181))
GPSX <- ifelse(GPS_obs$x_center > 180, GPS_obs$x_center-360, GPS_obs$x_center)
GPSY <- GPS_obs$y_center
#proj2 <- inla.mesh.projector(Mesh_GIA, projection = "mollweide", dims = c(361,181))
pdf(file = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/inla1Mean.pdf", width = 8.5, height = 11)
par(mfrow = c(3,1))
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_mpost)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- mean")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_mpred)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- mean")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(Mesh_GIA_sp@data$GIA_m)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "ICE6G GIA")
points(GPSX, GPSY, pch = 20)
dev.off()

pdf(file = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/inla1Var.pdf", width = 8.5, height =11)
par(mfrow = c(3,1))
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_spost)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginal standard error")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_spred)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginal standard error")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(diff1)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Difference between posterier marginal mean and ICE6G")
points(GPSX, GPSY, pch = 20)
dev.off()


par(mfrow = c(3,2))
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_mpost)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- mean")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_spost)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginal standard error")
points(GPSX, GPSY, pch = 20)



image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_mpred)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier prediction -- mean")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIA1_spred)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier prediction standard error")
points(GPSX, GPSY, pch = 20)


image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(Mesh_GIA_sp@data$GIA_m)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "ICE6G GIA")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(diff1)), col = topo.colors(20),
           xlab = "Longitude", ylab = "Latitude", main = "Difference between posterier prediction and ICE6G")
points(GPSX, GPSY, pch = 20)


#### 2 Test when GPS error ---> 0
## 2.1 produce sythetic GSPS data
GPS_sp <- SpatialPoints(data.frame(lon = GPS_obsU$x_center, lat = GPS_obsU$y_center), proj4string = CRS("+proj=longlat")) #This convert GIA_ice6g a SpatialPointDataFrame
GPS_data_new <- over(GPS_sp, GIA_ice6g_sp)$trend + rnorm(455)*0.01

## new redo previous by setting y = GPS_data_new


































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




