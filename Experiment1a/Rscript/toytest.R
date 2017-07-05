#### Test INLA and MCMC on a toy example 
library(INLA)
library(lattice)


#### 1 Test on the plane (example from JSS paper)
#### 1.1 Generate the True process and the data
## Simulate points -- uniform on [0,1] * [0,1]
m = 400
yloc <- matrix(runif(m*2), m, 2)
mesh <- inla.mesh.2d(loc = yloc, cutoff = 0.05, offset = c(0.1, 0.4), max.edge = c(0.05, 0.5))

## simulate x ~ GP(0, Sigma(0.4, 4))
range0 <- 0.4
sigma0 <- 1
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, 1),
                          B.kappa = cbind(log(kappa0), 0, -1), 
                          theta.prior.mean = c(0,0), theta.prior.prec = c(0.1, 1))
Q <- inla.spde.precision(spde, theta = c(log(2),0))
xSim <- as.numeric(inla.qsample(n=1, Q))

## Plot the simulated process
proj <- inla.mesh.projector(mesh, dims = c(100,100))
simData <- data.frame(x = proj$x, y = proj$y, z = as.numeric(inla.mesh.project(proj, xSim)))
levelplot(row.values = proj$x, column.values = proj$y, x = inla.mesh.project(proj, xSim), 
          contour = TRUE, aspect = "fill", labels = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))

## Use only a fraction of the data and discard some points in the middle
ff <- 0.7
indy <- sample(1:m, m*ff)
y <- y0[indy]
yyloc <- yloc[indy,]



## Find the projectiong matrix from mesh vertices to the locations
A <- inla.spde.make.A(mesh, loc = yyloc)

## Generate the data y = b + Ax +  e
errors <- abs(rnorm(m*ff)*0.01)
y <- as.vector(A %*% xSim) + rnorm(m*ff)

## Add some signals to the data
## 15% points nearer southwest are lower
n1 = round(m*ff*0.1)
id1 <- sample(which(yyloc[,1] <= 0.5 & yyloc[,2] <= 0.5), n1)
y[id1] <- y[id1] - abs(rnorm(n1, sd = 1))
levelplot(row.values = proj$x, column.values = proj$y, x = inla.mesh.project(proj, xSim), 
          panel=function(...){
            panel.levelplot(...)
            panel.points(yyloc[id1,], col = "red", pch = 19)
          },
          contour = TRUE, aspect = "fill", labels = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))

## 5% points nearer north are higher
n2 <- round(m*ff*0.10)
id2 <- sample(which(yyloc[,2] > 0.7), n2)
y[id2] <- y[id2] + abs(rnorm(n2, sd = 2))
levelplot(row.values = proj$x, column.values = proj$y, x = inla.mesh.project(proj, xSim), 
          panel=function(...){
            panel.levelplot(...)
            panel.points(yyloc[id1,], col = "red", pch = 19)
            panel.points(yyloc[id2,], col = "blue", pch = 19)
          },
          contour = TRUE, aspect = "fill", labels = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))


#### 1.2 Now start INLA analysis
ydata <- y - as.vector(A %*% xSim)

st.est <- inla.stack(data = list(y=ydata), A = list(A),
                     effects = list(GIA = 1:spde$n.spde), tag = "est")

## Predict at the y location, mesh grid and predict error location
xg <- seq(0, 1, length.out = 10)
yg <- seq(0, 1, length.out = 10)
xyg <- as.matrix(expand.grid(x = xg, y = yg))
A_pred <- inla.spde.make.A(mesh = mesh, loc = xyg)

st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA=1:spde$n.spde), tag = "pred")

stGIA <- inla.stack(st.est, st.pred)


hyper <- list(prec = list(fixed = TRUE, initial = log(1e3)))
formula = y ~ -1 + f(GIA, model = spde)
prec_scale <- c(1/errors, rep(1, 100))
res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = spde), family = "gaussian",
                 scale =prec_scale, 
                 control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))
summary(res_inla)

## Plot the posteriors of the parameters
result <- inla.spde2.result(res_inla, "GIA", spde)
par(mfrow= c(1,2))
plot(result[["marginals.range.nominal"]][[1]], type = "l",
     main = "Nominal range, posterior density")
plot(result[["marginals.variance.nominal"]][[1]], type = "l",
     main = "Nominal variance, posterior density")

## Plot the predicted GIA field mean and variance
pidx <- inla.stack.index(stGIA, tag = "pred")
GIA_mpost <- res_inla$summary.random$GIA$mean + xSim
GIA_spost <- res_inla$summary.random$GIA$sd

xyg_mpost <- res_inla$summary.linear.predictor$mean[pidx$data]
xyg_spost <- res_inla$summary.linear.predictor$sd[pidx$data]

y_mpost <- res_inla$summary.linear.predictor$mean[1:(m*ff)] + A%*%xSim
y_spost <- res_inla$summary.linear.predictor$sd[1:(m*ff)]


## Plot the variance
par(mfrow = c(2,1))
theta_mode <- result$summary.theta$mode
Q2 <- inla.spde2.precision(spde, theta = theta_mean)
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(sqrt(1/diag(Q2)))), col = topo.colors(40),
           xlim = c(0,1), ylim = c(0,1),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field")
points(yyloc[,1], yyloc[,2], cex = y_spost, pch = 1)

## The standard error
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_spost)), col = topo.colors(40),
           xlim = c(0,1), ylim = c(0,1),
           xlab = "Longitude", ylab = "Latitude", main = "Updated posterior error field")
points(yyloc[,1], yyloc[,2], cex = y_spost, pch = 1)


## Plot the error on mean
par(mfrow = c(2,1))
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(xSim)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field",
           xlim = c(0, 1), ylim = c(0,1))
points(yyloc[id1,], col = "blue", pch = 19)
points(yyloc[id2,], col = "red", pch = 19)
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_mpost)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field",
           xlim = c(0, 1), ylim = c(0,1))
points(xyg[,1], xyg[,2], cex = xyg_spost)
points(yyloc[,1], yyloc[,2], cex = y_spost, pch = 1, col = 2)
