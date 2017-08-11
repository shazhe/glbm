###################################################################
##          A Bayesian Gaussian Update (INLA)                    ##
## This script performs a full Bayesian Gaussian Update          ##
## of the GIA processes                                          ##
###################################################################
#### 0 setup codes
source("glbm/Experiment1a/Rscript/MVSTplus.R")

#### 1 Setup for spde parameters
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)
lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

#### 2 Build up the SPDE model
GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

#### 3 Link the process and observations
## Find the mapping between observations and processes basis
A_data <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)

## Create a stack of locations to plot the posterior (resolution = 5 degree apart)
s_lat <- seq(-89.5, 89.5, 10)
s_lon <- seq(-179.5, 179.5, 10)
sll <- expand.grid(x = s_lon, y = s_lat)
sll_loc <- do.call(cbind, Lll2xyz(lat = sll$y, lon = sll$x))
A_pred <- inla.spde.make.A(mesh = Mesh_GIA, loc = rbind(Mesh_GIA$loc, sll_loc))

## Create the estimation and stack
st.est <- inla.stack(data = list(y=GPS_data), A = list(A_data),
                     effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
## Create the prediction stack
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
## Combine the two stack
stGIA <- inla.stack(st.est, st.pred)

## Fix the GPS errors
hyper <- list(prec = list(fixed = TRUE, initial = 0))
formula = y ~ -1 +  f(GIA, model = GIA_spde)
prec_scale <- c(1/GPS_obs$std^2, rep(1, nrow(A_pred)))
res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                   control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))

summary(res_inla)

save(res_inla, GIA_spde, stGIA, file = paste0(outdir, exname, "inla.RData"))
