###################################################################
##          A Bayesian Gaussian Update                           ##
## This script performs a full Bayesian Gaussian Update          ##
## of the GIA processes                                          ##
## Given                                                         ##
## 1. prior assumptions of the parameters                        ##
## 2. subset of GPS obs                                          ##
## 3. GIA forward model priors                                   ##
###################################################################
## Pseudo Code of the "Slice within Gibbs"                       ##
##  1. Model                                                     ##
## Y | X, e ~ N(Ay %*% X, I * e)                                 ##
## X - Xm ~ GP (0, Sigma(rho, sigma))                            ##
## Xm = X_GIA -- offset                                          ##
## rho ~ U(a, b)                                                 ##
## sigma, e ~ IG (alpha, beta)                                   ##
##  2. Targets: X|Y, rho|Y, Sigma | y, e|Y                       ##
##  3.Algorithm (no coloring!)                                   ##
## r0, s0, e0 intial values, nsamples                            ##
## At Step i                                                     ##
## (1) Update X(i)|sigma(i-1), rho(i-1), e(i-1), y               ##
## (2) Update rho(i)|X(i),y                                      ##
## (3) Update sigma(i)| X(i), y                                  ##
## (4) Update e(i) | X(i), y                                     ##    
###################################################################
#### Setup

### load packages
library(GEOmap)
library(INLA)
if(runserver){
  INLA:::inla.dynload.workaround()
}#Use this on server with old gcc library
library(MVST)
library(fields)
library(slice)
library(actuar)
library(spam)
library(parallel)
library(coda)
library(lattice)
source("glbm/Experiment1a/Rscript/MVSTplus.R")

### set constants
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obsU$y_center, lon = GPS_obsU$x_center))
GPSX <- ifelse(GPS_obsU$x_center > 180, GPS_obsU$x_center-360, GPS_obsU$x_center)
GPSY <- GPS_obsU$y_center

CMat <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
nObs <- nrow(CMat)
nMesh <- ncol(CMat)
x_mu <- matrix(Mesh_GIA_sp@data$GIA_m, nrow = nMesh, ncol = 1)
igshape_new <- igshape0 + nObs/2

### Generate synthetic GPS data
y0 <- CMat %*% x_mu
ydata <- y0 + err

### Transform hyper parameters 
## GMRF Matern covariance parameters
## sigma ~ lognormal(mu_s, v_s)
## rho ~ lognormal(mu_r, v_r)
## Transform to log scale
## log(sigma) = log(sigma0) + theta_1 ~ normal(lmu_s, lv_s)
## log(rho) = log(rho0) + theta_2 ~ normal(lmu_r, lv_r)
tsigma <- Tlognorm(mu_s, v_s)
trho <- Tlognorm(mu_r, v_r)
lmu_s <- tsigma[1]
lv_s <- tsigma[2] # var(theta_1)
lmu_r <- trho[1]
lv_r <- trho[2] # var(theta_2)
## Priors for the SPDE parameters log(kappa) and log(tau):
## log(kappa) = log(kappa0) - theta_2 ~ normal(log(kappa0), lv_r)
## log(tau) = log(tau0) - theta1 + theta2 ~ norm(log(tau0), lv_s + lv_r)
## couple (rho, sigma) and (kappa, tau) through SPDE connection 
lkappa0 <- log(8)/2 - lmu_r
ltau0 <- 0.5*log(1/(4*pi)) - lmu_s - lkappa0
## Build the SPDE
GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/lv_s), sqrt(1/lv_r)))



mcmcGIA <- function(initial_vals){
  ### Initial values
  ## Observations
  sigma_e0 <- initial_vals[1]
  Q_obs0 <- Matrix(0, nrow = nObs, ncol = nObs, sparse = TRUE)
  diag(Q_obs0) <- 1/sigma_e0
  ## Latent GMRF
  theta0 <- initial_vals[2:3]
  Q_GIA0 <- inla.spde.precision(GIA_spde, theta=theta0)
  
  ### Allocate storage
  x_samp <- matrix(0,numsamples, nMesh)
  e_samp <- rep(0, numsamples)
  theta12_samp <- matrix(0, numsamples,2)
  
  
  ### MCMC alogrithms
  if(sampler == "MH_RW"){
    rwsd <- 1
  }
  
  if(sampler == "slice1"){
    slice_theta1 <- slice_theta2 <- slice(d = 1)
    nlearn <- 100
    lscale1 <- lscale2 <- rep(0, nlearn)
  }
  
  if(sampler == "slice2"){
    slice_theta <- slice(d = 2)
    nlearn <- 100 # learning step for the slice sampler
    lscale <- rep(0, nlearn)
  }
  
  ### Start MCMC simulation
  Q_obs <- Q_obs0
  Q_GIA <- Q_GIA0
  theta_t <- theta0
  mmchol <- summary(chol(as.spam.dgCMatrix(Q_GIA)))
  m <- mm <- 0
  t1 <- proc.time()
  pb <- txtProgressBar(min = 0, max = numsamples, style = 3)
  while (m  < numsamples){
    ### 1 Update the latent process
    Q_new <- as.spam.dgCMatrix(crossprod(CMat, Q_obs) %*% CMat + Q_GIA)
    bt <- crossprod(CMat, Q_obs)%*%ydata + Q_GIA %*% x_mu
    cholQ_new <- chol(Q_new,memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices))
    zm_new <- backsolve(cholQ_new, cbind(rnorm(nMesh), forwardsolve(cholQ_new, bt)))
    x_new <- rowSums(zm_new)
    
    ### 2 Update the measurement error
    res <- ydata - CMat%*%x_new
    igscale_new <- as.numeric(igscale0 + crossprod(res)/2)
    e_new <- rinvgamma(1,shape=igshape_new, scale = igscale_new)
    
    ### 3 Update the SPDE parameters
    z_GIA <- x_new - x_mu
    
    if(sampler == "MH_RW"){
      log_postcond_theta <- function(theta){
        Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=theta))
        ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
        as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + crossprod(theta/c(sqrt(lv_s), sqrt(lv_r)))) + ldetQ)
      }
      
      r <- 1
      a <- 0
      while(r > a){
        theta_p <- theta_t + rnorm(2, sd = 0.07)
        a <- exp(log_postcond_theta(theta_p) - log_postcond_theta(theta_t))
        r <- runif(1)
      }
      theta_t <- theta_p
    }
    
    ## bivariate slice sampler
    if(sampler == "slice2"){
      log_postcond_theta <- function(theta){
        Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=theta))
        ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
        as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + crossprod(theta/c(sqrt(lv_s), sqrt(lv_r)))) + ldetQ)
      }
      if(m <= nlearn){
        theta_p <- slice_theta(theta_t, log_postcond_theta, learn = TRUE)
        lscale[m] <- getscales(slice_theta)
      }else{
        theta_p <- slice_theta(theta_t, log_postcond_theta, learn = FALSE)
      }
      theta_t <- theta_p
    }
    
    ## univariate slice sampler
    if(sampler == "slice1"){
      log_postcond_theta1 <- function(theta1){
        Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=c(theta1, theta_t[2])))
        ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
        as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + theta1^2/lv_s) + ldetQ)
      }
      if(m <= nlearn){
        theta1_p <- slice_theta1(theta_t[1], log_postcond_theta1, learn = TRUE)
        lscale1[m] <- getscales(slice_theta1)
      }else{
        theta1_p <- slice_theta1(theta_t[1], log_postcond_theta1, learn = FALSE)
      }
      
      
      log_postcond_theta2 <- function(theta2){
        Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=c(theta1_p, theta2)))
        ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
        as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + theta2^2/lv_r) + ldetQ)
      }
      if(m <= nlearn){
        theta2_p <- slice_theta2(theta_t[2], log_postcond_theta2, learn = TRUE)
        lscale2[m] <- getscales(slice_theta2)
      }else{
        theta2_p <- slice_theta2(theta_t[2], log_postcond_theta2, learn = FALSE)
      }
      theta_t <- c(theta1_p, theta2_p)
    }
    
    ### 4 Store samples and new values
    if(mm >= burnin){
      if((mm-burnin)%%thinning == 0){
        m <- m +1
        x_samp[m,] <- x_new
        e_samp[m] <- e_new
        theta12_samp[m,] <- theta_t
      }
    }
    mm <- mm +1
    Q_GIA <- inla.spde.precision(GIA_spde, theta=theta_t)
    diag(Q_obs) <- 1/e_new
    
    ### Print the Progression Bar
    setTxtProgressBar(pb, m)
    
  }
  
  ### close the progress bar on finishing
  close(pb)
  
  t2 <- proc.time()
  ttime <- t2-t1
  
  return(list(time = ttime, samples = list(sigma_e = e_samp, theta1 = theta12_samp[,1],
                                           theta2 = theta12_samp[,2], xGIA = x_samp)))
  
}



res <- mclapply(ini_vals, mcmcGIA)

sigma_e <- mcmc.list(lapply(res, function(x) as.mcmc(x$samples$sigma_e, start = burnin+1, 
                                                     end = burnin+numsamples*thin, thin = thin)))
theta_1 <- mcmc.list(lapply(res, function(x) as.mcmc(x$samples$theta1, start = burnin+1, 
                                                     end = burnin+numsamples*thin, thin = thin)))
theta_2 <- mcmc.list(lapply(res, function(x) as.mcmc(x$samples$theta2, start = burnin+1, 
                                                     end = burnin+numsamples*thin, thin = thin)))
xGIA <- list()
for (i in c(1, 50, 800)){
  xGIA[[i]] <- mcmc.list(lapply(res, function(x) as.mcmc(x$samples$xGIA[,i], start = burnin+1, 
                                                         end = burnin+numsamples*thin, thin = thin)))
}

#### plot Sample analysis
## traceplot
pdf(file = paste0(wkdir, exname, "_MCMCanalysis.pdf"), width = 8, height = 6)
par(mfrow = c(2,3))
traceplot(sigma_e, main = expression(sigma[e]^2))
traceplot(theta_1, main = expression(theta[1]))
traceplot(theta_2, main = expression(theta[2]))

traceplot(xGIA[[1]], main = expression(x[1]))
traceplot(xGIA[[50]], main = expression(x[50]))
traceplot(xGIA[[800]], main = expression(x[800]))


#plot.new()
## sample correlations
#par(mfrow = c(2,3))
acfplot(sigma_e, main = expression(sigma[e]^2))
acfplot(theta_1, main = expression(theta[1]))
acfplot(theta_2, main = expression(theta[2]))

acfplot(xGIA[[1]], main = expression(x[1]))
acfplot(xGIA[[50]], main = expression(x[50]))
acfplot(xGIA[[800]], main = expression(x[800]))

#plot.new()
## density plot
#par(mfrow = c(2,3))
densityplot(sigma_e, expression(sigma[e]^2))
densityplot(theta_1, main = expression(theta[1]))
densityplot(theta_2, main = expression(theta[2]))

densityplot(xGIA[[1]], main = expression(x[1]))
densityplot(xGIA[[50]], main = expression(x[50]))
densityplot(xGIA[[800]], main = expression(x[800]))
dev.off()

#### The GIA field
## Posterior mean
GIApost_mean <- lapply(res, function(x) colMeans(x$samples$xGIA))
GIApost_sd <- lapply(res, function(x) apply( x$samples$xGIA, 2, sd))

proj1 <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(361,181))

for (i in (1:n.chains)){
pdf(file = paste0(wkdir, exname, "_GIAfield", i, ".pdf"), width = 12, height = 12)
par(mfrow = c(3,1))
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(Mesh_GIA_sp@data$GIA_m)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "The ICE6g")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIApost_mean[[i]])), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- mean")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIApost_sd[[i]])), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- sd")
points(GPSX, GPSY, pch = 20)
dev.off()
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

pdf(file = paste0(wkdir, exname, "_hyperpars.pdf"), width = 8, height = 4)
for (i in (1:n.chains)){
e_samp_sqrt <- sqrt(sigma_e[[i]])
sigma_e_mode <- getmode(e_samp_sqrt)
range_samp <- exp(lmu_r + theta_2[[i]])
rho_mode <- getmode(range_samp)
sigma_samp <- exp(lmu_s + theta_1[[i]])
sigma_mode <- getmode(sigma_samp)

par(mfrow=c(1,3))
plot(density(e_samp_sqrt), main = bquote(bold(sigma[e]("mode") == .(round(sigma_e_mode, 4)))))
plot(density(range_samp), main = bquote(bold(rho("mode") == .(round(rho_mode, 4)))))
plot(density(sigma_samp), main = bquote(bold(sigma("mode") == .(round(sigma_mode, 4)))))

}
dev.off()

save.image(file = paste0(wkdir, exname, ".RData"))





