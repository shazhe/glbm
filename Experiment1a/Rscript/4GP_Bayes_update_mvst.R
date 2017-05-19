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
#### INLA Approximation
library(GEOmap)
library(INLA)
INLA:::inla.dynload.workaround() 
library(MVST)
library(fields)
library(slice)
library(actuar)
library(spam)
load("experimentBHM/Mesh_GIAs.RData")
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
#load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/Mesh_GIA.RData")
#GPS_obs <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obsU$y_center, lon = GPS_obsU$x_center))

#### 1 Setup contants and initial values
### 1.1 Contants
## 1.1.1 data constants
CMat <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
y_obs <- GPS_obsU$trend
x_mu <- matrix(Mesh_GIA_sp@data$GIA_m, nrow = Mesh_GIA$n, ncol = 1)
#x_mu <- matrix(0, nrow = Mesh_GIA$n, ncol = 1)
nMesh <- length(x_mu)
nObs <- length(y_obs)
## 1.1.2 hyper-parameters for the priors
## measurment errors e ~ IG (alpha0, beta0), very vague
alpha0 <-  0.1
beta0 <- 1
alpha_new <- alpha0 + nObs/2
## theta1, theta2 for latent precision, very vague zero mean Gaussian
t_mu <- c(0,0)
t_sd <- c(1, 1)
## 1.2.3  MCMC parameters
numsamples = 5e4  # number of samples
#burnin = 500 
thinning = 25

### 1.2 Initial values for the parameters
## 1.2.1 measurement errors for obs
Q_obs <- Matrix(0, nrow = length(y_obs), ncol = length(y_obs), sparse = TRUE)
diag(Q_obs) <- (1/GPS_obsU$std)^2
#diag(Q_obs) <- 1/100
## 1.2.2 precision matrix for the latent field
## Use the inla parameterization -- why? inla use greate circle distance for S2
## assuming nu = 1, then alpha = 2 when d = 2.
rho0 <- 1000/6371
sigma0 <- 0.5
## theta20 = log(kappa0) = log(8*nu)/2 - log(rho0)
lkappa0 <- log(8)/2 - log(rho0)
## theta10 = log(tau0) = 1/2*log(gamma(nu)/(gamma(alpha)*(4*pi)^(d/2))) - log(sigma0) - nu*log(kappa0)
ltau0 <- 0.5*log(1/(4*pi)) - log(sigma0) - lkappa0
## after transformation the parameters becomes
## log(tau) = log(tau0) - theta1 + nu * theta2
## log(kappa) = log(kappa0) - theta2
## put vague gaussian prior on theta1 and theta2
theta_old <- c(0,0)
## parameterisation 1
## theta1 controls change in log(sigma), theta2 controls change in log(rho)
GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3))
Q_GIA <- inla.spde.precision(GIA_spde, theta=theta_old)

## parameterisation 2 same as one at the beginning
#GIA_spde2 <- inla.spde2.matern(Mesh_GIA)
#Q_GIA2 <- inla.spde.precision(GIA_spde2, theta=c(ltau0, lkappa0))



### 1.3 Initializing the MCMC sampler
## Create MCMC objects to store the samples
## The latent processes
x_samp <- matrix(0,nMesh,numsamples)  #post process samples in one matrix row = n.processes, col = n.samples
## The prameters
e_samp <- rep(0, numsamples)
theta12_samp <- matrix(0,2,numsamples)
#slice_theta <- slice(d = 2)
slice_theta1 <- slice_theta2 <- slice(d = 1)
nlearn <- 100 # learning step for the slice sampler
#lscale <- rep(0, nlearn)
lscale1 <- lscale2 <- rep(0, nlearn)
mmchol <- summary(chol(as.spam.dgCMatrix(Q_GIA)))
m <- mm <- 0
#### The Slice within Gibbs Sampler
### Setup a progress bar
set.seed(2)

t1 <- proc.time()
pb <- txtProgressBar(min = 0, max = numsamples, style = 3)
while (m  <= numsamples){
  ### 1 Update the latent process
  Q_new <- as.spam.dgCMatrix(crossprod(CMat, Q_obs) %*% CMat + Q_GIA)
  bt <- crossprod(CMat, Q_obs)%*%y_obs + Q_GIA %*% x_mu
  cholQ_new <- chol(Q_new,memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices))
  zm_new <- backsolve(cholQ_new, cbind(rnorm(nMesh), forwardsolve(cholQ_new, bt)))
  x_new <- rowSums(zm_new)
 
  ### 2 Update the measurement error
  res <- y_obs - CMat%*%x_new
  beta_new <- as.numeric(beta0 + crossprod(res)/2)
  e_new <- rinvgamma(1,shape=alpha_new, scale = beta_new)
  
  ### 3 Update the SPDE parameters
  z_GIA <- x_new - x_mu
  
  ## Rejection
  log_postcond_theta <- function(theta){
    Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=theta))
    ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
    as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + crossprod(theta/t_sd)) + ldetQ)
  }

  r <- 1
  a <- 0
  while(r > a){
    theta_p <- theta_old + rnorm(2, sd = 0.07)
    a <- exp(log_postcond_theta(theta_p) - log_postcond_theta(theta_old))
    r <- runif(1)
  }
  theta_new <- theta_p

  # ## slice sampler
  # log_postcond_theta <- function(theta){
  #  Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=theta))
  #  ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
  #  as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + crossprod(theta/t_sd)) + ldetQ)
  # }
  # if(m <= nlearn){
  #  theta_new <- slice_theta(theta_old, log_postcond_theta, learn = TRUE)
  #  lscale[m] <- getscales(slice_theta)
  # }else{
  #  theta_new <- slice_theta(theta_old, log_postcond_theta, learn = FALSE)
  # }
  # 
  
  
  # log_postcond_theta1 <- function(theta1){
  #   Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=c(theta1, theta_old[2])))
  #   ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
  #   as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + (theta1/t_sd[1])^2) + ldetQ)
  # }
  # if(m <= nlearn){
  #   theta1_new <- slice_theta1(theta_old[1], log_postcond_theta1, learn = TRUE)
  #   lscale1[m] <- getscales(slice_theta1)
  # }else{
  #   theta1_new <- slice_theta1(theta_old[1], log_postcond_theta1, learn = FALSE)
  # }
  # 
  # 
  # log_postcond_theta2 <- function(theta2){
  #   Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=c(theta1_new, theta2)))
  #  ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
  #   as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + (theta2/t_sd[2])^2) + ldetQ)
  # }
  # if(m <= nlearn){
  #   theta2_new <- slice_theta2(theta_old[2], log_postcond_theta2, learn = TRUE)
  #   lscale2[m] <- getscales(slice_theta2)
  # }else{
  #   theta2_new <- slice_theta2(theta_old[2], log_postcond_theta2, learn = FALSE)
  # }
  # theta_new <- c(theta1_new, theta2_new)

  
  ### 4 Store samples and new values
  if(mm%%thinning == 0){
    m <- m +1
    x_samp[,m] <- x_new
    e_samp[m] <- e_new
    theta12_samp[,m] <- theta_new
  }
  
  mm <- mm +1
  Q_GIA <- inla.spde.precision(GIA_spde, theta=theta_new)
  diag(Q_obs) <- 1/e_new
  theta_old <- theta_new
  
  ### Print the Progression Bar
  setTxtProgressBar(pb, m)
 
  }
 
### close the progress bar on finishing
close(pb)

t2 <- proc.time()
ttime <- t2-t1
save(rho0, sigma0, x_samp, e_samp, theta12_samp, ttime, file = "rejectLon2.RData")












#### Sample analysis
## traceplot
par(mfrow = c(2,3))
plot(e_samp, type = "l")
plot(theta12_samp[1,], type  ="l")
plot(theta12_samp[2,], type  ="l")

plot(x_samp[1,], type = "l")
plot(x_samp[500,], type = "l")
plot(x_samp[1000,], type = "l")


## sample correlations
par(mfrow = c(2,3))
acf(e_samp)
acf(theta12_samp[1,])
acf(theta12_samp[2,])

acf(x_samp[1,])
acf(x_samp[500,])
acf(x_samp[1000,])


## density plot
par(mfrow = c(2,3))
plot(density(e_samp))
plot(density(theta12_samp[1,]))
plot(density(theta12_samp[2,]))

plot(density(x_samp[1, ]))
plot(density(x_samp[500, ]))
plot(density(x_samp[1000, ]))


#### The GIA field
## Posterior mean
GIApost_mean <- rowMeans(x_samp[,])
GIApost_sd <- apply(x_samp[,], 1, sd)

proj1 <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(361,181))
GPSX <- ifelse(GPS_obsU$x_center > 180, GPS_obsU$x_center-360, GPS_obsU$x_center)
GPSY <- GPS_obsU$y_center

pdf(file = "reject1GIA.pdf", width = 8.5, height = 11)
par(mfrow = c(3,1))
image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(Mesh_GIA_sp@data$GIA_m)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "The ICE6g")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIApost_mean)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- mean")
points(GPSX, GPSY, pch = 20)

image.plot(proj1$x, proj1$y, inla.mesh.project(proj1, as.vector(GIApost_sd)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- sd")
points(GPSX, GPSY, pch = 20)
dev.off()

pdf(file = "reject1hyperpars.pdf", width = 8.5, height = 11)
par(mfrow = c(3,3))
plot(e_samp, type = "l")
plot(theta12_samp[1,], type  ="l")
plot(theta12_samp[2,], type  ="l")

plot(density(e_samp), main = "error, mean = 2.27")
plot(density(theta12_samp[1,]), main = "theta1")
plot(density(theta12_samp[2,]), main = "theta2")

range_samp <- exp(log(rho0) + theta12_samp[2,])
plot(density(range_samp), main = "range, range = 2176km")

sigma_samp <- exp(log(sigma0) + theta12_samp[1,])
plot(density(sigma_samp), main = "sigma, sigma = 0.6079")
dev.off()

mean(e_samp)
mean(theta12_samp[1,]) 
mean(theta12_samp[2,])
mean(range_samp)
mean(sigma_samp)

## Posterior sd