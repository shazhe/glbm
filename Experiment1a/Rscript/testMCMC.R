################################
####  Test MCMC samplers    ####
####                        ####
################################
#### Load mesh and GPS data
wkdir <- getwd()
runserver = TRUE
load("experimentBHM/Mesh_GIAs.RData")
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]
## generate error for sythetic GPS data
err <- rnorm(nrow(GPS_obsU), sd = 2)

#### MCMC parameters
set.seed(7)
numsamples = 100  
burnin = 0
thinning = 1
sampler = "slice1"
n.chains = 3
ini_vals <- list(vals1 = c(1, 0, 0), vals2 = c(0.5, 1, 1), vals3 = c(2, -1, -1))

#### Set priors for the hyper parameters
### Measurement error sigma_e
## sigma_e ~ InverseGamma(shape = 1, rate = 5e-5)
igshape0 <-  1
igscale0 <- 1

### SPDE parameters rho and sigma
## sigma ~ lognormal(mu_s, v_s)
## rho ~ lognormal(mu_r, v_r)
mu_s <- 20
v_s <- 40^2
mu_r <- 500/6371
v_r <- (1000/6371)^2

exname <- paste0("/experimentBHM/", sampler, "_", as.character(numsamples))
source("glbm/Experiment1a/Rscript/4GP_Bayes_update_mcmc.R")
source("glbm/Experiment1a/Rscript/4GP_Bayes_plot_mcmc.R")