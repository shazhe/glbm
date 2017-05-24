###############################
#### INLA test error       ####
#### plug-in approximation ####
###############################
wkdir <- getwd()
exname <- "/experimentBHM/LargeMesh_smallError_"
## Priors mean and variance for rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

#### Set up text experiment inputs
## 1 Test measurement error size with small Mesh
load("experimentBHM/Mesh_GIA.RData")
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]
err <- rnorm(nrow(GPS_obsU), sd = 0.1)

source("glbm/Experiment1a/Rscript/4GP_Bayes_update_inla.R")
