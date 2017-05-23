###############################
#### INLA test error       ####
#### plug-in approximation ####
###############################
wkdir <- getwd()
exname <- "/smallError_"
## Priors mean and variance for rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

#### Set up text experiment inputs
## 1 Test measurement error size with small Mesh
load("experimentBHM/Mesh_GIAs.RData")
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]




err1 <- rnorm(nrow(GPS_obsU), sd = 1.5)

