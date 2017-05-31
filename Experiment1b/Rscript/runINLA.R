####################################
#### Experiment 1b Run INLA     ####
#### plug-in approximation      ####
####################################
wkdir <- getwd()
## Priors mean and variance for rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

#### Set up experiment inputs
GPS_obs <- read.table("experimentBHM/GPS_combined_20170530_final.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]

## set seed
set.seed(6)


runserver = TRUE
meshSize <- c("s", "m", "l")
for (i in 1:3){
    exname <- paste0("/experimentBHM/1b",meshSize[i], "Mesh_")
    load(paste0("experimentBHM/Mesh_GIA", meshSize[i], ".RData"))
    source("glbm/Experiment1b/Rscript/2GP_Bayes_update_inla.R")
  }
