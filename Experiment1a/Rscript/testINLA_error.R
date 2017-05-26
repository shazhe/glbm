###############################
#### INLA test error       ####
#### plug-in approximation ####
###############################
wkdir <- getwd()
## Priors mean and variance for rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

#### Set up text experiment inputs
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]

## set seed
set.seed(6)


### 1 Test pointwise error size
runserver = TRUE
exname <- "/experimentBHM/PointErr_"
err <- rnorm(nrow(GPS_obsU), sd = 1)
load("experimentBHM/Mesh_GIAs.RData")
## Find a single and isolated point to test on with large and small error 
## cho0se by eye -- no. 417 and to illustrate the effect choose it to have the same sign 
## as the nearest points --no.128
err1 <- err2 <- err
err1[417] <- sign(err[128])*0.001
err2[417] <- sign(err[128])*5
source("glbm/Experiment1a/Rscript/4GP_inla_testErr.R")


### 2 Test general error size
runserver = TRUE
meshSize <- c("s", "m", "l")
errSize <- c("s", "L")
errors <- c(0.01, 2)

for(j in 1:2){
  err <- rnorm(nrow(GPS_obsU), sd = errors[j])
  for (i in 1:3){
    exname <- paste0("/experimentBHM/",meshSize[i], "Mesh_", errSize[j], "Err_")
    load(paste0("experimentBHM/Mesh_GIA", meshSize[i], ".RData"))
    source("glbm/Experiment1a/Rscript/4GP_Bayes_update_inla.R")
  }
}


### 3 Test effect of mesh size
## Assume has alreading run scripts for different mesh sizes in Test 1
runserver = TRUE
exname <- "/experimentBHM/MeshSize_"
source("glbm/Experiment1a/Rscript/4GP_inla_testMesh.R")