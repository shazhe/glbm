####################################
#### Experiment 1b Run INLA     ####
#### plug-in approximation      ####
####################################

#### 0 path and file setup
wkdir <- getwd()
.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.3/")
data_path <- "/./projects/GlobalMass/WP1-BHM/Experiment1b/Data_inputs/"
GIA_input <- "GIA_Pel-6-VM5"
GPS_input <- "GPSdataset_v03b_QC_20170809_comb"

GIA_file <- paste0(data_path, "GIA/", GIA_input, ".txt")
GPS_file <- paste0(data_path, "GPS/", GPS_input, ".txt")


#### 1 Prep Priors 
## Priors mean and variance for rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

## Prep GIA data and mesh




#### 3 Run INLA estimate and predict
set.seed(6)



#### 4 Plot the result



load("experimentBHM/mesh_reg.RData")
load("experimentBHM/GIA_sp_info.RData")

runserver = TRUE
meshSize <- c("s", "m", "l")
for (i in 1:3){
    exname <- paste0("/experimentBHM/1b",meshSize[i], "Mesh_")
    mesh_temp <- mesh_regs[[i]]
    source("glbm/Experiment1b/Rscript/2GP_Bayes_update_inla.R")
  }
