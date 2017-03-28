####################################################################### 
## Install R packages in a cluster environment with job scheduler    ##
## Run this script before runing any BHM scripts.                    ##
## Script copied from https://orfe.princeton.edu/help/r-packages     ##
#######################################################################

## Create the personal library if it doesn't exist. 
## Ignore a warning if the directory already exists.
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

## Install the packages.
#install.packages(c("rgdal", "rgeos", "dplyr"), Sys.getenv("R_LIBS_USER"), 
#repos = "https://www.stats.bris.ac.uk/R/" )

## Install a package that you have copied to the remote system.
install.packages("myPkgs/gpclib_1.5-5.tar.gz", Sys.getenv("R_LIBS_USER"))
#install.packages("myPkgs/MVST_1.0.1.tar.gz", Sys.getenv("R_LIBS_USER"))

## Install packages from other source
#install.packages("INLA", Sys.getenv("R_LIBS_USER"), 
#                 repos="https://www.math.ntnu.no/inla/R/stable")

## Submit the command R CMD BATCH install-Rpkg.R to the cluster's queue or job scheduler.

## Confirm installation by listing the contents of ~/R.

## Retrieve any error messages from install-Rpkg.Rout, 
## which is generated as a result of running install-Rpkg.R.