####################################################################### 
## Install R packages in a cluster environment with job scheduler    ##
## Run this script before runing any BHM scripts.                    ##
## Script copied from https://orfe.princeton.edu/help/r-packages     ##
## Note BlueCrystal does not allow installing from remote reop,      ##
## so all required packages are installed from source files under    ##
## ../myPkgs                                                         ##
#######################################################################

## Create the personal library if it doesn't exist. 
## Ignore a warning if the directory already exists.
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

## Packages need special care
## rgdal -- need to load module gdal and proj and user config path
install.packages('~/myPkgs/rgdal_1.2-5.tar.gz', type = "source", 
configure.args= c('--with-proj-include=/cm/shared/libraries/gnu_builds/proj-4.9.3/include',
'--with-proj-lib=/cm/shared/libraries/gnu_builds/proj-4.9.3/lib'))

install.packages('~/myPkgs/gpclib_1.5-5.tar.gz', type = "source")

## Some basic packages
pkgs1 <- c("R.utils", "magic", "rgeos")
install.packages(pkgs1, Sys.getenv("R_LIBS_USER"), repos = "http://www.stats.bris.ac.uk/R/", 
dependencies = TRUE)

## Packages MVST depends on
pkgs2 <- c("spam", "deldir", "SDMTools", "network", "fields", "matlab", "actuar", "akima", "geometry", "GEOmap")
install.packages(pkgs2, Sys.getenv("R_LIBS_USER"), repos = "http://www.stats.bris.ac.uk/R/", 
dependencies = TRUE)

## Finally install INLA and MVST
install.packages('INLA', Sys.getenv("R_LIBS_USER"),
 repos="https://www.math.ntnu.no/inla/R/stable")

install.packages('~/myPkgs/MVST_1.0.1.tar.gz', type = "source")

## Submit the command R CMD BATCH install-Rpkg.R to the cluster's queue or job scheduler.

## Confirm installation by listing the contents of ~/R.

## Retrieve any error messages from install-Rpkg.Rout, 
## which is generated as a result of running install-Rpkg.R.



#### Some other examples
## Install a package that you have copied to the remote system.
## Standard packages
#system("ls ~/myPkgs/std | grep 'tar.gz' > pkgnames.txt")
#pkgs <- paste0("~/myPkgs/std/", scan("pkgnames.txt", "char"))
#install.packages(pkgs, Sys.getenv("R_LIBS_USER"), repos = NULL)


#install.packages("~/myPkgs/gpclib_1.5-5.tar.gz", Sys.getenv("R_LIBS_USER"))
#install.packages("~/myPkgs/MVST_1.0.1.tar.gz", Sys.getenv("R_LIBS_USER"))

## install.packages("~/myPkgs/MVST_1.0.1.tar.gz", Sys.getenv("R_LIBS_USER"))


## Install packages from other source
## install.packages("INLA", Sys.getenv("R_LIBS_USER"), 
##                 repos="https://www.math.ntnu.no/inla/R/stable")

## Install Views
## library("ctv")
## install.views("Spatial")

## Install the packages.
## install.packages(c("rgdal", "rgeos", "dplyr"), Sys.getenv("R_LIBS_USER"),
## repos = "https://cran.r-project.org")


