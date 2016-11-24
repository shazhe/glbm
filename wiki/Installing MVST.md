# Instruction for installing MVST on R > 3.00

We document the instruction for installing MVST on R > 3.00 here. Instruction will be updated with the package, dependent packages and and new release of R. 

The MVST package can be installed from the gitHub repo https://github.com/andrewzm/MVST. Current version 1.0 (Last commit on 1st OCT, 2015).

### Step 0: If you are a Windows user, otherwise skip to Step 1
#### Install Rtools for building R packages under Windows
- Download the compatible version from https://cran.r-project.org/bin/windows/Rtools/ 
- Make sure you check the box for the path, details to be found at https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows
- Note, you can have multiple versions of Rtools for multiple versions of R, See https://stat.ethz.ch/pipermail/r-sig-windows/2016q4/000055.html


### Step 1: Install package "gcplib" from source

That's why you need Rtools :joy: The package is on CRAN but is deprecated. New release is going to remove the pacakge and replace it by "rgeos". Not sure where is the bit that MVST depends on gcplib. If known, we can probably try to update it.

- Download the package source from https://cran.r-project.org/web/packages/gpclib/index.html
- Run the followling line in R 
```R
install.packages("MyPathTo/gpclib_1.5-5.tar.gz", repos = NULL, type = "source")
```
where "MyPathTo" is the path to your downloaded source file.

### Step 2: Make sure you have the following list of packages installed

All available on CRAN and can be directly installed.

- devtools
- knitr
- igraph
- R.utils
- maps
- magic

You might need more spatial tools in R which can be found at the following two task views
- Spatial https://cran.r-project.org/web/views/Spatial.html
- Spatial Temporal https://cran.r-project.org/web/views/SpatioTemporal.html

You can install all packages under the task view by running 
```R
install.views("Spatial")

```
But this may take a while and not really necessary.

### Step 3: Install "MVST"!
Finally, we can install "MVST" smoothly by running
```R
library(devtools)
install_github("andrewzm/MVST",build_vignettes=F,dependencies=T)
```

### Linux and others
Some notes on installing "rgdal" for Linux user are given at the gitHub repo.
