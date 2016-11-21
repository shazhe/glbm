## Test script for generating global mesh using pakcage INLA

## If the packages is not installued use
## install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
library(INLA)

## A simplest way to generate coarse and uniform mesh on the globe
mesh.g <- inla.mesh.create(globe = 10)
plot(mesh.g)
## For a 3d interactive plot use rgl, if rgl not installed, run
## install.packages("rgl")
plot(mesh.g, rgl = TRUE)
