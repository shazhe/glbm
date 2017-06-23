#### This file contains code for generating different meshes on sphere


#### Generate equally distributed points around a sphere -- Fibonacci spiral
fiboSphere1 <- function(n = 1000L) {
  
  phi <- (sqrt(5) + 1) / 2  # golden ratio
  i <- seq(-(n-1), (n-1),2)
  theta <- 2 * pi * i / phi
  sphi <- i/n
  cphi <- sqrt((n+i) * (n-i))/n
  
  x <- cphi * sin(theta)
  y <- cphi * cos(theta)
  z <- sphi
  return(cbind(x,y,z))
}

#### Generate equally distributed points around a sphere -- Fibonacci spiral
fiboSphere2 <- function(n = 1000L) {
  phi <- (sqrt(5) + 1) / 2 - 1 # golden ratio
  ga <- phi * 2 * pi           # golden angle
  
  i <- 1L:n
  lon <- ga * i / (2 * pi)
  lon <- 2 * pi * (lon - floor(lon))
  lon <- ifelse(lon <= pi, lon, lon - 2 * pi)/pi*180
  lat <- asin(-1 + 2 * i / n)/pi*180
  
  cbind(lon = lon, lat = lat)
}
library(GEOmap)

aa <- fiboSphere2(1000)
bb <- do.call(cbind, Lll2xyz(lat = aa[,2], lon = aa[,1]))

library(INLA)
Mesh_b <- inla.mesh.2d(loc = aa, cutoff = 0.001, max.edge = 5) ## small
plot(Mes_b, rgl = T)

