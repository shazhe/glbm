#### In this experiment we do not consider land-ocean boundaries
#### Generate the ocean mesh over the entire sphere and assume the process is also global
#### Then predict NA for points on the land

## load useful functions
source("glbm/Experiment2/Exp2Fun.R")

### 0 load data
## Altimetry = GIA + GRACE + Steric

## 0.1 Altimetry
library(ncdf4)
alt_nc <- nc_open("Z:/WP3-Ocean/Data/Experiment2_inputs/trend_SSH_CCI_200501_201512.nc")
print(alt_nc)
lon <- ncvar_get(alt_nc, "lon")
lat <- ncvar_get(alt_nc, "lat")
trend_ssh <- ncvar_get(alt_nc, "trend_ssh") #note that there re NAs for land datat
err_ssh <- ncvar_get(alt_nc, "err")

alt_data <- data.frame(trend_ssh = as.numeric(trend_ssh), err_ssh = as.numeric(err_ssh))
alt_data$lon <- rep(lon, 180) 
alt_data$lat <- rep(lat, each = 360)
alt_data$lon <- ifelse( alt_data$lon > 180, alt_data$lon - 360, alt_data$lon)

## 0.2 GIA forward model solution
ice6g <- read.table("Z:/WP2-SolidEarth/GIAforwardModels/textfiles/GIA_Pel-6-VM5.txt", header = T)
ice6g$x_center <- ifelse(ice6g$x_center >= 180, ice6g$x_center - 360, ice6g$x_center)

## 0.3 GRACE
grace_trend <- read.table("Z:/WP2-SolidEarth/GRACE/GSFC_global_mascons_v01.1/GRACE_trends_GSFC_global_mascons_v01.1_version01.txt", header = T)
grace_location <- read.table("Z:/WP2-SolidEarth/GRACE/GSFC_global_mascons_v01.1/GRACE_locations_GSFC_global_mascons_v01.1_version01.txt", 
                             fill = TRUE, skip = 1)
names(grace_location) <- c("id", "area", "lonC", "latC", paste0(c("lon", "lat"), rep(1:19,each=2)))
mypoly <- function(x){
  lon <- na.omit(as.numeric(x[seq(5,42,2)]))
  lat <- na.omit(as.numeric(x[seq(6,42,2)]))
  coord <- cbind(lon, lat)
  return(Polygon(coord))
}

aa <- is.na(grace_location)
bb <- (42-rowSums(aa))/2
id <- which(bb>5)
grace_location2 <- grace_location[id,]


grace_polys <- list()
for (i in 1:nrow(grace_location2)){
  grace_polys[[i]] <- mypoly(grace_location2[i,])
}

## 0.4 Priors for the hyper parameters
## Priors mean and variance for the parameters: rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

## Transform the parameters for the SPDE_GMRF approximation
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

## 0.5 Load the coast line polygons
library(rgdal)
lands <- readOGR(dsn = "C:/ZSwork/glbm/Experiment2/ne_110m_land", layer = "ne_110m_land")
spLand <- SpatialPolygons(lands@polygons)
## Remove small islands -- area smaller than 20
landsareas <- sapply(spLand@polygons, function(x) slot(x, "area"))
Llands <- landsareas > 20
spLand2 <- spLand[Llands]
proj4string(spLand2) <- CRS("+proj=longlat")

### 1 generate mesh
## 1.1 Generate the "global Mesh" for steric
library(INLA)
Mesh_steric <- inla.mesh.create(globe= 64)

## 1.2 select the mesh points over the ocean
library(GEOmap)
M_coords <- Lxyz2ll(list(x=Mesh_steric$loc[,1], y = Mesh_steric$loc[,2], z = Mesh_steric$loc[,3]))
Mesh_sp <- SpatialPoints(data.frame(lon = M_coords$lon, lat = M_coords$lat), proj4string = CRS("+proj=longlat")) 
Mesh_land <- over(Mesh_sp, spLand2)
plot(M_coords$lon[is.na(Mesh_land)], M_coords$lat[is.na(Mesh_land)])

## 1.3 Find the Altimetry values for the mesh points
coordinates(alt_data) <- c("lon", "lat")
proj4string(alt_data) <- CRS("+proj=longlat")
gridded(alt_data) <- TRUE
mesh_alt <- over(Mesh_sp, alt_data)


## 1.3 Find the GIA values for the mesh points
coordinates(ice6g) <- c("x_center", "y_center")
proj4string(ice6g) <- CRS("+proj=longlat")
gridded(ice6g) <- TRUE
mesh_gia <- over(Mesh_sp, ice6g)[,c(2,3)]

## 1.4 Find the GRACE values for the mesh points
mesh_grace <- over(Mesh_sp, grace_sp)



### 2 INLA analysis
## Run INLA
## The SPDE model
lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0
GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

## Link the process to observations and predictions
A_data <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
A_pred <- inla.spde.make.A(mesh = Mesh_GIA, loc = rbind(GPS_loc, Mesh_GIA$loc))

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=GPS_data), A = list(A_data),
                     effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

## Fix the GPS errors
hyper <- list(prec = list(fixed = TRUE, initial = 0))
formula = y ~ -1 +  f(GIA, model = GIA_spde)
prec_scale <- c(1/GPS$std^2, rep(1, nrow(A_pred)))
res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))


### 3 prediction
INLA_pred <- res_inla$summary.linear.predictor

## Extract and project predictions
pred_idx <- inla.stack.index(stGIA, tag = "pred")$data
GPS_idx <- pred_idx[1:nrow(GPS)]
GIA_idx <- pred_idx[-c(1:nrow(GPS))]

## GPS 
GPS_u <- INLA_pred$sd[GPS_idx]
GPS_pred <- data.frame(lon = GPSX, lat = GPS$lat, u = GPS_u)

## GIA
GIA_diff <- INLA_pred$mean[GIA_idx] 
GIA_m <- GIA_diff + GIA_prior
GIA_u <- INLA_pred$sd[GIA_idx]
proj <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(720,362), xlim = c(0, 359), ylim = c(-90, 90))
GIA_grid <- expand.grid(proj$x, proj$y)
GIA_pred <- data.frame(lon = GIA_grid[,1], lat = GIA_grid[,2],
                       diff = as.vector(inla.mesh.project(proj, as.vector(GIA_diff))),
                       mean = as.vector(inla.mesh.project(proj, as.vector(GIA_m))),
                       u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))))


