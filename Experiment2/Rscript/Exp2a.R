
## ----loadlibs, message = FALSE, warning = FALSE--------------------------
library(rgdal); library(sp);library(GEOmap)
library(INLA)
source("~/glbm/BHM_sphere/functions.R")
source("~/glbm/BHM_sphere/functions-barriers-dt-models-march2017.R")

## ----alt_data------------------------------------------------------------
library(ncdf4)
alt_nc <- nc_open("/./projects/GlobalMass/WP3-Ocean/BHMinputs/trend_SSH_CCI_200501_201512.nc")
print(alt_nc)
lon <- ncvar_get(alt_nc, "lon")
lat <- ncvar_get(alt_nc, "lat")
trend_ssh <- ncvar_get(alt_nc, "trend_ssh") #note that there re NAs for land datat
err_ssh <- ncvar_get(alt_nc, "err")

alt_data <- data.frame(trend_ssh = as.numeric(trend_ssh), err_ssh = as.numeric(err_ssh))
alt_data$lon <- rep(lon, 180) 
alt_data$lat <- rep(lat, each = 360)
alt_data2 <- na.omit(alt_data)
## Find the xyz coords of the altimetry data
alt_loc <- do.call(cbind, Lll2xyz(lon = alt_data2$lon, lat = alt_data2$lat))


## ----prior---------------------------------------------------------------
## Priors mean and variance for the parameters: rho and sigma
mu_r <- 3000/6371
v_r <- 1
mu_s <- 5
v_s <- 10^2

## ----ssh_mesh------------------------------------------------------------
## Genereate Fibonacci points on the sphere
fibo_points <- fiboSphere(N = 12960, L0 = TRUE)
fibo_points_xyz <- do.call(cbind, Lll2xyz(lat = fibo_points[,2], lon = fibo_points[,1]))
mesh0 <- inla.mesh.2d(loc = fibo_points_xyz, cutoff = 0.01, max.edge = 1)
## Make this "smoother"
mesh0 <- inla.mesh.2d(loc = mesh0$loc, cutoff = 0.01, max.edge = 1)

## ----ssh_mesh2-----------------------------------------------------------
## Load the Ocean polygon
Ocean <- readOGR(dsn = "/./projects/GlobalMass/WP1-BHM/maps/ne_110m_ocean", layer = "ne_110m_ocean")

## Remove mesh not in the ocean
mesh0 <- dt.mesh.addon.posTri(mesh = mesh0, globe = TRUE) 
Tlonlat <- Lxyz2ll(list(x = mesh0$posTri[,1], y = mesh0$posTri[,2], z = mesh0$posTri[,3]))
mesh0$Trill <- cbind(lon = Tlonlat$lon, lat =Tlonlat$lat)
TinOcean <- unlist(over(Ocean, SpatialPoints(coords=mesh0$Trill, proj4string = CRS(proj4string(Ocean))), returnList=T))
TAll <- 1:mesh0$t
ToutOcean <- TAll[-TinOcean]
Omega = dt.Omega(list(TinOcean, 1:mesh0$t), mesh0)

mesh_ssh <- mesh.sub(mesh0, Omega, 1)
summary(mesh_ssh)

## ----spde----------------------------------------------------------------
## Transform the parameters for the SPDE_GMRF approximation
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

## Build the SPDE model with the prior
lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

SSH_spde <- inla.spde2.matern(mesh_ssh, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))


## ----link_data-----------------------------------------------------------
## Link the process to observations and predictions
A_data <- inla.spde.make.A(mesh = mesh_ssh, loc = alt_loc)
A_pred <- inla.spde.make.A(mesh = mesh_ssh, loc = rbind(mesh_ssh$loc))

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=alt_data2$trend_ssh), A = list(A_data),
                     effects = list(SSH = 1:SSH_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(SSH=1:SSH_spde$n.spde), tag = "pred")
stSSH <- inla.stack(st.est, st.pred)

## ----inla_run, include = TRUE, eval = FALSE------------------------------
## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))
prec_scale <- c(1/alt_data2$err_ssh^2, rep(1, nrow(A_pred)))

## The formular for modelling the SSH mean
formula = y ~ -1 +  f(SSH, model = SSH_spde)

## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stSSH, spde = SSH_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stSSH), compute =TRUE))

## ----inla_res, include = TRUE, eval = FALSE------------------------------
INLA_pred <- res_inla$summary.linear.predictor
## Extract and project predictions
pred_idx <- inla.stack.index(stSSH, tag = "pred")$data

## SSH
SSH_m <- INLA_pred$mean[pred_idx]
SSH_u <- INLA_pred$sd[pred_idx]
proj <- inla.mesh.projector(mesh_ssh, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-89.5, 89.5))
SSH_grid <- expand.grid(proj$x, proj$y)
SSH_pred <- data.frame(lon = SSH_grid[,1], lat = SSH_grid[,2],
                       mean = as.vector(inla.mesh.project(proj, as.vector(SSH_m))),
                       u = as.vector(inla.mesh.project(proj, as.vector(SSH_u))))

res_SSH <- list(res_inla = res_inla, spde = SSH_spde, st = stSSH,
            mesh = mesh_ssh,  SSH_pred = SSH_pred)

save(res_SSH, file ="/./projects/GlobalMass/WP1-BHM/Experiment2a/exp2a_ssh.RData")

## ----grace_data----------------------------------------------------------
grace_data <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_trends_v03.txt", header = T)
grace_loc <-  read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_loc_v03.txt", skip = 1)
n_grace <- ncol(grace_loc)
n_ll <- (n_grace-4)/2
names(grace_loc) <- c("id", "area", "lon_c", "lat_c", paste0(c("lon", "lat"), rep(1:n_ll, each = 2)))

## Create spatial polygons data frame
Polygon_list <- list()
for(i in 1:nrow(grace_loc)){
  lons <- na.omit(as.numeric(c(grace_loc[i, seq(5,ncol(grace_loc), 2)], grace_loc[i, 5])))
  lats <- na.omit(as.numeric(c(grace_loc[i, seq(6,ncol(grace_loc), 2)], grace_loc[i, 6])))
  Polygon_list[[i]] <- Polygon(cbind(lons, lats))
}

Polygons_list <- lapply(1:length(Polygon_list), function(x) Polygons(list(Polygon_list[[x]]), x))
SpPolygon <- SpatialPolygons(Polygons_list, proj4string = CRS("+proj=longlat"))

grace_sp <- SpatialPolygonsDataFrame(SpPolygon,grace_data)

## Priors mean and variance for the parameters: rho and sigma
mu_r <- 2000/6371
v_r <- 1
mu_s <- 16
v_s <- 30^2

## Transform the parameters for the SPDE_GMRF approximation
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

## Build the SPDE model with the prior
lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

mass_spde <- inla.spde2.matern(mesh0, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                               theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))


## ----grace_link----------------------------------------------------------
## For each polygon observation we generate the regular spaced grid and the number of grid cell is proportional to the area of the polygon
grace_area <- geosphere::areaPolygon(grace_sp)/(1000^2)
plot(grace_loc$area, grace_area)
abline(a = 0, b =1)
grace_sp$area <-grace_area
area_mean <- mean(grace_area)

## Generate the integration grid for each polygons
poly_block <- function(i, dis = 10){
  sp_i <- SpatialPolygons(list(grace_sp@polygons[[i]]), proj4string=CRS("+proj=longlat"))
  area_i <- grace_sp$area[i]
  grid_i <- spsample(sp_i, n = round(area_i/dis^2), type = "regular", offset=c(0.5, 0.5))
  ngrid_i <- length(grid_i)
  grid_xyz <- do.call(cbind, Lll2xyz(lat = grid_i@coords[,2], lon = grid_i@coords[,1]))
  block_i <- rep(i, ngrid_i) 
  weights <- rep(area_i/1e4/ngrid_i, ngrid_i)
  return(list(grid_xyz = grid_xyz, block = block_i, weights = weights, ngrid = ngrid_i))
}

grace_block <- lapply(1:nrow(grace_sp), poly_block, dis = 10)

grid_xyz <- do.call(rbind, lapply(grace_block, "[[", "grid_xyz"))
grid_block <- do.call(c, lapply(grace_block, "[[", "block"))
weights <- do.call(c, lapply(grace_block, "[[", "weights"))

A_GRACE_data <- inla.spde.make.A(mesh = mesh0, loc = grid_xyz, block = grid_block,  weights = weights)

## Same for prediction on a 1 degree resolution grid, we need to know the area of the grid for integration
## generate the prediction grid
gx <- seq(0, 359, 1)
gy <- seq(-89.5, 89.5)
grid_ll <- expand.grid(gx, gy)
pred_data <- data.frame(lon = grid_ll[,1], lat = grid_ll[,2])
coordinates(pred_data) <-c("lon", "lat")
gridded(pred_data) <- TRUE
pred_data <- as(pred_data, "SpatialPolygons")
proj4string(pred_data) <- CRS("+proj=longlat")
areas <- geosphere::areaPolygon(pred_data)/(1000^2)
grid_pred <- do.call(cbind,Lll2xyz(lat = grid_ll[,2], lon = grid_ll[,1]))
A_mass_pred <- inla.spde.make.A(mesh = mesh0, loc = grid_pred, weights = areas/1e4)

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=grace_sp$mmweq), A = list(A_GRACE_data),
                     effects = list(mass = 1:mass_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(rbind(A_GRACE_data, A_mass_pred)),
                      effects = list(mass=1:mass_spde$n.spde), tag = "pred")
stmass <- inla.stack(st.est, st.pred)

## ----inla_run_grace, include = TRUE, eval = FALSE------------------------
## ## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))
prec_scale <- c(1/grace_sp$std^2, rep(1, nrow(A_mass_pred) + nrow(A_GRACE_data)))

## The formular for modelling the SSH mean
formula = y ~ -1 +  f(mass, model = mass_spde)

## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stmass, spde = mass_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stmass), compute =TRUE))

## ----inla_res_grace, include = TRUE, eval = FALSE------------------------
INLA_pred <- res_inla$summary.linear.predictor
## Extract and project predictions
pred_idx <- inla.stack.index(stmass, tag = "pred")$data
idx_grace <- pred_idx[1:nrow(A_GRACE_data)]
idx_grid <- pred_idx[-(1:nrow(A_GRACE_data))]

## SSH
mass_m <- INLA_pred$mean[idx_grid]
mass_u <- INLA_pred$sd[idx_grid]
proj <- inla.mesh.projector(mesh0, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-89.5, 89.5))
mass_grid <- expand.grid(proj$x, proj$y)
mass_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
                       mean = mass_m,
                       u = mass_u)

res_mass <- list(res_inla = res_inla, spde = mass_spde, st = stmass,
            mesh = mesh0,  mass_pred = mass_pred)

grace_m <- INLA_pred$mean[idx_grace]
grace_u <- INLA_pred$sd[idx_grace]
grace_sp$pred_mean <- grace_m
grace_sp$pred_u <- grace_u


save(res_mass, grace_sp , file ="/./projects/GlobalMass/WP1-BHM/Experiment2a/exp2a_mass.RData")

## ----gia_data------------------------------------------------------------
ice6g <- read.table("/./projects/WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt", header = T)
ice6g2<- ice6g[order(ice6g$y_center),]
ice6g2 <- ice6g2[order(ice6g$x_center),]

## ----setric--------------------------------------------------------------
steric_m <- res_SSH$SSH_pred$mean - res_mass$mass_pred$mean - ice6g2$trend
steric_u <- res_SSH$SSH_pred$u - res_mass$mass_pred$u
steric_data <- data.frame(lon = ice6g2$x_center, lat = ice6g2$y_center, mean = steric_m, u = steric_u)
save(res_SSH, res_mass, grace_sp, steric_data, "/./projects/WP1-BHM/Experiment2a/exp2a_all.RData")

