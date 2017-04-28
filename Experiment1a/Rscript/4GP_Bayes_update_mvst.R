###################################################################
##          A Bayesian Gaussian Update                           ##
## This script performs a full Bayesian Gaussian Update          ##
## of the GIA processes                                          ##
## Given                                                         ##
## 1. prior assumptions of the parameters                        ##
## 2. subset of GPS obs                                          ##
## 3. GIA forward model priors                                   ##
###################################################################
## Pseudo Code of the "Slice within Gibbs"                       ##
##  1. Model                                                     ##
## Y | X, e ~ N(Ay %*% X, I * e)                                 ##
## X - Xm ~ GP (0, Sigma(rho, sigma))                            ##
## Xm = X_GIA -- offset                                          ##
## rho ~ U(a, b)                                                 ##
## sigma, e ~ IG (alpha, beta)                                   ##
##  2. Targets: X|Y, rho|Y, Sigma | y, e|Y                       ##
##  3.Algorithm                                                  ##
## r0, s0, e0 intial values, nsamples                            ##
## At Step i                                                     ##
## (1) Update X(i)|sigma(i-1), rho(i-1), e(i-1), y               ##
## (2) Update rho(i)|X(i),y                                      ##
## (3) Update sigma(i)| X(i), y                                  ##
## (4) Update e(i) | X(i), y                                     ##    
###################################################################
#### INLA Approximation
library(rgl)
library(GEOmap)
library(INLA)
#INLA:::inla.dynload.workaround() 
library(MVST)
library(fields)
#load("experimentBHM/Mesh_GIA.RData")
#GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/Mesh_GIA.RData")
GPS_obs <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obsU$y_center, lon = GPS_obsU$x_center))

## Find the mapping between observations and processes basis
CMat <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
GIA_spde <- inla.spde2.matern(Mesh_GIA)
y <- as.vector(GPS_obsU$trend)


#### Fixed parameters for the precision matrix
## 1 Use Priors only
GIA_GMRF<- GMRF(mu = GIA_mu, Q = GIA_Q0) 
GPS_obsDF <- Obs(df=data.frame(x = GPS_loc[,1], y = GPS_loc[,2], w = GPS_loc[,3], z = GPS_obs$trend, 
                               std = GPS_obs$std, t = rep(0, nrow(GPS_loc))))
L <- link(GIA_GMRF, GPS_obsDF, Cmat = CMat)
e <- new("link_list",list(L))
v <- new("block_list",list(G1 = GIA_GMRF, O = GPS_obsDF))
G <- new("Graph",e = e,v = v)
G_reduced <- compress(G)

Olist <- .extractClass(G_reduced@v,"Obs")
Glist <- .extractClass(G_reduced@v,"process")
Obs <- Olist[[1]]
Field <- Glist[[1]]
y_tot <- getDf(Obs)
Qobs <- getPrecision(Obs)
Q_full <- getPrecision(Field)
x_prior <- getMean(Field)
C_full <- G_reduced@e[[1]]@Cmat
C_full_alt_only <- G_reduced_alt_only@e[[1]]@Cmat
x_des <- Field@rep
x_des$n <- 1:nrow(x_des)
x_des$prior <- 0
x_des$mean <- x_des$prior


##########################################################################################################
# GRAPH COLOURING
x_des$block <- x_des$class
x_des$wave <- x_des$colour
bl_list <- sort(unique(x_des$block)) 
numblocks <- length(bl_list)
ib <- X <- L <- Q <- vector("list",numblocks)

for (i in 1:numblocks) ib[[i]] <-  which(x_des$block==bl_list[i])
block_wave_map <- as.data.frame(unique(cbind(x_des$block,x_des$wave)))
names(block_wave_map) <- c("block","wave")
block_wave_map <- block_wave_map[order(block_wave_map$block),]

# Send x_des,C_full,y_tot,block_wave_map, Qobs,Q_full and relevant functions
# save(ib,x_des,C_full,Qobs,Q_full,y_tot,block_wave_map,file="data_for_Gibbs_sampler.rda")

##########################################################################################################

# From here on run on the HPC
numsamples = 2000
dither = 0
burnin = 500
x_samp <- matrix(0,nrow(x_des),numsamples)
x_samp[,1] <- Results$Post_GMRF@mu
in_AP <-  which(x_des$x < -2000 & x_des $ y > 750 & !(x_des$name == "GIA"))
in_cx <- which((x_des$x^2 + x_des$y^2) < 425^2)
x_samp[in_AP,1] <- 0
x_samp[in_cx,1] <- 0

alpha_samp <- rep(alpha_init,numsamples)
gammaG_samp <- rep(1,numsamples)
# eta_samp <- matrix(c(10,2,10,2),4,numsamples)
# eta_accept = 0
# eta_propose <- rep(0,4)
basin_list <- unique(x_des$basins_rignot11)
#gammaI_samp <- lapply(basin_list,function(x) return(rep(100,nsamples)))

# Now get C matrix from ICESAT/ENVISAT
e_IE_only <- new("link_list",list(L1,L2,L5,L6,L10,L11))
v_IE_only <- new("block_list",list(G1=GIA,G2=SURF_FIRN,G3=ICE,O1=ICESAT,O2=ENVISAT))
G_IE_only <- Graph(e=e_IE_only,v=v_IE_only)
G_reduced_IE_only <- compress(G_IE_only)
C_full_IE <- getC(G_reduced_IE_only)
y_tot_IE <- getDf(G_reduced_IE_only)

hyp_eta_mu <- c(-12,2)
hyp_eta_sigma <- c(1,1)
eta_samp <- matrix(c(-12,2),2,numsamples)
var_delta <- exp(eta_samp[1,1] + eta_samp[2,1]*y_tot_IE$fine_scale)
Qobs@x[which(y_tot$obs_name %in% c("ICESAT","ENVISAT"))] <- 1/(y_tot_IE$std^2 + var_delta)


# GRACE priors
# GRACE we expect an order of 1-20 times the reported error
hyp_GRACE = optim(par=c(4,4e-4),Findalphabeta_invgamma, p5=1, p95=20)
hyp_GRACEprec_shape <- hyp_GRACE$par[1]
hyp_GRACEprec_scale <- hyp_GRACE$par[2]
gammaG_bp = hyp_GRACEprec_scale # rate of gamma is scale of inverse gamma
gammaG_ap = hyp_GRACEprec_shape
## Qalpha_p <- 1/(0.1^2)
## alpha_p = 0.07
hyp_alpha = optim(par=c(2,6),Findalphabeta_beta, p95=0.12)
hyp_alpha_shape1 <- hyp_alpha$par[1]
hyp_alpha_shape2 <- hyp_alpha$par[2]

# GRACE maps
mGRACE <- nrow(subset(GRACE,t==1))
XX <- as.matrix(dist(subset(getDf(GRACE),t==1,select=c("x","y"))))
X2 <- XX*(XX<450)
GRACE_neighb <- neighb_from_prec(X2)
GRACE["est_smooth"] <- 0
alpha0=0.1
Smooth_mat <- matrix(0,mGRACE,mGRACE)
# Numerical diffusion from http://pauli.uni-muenster.de/tp/fileadmin/lehre/NumMethoden/WS0910/ScriptPDE/Heat.pdf
for(i in 1:mGRACE) {
  nn <- length(getDf(GRACE)[GRACE_neighb[[i]],]$x)
  Smooth_mat[i,i] <- (1-nn*alpha0)
  Smooth_mat[i,GRACE_neighb[[i]]] <- alpha0
}
D2 <- Smooth_mat
D2[D2 > 0] <- 1
D2 <- D2 - diag(diag(D2))
D1 <- diag(rowSums(D2))

D1_big <- Reduce("bdiag",lapply(t_axis,function(x) {return(D1)}))
D2_big <- Reduce("bdiag",lapply(t_axis,function(x) {return(D2)}))

# Which GRACE elements to update?
intGRACE <- which(y_tot$obs_name == "GRACE")
C_full2 <- as(C_full,"dgTMatrix")
GRACE_in_mat_positions <- which((C_full2@i+1) %in%  intGRACE)
rm(C_full2)

# Which ICESAT elements to update?
intICE <- which(x_des$name %in% c("ICE"))
intBETA <- which(x_des$name %in% c("ICE1","ICE2"))
intICEALL <- which(x_des$name %in% c("ICE","ICE1","ICE2"))
Qtemp <- as(Q_full,"dgTMatrix")
ICE_in_mat_positions <- which(((Qtemp@i+1) %in%  intICEALL) & ((Qtemp@j+1) %in%  intICEALL))
rm(Qtemp)

slice_fn1 <- slice(0.5)
slice_fn2 <- slice(0.5)
slice_fn3 <- slice(0.01)
learn=T 

# transformation of eta
#A <- matrix(c(0.3*0.4,-0.0333,0.4,0.01),2,2)
#b <- matrix(c(0.4*1.5,-0.42),2,1)

load("~/cache/fine_scale_transform.rda")


############
# THE SAMPLER
#############
m = 1
mm = 1
count = 1

while (m < numsamples) {
  m = m+1
  x_des$current_samp <- x_samp[,(m-1)]
  wave_nums <- sort(unique(block_wave_map$wave))
  for (j in wave_nums) {
    # These can be done in parallel
    for(i in subset(block_wave_map,wave==j)$block) {
      Q[[i]] <- t(C_full[,ib[[i]]])%*%Qobs%*%C_full[,ib[[i]]] + Q_full[ib[[i]],ib[[i]]]
      X[[i]] <- cholPermute(Q[[i]])
      
      ib_comp <- unlist(ib[-i])
      ybar <- t(C_full[,ib[[i]]])%*%Qobs%*%(y_tot$z - C_full[,ib_comp]%*%x_des$current_samp[ib_comp]) + Q_full[ib[[i]],ib[[i]]]%*%x_des$prior[ib[[i]]] - 
        Q_full[ib[[i]],ib_comp] %*%(x_des$current_samp[ib_comp] - x_des$prior[ib_comp])
      this_mean <- as.matrix(cholsolve(Q[[i]],ybar,perm=T,cholQp = X[[i]]$Qpermchol, P = X[[i]]$P))
      x_des$current_samp[ib[[i]]] <-   as.vector(this_mean + X[[i]]$P %*% solve(t(X[[i]]$Qpermchol),t(X[[i]]$P)%*%rnorm(length(ib[[i]]))))
      if(!(class(x_des$current_samp) == "numeric")) break
      #print(i)
    }
  }
  x_samp[,m] =x_des$current_samp 
  
  # Slice sampler
  log_cond_posterior1 <- function(tau1) {
    #eta <- solve(A) %*% (rbind(tau1,tau2) - b)
    eta <-  L%*%rbind(tau1,tau2) + b
    var_delta <- exp(eta[1] + eta[2]*y_tot_IE$fine_scale)
    var_obs_delta <- y_tot_IE$std^2 + var_delta
    return(-0.5 * sum(log(var_obs_delta)) - 0.5*sum(residuals^2/var_obs_delta) -
             (eta[1] - hyp_eta_mu[1])^2/(2*hyp_eta_sigma[1]^2) -
             (eta[2] - hyp_eta_mu[2])^2/(2*hyp_eta_sigma[2]^2))
  }
  
  log_cond_posterior2 <- function(tau2) {
    #eta <- solve(A) %*% (rbind(tau1,tau2) - b)
    eta <-  L%*%rbind(tau1,tau2) + b
    var_delta <- exp(eta[1] + eta[2]*y_tot_IE$fine_scale)
    var_obs_delta <- y_tot_IE$std^2 + var_delta
    return(-0.5 * sum(log(var_obs_delta)) - 0.5*sum(residuals^2/var_obs_delta) -
             (eta[1] - hyp_eta_mu[1])^2/(2*hyp_eta_sigma[1]^2) -
             (eta[2] - hyp_eta_mu[2])^2/(2*hyp_eta_sigma[2]^2))
  }
  residuals <- y_tot_IE$z - C_full_IE %*% x_samp[,m] # This is not dependent on eta
  
  if (m > 50) learn <- F
  
  # Sample eta1
  tau_orig <- solve(L) %*% (eta_samp[,(m-1)] - b)
  tau1 <- tau_orig[1]; tau2 <- tau_orig[2];
  tau1_samp <- slice_fn1(tau1,log_cond_posterior1,learn=learn)
  eta_samp[1,m] <- (L %*% rbind(tau1_samp,tau2) + b)[1]
  var_delta <- exp(eta_samp[1,m] + eta_samp[2,(m-1)]*y_tot_IE$fine_scale)
  Qobs@x[which(y_tot$obs_name %in% c("ICESAT","ENVISAT"))] <- 1/(y_tot_IE$std^2 + var_delta)
  
  # Sample eta2
  tau1 <- tau1_samp
  tau2_samp <- slice_fn2(tau2,log_cond_posterior2,learn=learn)
  eta_samp[2,m] <- (L %*% rbind(tau1,tau2_samp) + b)[2]
  var_delta <- exp(eta_samp[1,m] + eta_samp[2,m]*y_tot_IE$fine_scale)
  Qobs@x[which(y_tot$obs_name %in% c("ICESAT","ENVISAT"))] <- 1/(y_tot_IE$std^2 + var_delta)
  
  #  # Sample alpha         
  #   Q_alpha <- as.numeric(t(x_samp[,m])%*% t(C_GRACE) %*% t(D2_big - D1_big) %*% Qobs[intGRACE,intGRACE] %*%  (D2_big - D1_big) %*% C_GRACE %*% x_samp[,m] + Qalpha_p)
  #   mu_alpha <- as.numeric(1/Q_alpha * (t(x_samp[,m]) %*% t(C_GRACE) %*% t(D2_big - D1_big) %*% Qobs[intGRACE,intGRACE] %*% (y_tot$z[intGRACE] - C_GRACE %*% x_samp[,m]) +  Qalpha_p*alpha_p))
  #   alpha_samp[m] <- mu_alpha + sqrt(1/Q_alpha)*rnorm(1)
  #   cat(alpha_samp[m],sep="\n")
  #   
  #   # Update C
  #   GRACE <- setalpha(GRACE,alpha_samp[m],450)
  #   L1 <- link(ICE,GRACE,n_grid = 400,mul_factor=consts$rho_ICE_MT_per_mkm2,mask = "in_land")
  #   L2 <- link(SURF_FIRN,GRACE,n_grid = 40,mul_factor=c(1,0),mulfun = rho_cont_land) # mulfun acts like a mask in this case
  #   L3 <- link(GIA,GRACE,n_grid = 400,mul_factor=consts$rho_ROCK_MT_per_mkm2)
  #   e <- new("link_list",list(L1,L2,L3))
  #   v <- new("block_list",list(G1=GIA,G2=SURF_FIRN,G3=ICE,O1=GRACE))
  #   G <- new("Graph",e=e,v=v)
  #   G_reduced <- compress(G)
  #   #C_full[intGRACE,] <- G_reduced@e[[1]]@Cmat
  #   C_new_elements <- G_reduced@e[[1]]@Cmat@x
  #   #C_full@x[GRACE_in_mat_positions] <- C_new_elements[!(C_new_elements == 0)]
  #   C_full@x[GRACE_in_mat_positions] <- C_new_elements
  
  
  
  ## Slice sampler
  log_cond_posterior3 <- function(alpha) {
    if(alpha <= 0 | alpha >= 0.14) {
      -Inf
    } else {
      P <- Find_Smooth_mat(df=subset(GRACE@df,t==1),alpha0=alpha,av_dist=450)
      P_big <- Reduce("bdiag",lapply(t_axis[-1],function(x) P))
      C_GRACE_new <- P_big %*% C_GRACE
      residuals <- y_GRACE - C_GRACE_new %*% x_samp[,m]
      as.vector(-0.5*sum(residuals^2 * q_obs_GRACE) )
      #+(hyp_alpha_shape1 - 1)*log(alpha) + (hyp_alpha_shape2 -1)*log(1 - alpha))
      #- 0.5*(alpha-alpha_p)^2*Qalpha_p)
    }
  }
  y_GRACE <- subset(y_tot,obs_name=="GRACE")$z
  q_obs_GRACE <- Qobs@x[intGRACE]
  alpha_samp[m] <- slice_fn3(alpha_samp[m-1],log_cond_posterior3,learn=learn)
  P <- Find_Smooth_mat(df=subset(GRACE@df,t==1),alpha0=alpha_samp[m],av_dist=450)
  P_big <- Reduce("bdiag",lapply(t_axis[-1],function(x) P))
  C_new_elements <- (P_big %*% C_GRACE)@x
  C_full@x[GRACE_in_mat_positions] <- C_new_elements
  cat(paste0("alpha = ", alpha_samp[m],"\n"))
  
  # Sample gammaI
  #residuals <- x_samp[intICE,m] - Big_B %*% x_samp[intBETA,m]
  #gammaI_beta <- as.numeric(gammaI_bp + 0.5*t(residuals) %*% Q_ICE_no_gamma %*% residuals)
  #gammaI_alpha <- as.numeric(gammaI_ap + length(intICE)/2)
  #gammaI_samp[m] <- rgamma(shape=gammaI_alpha,rate=gammaI_beta,n=1) 
  #cat(paste("gammaI_samp ",gammaI_samp[m]),sep="\n")
  #Q_ICE_with_gamma <- function(k) {  
  #  return(gammaI_samp[m]*sparsediag(prec_vs_vel))
  #}
  #ICE_VAR <- new("VAR_Gauss", mu = mu_ICE,A=A_ICE, B=B_ICE, Qw = Q_ICE_with_gamma, Qb = bdiag(Q_beta0,Q_beta1),t_axis = t_axis,name="ICE")
  #ICE <- GMRF_basis(G = ICE_VAR, Basis = Meshes$ICE)
  ##Q_full[intICEALL,intICEALL] <- ICE_VAR@Q
  #Q_new_elements <- ICE_VAR@Q@x
  #Q_full@x[ICE_in_mat_positions] <- Q_new_elements[!(Q_new_elements == 0)]
  
  
  # Sample gammaG
  attributes(Qobs)$x[intGRACE] <- 1/(((GRACE["std"])^2))
  gamma_beta <- as.numeric(gammaG_bp + 0.5*t(y_tot$z[intGRACE] - C_full[intGRACE,]%*%x_samp[,m])%*%Qobs[intGRACE,intGRACE]%*%(y_tot$z[intGRACE] - C_full[intGRACE,]%*%x_samp[,m]))
  gamma_alpha <- as.numeric(gammaG_ap + nrow(GRACE)/2)
  gammaG_samp[m] <- actuar::rinvgamma(shape=gamma_alpha,scale=gamma_beta,n=1)
  cat(paste0("gamma = ", gammaG_samp[m],"\n"))
  # Update Qobs
  attributes(Qobs)$x[intGRACE] <- 1/(gammaG_samp[m]*(GRACE["std"])^2)   
  
  #x_des$mean <- apply(x_samp[,burnin:(m-1)],1,mean)
  #x_des$mean <- x_des$current_samp
  #x_des$mean <- x_samp[,20]
  #Meshes$ICE["plot"] <- subset(x_des,name=="ICE1")$mean* Meshes$ICE["in_land"]
  #Meshes$ICE["plot"] <- subset(x_des,name=="ICE" & t == 5)$mean * Meshes$ICE["in_land"]
  #Meshes$SURF["plot"] <- head(subset(x_des,name=="SURF_FIRN" & t==3)$mean,1984)
  #Meshes$GIA["plot"] <- subset(x_des,name=="GIA")$mean
  #PlotAntarctica2(plot(Meshes$ICE,"plot",min=-0.5,max=0.5),shapefiles)
  #plot(Meshes$SURF,"plot",min=-0.5,max=0.5)
  #PlotAntarctica2(plot(Meshes$GIA,"plot",min=-0.01,max=0.01),shapefiles)
  
  
  # Dithering
  if (m > mm + dither) {
    x_samp[,mm] =x_samp[,m]
    mm = mm + 1
    m = mm
  }
  print(m)
  if(m %% 100 == 0)  save(gammaG_samp,alpha_samp,eta_samp,x_samp,x_des,burnin,m,Meshes,file="~/cache/MCMCresults.rda")
  if(m %% 2000 == 0) {
    save(gammaG_samp,alpha_samp,eta_samp,x_samp,x_des,burnin,file=paste0("~/cache/MCMCResultsRound",count,".rda"))      
    gammaG_samp[1] <- gammaG_samp[m]
    alpha_samp[1] <- alpha_samp[m]
    eta_samp[,1] <- eta_samp[,m]
    x_samp[,1] <- x_samp[,m]
    m <- 1            
    count <- count+ 1
  }
}