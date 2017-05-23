##### Functions for extending MVST and INLA


## Prior setup
## Functions for transforming priors from original space to log-normal space
Tlognorm <- function(mu, v){
  logv <- log(1 + v/mu^2)
  logmu <- log(mu^2) - 0.5*log(mu^2 + v)
  return(c(logmu, logv))
}

## Functions for transforming priors from log-normal space back to original space
lognormT <- function(logmu, logv){
  mu <- exp(logmu + logv/2)
  v <- (exp(logv)-1)*exp(2*logmu + logv)
  return(c(mu, v))
}


## for printing plots
my.dev.print <- function(stem) {
  stem <- gsub("\\..+$", "", stem)
  fname <- paste0("inla_", stem, ".png") 
  dev.print(png, fname,
            width = par("din")[1], height = par("din")[2],
            units = "in", pointsize = 12, res = 720,
            bg = "transparent")
  fname <- paste0("jgi2_", stem, ".pdf") 
  dev.print(pdf, fname)
}
