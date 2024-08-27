library(R2WinBUGS)
library(boa)

# set experiment date
timestamp = strftime(Sys.time(), "%Y%m%d-%H")
winBUGS.path = "D:/pkgs/WinBUGS14/"
datapath = paste0(getwd(),'/data')
print(datapath)


# data size
iter = 10
N = 500
P = 10
Nlam = 7
Neta = 1 
Nxi = 2
Ngam = 2

# containers for Bayesian estimates and standard errors
uby.E <- matrix(data = NA, nrow = iter, ncol = P)
uby.SE <- matrix(data = NA, nrow = iter, ncol = P)
lam.E <- matrix(data = NA, nrow = iter, ncol = Nlam)
lam.SE <- matrix(data = NA, nrow = iter, ncol = Nlam)
gam.E <- matrix(data = NA, nrow = iter, ncol = Ngam)
gam.SE <- matrix(data = NA, nrow = iter, ncol = Ngam)
phx.E <- matrix(data = NA, nrow = iter, ncol = Nxi^2)
phx.SE <- matrix(data = NA, nrow = iter, ncol = Nxi^2)
ubeta.E <- numeric(iter)
ubeta.SE <- numeric(iter)
# only continuous part of y's variances need to estimate
sgm.E <- matrix(data = NA, nrow = iter, ncol = 4) 
sgm.SE <- matrix(data = NA, nrow = iter, ncol = 4)
sgd.E <- numeric(iter)
sgd.SE <- numeric(iter)

# containers for HPD (Highest Probability Density) intervals
uby.hpd <- array(NA, c(iter, P, 2))
ubeta.hpd <- array(NA, c(iter, 2))
lam.hpd <- array(NA, c(iter, Nlam, 2))
gam.hpd <- array(NA, c(iter, Ngam, 2))
phx.hpd <- array(NA, c(iter, Nxi^2, 2))
sgm.hpd <- array(NA, c(iter, 4, 2))
sgd.hpd <- array(NA, c(iter, 2))

# container for DIC values
DIC = numeric(iter)

# parameters to be estimated
parameters <- c("uby", "ubeta", "lam", "gam", "phx", "sgm", "sgd")

# set important prior param R_0
R0 = matrix(c(7.0, 2.1, 2.1, 7.0), nrow = 2)

# initial values for MCMC in WinBUGS
init1 <- list(
  uby = rep(0.5, P),
  ubeta = 0.5,
  lam = rep(0.5, Nlam),
  gam = rep(0.5, Ngam),
  phi = matrix(c(1, 0.5, 0.5, 1), nrow = Nxi),
  psi = rep(1, 4),
  psd = 1,
  xi = matrix(data = rep(0.3, N * Nxi), ncol = Nxi)
)

init2 <- list(
  uby = rep(0, P),
  ubeta = 0,
  lam = rep(0, Nlam),
  gam = rep(0, Ngam),
  phi = matrix(c(2, 0, 0, 2), nrow = Nxi),
  psi = rep(2, 4),
  psd = 2,
  xi = matrix(data = rep(0.3, N * Nxi), ncol = Nxi)
)

init3 <- list(
  uby = rep(0, P),
  ubeta = 0.3,
  lam = rep(0.8, Nlam),
  gam = rep(0.45, Ngam),
  phi = matrix(c(1, 0.2, 0.2, 1), nrow = Nxi),
  psi = rep(0.3, 4),
  psd = 0.4,
  xi = matrix(data = rep(0.3, N * Nxi), ncol = Nxi)
)

init4 <- list(
  uby = rep(0, P),
  ubeta = 0.4,
  lam = rep(0.9, Nlam),
  gam = rep(0.5, Ngam),
  phi = matrix(c(1, 0.3, 0.3, 1), nrow = Nxi),
  psi = rep(0.4, 4),
  psd = 0.3,
  xi = matrix(data = rep(0.3, N * Nxi), ncol = Nxi)
)

inits <- list(init1, init3)


# Do simulations based on 10 replications
for (t in 1:iter) {
  iterpath = paste0(getwd(),"/Q1.2_", t)
  dir.create(iterpath, showWarnings = FALSE, recursive = TRUE)
  
  Y <- as.matrix(read.table(paste(datapath, "/Y-", t, ".txt", sep = "")))
  D <- read.table(paste(datapath, "/D-", t, ".txt", sep = ""))$x

  data <- list(N = N, P = P, R = R0,
              d = D, z = Y,
              low = c(-2000, 0),
              high = c(0, 2000))
  
  model = bugs(data, inits, parameters, 
               model.file = paste0(getwd(),"/../model1.txt"), 
               n.chains = 2, 
               n.iter = 25000, 
               n.burnin = 15000, 
               n.thin = 1, 
               DIC = TRUE,
               bugs.directory = winBUGS.path,
               working.directory = iterpath, 
               debug = FALSE)
  
  # save estimates and standard errors
  uby.E[t, ] = model$mean$uby
  uby.SE[t, ] = model$sd$uby
  ubeta.E[t] = model$mean$ubeta
  ubeta.SE[t] = model$sd$ubeta
  lam.E[t, ] = model$mean$lam
  lam.SE[t, ] = model$sd$lam
  gam.E[t, ] = model$mean$gam
  gam.SE[t, ] = model$sd$gam
  phx.E[t, ] = c(model$mean$phx)
  phx.SE[t, ] = c(model$sd$phx)
  sgm.E[t, ] = model$mean$sgm
  sgm.SE[t, ] = model$sd$sgm
  sgd.E[t] = model$mean$sgd
  sgd.SE[t] = model$sd$sgd
  
  # save HPD intervals
  for (k in 1:P) {
    temp = model$sims.array[ , 1, k]
    uby.hpd[t, k, ] = boa.hpd(temp, 0.05)
  }
  temp = model$sims.array[ , 1, P + 1]
  ubeta.hpd[t, ] = boa.hpd(temp, 0.05)
  for (k in 1:Nlam) {
    temp = model$sims.array[ , 1, P + 1 + k]
    lam.hpd[t, k, ] = boa.hpd(temp, 0.05)
  }
  for (k in 1:Ngam) {
    temp = model$sims.array[ , 1, P + 1 + Nlam + k]
    gam.hpd[t, k, ] = boa.hpd(temp, 0.05)
  }
  for (k in 1:Nxi^2) {
    temp = model$sims.array[ , 1, P + 1 + Nlam + Ngam + k]
    phx.hpd[t, k, ] = boa.hpd(temp, 0.05)
  }
  for (k in 1:4) {
    temp = model$sims.array[ , 1, P + 1 + Nlam + Ngam + Nxi^2 + k]
    sgm.hpd[t, k, ] = boa.hpd(temp, 0.05)
  }
  temp = model$sims.array[ , 1, P + 1 + Nlam + Ngam + Nxi^2 + 4 + 1]
  sgd.hpd[t, ] = boa.hpd(temp, 0.05)
  
  # save DIC values
  DIC[t] = model$DIC
  
  print(model$summary)
}

metr_params = list(
  E = list(
    uby = uby.E, ubeta = ubeta.E, lam = lam.E, gam = gam.E, 
    phx = phx.E, sgm = sgm.E, sgd = sgd.E
  ),
  SE = list(
    uby = uby.SE, ubeta = ubeta.SE, lam = lam.SE, gam = gam.SE, 
    phx = phx.SE, sgm = sgm.SE, sgd = sgd.SE
  ),
  HPD = list(
    uby = uby.hpd, ubeta = ubeta.hpd, lam = lam.hpd, gam = gam.hpd, 
    phx = phx.hpd, sgm = sgm.hpd, sgd = sgd.hpd
  )
)

save(metr_params, file = paste0(getwd(), "/Q1.2_metrparams_", timestamp, ".RData"))
