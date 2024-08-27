library(mvtnorm)

set.seed(5020)

# set data repository
datapath = paste0(getwd(),'/data')
dir.create(datapath, showWarnings = FALSE, recursive = TRUE)

# data size
iter = 10
N = 500
P = 10
Neta = 1 # not used here
Nxi = 2 #
Ngam = 2 #

# set the true values of parameters
uby <- rep(0, P)
lam <- c(0.8, 0.8, 0.7, 0.9, 0.7, 0.9, 0.8)
sgm <- c(1, 1, 1, 0.3, 0.3, 0.25, 0.25)
ubeta <- 0.3
gam <- c(0.4, 0.5)
phx <- matrix(data = c(1, 0.2, 0.2, 0.81), ncol = 2)
sgd <- 0.36

# set important prior param R_0
# R0 <- matrix(c(7.0, 2.1, 2.1, 7.0), nrow = 2)


# containers for generated data
Y <- matrix(data = NA, nrow = N, ncol = P)
D <- numeric(N)
p <- numeric(P)
v <- numeric(P)


# generate data
for (t in 1:iter) {
  for (i in 1:N) {
    # BD[i] = rt(1, 5)
    # BC[i] = rt(1, 5)
    
    # generate the fixed covariates in SE (from Bernoulli(0.7))
    d <- rbinom(1, 1, 0.7)
    D[i] <- d
    
    # generate xi
    xi <- rmvnorm(1, c(0, 0), phx)
    
    # generate error term in SE
    del <- rnorm(1, 0, sqrt(sgd))
    
    # generate eta according to the SE
    eta <- ubeta * d + gam[1] * xi[1] + gam[2] * xi[2] + del
    
    # generate error term in ME
    eps <- numeric(7)
    for (k in 1:7) { eps[k] <- rnorm(1, 0, sgm[k]) }
    
    # generate theta in ME
    v[1] <- uby[1] + eta + eps[1]
    v[2] <- uby[2] + lam[1] * eta + eps[2]
    v[3] <- uby[3] + lam[2] * eta + eps[3]
    Y[i, 4] <- uby[4] + xi[1] + eps[4]
    Y[i, 5] <- uby[5] + lam[3] * xi[1] + eps[5]
    Y[i, 6] <- uby[6] + lam[4] * xi[1] + eps[6]
    Y[i, 7] <- uby[7] + lam[5] * xi[1] + eps[7]
    v[8] <- uby[8] + xi[2]
    v[9] <- uby[9] + lam[6] * xi[2]
    v[10] <- uby[10] + lam[7] * xi[2]
    
    # transform theta to ordinal variables
    for (j in 1:3) {
      if (v[j] > 0) Y[i, j] <- 1
      else Y[i, j] <- 0
    }

    # transform theta to binary variables
    for (j in 8:10) {
      p[j] <- exp(v[j]) / (1 + exp(v[j]))
      Y[i, j] <- rbinom(1, 1, p[j])
    }
  }
  
  # save data matrix
  write.table(Y, paste(datapath, "/Y-", t, ".txt", sep = ""))
  write.table(D, paste(datapath, "/D-", t, ".txt", sep = ""))
}

true_params = list(
  lam = lam,
  uby = uby,
  sgm = sgm,
  ubeta = ubeta,
  gam = gam,
  phx = phx,
  sgd = sgd
)

save(true_params, file = paste0(datapath, "/trueparams.RData"))

