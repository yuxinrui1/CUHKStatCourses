# Simulation for a log linear model with count outcome
# Set the number of replication, sample size, and true model parameters
nrep = 10
n <- 300
beta <- c(0.5, 1, -1)
Ebeta = matrix(NA, nrow = nrep, ncol = 3)

for (irep in 1:nrep) {
  #generate data
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x <- rbind(rep(1,n), x1, x2)
  mu <- exp(beta %*% x)
  y = rpois(n, mu)
  #fit the model
  glm_fit = glm(y ~ x1 + x2, family = poisson)
  # store the estimates
  Ebeta[irep,] = glm_fit$coefficients
}

# Evaluate model performance through average bias and root mean squre error
Avg_bias = colMeans(Ebeta) - beta
RMSE = sqrt(rowMeans(apply(Ebeta, 1, function(x){return ((x - beta)^2)})))