library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rm(list = ls())
# question 1 --------------------------------------------------------------

R <- 10
n <- 500

b <- 0.3
vec_lambda <- c(1, 0.8, 0.8, 1, 0.7, 0.9, 0.7, 1, 0.8)
vec_gamma <- c(0.2, 0.5)
mat_phi <- matrix(data = c(1, 0.2, 0.2, 0.81), ncol = 2)
vec_psi_eps <- c(1, 1, 1, 0.3, 0.3, 0.4, 0.4)
psi_delta <- 0.36


# data generation ---------------------------------------------------------

set.seed(128)
mat_xi <- matrix(data = mvtnorm::rmvnorm(n, mean = rep(0, 2), sigma = mat_phi), 
                 ncol = 2)
vec_d <- rbinom(n, size = 1, prob = .7)

vec_eta <- b * vec_d + mat_xi %*% vec_gamma + 
  rnorm(n, mean = 0, sd = sqrt(psi_delta))

# mat_y <- matrix(data = NA, nrow = n, ncol = 9)

mat_ys <- matrix(rep(vec_lambda[1:3], each = n), ncol = 3) *
  matrix(rep(vec_eta, 3), ncol = 3) + 
  sapply(vec_psi_eps[1:3], 
         function (var) rnorm(n, mean = 0, sd = sqrt(var)))

mat_y1 <- sign(sign(mat_ys) + 1)

mat_y2 <- matrix(rep(vec_lambda[4:7], each = n), ncol = 4) *
  matrix(rep(mat_xi[, 1], 4), ncol = 4) + 
  sapply(vec_psi_eps[4:7], 
         function (var) rnorm(n, mean = 0, sd = sqrt(var)))

mat_vart <- matrix(rep(vec_lambda[8:9], each = n), ncol = 2) *
  matrix(rep(mat_xi[, 2], 2), ncol = 2)

mat_y3 <- mat_vart %>% faraway::ilogit() %>% 
  apply(c(1, 2), function(x) rbinom(1, 1, x))

ls_dat1 <- list(N = n,
                # y1 = mat_y1[, 1],
                # y2 = mat_y1[, 2],
                # y3 = mat_y1[, 3],
                n1 = sum(mat_y1[, 1] == 1),
                n2 = sum(mat_y1[, 2] == 1),
                n3 = sum(mat_y1[, 3] == 1),
                sub1 = which(mat_y1[, 1] == 1),
                sub1c = which(mat_y1[, 1] == 0),
                sub2 = which(mat_y1[, 2] == 1),
                sub2c = which(mat_y1[, 2] == 0),
                sub3 = which(mat_y1[, 3] == 1),
                sub3c = which(mat_y1[, 3] == 0),
                y4 = mat_y2[, 1], 
                y5 = mat_y2[, 2],
                y6 = mat_y2[, 3],
                y7 = mat_y2[, 4], 
                y8 = mat_y3[, 1],
                y9 = mat_y3[, 2],
                d = vec_d,
                lw_bd = c(-Inf, 0),
                up_bd = c(0, Inf)) 

fit_m1 <- stan(
  file = "Asg2/stan_file/q1.stan",
  data = ls_dat1,                   # named list of data
  chains = 4,                       # number of Markov chains
  warmup = 1000,                    # number of warmup iterations per chain
  iter = 2000                       # total number of iterations per chain
)

# save.image("Asg2/output/q1.RData")
# load("Asg2/output/q1.RData")

print(
  fit_m1, 
  pars = c('lp__', 'gam', 'b',
         'lam2', 'lam3', 'lam5', 'lam6', 'lam7', 'lam9', 
         'psi_e', 'psi_d', 'mat_phi'), 
  probs=c(.05, .95)
)

traceplot(
  fit_m1, 
  pars = c('lp__', 'gam', 'b',
           'lam2', 'lam3', 'lam5', 'lam6', 'lam7', 'lam9', 
           'psi_e', 'psi_d', 'mat_phi')
)

rstan::extract(
  fit_m1, 
  pars = c('gam', 'b',
           'lam2', 'lam3', 'lam5', 'lam6', 'lam7', 'lam9', 
           'psi_e', 'psi_d', 'mat_phi')
) %>% as.data.frame() %>% sapply(mean)




rm(list = ls())

## set parameters

R <- 10
n <- 500

b <- 0.3
vec_lambda <- c(1, 0.8, 0.8, 1, 0.7, 0.9, 0.7, 1, 0.8)
vec_gamma <- c(0.2, 0.5)
mat_phi <- matrix(data = c(1, 0.2, 0.2, 0.81), ncol = 2)
vec_psi_eps <- c(1, 1, 1, 0.3, 0.3, 0.4, 0.4)
psi_delta <- 0.36

fun_sim <- function() {
  mat_xi <- matrix(data = mvtnorm::rmvnorm(n, mean = rep(0, 2), sigma = mat_phi), 
                   ncol = 2)
  vec_d <- rbinom(n, size = 1, prob = .7)
  
  vec_eta <- b * vec_d + mat_xi %*% vec_gamma + 
    rnorm(n, mean = 0, sd = sqrt(psi_delta))
  
  mat_ys <- matrix(rep(vec_lambda[1:3], each = n), ncol = 3) *
    matrix(rep(vec_eta, 3), ncol = 3) + 
    sapply(vec_psi_eps[1:3], 
           function (var) rnorm(n, mean = 0, sd = sqrt(var)))
  
  mat_y1 <- sign(sign(mat_ys) + 1)
  
  mat_y2 <- matrix(rep(vec_lambda[4:7], each = n), ncol = 4) *
    matrix(rep(mat_xi[, 1], 4), ncol = 4) + 
    sapply(vec_psi_eps[4:7], 
           function (var) rnorm(n, mean = 0, sd = sqrt(var)))
  
  mat_vart <- matrix(rep(vec_lambda[8:9], each = n), ncol = 2) *
    matrix(rep(mat_xi[, 2], 2), ncol = 2)
  
  mat_y3 <- mat_vart %>% faraway::ilogit() %>% 
    apply(c(1, 2), function(x) rbinom(1, 1, x))
  
  ls_dat1 <- list(N = n,
                  n1 = sum(mat_y1[, 1] == 1),
                  n2 = sum(mat_y1[, 2] == 1),
                  n3 = sum(mat_y1[, 3] == 1),
                  sub1 = which(mat_y1[, 1] == 1),
                  sub1c = which(mat_y1[, 1] == 0),
                  sub2 = which(mat_y1[, 2] == 1),
                  sub2c = which(mat_y1[, 2] == 0),
                  sub3 = which(mat_y1[, 3] == 1),
                  sub3c = which(mat_y1[, 3] == 0),
                  y4 = mat_y2[, 1], 
                  y5 = mat_y2[, 2],
                  y6 = mat_y2[, 3],
                  y7 = mat_y2[, 4], 
                  y8 = mat_y3[, 1],
                  y9 = mat_y3[, 2],
                  d = vec_d,
                  lw_bd = c(-Inf, 0),
                  up_bd = c(0, Inf)) 
  
  fit_m1 <- stan(
    file = "Asg2/stan_file/q1.stan",
    data = ls_dat1,
    chains = 4,
    warmup = 1000,
    iter = 2000
  )
  rstan::extract(
    fit_m1, 
    pars = c('gam', 'b', 'mu',
             'lam2', 'lam3', 'lam5', 'lam6', 'lam7', 'lam9', 
             'psi_e', 'psi_d', 'mat_phi')
  ) %>% as.data.frame() %>% sapply(mean)
}

set.seed(128)
mat_result <- replicate(10, fun_sim())

mat_result

save(mat_result, file = "Asg2/output/result.RData")
load(file = "Asg2/output/result.RData")

vec_par <- c(vec_gamma, b, rep(0, 9), vec_lambda[c(-1, -4, -8)], 
             vec_psi_eps[-3:-1], psi_delta, as.vector(mat_phi))

mat_bias_mse <- 
  (mat_result - vec_par) %>% 
  apply(1, function (vec) {
    c(BIAS = mean(vec), 
      MSE = sqrt(mean(vec^2)))
  })

mat_bias_mse %>% round(digits = 5) %>% write.csv("Asg2/output/q1_BM.csv")







# test --------------------------------------------------------------------

stan(
  data = list(p = 2, 
              a = c(1, 2)),
  file = "Asg2/stan_file/test.stan",
  chains = 1,
  warmup = 10,
  iter = 20
)




