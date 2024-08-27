library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rm(list = ls())

# Parameters --------------------------------------------------------------

R <- 10  ## replicate
N <- 500 ## sample size
p <- 9

## measurement equation
vec_a <- rep(c(0.5, -0.5, 0.8), each = 3)
lam_21 <- lam_52 <- lam_83 <- 0.8
lam_31 <- lam_62 <- lam_93 <- 0.5
vec_psi_epsilon <- rep(c(0.3, 0.5, 0.4), each = 3)

## structural equation
b <- -1
vec_gamma <- rep(c(0.5, 0.3), each = 2)
phi <- 1  ## variance of xi
psi_delta <- .36  ## variance of delta

mat_loading <- matrix(data = 0, nrow = p, ncol = 3)
mat_loading[1:3, 1] <- c(1, lam_21, lam_31)
mat_loading[4:6, 2] <- c(1, lam_52, lam_62)
mat_loading[7:9, 3] <- c(1, lam_83, lam_93)
mat_loading

# Data Generation ---------------------------------------------------------

set.seed(128)

vec_xi1 <- rnorm(N, sd = sqrt(phi))
vec_xi2 <- rnorm(N, sd = sqrt(phi))
vec_d <- rnorm(N, sd = 1)
vec_c <- rnorm(N, sd = 1)

vec_eta <- b * vec_d + vec_gamma[1] * vec_xi1 + vec_gamma[2] * vec_xi2 + 
  vec_gamma[3] * vec_xi1 * vec_xi2 + vec_gamma[4] * vec_xi1^2 + 
  rnorm(N, sd = sqrt(psi_delta))

mat_y <- vec_a %*% t(vec_c) + 
  mat_loading %*% t(cbind(vec_eta, vec_xi1, vec_xi2)) +
  t(sapply(vec_psi_epsilon, 
           function(var) rnorm(n = N, sd = sqrt(var))))


# Model Fitting -----------------------------------------------------------

dat_mod1 <- list(N = N, p = p,
                 Y = mat_y, c = vec_c, d = vec_d)

mod_NL1 <- stan(file = "Asg1/stan_file/NonLinearSEM1.stan",
                data = dat_mod1)

rstan::extract(mod_NL1, pars = c('mu', 'a', 'b', 'gam', 
                                 'lam_21', 'lam_52', 'lam_83', 
                                 'lam_31', 'lam_62', 'lam_93', 
                                 'psi_eps', 'phi', 'psi_delta')) %>% 
  as.data.frame() %>% 
  sapply(mean)

rstan::extract(mod_NL1, pars = c('mu', 'a', 'b', 'gam', 
                                 'lam_21', 'lam_52', 'lam_83', 
                                 'lam_31', 'lam_62', 'lam_93', 
                                 'psi_eps', 'phi', 'psi_delta')) %>% 
  as.data.frame() %>% 
  sapply(function(x) c("mean" = mean(x), 
                       # "sd" = sd(x),
                       quantile(x, probs = .025),
                       quantile(x, probs = .975))) %>% 
  as.data.frame() %>% rownames_to_column(var = "value") %>% 
  tibble() %>% print(width = Inf)

rstan::extract(mod_NL1, pars = c('gam', 'b', 'lam_21', 'lam_52', 'lam_83', 
                                 'lam_31', 'lam_62', 'lam_93', 'phi', 'lp__')) %>% 
  as.data.frame() %>% tibble() %>% 
  mutate(iter = 1:n()) %>% 
  gather(-iter, key = "par", value = "value") %>% 
  mutate(par = factor(par, levels = unique(par))) %>% 
  ggplot() +
  geom_point(mapping = aes(x = iter, y = value, color = par),
             show.legend = FALSE) +
  facet_wrap(par ~., scales = "free_y", nrow = 3) + 
  labs(title = "Convergence Scenario") +
  theme(plot.title = element_text(hjust = .5))

pdf("Asg1/output/conv.pdf", width = 18, height = 12)
dev.off()

# function for simulation -------------------------------------------------

fun_sim1 <- function() {
  vec_xi1 <- rnorm(N, sd = sqrt(phi))
  vec_xi2 <- rnorm(N, sd = sqrt(phi))
  vec_d <- rnorm(N, sd = 1)
  vec_c <- rnorm(N, sd = 1)
  
  vec_eta <- b * vec_d + vec_gamma[1] * vec_xi1 + vec_gamma[2] * vec_xi2 + 
    vec_gamma[3] * vec_xi1 * vec_xi2 + vec_gamma[4] * vec_xi1^2 + 
    rnorm(N, sd = sqrt(psi_delta))
  
  mat_y <- vec_a %*% t(vec_c) + 
    mat_loading %*% t(cbind(vec_eta, vec_xi1, vec_xi2)) +
    t(sapply(vec_psi_epsilon, 
             function(var) rnorm(n = N, sd = sqrt(var))))
  
  dat_mod1 <- list(N = N, p = p,
                   Y = mat_y, c = vec_c, d = vec_d)
  
  mod_NL1 <- stan(file = "Asg1/stan_file/NonLinearSEM1.stan",
                  data = dat_mod1)
  
  rstan::extract(mod_NL1, pars = c('mu', 'a', 'b', 'gam', 
                                   'lam_21', 'lam_52', 'lam_83', 
                                   'lam_31', 'lam_62', 'lam_93', 
                                   'psi_eps', 'phi', 'psi_delta')) %>% 
    as.data.frame() %>% 
    sapply(mean)
}

dat_sim1 <- replicate(R, fun_sim1())

vec_theta0 <- c(rep(0, 9), vec_a, b,  ## for mu, a, b
                vec_gamma, 
                lam_21, lam_52, lam_83,
                lam_31, lam_62, lam_93, 
                vec_psi_epsilon, 
                phi, phi, psi_delta)

apply(dat_sim1 - vec_theta0, 1, 
      function (x) c("BIAS" = mean(x),
                     "RMS" = sqrt(mean(x^2)))) %>% 
  round(digits = 5) %>% write.csv("Asg1/output/q2_bias_rms.csv")

vec_bias1 <- vec_mean1 - vec_theta0

vec_mse <- vec_bias1


# 2b ----------------------------------------------------------------------

dat_mod1 <- list(N = N, p = p,
                 Y = mat_y, c = vec_c, d = vec_d)

mod_NL1 <- stan(file = "Asg1/stan_file/NonLinearSEM1.stan",
                data = dat_mod1)

mod_NL1s <- stan(file = "Asg1/stan_file/NonLinearSEM1_prior.stan",
                data = dat_mod1)

rbind(rstan::extract(mod_NL1, pars = c('mu', 'a', 'b', 'gam', 
                                       'lam_21', 'lam_52', 'lam_83', 
                                       'lam_31', 'lam_62', 'lam_93', 
                                       'psi_eps', 'phi', 'psi_delta')) %>% 
        as.data.frame() %>% sapply(mean),
      rstan::extract(mod_NL1s, pars = c('mu', 'a', 'b', 'gam', 
                                        'lam_21', 'lam_52', 'lam_83', 
                                        'lam_31', 'lam_62', 'lam_93', 
                                        'psi_eps', 'phi', 'psi_delta')) %>% 
        as.data.frame() %>% sapply(mean)) %>% 
  round(digits = 4) %>% write.csv("Asg1/output/q2_prior.csv")


# Bayes Factor ------------------------------------------------------------

dat_mod2 <- list(N = N, p = p,
                 Y = mat_y, c = vec_c, d = vec_d, t = 1)

mod_NL2 <- stan(file = "Asg1/stan_file/NonLinearSEM2.stan",
     data = dat_mod2, chains = 1)

rstan::extract(mod_NL2, pars = c("U")) %>% 
  as.data.frame() %>% sapply(mean)

rstan::extract(mod_NL2, pars = c('gam', 'lp__')) %>% 
  as.data.frame() %>% tibble() %>% 
  mutate(iter = 1:n()) %>% 
  gather(-iter, key = "par", value = "value") %>% 
  mutate(par = factor(par, levels = unique(par))) %>% 
  ggplot() +
  geom_point(mapping = aes(x = iter, y = value, color = par)) +
  facet_wrap(par ~., scales = "free_y", nrow = 2) + 
  labs(title = "Convergence Scenario of Model 1") +
  theme(plot.title = element_text(hjust = .5))

vec_t <- seq(from = 0, to = 1, by = 0.1) 
s <- length(vec_t)
vec_u <- rep(0, t)
for (k in 1:s) {
  t <- vec_t[k]
  dat_mod2 <- list(N = N, p = p,
                   Y = mat_y, c = vec_c, d = vec_d, t = t)
  mod_NL2 <- stan(file = "Asg1/stan_file/NonLinearSEM2.stan",
                  data = dat_mod2, chains = 1)
  vec_u[k] <- rstan::extract(mod_NL2, pars = c("U")) %>% 
    as.data.frame() %>% sapply(mean)
}

logBmean(vec_u[1:(s-1)] + vec_u[2:s])/2





# Bayes Factor2 -----------------------------------------------------------

dat_mod2 <- list(N = N, p = p,
                 Y = mat_y, c = vec_c, d = vec_d, t = 1)

mod_L <- stan(file = "Asg1/stan_file/LinearSEM.stan",
                data = dat_mod2, chains = 1)

rstan::extract(mod_L, pars = c("U")) %>% 
  as.data.frame() %>% sapply(mean)

rstan::extract(mod_L, pars = c('gam', 'lp__')) %>% 
  as.data.frame() %>% tibble() %>% 
  mutate(iter = 1:n()) %>% 
  gather(-iter, key = "par", value = "value") %>% 
  mutate(par = factor(par, levels = unique(par))) %>% 
  ggplot() +
  geom_point(mapping = aes(x = iter, y = value, color = par)) +
  facet_wrap(par ~., scales = "free_y", nrow = 2) + 
  labs(title = "Convergence Scenario of Model 1") +
  theme(plot.title = element_text(hjust = .5))

vec_t <- seq(from = 0, to = 1, by = 0.1) 
s <- length(vec_t)
vec_u <- rep(0, t)
for (k in 1:s) {
  t <- vec_t[k]
  dat_mod2 <- list(N = N, p = p,
                   Y = mat_y, c = vec_c, d = vec_d, t = t)
  mod_L <- stan(file = "Asg1/stan_file/LinearSEM.stan",
                data = dat_mod2, chains = 1)
  vec_u[k] <- rstan::extract(mod_L, pars = c("U")) %>% 
    as.data.frame() %>% sapply(mean)
}

mean(vec_u[1:(s-1)] + vec_u[2:s])/2



