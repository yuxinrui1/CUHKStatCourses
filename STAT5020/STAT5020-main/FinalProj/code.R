
# loading -----------------------------------------------------------------
library(tidyverse)

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rm(list = ls())
# data preliminary --------------------------------------------------------

data.raw <- read_csv("data/ObesityDataSet_raw_and_data_sinthetic.csv")

data.raw %>% summary()

data_ <- data.raw %>%
  rename(FHO = family_history_with_overweight) %>% 
  mutate(MBI = Weight / (Height * Height),
         FHO = factor(FHO, levels = c("no", "yes"), labels = 1:2) %>% 
           as.integer()) %>% 
  select(FHO, MBI, Weight, NCP, CH2O, FCVC, FAF, TUE)

N <- nrow(data_)

# model fitting -----------------------------------------------------------

data.norm <- data_ %>% 
  rename(y1 = MBI, 
         y2 = Weight,
         y3 = NCP,
         y4 = CH2O,
         y5 = FCVC,
         y6 = FAF,
         y7 = TUE,
         g = FHO)

data.norm
data.norm %>% summary()


## model 0 one-sample with non-info prior -----------------------------------

mod0e <- stan(
  file = "stan_file/model0e.stan",
  data = list(
    N = N,
    y1 = data.norm$y1,
    y2 = data.norm$y2,
    y3 = data.norm$y3,
    y4 = data.norm$y4,
    y5 = data.norm$y5,
    y6 = data.norm$y6,
    y7 = data.norm$y7,
    mu0 = sapply(data.norm[paste0("y", 1:7)], mean) %>% as.vector()
  )
)

# save(mod0e, file = "output/mod0e.RData")
load(file = "output/mod0e.RData")

traceplot(
  mod0e,
  pars = c(
    'lp__', 
    'mu', 
    'lam',
    'gam')
)

traceplot(
  mod0e,
  pars = c(
    'psi_e',
    'psi_d',
    'mat_phi'
  )
)

print(
  mod0e,
  probs = c(.05, .5, .95),
  pars = c(
    'mu',
    'lam',
    'gam',
    'psi_e',
    'psi_d',
    'mat_phi'
  )
)


## model 0 -----------------------------------------------------------------

mod0 <- stan(
  file = "stan_file/model0.stan",
  data = list(
    N = N,
    y1 = data.norm$y1,
    y2 = data.norm$y2,
    y3 = data.norm$y3,
    y4 = data.norm$y4,
    y5 = data.norm$y5,
    y6 = data.norm$y6,
    y7 = data.norm$y7,
    mu0 = sapply(data.norm[paste0("y", 1:7)], mean) %>% as.vector()
  )
)

save(mod0, file = "output/mod0.RData")
load(file = "output/mod0.RData")

traceplot(
  mod0,
  pars = c(
    'lp__', 
    'mu', 
    'lam',
    'gam')
)

traceplot(
  mod0,
  pars = c(
    'psi_e',
    'psi_d',
    'mat_phi'
  )
)

print(
  mod0,
  probs = c(.05, .5, .95),
  pars = c(
    'mu',
    'lam',
    'gam',
    'psi_e',
    'psi_d',
    'mat_phi'
  )
)

pdf(file = "output/mod0.pdf",
    width = 10,
    height = 10)

traceplot(
  mod0,
  pars = c(
    'lp__', 
    'mu', 
    'lam',
    'gam',
    'psi_e',
    'psi_d',
    'mat_phi'
  ), 
  nrow = 7
)

dev.off()

est_mod0 <- rstan::extract(
  mod0,
  pars = c(
    'mu',
    'lam',
    'gam',
    'psi_e',
    'psi_d',
    'mat_phi'
  )
) %>% as.data.frame() %>% sapply(function (vec) {
  c(mean = mean(vec),
    sd = sd(vec))
}) 

est_mod0 %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "par") %>% # tibble() %>% 
  format(digits = 2) %>% 
  write_csv("output/mod0.csv")

## model 1 multi-sample ----------------------------------------------------

mod1 <- stan(
  file = "stan_file/model1.stan",
  data = list(
    N = N,
    g = data.norm$g,
    y1 = data.norm$y1,
    y2 = data.norm$y2,
    y3 = data.norm$y3,
    y4 = data.norm$y4,
    y5 = data.norm$y5,
    y6 = data.norm$y6,
    y7 = data.norm$y7,
    mu0 = sapply(data.norm[paste0("y", 1:7)], mean) %>% as.vector()
  )
)

# save(mod1, file = "output/mod1.RData")
load(file = "output/mod1.RData")

traceplot(
  mod1,
  pars = c(
    'lp__', 
    'mu', 
    'lam',
    'gam')
)

traceplot(
  mod1,
  pars = c(
    'psi_e',
    'psi_d',
    'mat_phi'
  )
)

print(
  mod1,
  probs = c(.05, .5, .95),
  pars = c(
    'mu',
    'lam',
    'gam',
    'psi_e',
    'psi_d',
    'mat_phi'
  )
)

pdf(file = "output/mod1.pdf",
    width = 15,
    height = 18)

traceplot(
  mod1,
  pars = c(
    'lp__',
    'mu',
    'lam',
    'gam',
    'psi_e',
    'psi_d',
    'mat_phi'
  ),
  nrow = 9
)

dev.off()


est_mod1 <- rstan::extract(
  mod1,
  pars = c(
    'mu',
    'lam',
    'gam',
    'psi_e',
    'psi_d',
    'mat_phi'
  )
) %>% as.data.frame() %>% sapply(mean)

tibble(par = names(est_mod1),
       value = est_mod1) %>% 
  mutate(gp = str_extract(par, "\\d$"),
         gp = paste0("group", gp),
         par = str_remove(par, ".\\d$")) %>% 
  spread(key = gp, value = value) %>% 
  as.data.frame() %>% 
  format(digits = 2) %>% 
  write_csv("output/mod1.csv")


## model 1a BF --------------------------------------------------------------

vec_t <- seq(from = 0, to = 1, by = 0.05)
J <- length(vec_t)


modt1a <- stan(
  file = "stan_file/modelt1a.stan",
  data = list(
    N = N,
    g = data.norm$g,
    y1 = data.norm$y1,
    y2 = data.norm$y2,
    y3 = data.norm$y3,
    y4 = data.norm$y4,
    y5 = data.norm$y5,
    y6 = data.norm$y6,
    y7 = data.norm$y7,
    mu0 = sapply(data.norm[paste0("y", 1:7)], mean) %>% as.vector(),
    t = 0
  ),
  chains = 1,
  warmup = 500,
  iter = 1500
)

rstan::extract(modt1a, pars = c("U")) %>% as.data.frame() %>% sapply(mean)
## 66.71792 (*2)


rep_U <- replicate(
  n = 5,
  expr = {
    modt1a <- stan(
      file = "stan_file/modelt1a.stan",
      data = list(
        N = N,
        g = data.norm$g,
        y1 = data.norm$y1,
        y2 = data.norm$y2,
        y3 = data.norm$y3,
        y4 = data.norm$y4,
        y5 = data.norm$y5,
        y6 = data.norm$y6,
        y7 = data.norm$y7,
        mu0 = sapply(data.norm[paste0("y", 1:7)], mean) %>% as.vector(),
        t = .5
      ),
      chains = 1,
      warmup = 500,
      iter = 1500
    )
    u_mean <- rstan::extract(modt1a, pars = c("U")) %>% as.data.frame() %>% sapply(mean)
    return(u_mean)
  }
)
# 65.82619 66.57391 67.88888 66.22054 62.90564 (*2)

vec_U <- sapply(
  X = vec_t, 
  FUN = function (t) {
    print(paste0("=====", t, "====="))
    modt1a <- stan(
      file = "stan_file/modelt1a.stan",
      data = list(
        N = N,
        g = data.norm$g,
        y1 = data.norm$y1,
        y2 = data.norm$y2,
        y3 = data.norm$y3,
        y4 = data.norm$y4,
        y5 = data.norm$y5,
        y6 = data.norm$y6,
        y7 = data.norm$y7,
        mu0 = sapply(data.norm[paste0("y", 1:7)], mean) %>% as.vector(),
        t = t
      ),
      chains = 1,
      warmup = 500,
      iter = 1500
    )
    u_mean <- rstan::extract(modt1a, pars = c("U")) %>% as.data.frame() %>% sapply(mean)
    return(u_mean)
  }
)

# *t
# U         U         U         U         U         U         U         U         U 
# 0.000000  3.431892 11.911836 27.281119 51.078600 66.100617 73.014773 76.164708 74.634696 
# U         U         U         U         U         U         U         U         U 
# 70.139466 65.035462 62.421820 55.830435 50.918417 42.723004 40.207634 38.325897 38.231044 
# U         U         U 
# 29.153170 25.316281 22.161019 

load(file = "output/Ut1a.RData"); vec_U
BF1a <- .5 * sum(diff(vec_t) * (vec_U[-1] + vec_U[-length(vec_U)])); BF1a
# 120.2323






# model 2a BF -------------------------------------------------------------

modt2a <- stan(
  file = "stan_file/modelt2a.stan",
  data = list(
    N = N,
    g = data.norm$g,
    y1 = data.norm$y1,
    y2 = data.norm$y2,
    y3 = data.norm$y3,
    y4 = data.norm$y4,
    y5 = data.norm$y5,
    y6 = data.norm$y6,
    y7 = data.norm$y7,
    mu0 = sapply(data.norm[paste0("y", 1:7)], mean) %>% as.vector(),
    t = .5
  ),
  chains = 1,
  warmup = 500,
  iter = 1500
)

rstan::extract(modt2a, pars = c("U")) %>% as.data.frame() %>% sapply(mean)
# 528.696 

vec_U <- sapply(
  X = vec_t, 
  FUN = function (t) {
    modt2a <- stan(
      file = "stan_file/modelt2a.stan",
      data = list(
        N = N,
        g = data.norm$g,
        y1 = data.norm$y1,
        y2 = data.norm$y2,
        y3 = data.norm$y3,
        y4 = data.norm$y4,
        y5 = data.norm$y5,
        y6 = data.norm$y6,
        y7 = data.norm$y7,
        mu0 = sapply(data.norm[paste0("y", 1:7)], mean) %>% as.vector(),
        t = t
      ),
      chains = 1,
      warmup = 500,
      iter = 1500
    )
    u_mean <- rstan::extract(modt2a, pars = c("U")) %>% as.data.frame() %>% sapply(mean)
    return(u_mean)
  }
)

load(file = "output/Ut2a.RData"); vec_U

BF2a <- .5 * sum(diff(vec_t) * (vec_U[-1] + vec_U[-length(vec_U)])); BF2a
# 119.2417


## model comparison --------------------------------------------------------

load(file = "output/Ut1a.RData"); vec_U
BF1a <- .5 * sum(diff(vec_t) * (vec_U[-1] + vec_U[-length(vec_U)])); BF1a
# 120.2323

load(file = "output/Ut2a.RData"); vec_U
BF2a <- .5 * sum(diff(vec_t) * (vec_U[-1] + vec_U[-length(vec_U)])); BF2a
# 119.2417

BF12 <- BF1a - BF2a; 2 * BF12
# 1.981227


# testing -----------------------------------------------------------------

stan(
  data = list(N = 5, 
              y = c(1, 1, 0, 1, 0)),
  file = "stan_file/test1.stan",
  chains = 4,
  warmup = 10,
  iter = 2000
)

stan(
  data = list(),
  file = "stan_file/test2.stan",
  chains = 4,
  warmup = 10,
  iter = 20
)



