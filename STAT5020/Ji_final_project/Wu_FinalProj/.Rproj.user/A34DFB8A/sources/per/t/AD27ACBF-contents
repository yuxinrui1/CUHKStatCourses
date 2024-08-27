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

vec_t <- seq(from = 0, to = 1, by = 0.05)
J <- length(vec_t)

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
    print(paste0("=====", t, "=====", u_mean, "====="))
    return(u_mean)
  }
)

save(vec_U, file = "output/Ut2a.RData")