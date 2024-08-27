## -------------------------------------------------------------------------------------------------------------------------
library(rstan)
library(tidyverse)
data_raw <- read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/Course/STAT5020/Ji_final_project/diabetes_012_health_indicators_BRFSS2015.csv")
data_raw = data_raw[1:1000, ]

data <- data_raw %>%
  rename(y1=Diabetes_012, y2=HighBP, y3=HighChol, y4=HvyAlcoholConsump, y5=Fruits, y6=Veggies, y7=PhysActivity, y8=Smoker, y9=DiffWalk, y10=BMI, y11=Stroke, y12=HeartDiseaseorAttack) %>% mutate(g = if_else(Income < 5, 1, 2)) %>% select(starts_with("y"), g) %>%  mutate(across(everything(), as.integer))

N = nrow(data)

data


## -------------------------------------------------------------------------------------------------------------------------
m0 <- stan(
  file = "m0.stan",
  data = list(
    N = N,
    y1 = data$y1,
    y2 = data$y2,
    y3 = data$y3,
    y4 = data$y4,
    y5 = data$y5,
    y6 = data$y6,
    y7 = data$y7,
    y8 = data$y8,
    y9 = data$y9,
    y10 = data$y10,
    y11 = data$y11,
    y12 = data$y12,
    g = data$g,
    mu0 = sapply(data[paste0("y", 1:12)], mean) %>% as.vector()
  )
)
save(m0, file = "./m0.RData")

## -------------------------------------------------------------------------------------------------------------------------
load(file = "./m0.RData")

traceplot(
  m0,
  pars = c(
    'lp__', 
    'mu', 
    'lam',
    'gam')
)

traceplot(
  m0,
  pars = c(
    'psi_e',
    'psi_d',
    'mat_phi'
  )
)

print(
  m0,
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

pdf(file = "./m0.pdf",
    width = 20,
    height = 20)

traceplot(
  m0,
  pars = c(
    'lp__', 
    'mu', 
    'lam',
    'gam',
    'psi_e',
    'psi_d',
    'mat_phi'
  ), 
  nrow = 12
)

dev.off()

est_mod0 <- rstan::extract(
  m0,
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
  write_csv("m0.csv")


## -------------------------------------------------------------------------------------------------------------------------
m1 <- stan(
  file = "m1.stan",
  data = list(
    N = N,
    y1 = data$y1,
    y2 = data$y2,
    y3 = data$y3,
    y4 = data$y4,
    y5 = data$y5,
    y6 = data$y6,
    y7 = data$y7,
    y8 = data$y8,
    y9 = data$y9,
    y10 = data$y10,
    y11 = data$y11,
    y12 = data$y12,
    g = data$g,
    mu0 = sapply(data[paste0("y", 1:12)], mean) %>% as.vector()
  )
)
save(m1, file = "m1.RData")


## -------------------------------------------------------------------------------------------------------------------------
load(file = "m1.RData")

traceplot(
  m1,
  pars = c(
    'lp__', 
    'mu', 
    'lam',
    'gam')
)

traceplot(
  m1,
  pars = c(
    'psi_e',
    'psi_d',
    'mat_phi'
  )
)

print(
  m1,
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

pdf(file = "m1.pdf",
    width = 30,
    height = 36)

traceplot(
  m1,
  pars = c(
    'lp__',
    'mu',
    'lam',
    'gam',
    'psi_e',
    'psi_d',
    'mat_phi'
  ),
  nrow = 15
)

dev.off()


est_mod1 <- rstan::extract(
  m1,
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
  write_csv("m1.csv")


## -------------------------------------------------------------------------------------------------------------------------
vec_t <- seq(from = 0, to = 1, by = 0.05)
J <- length(vec_t)


modt1a <- stan(
  file = "mt1a.stan",
  data = list(
    N = N,
    g = data$g,
    y1 = data$y1,
    y2 = data$y2,
    y3 = data$y3,
    y4 = data$y4,
    y5 = data$y5,
    y6 = data$y6,
    y7 = data$y7,
    y8 = data$y8,
    y9 = data$y9,
    y10 = data$y10,
    y11 = data$y11,
    y12 = data$y12,
    mu0 = sapply(data[paste0("y", 1:12)], mean) %>% as.vector(),
    t = 0
  ),
  chains = 1,
  warmup = 500,
  iter = 1500
)

rstan::extract(modt1a, pars = c("U")) %>% as.data.frame() %>% sapply(mean)

rep_U <- replicate(
  n = 5,
  expr = {
    modt1a <- stan(
      file = "mt1a.stan",
      data = list(
        N = N,
        g = data$g,
        y1 = data$y1,
        y2 = data$y2,
        y3 = data$y3,
        y4 = data$y4,
        y5 = data$y5,
        y6 = data$y6,
        y7 = data$y7,
        y8 = data$y8,
        y9 = data$y9,
        y10 = data$y10,
        y11 = data$y11,
        y12 = data$y12,
        mu0 = sapply(data[paste0("y", 1:12)], mean) %>% as.vector(),
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

vec_U <- sapply(
  X = vec_t, 
  FUN = function (t) {
    print(paste0("=====", t, "====="))
    modt1a <- stan(
      file = "mt1a.stan",
      data = list(
        N = N,
        g = data$g,
        y1 = data$y1,
        y2 = data$y2,
        y3 = data$y3,
        y4 = data$y4,
        y5 = data$y5,
        y6 = data$y6,
        y7 = data$y7,
        y8 = data$y8,
        y9 = data$y9,
        y10 = data$y10,
        y11 = data$y11,
        y12 = data$y12,
        mu0 = sapply(data[paste0("y", 1:12)], mean) %>% as.vector(),
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

BF1a <- .5 * sum(diff(vec_t) * (vec_U[-1] + vec_U[-length(vec_U)])); BF1a


## -------------------------------------------------------------------------------------------------------------------------
modt2a <- stan(
  file = "mt2a.stan",
  data = list(
    N = N,
    g = data$g,
    y1 = data$y1,
    y2 = data$y2,
    y3 = data$y3,
    y4 = data$y4,
    y5 = data$y5,
    y6 = data$y6,
    y7 = data$y7,
    y8 = data$y8,
    y9 = data$y9,
    y10 = data$y10,
    y11 = data$y11,
    y12 = data$y12,
    mu0 = sapply(data[paste0("y", 1:12)], mean) %>% as.vector(),
    t = .5
  ),
  chains = 1,
  warmup = 500,
  iter = 1500
)

rstan::extract(modt2a, pars = c("U")) %>% as.data.frame() %>% sapply(mean)

vec_U <- sapply(
  X = vec_t, 
  FUN = function (t) {
    modt2a <- stan(
      file = "mt2a.stan",
      data = list(
        N = N,
        g = data$g,
        y1 = data$y1,
        y2 = data$y2,
        y3 = data$y3,
        y4 = data$y4,
        y5 = data$y5,
        y6 = data$y6,
        y7 = data$y7,
        y8 = data$y8,
        y9 = data$y9,
        y10 = data$y10,
        y11 = data$y11,
        y12 = data$y12,
        mu0 = sapply(data[paste0("y", 1:12)], mean) %>% as.vector(),
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

vec_U

BF2a <- .5 * sum(diff(vec_t) * (vec_U[-1] + vec_U[-length(vec_U)])); BF2a



## -------------------------------------------------------------------------------------------------------------------------
BF1a <- .5 * sum(diff(vec_t) * (vec_U[-1] + vec_U[-length(vec_U)])); BF1a
BF2a <- .5 * sum(diff(vec_t) * (vec_U[-1] + vec_U[-length(vec_U)])); BF2a
BF12 <- BF1a - BF2a; 2 * BF12

