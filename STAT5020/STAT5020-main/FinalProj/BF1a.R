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
