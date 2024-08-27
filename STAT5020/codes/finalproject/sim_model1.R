library(R2WinBUGS)

setwd("D:/OneDrive/OneDrive - sjtu.edu.cn/Records/CUHK/Course/STAT5020/codes/finalproj")
timestamp = strftime(Sys.time(), "%Y%m%d-%H")
winBUGS.path = "D:/pkgs/winbugs14_full_patched/WinBUGS14"
set.seed(5020)

for (t in 1:5) {
  # prepare data
  df.full = readRDS(paste0("simdata/simZ",t,".Rds"))
  df = df.full
  Z1 = as.matrix(df[df$g=="M",paste0("z",1:17)])
  Z2 = as.matrix(df[df$g=="F",paste0("z",1:17)])
  N1 = dim(Z1)[1]
  N2 = dim(Z2)[1]
  P = dim(Z1)[2]
  
  thdpsych = as.matrix(readRDS(paste0("simdata/thda",t,".Rds")))
  thdsleep = as.matrix(readRDS(paste0("simdata/thdb",t,".Rds")))
  
  R0 = matrix(c(
    3, 1, 1,
    1, 3, 1,
    1, 1, 3
  ), nrow = 3)
  
  # parameter to be estimate
  parameters = c("mu.y1", "mu.y2", "lam1", "lam2", 
                 "sgm1", "sgm2", "sgd1", "sgd2", 
                 "gam1", "gam2", "phx1", "phx2")
  
  
  # initial values for MCMC in WinBUGS
  init1 <- list(
    mu.y1 = rep(-0.2, 17),
    mu.y2 = rep(-0.2, 17),
    lam1 = c(rep(0.9, 8), rep(0.3, 5)),
    lam2 = c(rep(0.9, 8), rep(0.3, 5)),
    psi1 = c(rep(0.9, 11), rep(6, 6)),
    psi2 = c(rep(0.9, 11), rep(6, 6)),
    psd1 = 1,
    psd2 = 1,
    gam1 = rep(0.3, 3),
    gam2 = rep(0.3, 3),
    phi1 = matrix(c(1, 0, 0,
                    0, 1, 0,
                    0, 0, 1), nrow = 3),
    phi2 = matrix(c(1, 0, 0,
                    0, 1, 0,
                    0, 0, 1), nrow = 3)
  )
  
  init2 <- list(
    mu.y1 = rep(0.2, 17),
    mu.y2 = rep(0.2, 17),
    lam1 = c(rep(1.2, 8), rep(0.5, 5)),
    lam2 = c(rep(0.9, 8), rep(0.5, 5)),
    psi1 = c(rep(1.1, 11), rep(4, 6)),
    psi2 = c(rep(1.1, 11), rep(4, 6)),
    psd1 = 3,
    psd2 = 3,
    gam1 = rep(0.7, 3),
    gam2 = rep(0.7, 3),
    phi1 = matrix(c(2, 0, 0,
                    0, 2, 0,
                    0, 0, 2), nrow = 3),
    phi2 = matrix(c(2, 0, 0,
                    0, 2, 0,
                    0, 0, 2), nrow = 3)
  )
  
  inits = list(init1, init2)
  
  
  # set experiment path
  iterpath = paste0(getwd(),"/sim_model1_", t)
  dir.create(iterpath, showWarnings = FALSE, recursive = TRUE)
  
  data <- list(N1 = N1, N2 = N2,
               R = R0,
               z1 = Z1, z2 = Z2,
               thdpsych = thdpsych,
               thdsleep = thdsleep)
  
  model = bugs(data, inits, parameters, 
               model.file = paste0(getwd(),"/../model1.txt"), 
               n.chains = 2, 
               n.iter = 5000, 
               n.burnin = 2500, 
               n.thin = 1, 
               DIC = TRUE,
               bugs.directory = winBUGS.path,
               working.directory = iterpath, 
               debug = FALSE)
  
  save(model, file = paste0(iterpath, "/model.RData"))
}


