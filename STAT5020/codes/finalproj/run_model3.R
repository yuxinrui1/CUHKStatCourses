library(R2WinBUGS)

setwd("D:/OneDrive/OneDrive - sjtu.edu.cn/Records/CUHK/Course/STAT5020/codes/finalproj")
timestamp = strftime(Sys.time(), "%Y%m%d-%H")
winBUGS.path = "D:/pkgs/winbugs14_full_patched/WinBUGS14"
set.seed(5020)

# prepare data
df.full = readRDS("abcd_sel.cleaned.onehot.Rds")
df = df.full
Z1 = as.matrix(df[df$g=="M",paste0("z",1:17)])
Z2 = as.matrix(df[df$g=="F",paste0("z",1:17)])
N1 = dim(Z1)[1]
N2 = dim(Z2)[1]
P = dim(Z1)[2]

thdpsych = matrix(rep(0,6*4), nrow = 6)
thdsleep = matrix(rep(0,5*6), nrow = 5)
for (j in 1:6) {
  thdpsych[j,1] = -200
  thdpsych[j,2] = min(200, qnorm(mean(Z1[,j]<2)))
  thdpsych[j,3] = min(200, qnorm(mean(Z1[,j]<3)))
  thdpsych[j,4] = 200
}
for (j in 1:5) {
  thdsleep[j,1] = -200
  thdsleep[j,2] = min(200, qnorm(mean(Z1[,j+6]<2)))
  thdsleep[j,3] = min(200, qnorm(mean(Z1[,j+6]<3)))
  thdsleep[j,4] = min(200, qnorm(mean(Z1[,j+6]<4)))
  thdsleep[j,5] = min(200, qnorm(mean(Z1[,j+6]<5)))
  thdsleep[j,6] = 200
}

R0 = matrix(c(
  5, 2, 2,
  2, 5, 2,
  2, 2, 5
), nrow = 3)

# parameter to be estimate
parameters = c("mu.y1", "mu.y2", "lam", 
               "sgm1", "sgm2", "sgd1", "sgd2", 
               "gam1", "gam2", "phx")


# initial values for MCMC in WinBUGS
init1 <- list(
  mu.y1 = rep(0.5, 17),
  mu.y2 = rep(0.5, 17),
  lam = rep(0.5, 13),
  psi1 = rep(0.3, 17),
  psi2 = rep(0.3, 17),
  psd1 = 0.5,
  psd2 = 0.5,
  gam1 = rep(0.5, 3),
  gam2 = rep(0.5, 3),
  phi = matrix(c(1, 0, 0,
                  0, 1, 0,
                  0, 0, 1), nrow = 3),
  xi1 = matrix(rep(0.3, N1 * 3), ncol = 3),
  xi2 = matrix(rep(0.3, N2 * 3), ncol = 3)
)

inits = list(init1)


# set experiment path
iterpath = paste0(getwd(),"/model3_", timestamp)
dir.create(iterpath, showWarnings = FALSE, recursive = TRUE)

data <- list(N1 = N1, N2 = N2,
             R = R0,
             z1 = Z1, z2 = Z2,
             thdpsych = thdpsych,
             thdsleep = thdsleep)

model = bugs(data, inits, parameters, 
             model.file = paste0(getwd(),"/../model3.txt"), 
             n.chains = 1, 
             n.iter = 10000, 
             n.burnin = 5000, 
             n.thin = 1, 
             DIC = TRUE,
             bugs.directory = winBUGS.path,
             working.directory = iterpath, 
             debug = FALSE)

save(model, file = paste0(iterpath, "/model3.RData"))

