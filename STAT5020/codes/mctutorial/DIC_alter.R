library(mvtnorm)
library(R2WinBUGS)

timestamp = strftime(Sys.time(), "%Y%m%d-%H")
winBUGS.path = "D:/pkgs/WinBUGS14/"
datapath = paste0(getwd(),'/data')


iter = 10
NY = 9  # dimension of Y
Neta = 1  # dimension of eta
Nxi = 2  # dimension of xi
Ngam = 5  # dimension of gamma
R = matrix(c(1, 0.3, 0.3, 1), nrow = 2)


dic = numeric(iter)


parameters = c("u", "lam", "b", "a", "gam", "sgm", "sgd", "phx")

init1 = list(u = rep(0, NY), lam = rep(0, NY - Neta - Nxi), b = 0, 
             a = rep(0, NY), gam = rep(0, Ngam), psi = rep(1, NY), 
             psd = 1, phi = matrix(c(1, 0, 0, 1), nrow = 2))

inits = list(init1)

for (r in 1:iter) {
  iterpath = paste0(getwd(),"/dicalter",r)
  dir.create(iterpath, showWarnings = FALSE, recursive = TRUE)
  
  # load previous dataset
  Y = as.matrix(read.table(paste0(datapath, "/Y-", r, ".txt")))
  BD = read.table(paste0(datapath, "/BD-", r, ".txt"))$x
  BC = read.table(paste0(datapath, "/BC-", r, ".txt"))$x
  
  # Run WINBUGS
  data = list(N = 500, zero = c(0, 0), d = BD, c = BC, R = R, y = Y)
  
  model = bugs(data, inits, parameters, model.file = paste0(getwd(),"/../model_alter.txt"), 
               n.chains = 1, n.iter = 3000, n.burnin = 1000, 
               n.thin = 1, bugs.directory = winBUGS.path,
               working.directory = iterpath, debug = FALSE)
  
  dic[r] = model$DIC
}

resultlst = list(
  dic = dic, dic.mean = mean(dic), dic.sd = sd(dic)
)

print(resultlst)

save(resultlst, file = paste0(getwd(),'/dicalter-', timestamp, ".RData"))


