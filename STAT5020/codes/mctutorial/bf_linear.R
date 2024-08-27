library(R2WinBUGS)  #Load R2WinBUGS package

timestamp = strftime(Sys.time(), "%Y%m%d-%H")
winBUGS.path = "D:/pkgs/WinBUGS14/"
datapath = paste0(getwd(),'/data')

iter = 10
cut = 20

NY = 9  # dimension of Y
Neta = 1  # dimension of eta
Nxi = 2  # dimension of xi
Ngam = 4  # dimension of gamma
R = matrix(c(1, 0.3, 0.3, 1), nrow = 2)


lbf = numeric(iter)


parameters = c("ubar")

init1 = list(u = rep(0, NY), lam = rep(0, NY - Neta - Nxi), b = 0, 
             a = rep(0, NY), gam = rep(0, Ngam), psi = rep(1, NY), 
             psd = 1, phi = matrix(c(1, 0, 0, 1), nrow = 2))

inits = list(init1)

# Path sampling
for (r in 1:10) {
  iterpath = paste0(getwd(), "/bflinear", r)
  dir.create(iterpath, showWarnings = FALSE, recursive = TRUE)
  
  # load previous dataset
  Y = as.matrix(read.table(paste0(datapath, "/Y-", r, ".txt")))
  BD = read.table(paste0(datapath, "/BD-", r, ".txt"))$x
  BC = read.table(paste0(datapath, "/BC-", r, ".txt"))$x
  
  data = list(N = 500, zero = c(0, 0), d = BD, c = BC, R = R, y = Y, t = NA)
  
  u = numeric(cut)
  for (i in 1:cut) {
    data$t <- (i - 1)/(cut - 1)
    
    model = bugs(data, inits, parameters, 
                 model.file = paste0(iterpath, "/../model_BF_linear.txt"), 
                 n.chains = 1, n.iter = 3000,
                 n.burnin = 1000, n.thin = 1, bugs.directory = winBUGS.path,
                 working.directory = iterpath)
    
    u[i] <- model$mean$ubar
    
  }
  
  # Caluate log Bayes factor
  logBF = 0
  for (i in 1:(cut - 1)) {
    logBF = logBF + (u[i + 1] + u[i])/(2 * (cut - 1))
  }
  
  lbf[r] = logBF
}

resultlst = list(
  lbf = lbf, lbf.mean = mean(lbf), lbf.sd = sd(lbf)
)

save(resultlst, file = paste0(getwd(), '/bflinear-', timestamp, ".RData"))