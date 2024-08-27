library(mvtnorm)


work.path = "D:/OneDrive/OneDrive - sjtu.edu.cn/Records/CUHK/Course/STAT5020/codes/finalproj"
data.path = paste0(work.path, "/simdata")
dir.create(data.path, showWarnings = FALSE, recursive = TRUE)

set.seed(5020)

# simulation setting
iter = 5
Ng = 1000

mu.y1 = rep(0.1, 17)
mu.y2 = rep(0.0, 17)
lam1 = c(0.9, 0.8, 0.8, 0.7,
         1.1, 1.2, 1.3, 1.2,
         rep(0.3, 5))
lam2 = c(1.1, 1.2, 1.3, 1.2,
         0.9, 0.8, 0.8, 0.7,
         rep(0.5, 5))
sgm1 = c(rep(1, 11), rep(0.2, 6))
sgm2 = c(rep(0.7, 11), rep(0.2, 6))
sgd1 = 0.5
sgd2 = 0.5
gam1 = c(0.5, 0.2, 0.6)
gam2 = c(0.5, 0.2, 0.6)
phx1 = matrix(c(1, -0.2, 0,
                -0.2, 1, -0.5,
                0, -0.5, 1), ncol = 3)
phx2 = matrix(c(1, -0.2, 0,
               -0.2, 1, -0.5,
               0, -0.5, 1), ncol = 3)
rata = c(7, 2, 1)
ratb = c(6, 2, 1 ,0.5, 0.5)
cumpa = c(0,cumsum(rata)/sum(rata))
cumpb = c(0,cumsum(ratb)/sum(ratb))

Z1 = matrix(data = NA, nrow = Ng, ncol = 17)
Z2 = matrix(data = NA, nrow = Ng, ncol = 17)
V1 = matrix(data = NA, nrow = Ng, ncol = 17)
V2 = matrix(data = NA, nrow = Ng, ncol = 17)


for (t in 1:iter) {
  for (i in 1:Ng) {
    # generate xi
    xi1 = rmvnorm(1, c(0,0,0), phx1)
    xi2 = rmvnorm(1, c(0,0,0), phx2)
    # generate error term in SE
    del1 = rnorm(1, 0, sqrt(sgd1))
    del2 = rnorm(1, 0, sqrt(sgd2))
    # generate eta according to the SE
    eta1 = gam1[1]*xi1[1] + gam1[2]*xi1[2] + gam1[3]*xi1[3] + del1
    eta2 = gam2[1]*xi2[1] + gam2[2]*xi2[2] + gam2[3]*xi2[3] + del2
    # generate error term in ME
    eps1 = numeric(17)
    eps2 = numeric(17)
    for (k in 1:17) {
      eps1[k] = rnorm(1, 0, sqrt(sgm1[k]))
      eps2[k] = rnorm(1, 0, sqrt(sgm2[k]))
    }
    # generate theta in ME
    V1[i,1] = mu.y1[1] + eta1 + eps1[1]
    V1[i,2] <- mu.y1[2] + xi1[1] + eps1[2]
    V1[i,3] <- mu.y1[3] + lam1[1] * xi1[1] + eps1[3]
    V1[i,4] <- mu.y1[4] + lam1[2] * xi1[1] + eps1[4]
    V1[i,5] <- mu.y1[5] + lam1[3] * xi1[1] + eps1[5]
    V1[i,6] <- mu.y1[6] + lam1[4] * xi1[1] + eps1[6]
    V1[i,7] <- mu.y1[7] + xi1[2] + eps1[7]
    V1[i,8] <- mu.y1[8] + lam1[5] * xi1[2] + eps1[8]
    V1[i,9] <- mu.y1[9] + lam1[6] * xi1[2] + eps1[9]
    V1[i,10] <- mu.y1[10] + lam1[7] * xi1[2] + eps1[10]
    V1[i,11] <- mu.y1[11] + lam1[8] * xi1[2] + eps1[11]
    V1[i,12] <- mu.y1[12] + xi1[3] + eps1[12]
    V1[i,13] <- mu.y1[13] + lam1[9] * xi1[3] + eps1[13]
    V1[i,14] <- mu.y1[14] + lam1[10] * xi1[3] + eps1[14]
    V1[i,15] <- mu.y1[15] + lam1[11] * xi1[3] + eps1[15]
    V1[i,16] <- mu.y1[16] + lam1[12] * xi1[3] + eps1[16]
    V1[i,17] <- mu.y1[17] + lam1[13] * xi1[3] + eps1[17]
    V2[i,1] = mu.y2[1] + eta2 + eps2[1]
    V2[i,2] <- mu.y2[2] + xi2[1] + eps2[2]
    V2[i,3] <- mu.y2[3] + lam2[1] * xi2[1] + eps2[3]
    V2[i,4] <- mu.y2[4] + lam2[2] * xi2[1] + eps2[4]
    V2[i,5] <- mu.y2[5] + lam2[3] * xi2[1] + eps2[5]
    V2[i,6] <- mu.y2[6] + lam2[4] * xi2[1] + eps2[6]
    V2[i,7] <- mu.y2[7] + xi2[2] + eps2[7]
    V2[i,8] <- mu.y2[8] + lam2[5] * xi2[2] + eps2[8]
    V2[i,9] <- mu.y2[9] + lam2[6] * xi2[2] + eps2[9]
    V2[i,10] <- mu.y2[10] + lam2[7] * xi2[2] + eps2[10]
    V2[i,11] <- mu.y2[11] + lam2[8] * xi2[2] + eps2[11]
    V2[i,12] <- mu.y2[12] + xi2[3] + eps2[12]
    V2[i,13] <- mu.y2[13] + lam2[9] * xi2[3] + eps2[13]
    V2[i,14] <- mu.y2[14] + lam2[10] * xi2[3] + eps2[14]
    V2[i,15] <- mu.y2[15] + lam2[11] * xi2[3] + eps2[15]
    V2[i,16] <- mu.y2[16] + lam2[12] * xi2[3] + eps2[16]
    V2[i,17] <- mu.y2[17] + lam2[13] * xi2[3] + eps2[17]
  }
  # transform theta to ordinal variables
  avg1 = apply(V1[,1:11],2,mean)
  sd1 = apply(V1[,1:11],2,sd)
  thda = t(apply(cbind(avg1, sd1)[1:6,], 1, function(x){
    return(qnorm(cumpa, mean = x[1], sd = x[2]))
  }))
  thdb = t(apply(cbind(avg1, sd1)[7:11,], 1, function(x){
    return(qnorm(cumpb, mean = x[1], sd = x[2]))
  }))
  thda[,1] = -200
  thdb[,1] = -200
  thda[,ncol(thda)] = 200
  thdb[,ncol(thdb)] = 200
  for (j in 1:6) {
    Z1[,j] = (V1[,j] > thda[j,1]) + (V1[,j] > thda[j,2]) + (V1[,j] > thda[j,3])
    Z2[,j] = (V2[,j] > thda[j,1]) + (V2[,j] > thda[j,2]) + (V2[,j] > thda[j,3])
  }
  for (j in 7:11) {
    Z1[,j] = (V1[,j] > thdb[j-6,1]) + (V1[,j] > thdb[j-6,2]) + 
      (V1[,j] > thdb[j-6,3]) + (V1[,j] > thdb[j-6,4]) + (V1[,j] > thdb[j-6,5])
    Z2[,j] = (V2[,j] > thdb[j-6,1]) + (V2[,j] > thdb[j-6,2]) + 
      (V2[,j] > thdb[j-6,3]) + (V2[,j] > thdb[j-6,4]) + (V2[,j] > thdb[j-6,5])
  }
  for (j in 12:17) {
    Z1[,j] = V1[,j]
    Z2[,j] = V2[,j]
  }
  g = "M"
  colnames(Z1) = paste0("z",1:17)
  Z1.cb = cbind(g, data.frame(Z1))
  g = "F"
  colnames(Z2) = paste0("z",1:17)
  Z2.cb = cbind(g, data.frame(Z2))
  Z = rbind(Z1.cb, Z2.cb)
  saveRDS(Z, file = paste0(data.path,"/simZ",t,".Rds"))
  saveRDS(thda, file = paste0(data.path,"/thda",t,".Rds"))
  saveRDS(thdb, file = paste0(data.path,"/thdb",t,".Rds"))
}

true_params = list(
  mu.y1 = mu.y1,
  mu.y2 = mu.y2,
  lam1 = lam1,
  lam2 = lam2,
  sgm1 = sgm1,
  sgm2 = sgm2,
  sgd1 = sgd1,
  sgd2 = sgd2,
  gam1 = gam1,
  gam2 = gam2,
  phx1 = phx1,
  phx2 = phx2,
  rata = rata,
  ratb = ratb
  
)
saveRDS(true_params, file = paste0(data.path,"/true_params.Rds"))

