
setwd("D:/OneDrive/OneDrive - sjtu.edu.cn/Records/CUHK/Course/STAT5020/codes/finalproj")
datestamp = "20230517-12"
load(paste0("model1_",datestamp,"/model.RData")); model1 = model
load(paste0("model2_",datestamp,"/model2.RData")); model2 = model
load(paste0("model3_",datestamp,"/model3.Rdata")); model3 = model
models = list(model1, model2, model3)


for (m in 1:3) {
  print(paste("model", m, "DIC and DICbyR"))
  print(models[[m]]$DIC)
}

est = model1$summary
write.csv(est, file = "realdatamodel1est.csv")
chain = model1$sims.matrix

rho12g1 = chain[,"phx1[1,2]"] / sqrt(chain[,"phx1[1,1]"] * chain[,"phx1[2,2]"])
rho13g1 = chain[,"phx1[1,3]"] / sqrt(chain[,"phx1[1,1]"] * chain[,"phx1[3,3]"])
rho23g1 = chain[,"phx1[2,3]"] / sqrt(chain[,"phx1[2,2]"] * chain[,"phx1[3,3]"])
rho12g2 = chain[,"phx2[1,2]"] / sqrt(chain[,"phx2[1,1]"] * chain[,"phx2[2,2]"])
rho13g2 = chain[,"phx2[1,3]"] / sqrt(chain[,"phx2[1,1]"] * chain[,"phx2[3,3]"])
rho23g2 = chain[,"phx2[2,3]"] / sqrt(chain[,"phx2[2,2]"] * chain[,"phx2[3,3]"])
rho.df = data.frame(cbind(rho12g1, rho12g2, rho13g1, rho13g2, rho23g1, rho23g2))
apply(rho.df, 2, mean)
apply(rho.df, 2, sd)

true_params = readRDS("simdata/true_params.Rds")
params = names(true_params)[1:(length(true_params)-2)]
true_arr = c()
name_arr = c()
for (i in 1:length(params)) {
  arr = c(true_params[[params[i]]])
  nmm = paste0(params[i], "[", 1:length(arr), "]")
  if (params[i] == "phx1" | params[i] == "phx2"){
    tmp = c("1,1","1,2","1,3","2,1","2,2","2,3","3,1","3,2","3,3")
    nmm = paste0(params[i], "[", tmp, "]")
  } else if (params[i] == "sgd1" | params[i] == "sgd2") {
    nmm = params[i]
  }
  true_arr = c(true_arr, arr)
  name_arr = c(name_arr, nmm)
}

summary.df = data.frame(true = true_arr, row.names = name_arr)
for (t in 1:5) {
  load(paste0("sim_model1_",t,"/model.RData"))
  simmodel = model
  tmp.df = data.frame(simmodel$summary)
  colsel = tmp.df[rownames(summary.df),"mean"]
  summary.df[,paste0("rep",t)] = colsel
}

saveRDS(summary.df, "simdata.summary.Rds")

summary.short = t(apply(summary.df, 1, function(x){
  tval = x[1]
  ests = x[2:6]
  return(c(tval, mean(ests)-tval, sqrt(mean((ests-tval)^2))))
}))

colnames(summary.short) = c("true", "bias", "RMSE")

write.csv(summary.short, file = "simdata.summary.csv")
