load(paste0(getwd(),'/data/trueparams.RData'))
load(paste0(getwd(),'/q1_metrparams.RData'))
load(paste0(getwd(),'/Q1.2_metrparams.RData'))
load(paste0(getwd(),'/Q1.2_metrparams_20230503-14.RData'))



parameters = c("uby", "lam", "gam", "phx", "sgm")
other = c("ubeta", "sgd")

reportq13 <- function(est, tru) {
  if (length(tru)==1){
    return(list(
      bias = mean(est - tru),
      rmse = sqrt(mean((est - tru)^2))
    ))
    
  }else{
    return(list(
      bias = apply(sweep(est, 2, tru), 2, mean),
      rmse = sqrt(apply(sweep(est, 2, tru)^2, 2, mean))
    ))
  }
}

for (v in c(parameters,other)) {
  est = metr_params[['E']][[v]]
  tru = c(true_params[[v]])
  if (v == 'sgm') {
    tru = tru[4:7]
  }
  # Avg_bias = colMeans(est) - tru
  # RMSE = sqrt(colMeans(apply(est, 1, function(x){return ((x - tru)^2)})))
  print(v)
  res = reportq13(est, tru)
  print(res)
  # print('bias')
  # print(Avg_bias)
  # print('RMSE')
  # print(RMSE)
}

p = 0.5
n = 20
x = 0:20
b = dbinom(x,n,p)
xx = 0:200/10
a = dnorm(xx,n*p,sqrt(n*p*(1-p)))

plot(x,b,col='red', xlab = "x", ylab = "P")
lines(xx,a,col='blue',type='l')
legend(0, 0.16, legend=c("20 Bernoullis", "CLT approximation"),
       col=c("red", "blue"), lty = 1)
