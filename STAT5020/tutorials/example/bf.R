library(R2WinBUGS) #Load R2WinBUGS package

NY=9 # dimension of Y
Neta=1 # dimension of eta
Nxi=2 # dimension of xi
Ngam=3 # dimension of gamma

Y=matrix(as.matrix(read.table("Y-1.txt")),ncol = NY)
BD = read.table("BD-1.txt")[,1]
BC = read.table("BC-1.txt")[,1]

data=list(N=500,zero=c(0,0),d=BD,c=BC,R=matrix(c(1, 0, 0, 1), nrow=2),y=Y,t=NA)

#data=list(N=300, P=9, ux=c(0,0), R=matrix(c(1, 0, 0, 1), nrow=2), y=Y, t=NA) #Data

init1=list(u=rep(0,NY),lam=rep(0,NY-Neta-Nxi),b=0,a=rep(0,NY),gam=rep(0,Ngam),psi=rep(1,NY),psd=1,phi=matrix(c(1,0,0,1),nrow=2))

inits=list(init1)

parameters=c("ubar")

cut = 5
u=numeric(cut)

#Path sampling
for (i in 1:cut) {
  data$t<-(i-1)/(cut-1)
  
  model=bugs(data,inits,parameters,model.file="model_BF.txt",n.chains=1,n.iter = 3000,n.burnin = 1000,n.thin=1,bugs.directory = "D:/winbugs14_full_patched/WinBUGS14",working.directory = getwd())
  
  u[i]<-model$mean$ubar
}

#Caluate log Bayes factor
logBF=0
for (i in 1:(cut-1)) { logBF=logBF+(u[i+1]+u[i])/(2*(cut-1)) }
print(logBF)