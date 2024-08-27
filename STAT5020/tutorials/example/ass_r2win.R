library(mvtnorm)
library(R2WinBUGS)

iter = 1
NY=9 # dimension of Y
Neta=1 # dimension of eta
Nxi=2 # dimension of xi
Ngam=3 # dimension of gamma

N = 500
BD=numeric(N)
BC=numeric(N)
XI=matrix(NA,nrow=N,ncol=2)
Eta=numeric(N)
Y=matrix(NA,nrow=N,ncol=NY)

# The covariance matrix of xi
phi=matrix(c(1,0.3,0.3,1),nrow=2)

# Estimates and standard error estimates
Eu=matrix(NA,nrow=iter,ncol=NY)
SEu=matrix(NA,nrow=iter,ncol=NY)
Elam=matrix(NA,nrow=iter,ncol=NY-Neta-Nxi)
SElam=matrix(NA,nrow=iter,ncol=NY-Neta-Nxi)
Eb=numeric(iter)
SEb=numeric(iter)
Ea=matrix(NA,nrow=iter,ncol=NY)
SEa=matrix(NA,nrow=iter,ncol=NY)
Egam=matrix(NA,nrow=iter,ncol=Ngam)
SEgam=matrix(NA,nrow=iter,ncol=Ngam)
Esgm=matrix(NA,nrow=iter,ncol=NY)
SEsgm=matrix(NA,nrow=iter,ncol=NY)
Esgd=numeric(iter)
SEsgd=numeric(iter)
Ephx=matrix(NA,nrow=iter,ncol=3)
SEphx=matrix(NA,nrow=iter,ncol=3)

R=matrix(c(1,0.3,0.3,1),nrow=2)

parameters=c("u","lam","b","a","gam","sgm","sgd","phx")

init1=list(u=rep(0,NY),lam=rep(0,NY-Neta-Nxi),b=0,a=rep(0,NY),gam=rep(0,Ngam),psi=rep(1,NY),psd=1,phi=matrix(c(1,0,0,1),nrow=2))

init2=list(u=rep(1,NY),lam=rep(1,NY-Neta-Nxi),b=1,a=rep(1,NY),gam=rep(1,Ngam),psi=rep(2,NY),psd=2,phi=matrix(c(2,0,0,2),nrow=2))

inits=list(init1,init2)

eps=numeric(NY)

for(t in 1:iter){
  # generate data
  for(i in 1:N){
    BD[i]=rt(1,5)
    BC[i]=rt(1,5)
    
    XI[i,]=rmvnorm(1,c(0,0),phi)
    
    delta=rnorm(1,0,sqrt(0.36))
    
    Eta[i]=0.5*BD[i]+0.4*XI[i,1]+0.3*XI[i,2]+0.5*XI[i,1]*XI[i,2]+delta
    
    eps[1:3]=rnorm(3,0,sqrt(0.3))
    eps[4:6]=rnorm(3,0,sqrt(0.4))
    eps[7:9]=rnorm(3,0,sqrt(0.5))
    
    Y[i,1]=0.2*BC[i]+Eta[i]+eps[1]
    Y[i,2]=-0.2*BC[i]+0.9*Eta[i]+eps[2]
    Y[i,3]=0.4*BC[i]+0.6*Eta[i]+eps[3]
    Y[i,4]=0.3*BC[i]+XI[i,1]+eps[4]
    Y[i,5]=-0.2*BC[i]+0.7*XI[i,1]+eps[5]
    Y[i,6]=0.4*BC[i]+0.9*XI[i,1]+eps[6]
    Y[i,7]=0.5*BC[i]+XI[i,2]+eps[7]
    Y[i,8]=-0.4*BC[i]+0.8*XI[i,2]+eps[8]
    Y[i,9]=0.3*BC[i]+0.6*XI[i,2]+eps[9]
    
    
  }
  
  # Run WINBUGS
  data=list(N=500,zero=c(0,0),d=BD,c=BC,R=R,y=Y)
  
  write.table(Y,paste("Y-",t,".txt",sep=""))
  write.table(BD,paste("BD-",t,".txt",sep=""))
  write.table(BC,paste("BC-",t,".txt",sep=""))
  
  model=bugs(data,inits,parameters,model.file="model.txt",n.chains=2,n.iter = 3000,n.burnin = 1000,n.thin=1,
             bugs.directory = "D:/winbugs14_full_patched/WinBUGS14",working.directory = getwd(), debug = TRUE)
  
  # save estimates
  Eu[t,]=model$mean$u
  SEu[t,]=model$sd$u
  Elam[t,]=model$mean$lam
  SElam[t,]=model$sd$lam
  Eb[t]=model$mean$b
  SEb[t]=model$sd$b
  Ea[t,]=model$mean$a
  SEa[t,]=model$sd$a
  Egam[t,]=model$mean$gam
  SEgam[t,]=model$sd$gam
  Esgm[t,]=model$mean$sgm
  SEsgm[t,]=model$sd$sgm
  Esgd[t]=model$mean$sgd
  SEsgd[t]=model$sd$sgd
  Ephx[t,1]=model$mean$phx[1,1]
  SEphx[t,1]=model$sd$phx[1,1]
  Ephx[t,2]=model$mean$phx[1,2]
  SEphx[t,2]=model$sd$phx[1,2]
  Ephx[t,3]=model$mean$phx[2,2]
  SEphx[t,3]=model$sd$phx[2,2]
}