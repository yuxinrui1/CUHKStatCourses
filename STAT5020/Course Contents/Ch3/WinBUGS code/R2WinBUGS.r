library(mvtnorm)   #Load mvtnorm package
library(R2WinBUGS) #Load R2WinBUGS package

N=500                         #Sample size
BZ=numeric(N)                 #Fixed covariate in structural equation
XI=matrix(NA, nrow=N, ncol=2) #Explanatory latent variables
Eta=numeric(N)                #Outcome latent variables
Y=matrix(NA, nrow=N, ncol=10)  #Observed variables

#The covariance matrix of xi
phi=matrix(c(1, 0.3, 0.3, 1), nrow=2)

#Estimates and standard error estimates
Eu=matrix(NA, nrow=100, ncol=10);    SEu=matrix(NA, nrow=100, ncol=10)
Elam=matrix(NA, nrow=100, ncol=7);   SElam=matrix(NA, nrow=100, ncol=7)
Eb=numeric(100);                     SEb=numeric(100)
Egam=matrix(NA, nrow=100, ncol=5);   SEgam=matrix(NA, nrow=100, ncol=5)
Esgm=matrix(NA, nrow=100, ncol=10);  SEsgm=matrix(NA, nrow=100, ncol=10)
Esgd=numeric(100);                   SEsgd=numeric(100)
Ephx=matrix(NA, nrow=100, ncol=3);   SEphx=matrix(NA, nrow=100, ncol=3)

R=matrix(c(1.0, 0.3, 0.3, 1.0), nrow=2)

parameters=c("u", "lam", "b", "gam", "sgm", "sgd", "phx")

init1=list(u=rep(0,10), lam=rep(0,7), b=0, gam=rep(0,5), psi=rep(1,10),
           psd=1, phi=matrix(c(1, 0, 0, 1), nrow=2))

init2=list(u=rep(1,10), lam=rep(1,7), b=1, gam=rep(1,5), psi=rep(2,10),
           psd=2, phi=matrix(c(2, 0, 0, 2), nrow=2))

inits=list(init1, init2)

eps=numeric(10)

for (t in 1:100) {
    #Generate Data
    for (i in 1:N) {
        BZ[i]=rt(1, 5)

        XI[i,]=rmvnorm(1, c(0,0), phi)

        delta=rnorm(1, 0, sqrt(0.36))
        Eta[i]=0.5*BZ[i]+0.4*XI[i,1]+0.4*XI[i,2]+0.3*XI[i,1]*XI[i,2]
               +0.2*XI[i,1]*XI[i,1]+0.5*XI[i,2]*XI[i,2]+delta

        eps[1:3]=rnorm(3, 0, sqrt(0.3))
        eps[4:7]=rnorm(4, 0, sqrt(0.5))
        eps[8:10]=rnorm(3, 0, sqrt(0.4))
        Y[i,1]=Eta[i]+eps[1]
        Y[i,2]=0.9*Eta[i]+eps[2]
        Y[i,3]=0.7*Eta[i]+eps[3]
        Y[i,4]=XI[i,1]+eps[4]
        Y[i,5]=0.9*XI[i,1]+eps[5]
        Y[i,6]=0.7*XI[i,1]+eps[6]
        Y[i,7]=0.5*XI[i,1]+eps[7]
        Y[i,8]=XI[i,2]+eps[8]
        Y[i,9]=0.9*XI[i,2]+eps[9]
        Y[i,10]=0.7*XI[i,2]+eps[10]
    }

    #Run WinBUGS
    data=list(N=500, zero=c(0,0), z=BZ, R=R, y=Y)

    model<-bugs(data,inits,parameters,
                model.file="D:/model.txt",
                n.chains=2,n.iter=10000,n.burnin=4000,n.thin=1,
                bugs.directory="D:/WinBUGS14")

    #Save Estimates
    Eu[t,]=model$mean$u;               SEu[t,]=model$sd$u
    Elam[t,]=model$mean$lam;           SElam[t,]=model$sd$lam
    Eb[t]=model$mean$b;                SEb[t]=model$sd$b
    Egam[t,]=model$mean$gam;           SEgam[t,]=model$sd$gam
    Esgm[t,]=model$mean$sgm;           SEsgm[t,]=model$sd$sgm
    Esgd[t]=model$mean$sgd;            SEsgd[t]=model$sd$sgd
    Ephx[t,1]=model$mean$phx[1,1];     SEphx[t,1]=model$sd$phx[1,1]
    Ephx[t,2]=model$mean$phx[1,2];     SEphx[t,2]=model$sd$phx[1,2]
    Ephx[t,3]=model$mean$phx[2,2];     SEphx[t,3]=model$sd$phx[2,2]
}
