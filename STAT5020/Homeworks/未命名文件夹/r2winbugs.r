library(mvtnorm)   #Load mvtnorm package
library(R2WinBUGS) #Load R2WinBUGS package

rep<-10 #replications

N=500       

d = numeric(N)
z = matrix(NA, nrow=N, ncol=9)  
y = matrix(NA, nrow=N, ncol=9)  
eta = numeric(N)   
xi = matrix(NA, nrow=N, ncol=2) 

#Estimates and standard error estimates
Eu=matrix(NA, nrow=rep, ncol=9)  
SEu=matrix(NA, nrow=rep, ncol=9)
Elam=matrix(NA, nrow=rep, ncol=6)  
SElam=matrix(NA, nrow=rep, ncol=6)
Eb=numeric(rep)                     
SEb=numeric(rep)
Egam=matrix(NA, nrow=rep, ncol=2)   
SEgam=matrix(NA, nrow=rep, ncol=2)
Esgm=matrix(NA, nrow=rep, ncol=9)  
SEsgm=matrix(NA, nrow=rep, ncol=9)
Esgd=numeric(rep)                  
SEsgd=numeric(rep)
Ephx=matrix(NA, nrow=rep, ncol=3)  
SEphx=matrix(NA, nrow=rep, ncol=3)

parameters=c("uby", "lam", "ubeta", "gam", "sgm", "sgd", "phx")

init1=list(uby=rep(0.,9), lam=rep(0.5,6), ubeta=0.3, gam=rep(0.3,2), psi=rep(1.,4), # initial value for parameters of interest
           psd=1., phi=matrix(c(1., 0., 0., 1.), nrow=2), xi=matrix(data=rep(0.0,1000),ncol=2))

init2=list(uby=rep(0.5,9), lam=rep(0.5,6), ubeta=2, gam=rep(0.5,2), psi=rep(3,4),
           psd=2, phi=matrix(c(1.5, 0, 0, 1.5), nrow=2), xi=matrix(data=rep(0.0,1000),ncol=2))


# init1<-list(uby=rep(1.0,9),ubeta=0.3,lam=c(0.8,0.8,0.8,0.8,0.8,0.8),
#             gam=c(0.2,0.5),psd=3.33,phi=matrix(data=c(1.0989,-0.3267,-0.3267,
#                                                           1.0989),ncol=2,byrow=TRUE), xi=matrix(data=rep(0.0,1000),ncol=2),psi=c(3.33,3.33,2.5,2.5))
# init2<-list(uby=rep(1.0,9),ubeta=0.4,lam=rep(1.0,6),gam=c(0.1,0.4),psd=3.0,phi=matrix(data=c(1.0,0.0,0.0,1.0),ncol=2,byrow=TRUE),
#             xi=matrix(data=rep(0.0,1000),ncol=2),psi=c(3,3,2.5,2.5))
inits=list(init1, init2)
R=matrix(c(7.0, 2.1, 2.1, 7), nrow=2)

phi<-matrix(data=c(1.0,0.2,0.2,0.81),ncol=2)
# what we observe is denoted as z, the underlying is denoted as y
for (t in 1:10) {
    #Generate Data
    for (i in 1:N) {    
        # structural equation
        xi[i,]=rmvnorm(1, c(0,0), phi)
        d[i]=rbinom(1,1, 0.7)

        delta=rnorm(1, sd = 0.6)
        eta[i]=0.3*d[i]+0.2*xi[i,1]+0.5*xi[i,2] + delta

        #  measurement equation
        
        y[i,1]=eta[i] + rnorm(1, sd = 1)
        y[i,2]=0.8*eta[i]+ rnorm(1, sd = 1)
        y[i,3]=0.8*eta[i]+ rnorm(1, sd = 1)
        y[i,4]=xi[i,1]+ rnorm(1, sd = sqrt(0.3))
        y[i,5]=0.7*xi[i,1]+ rnorm(1, sd = sqrt(0.3))
        y[i,6]=0.9*xi[i,1]+ rnorm(1, sd = sqrt(0.4))
        y[i,7]=0.7*xi[i,1]+ rnorm(1, sd = sqrt(0.4))
        y[i,8]=xi[i,2]
        y[i,9]=0.8*xi[i,2]

        for (j in 1:3) { if (y[i,j]>0) z[i,j]<-1 else z[i,j]<-0 }
        for (j in 4:7) { z[i,j] <- y[i,j] }
        for (j in 8:9) {
            pij <- exp(y[i,j])/(1+exp(y[i,j]))
            z[i,j]<-rbinom(1,1,pij)
        }
    }
    #Run WinBUGS
    data=list(N=500, P=9, d=d, R=R, z=z, low=c(-2000,0),high=c(0,2000))

    model<-bugs(data,inits,parameters,
                model.file="C:/Users/s1155165569/OneDrive - The Chinese University of Hong Kong/courses/Structural Equation Models/hw2/model.txt",
                n.chains=2,n.iter=10000,n.burnin=3000,n.thin=1,
                bugs.directory="D:/stat/winbugs14_full_patched/WinBUGS14", debug = T)

    #Save Estimates
    Eu[t,]=model$mean$uby            
    SEu[t,]=model$sd$uby
    Elam[t,]=model$mean$lam        
    SElam[t,]=model$sd$lam
    Eb[t]=model$mean$b;                
    SEb[t]=model$sd$b
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