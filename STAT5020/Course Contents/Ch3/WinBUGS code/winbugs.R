##################winbugs generater data#################
library(mvtnorm)   #Load mvtnorm package

N=300                         #Sample size
XI=matrix(NA, nrow=N, ncol=2) #Explanatory latent variables
Eta=numeric(N)                #Outcome latent variables
Y=matrix(NA, nrow=N, ncol=9) #Observed variables

#The covariance matrix of xi
phi=matrix(c(1, 0.3, 0.3, 1), nrow=2)

eps=numeric(10)

  #Generate Data
  for (i in 1:N) {
    XI[i,]=rmvnorm(1, c(0,0), phi)
    delta=rnorm(1, 0, sqrt(0.5))
    Eta[i]=0.6*XI[i,1]+0.4*XI[i,2]+delta
    
    eps[1:3]=rnorm(3, -1, sqrt(0.4))
    eps[4:6]=rnorm(3, 0, sqrt(0.4))
    eps[7:9]=rnorm(3, 1, sqrt(0.4))
    
    Y[i,1]=Eta[i]+eps[1]
    Y[i,2]=0.7*Eta[i]+eps[2]
    Y[i,3]=0.5*Eta[i]+eps[3]
    Y[i,4]=XI[i,1]+eps[4]
    Y[i,5]=0.8*XI[i,1]+eps[5]
    Y[i,6]=0.8*XI[i,1]+eps[6]
    Y[i,7]=XI[i,2]+eps[7]
    Y[i,8]=0.8*XI[i,2]+eps[8]
    Y[i,9]=-0.8*XI[i,2]+eps[9]
  }
  
options(scipen = 200)
write(t(Y),file="Y.txt",sep=",",ncol=N*9)



