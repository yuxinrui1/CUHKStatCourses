
state.x77
states <- as.data.frame(state.x77[,c("Murder", "Population","Illiteracy", "Income", "Frost")])
cor(states)
library(car)
scatterplotMatrix(states)
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost,data=states)
summary(fit)
plot(fit)
qqPlot(fit)
durbinWatsonTest(fit)
crPlots(fit)
ncvTest(fit)
library(gvlma)
gvlma(fit)
PSr=((100-state.x77[,"Illiteracy"])+rnorm(50,0,0.04))/100
states1 <- as.data.frame(cbind(state.x77[,c("Murder", "Population","Illiteracy", "Income", "Frost")],PSr))
fit1 <- lm(Murder ~ Population + Illiteracy + Income + Frost + PSr,data=states1)
vif(fit)
vif(fit1)


USJudgeRatings
cor(USJudgeRatings[,-1])
cov(scale(USJudgeRatings[,-1]))
eigen(cor(USJudgeRatings[,-1]))
pc1=princomp(USJudgeRatings[,-1],cor=T)
summary(pc1)
screeplot(pc1,type="lines")
lines(seq(1,10,0.01),rep(1,901))
biplot(pc1)
predict(princomp(USJudgeRatings[,-1],cor=T))

library(psych)
pc1 <- principal(Harman23.cor$cov, nfactors=3,rotate="none")
pc1
plot(pc1$Structure[,c(1:2)],ylim=c(-1,1),xlim=c(0.6,0.9))
text(pc1$Structure[,c(1,2)],row.names(pc1$Structure),cex=0.8,pos=1)

factor.pa<-function(S, m){
  p<-nrow(S)
  diag_S<-diag(S) 
  sum_rank<-sum(diag_S)
  rowname<-paste("X", 1:p, sep="")
  colname<-paste("Factor", 1:m, sep="")
  A<-matrix(0, nrow=p, ncol=m, dimnames=list(rowname, colname))
  eig<-eigen(S)
  for (i in 1:m)
    A[,i]<-sqrt(eig$values[i])*eig$vectors[,i]
  h<-diag(A%*%t(A))
  
  rowname<-c("SS loadings", "Proportion Var", "Cumulative Var")
  B<-matrix(0, nrow=3, ncol=m, dimnames=list(rowname, colname))
  for (i in 1:m){
    B[1,i]<-sum(A[,i]^2)
    B[2,i]<-B[1,i]/sum_rank
    B[3,i]<-sum(B[1,1:i])/sum_rank
  }
  method<-c("Principal Component Method")
  list(method=method, loadings=A, 
       var=cbind(common=h, spcific=diag_S-h), B=B) 
}


factor.pa(Harman23.cor$cov,4)

principal(Harman23.cor$cov, nfactors=1,rotate="none")
principal(Harman23.cor$cov, nfactors=2,rotate="none")
principal(Harman23.cor$cov, nfactors=3,rotate="none")

library(GPArotation)
principal(Harman23.cor$cov, nfactors=3,rotate="varimax")
principal(Harman23.cor$cov, nfactors=3,rotate="promax")


factor.pa<-function(S, m){
  p<-nrow(S); diag_S<-diag(S); sum_rank<-sum(diag_S)
  rowname<-paste("X", 1:p, sep="")
  colname<-paste("Factor", 1:m, sep="")
  A<-matrix(0, nrow=p, ncol=m, 
            dimnames=list(rowname, colname))
  
  d = 1-diag(1/solve(S))
  kmax=5000; k<-1; h <- diag_S-d
  repeat{
    diag(S)<- h; h1<-h; eig<-eigen(S)
    for (i in 1:m)
      A[,i]<-sqrt(eig$values[i])*eig$vectors[,i]
    h<-diag(A %*% t(A))
    if ((sqrt(sum((h-h1)^2))<1e-4)|k==kmax) break
    k<-k+1
  }
  
  rowname<-c("SS loadings", "Proportion Var", "Cumulative Var")
  B<-matrix(0, nrow=3, ncol=m, 
            dimnames=list(rowname, colname))
  for (i in 1:m){
    B[1,i]<-sum(A[,i]^2)
    B[2,i]<-B[1,i]/sum_rank
    B[3,i]<-sum(B[1,1:i])/sum_rank
  }
  method<-c("Principal Axis Method")
  list(method=method, loadings=A, 
       var=cbind(common=h, spcific=diag_S-h), B=B, iterative=k) 
}   


ability.cov
R=cov2cor(ability.cov$cov)
pc=principal(R, nfactors=2,rotate="none",n.obs=112)
pa1=fa(R,nfactors=2,fm="pa",rotate = "none",n.obs=112)
pa2=factor.pa(R,2)
pa3=fa(R,nfactors=2,fm="pa",rotate = "varimax",n.obs=112)
pa4=fa(R,nfactors=2,fm="pa",rotate = "promax",n.obs=112)
fa.diagram(pa1,simple=F)
fa.diagram(pc,simple=F)
fa.diagram(pa3,simple=F)
fa.diagram(pa4,simple=F)



install.packages("lavaan", dependencies = TRUE)
library("lavaan")

WiscIV.cor <- lav_matrix_lower2full(c(1.000,
                                      0.72, 1.000,
                                      0.64, 0.63, 1.000,
                                      0.51, 0.48, 0.37, 1.000,
                                      0.37, 0.38, 0.38, 0.38, 1.000))

colnames(WiscIV.cor)<-rownames(WiscIV.cor)<-c("Info", "Sim", "Word", "Matrix", "Pict")
model <- '
# Measurement model
g=~ a*Info + b*Sim + c*Word + d*Matrix + e*Pict
'
WiscIV.fit <- sem(model,sample.cov= WiscIV.cor,sample.nobs =550)
summary(WiscIV.fit)
WiscIV.fit1 <- sem(model,sample.cov= WiscIV.cor,sample.nobs =550,std.lv=T)
summary(WiscIV.fit1,standardize=TRUE)

model2 <- '
g=~ NA*Info + b*Sim + c*Word + d*Matrix + e*Pict
g ~~ 1*g
'
WiscIV.fit2 <- sem(model2,sample.cov= WiscIV.cor,sample.nobs =550)
summary(WiscIV.fit2)

model3 <- '
# Measurement model
Verbal=~ a*Info + b*Sim + c*Word 
Fluid =~ d*Matrix + e*Pict
# error Variance and Covariance (psi)
Verbal ~~ Fluid
'
WiscIV.fit3 <- sem(model3,sample.cov= WiscIV.cor,sample.nobs =550)
summary(WiscIV.fit3)

model31 <- '
# Measurement model
Verbal=~ NA*Info + b*Sim + c*Word 
Fluid =~ NA*Matrix + e*Pict
# error Variance and Covariance (psi)
Verbal ~~ 1*Verbal
Fluid ~~ 1*Fluid
Verbal ~~ Fluid
'
WiscIV.fit31 <- sem(model31,sample.cov= WiscIV.cor,sample.nobs =550)
summary(WiscIV.fit31)

model32 <- '
# Measurement model
Verbal=~ NA*Info + b*Sim + c*Word 
Fluid =~ NA*Matrix + e*Pict
# error Variance and Covariance (psi)
Verbal ~~ 1*Verbal
Fluid ~~ 1*Fluid
Verbal ~~ Fluid
# error covariances
Word ~~ Matrix
Word ~~ Pict
'
WiscIV.fit32 <- sem(model32,sample.cov= WiscIV.cor,sample.nobs =550)
summary(WiscIV.fit32)

model322 <- '
# Measurement model
Verbal=~ NA*Info + b*Sim + c*Word 
Fluid =~ NA*Matrix + e*Pict
# error Variance and Covariance (psi)
Verbal ~~ 1*Verbal
Fluid ~~ 1*Fluid
Verbal ~~ Fluid
# error covariances
Info ~~ Matrix
Word ~~ Matrix
Word ~~ Pict
'
WiscIV.fit322 <- sem(model322,sample.cov= WiscIV.cor,sample.nobs =550)

model33 <- '
# Measurement model
Verbal=~ NA*Info + b*Sim + c*Word 
Fluid =~ NA*Matrix + e*Pict
# error Variance and Covariance (psi)
Verbal ~~ 1*Verbal
Fluid ~~ 1*Fluid
Verbal ~~ Fluid
# error covariances   [model is not identified]
Info ~~ Sim
Info ~~ Word
Word ~~ Pict
Info ~~ Pict
Matrix ~~ Pict
'
WiscIV.fit33 <- sem(model33,sample.cov= WiscIV.cor,sample.nobs =550) 

model34 <- '
# Measurement model
Verbal=~ NA*Info + b*Sim + c*Word 
Fluid =~ NA*Matrix + e*Pict + f*Word
# error Variance and Covariance (psi)
Verbal ~~ 1*Verbal
Fluid ~~ 1*Fluid
Verbal ~~ Fluid
'
WiscIV.fit34 <- sem(model34,sample.cov= WiscIV.cor,sample.nobs =550) 
summary(WiscIV.fit34)

fitMeasures(WiscIV.fit34)
anova(WiscIV.fit31,WiscIV.fit32)
modificationIndices(WiscIV.fit31,sort.=TRUE)
anova(WiscIV.fit31,WiscIV.fit34)


police.cor <- lav_matrix_upper2full(
  c(1.00, 0.50, 0.41, 0.33, 0.28, 0.30, -0.24,-0.23,-0.20,
          1.00, 0.35, 0.29, 0.26, 0.27, -0.19,-0.19,-0.18,
                1.00, 0.30, 0.27, 0.29, -0.17,-0.16,-0.14,
                      1.00, 0.52, 0.48, -0.13,-0.11,-0.15,
                            1.00, 0.44, -0.11,-0.09,-0.10,
                                  1.00, -0.15,-0.13,-0.13,
                                         1.00, 0.58, 0.47,
                                               1.00, 0.42,
                                                      1.00))
colnames(police.cor)<-rownames(police.cor)<-c("PS" , "RE" , "RT" , "HO" , "CO" , "ET" , "BU" , "VA" , "RO")


modelp1 <- '
# Measurement model
F1 =~ NA*PS + b*RE + c*RT
F2 =~ NA*HO + e*CO + f*ET
F3 =~ NA*BU + h*VA + i*RO
# error Variance and Covariance (psi)
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3
F1 ~~ j*F2
F1 ~~ k*F3
F2 ~~ l*F3
'
p1.fit=sem(modelp1,sample.cov= police.cor,sample.nobs =11000) 
summary(p1.fit)

modelp1.1 <- '
# Measurement model
F1 =~ NA*PS + b*RE + c*RT + d*HO + e*CO + f*ET
F3 =~ NA*BU + h*VA + i*RO
# error Variance and Covariance (psi)
F1 ~~ 1*F1
F3 ~~ 1*F3
F1 ~~ k*F3
'
p1.1.fit=sem(modelp1.1,sample.cov= police.cor,sample.nobs =11000) 
summary(p1.1.fit)
anova(p1.fit,p1.1.fit)

residuals(p1.fit)
residuals(p1.1.fit)

fitmeasures(p1.fit)
fitmeasures(p1.1.fit)

summary(p1.1.fit)
(parameterestimates(p1.fit)$est[1:9])^2

modificationindices(p1.fit,sort=T)

modelp2 <- '
# Measurement model
F1 =~ a*PS + b*RE + c*RT
F2 =~ d*HO + e*CO + f*ET + RT
F3 =~ g*BU + h*VA + i*RO
# error Variance and Covariance (psi)
F1 ~~ j*F2
F1 ~~ k*F3
F2 ~~ l*F3
'
p2.fit=sem(modelp2,sample.cov= police.cor,sample.nobs =11000,std.lv=T) 

summary(p2.fit)
anova(p1.fit,p2.fit)
modificationindices(p2.fit,sort=T)


modelp2.1 <- '
# Measurement model
F1 =~ a*PS + b*RE + c*RT
F2 =~ d*HO + e*CO + f*ET + RT
F3 =~ g*BU + h*VA + i*RO
# error Variance and Covariance (psi)
F1 ~~ j*F2
F1 ~~ k*F3
F2 ~~ l*F3
#covariance
HO ~~ RO
'
p2.1.fit=sem(modelp2.1,sample.cov= police.cor,sample.nobs =11000,std.lv=T) 
summary(p2.1.fit)
anova(p2.1.fit,p2.fit)


modelp3 <- '
# Measurement model
F1 =~ a*PS + b*RE + c*RT
F2 =~ d*HO + e*CO + f*ET 
F3 =~ g*BU + h*VA + i*RO
# error Variance and Covariance (psi)
F1 ~~ j*F2
#Structural model
F3 ~ k*F1+l*F2
# error Variance
F3 ~~ m*F3
' 
p3.fit=sem(modelp3,sample.cov= police.cor,sample.nobs =11000) 
summary(p3.fit)


modelp3.1 <- '
# Measurement model
F1 =~ a*PS + b*RE + c*RT
F2 =~ d*HO + e*CO + f*ET + RT
F3 =~ g*BU + h*VA + i*RO
#Structural model
F1 ~ j*F2
F3 ~ k*F1+l*F2
#define a new parameter 
ind := j*k
# error Variance
F3 ~~ m*F3
F1 ~~ n*F1
' 
p3.1.fit=sem(modelp3.1,sample.cov= police.cor,sample.nobs =11000) 
summary(p3.1.fit)


modelp3.2 <- '
# Measurement model
F1 =~ a*PS + b*RE + c*RT
F2 =~ d*HO + e*CO + f*ET + RT
F3 =~ g*BU + h*VA + i*RO
#Structural model
F1 ~ j*F2
F3 ~ k*F1
# error Variance
F3 ~~ m*F3
F1 ~~ n*F1
' 
p3.2.fit=sem(modelp3.2,sample.cov= police.cor,sample.nobs =11000) 
anova(p3.1.fit,p3.2.fit)
fitmeasures(p3.1.fit)
fitmeasures(p3.2.fit)





