install.packages("lavaan", dependencies = TRUE)
library("lavaan")

x<-lav_matrix_lower2full(
  c(1.000,
     0.20, 1.000,
     0.24, 0.30, 1.000,
     0.70, 0.80, 0.30, 1.000))
colnames(x) <- rownames(x) <- c("X1", "X2", "X3","Y")

regression.model<-'
# structural model for Y
Y ~ a*X1 + b*X2 + c*X3
# label the residual variance of Y
Y ~~ z*Y
'
regression.fit <- sem(regression.model, sample.cov = x, sample.nobs = 1000)
summary(regression.fit, rsquare = TRUE)


regression.model<-'
# structural model for Y
Y ~ a*X1 + b*X2 + c*X3

# error Variance and Covariance (psi)
Y ~~ Y
'

x<-lav_matrix_lower2full(
  c(0.594,
    0.483, 0.754,
    3.993, 3.626, 47.457,
    0.426, 1.757, 4.100, 10.267,
    0.500, 0.722, 6.394, 0.525, 2.675))
colnames(x) <- rownames(x) <- c("GPA_R", "GPA_E", "SAT","IQ","Motiv")

college.model <- "
# structural model for GPA
GPA_R ~ a*SAT
GPA_E ~ b*SAT + c*IQ + d*Motiv
# error Variance and Covariance (psi)
GPA_R ~~ e*GPA_R
GPA_E ~~ f*GPA_E
GPA_R ~~ g*GPA_E
"
# Fit Model 1 to the data
fit1 <- sem(college.model, sample.cov=x, sample.nobs=150)
summary(fit1)

xx=cov2cor(x)
fit2 <- sem(college.model, sample.cov=xx, sample.nobs=150)
summary(fit2)


x <- lav_matrix_lower2full(c(648.07, 30.05, 8.64, 140.18, 25.57, 233.21))
colnames(x) <- rownames(x) <- c("salary", "school", "iq")

ind.model <- '
# structural model
salary ~ a*school + c*iq
iq ~ b*school
#define a new parameter 
ind := b*c
# error Variance and Covariance (psi)
salary ~~ d* salary
iq ~~ e*iq
'
fit3 <- sem(ind.model, sample.cov=x, sample.nobs=1000)
summary(fit3,rsquare=T)
xx=cov2cor(x)
fit4 <- sem(ind.model, sample.cov=xx, sample.nobs=1000)
summary(fit4,rsquare=T)

wo.model <- '
# structural model
salary ~  c*iq
iq ~ b*school
#define a new parameter 
ind := b*c
# error Variance and Covariance (psi)
salary ~~ d* salary
iq ~~ e*iq
'
fit5 <- sem(wo.model, sample.cov=xx, sample.nobs=1000)
summary(fit5,rsquare=T)


wo1.model <- '
# structural model
salary ~ a*school
# error Variance and Covariance (psi)
salary ~~ d* salary
'
fit6 <- sem(wo1.model, sample.cov=xx, sample.nobs=1000)
summary(fit6,rsquare=T)

privSchool.cor <- c(1, 
                    0.178, 1,
                    0.230, 0.327, 1,
                    0.106, 0.245, 0.183, 1,
                    0.195, 0.356, 0.721, 0.178, 1)
privSchool.cor <- lav_matrix_lower2full(privSchool.cor)
colnames(privSchool.cor) <- rownames(privSchool.cor) <- c("Race", "SES", "CogAbil", "SchTyp", "AcadAch")
full.model <- '
# structural model
AcadAch ~ j*SchTyp + g*SES + d*Race + i*CogAbil
SchTyp ~ f*SES + b*Race + h*CogAbil
CogAbil ~ e*SES + c*Race
SES ~ a*Race
# error Variance and Covariance (psi)
AcadAch ~~  AcadAch
SchTyp ~~  SchTyp
CogAbil ~~  CogAbil
SES ~~ SES
'
full.fit <- sem(full.model, sample.cov=privSchool.cor, sample.nobs=18058)
summary(full.fit)



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

model321 <- '
# Measurement model
Verbal=~ NA*Info + b*Sim + c*Word 
Fluid =~ NA*Matrix + e*Pict
# error Variance and Covariance (psi)
Verbal ~~ 1*Verbal
Fluid ~~ 1*Fluid
Verbal ~~ Fluid
# error covariances
Word ~~ Matrix
'
WiscIV.fit321 <- sem(model321,sample.cov= WiscIV.cor,sample.nobs =550)
summary(WiscIV.fit321)
anova(WiscIV.fit31,WiscIV.fit321)

fitMeasures(WiscIV.fit34)
residuals(WiscIV.fit34)$cov
residuals(WiscIV.fit34,type="cor")$cov
sqrt(sum(residuals(WiscIV.fit34)$cov^2)/30)

fitMeasures(WiscIV.fit34,fit.measures = "GFI")
1-12*(12+1)/(2*15)*(1-fitMeasures(WiscIV.fit34,fit.measures = "GFI"))
fitMeasures(WiscIV.fit34,fit.measures = "AGFI")


modelbaseline <- '
# Measurement model
Info ~~ Info
Sim ~~ Sim
Word ~~ Word 
Matrix ~~ Matrix
Pict ~~ Pict
'
WiscIVbaseline <- sem(modelbaseline,sample.cov= WiscIV.cor,sample.nobs =550) 
fitMeasures(WiscIVbaseline)

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


model <- '
  # Measurement Models
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # Structural Models
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
  # error Variance and Covariance
    dem60 ~~ dem60
    dem65 ~~ dem65
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
fit <- sem(model, data = PoliticalDemocracy)
summary(fit)
fit <- sem(model, sample.cov=cov(PoliticalDemocracy),sample.nobs=75)
summary(fit)
fit <- sem(model, data = scale(PoliticalDemocracy,scale=F))
summary(fit,standard=T)

cov(scale(PoliticalDemocracy))
cor(PoliticalDemocracy)
fitsd <- sem(model, data = scale(PoliticalDemocracy))
summary(fitsd,standard=T)
fitsd <- sem(model, sample.cov=cor(PoliticalDemocracy),sample.nobs=75)


modelint <- '
  # Measurement Models
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # Structural Models
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
  # error Variance and Covariance
    dem60 ~~ dem60
    dem65 ~~ dem65
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
    #intercept
    x1 ~ 1
    x2 ~ 1
    x3 ~ 1
    y1 ~ 1
    y2 ~ 1
    y3 ~ 1
    y4 ~ 1
    y5 ~ 1
    y6 ~ 1
    y7 ~ 1
    y8 ~ 1 
'
fitint <- sem(modelint, data = PoliticalDemocracy)
summary(fitint)
fitint1 <- sem(model, data = PoliticalDemocracy,meanstructure=T)
summary(sem(modelint, sample.cov=cor(PoliticalDemocracy),sample.nobs=75))



sma_male.cov <- lav_matrix_lower2full(
  c(63,
    70, 110,
    41, 52, 60,
    30, 37, 36, 32))

colnames(sma_male.cov) <- rownames(sma_male.cov) <- c("sales1", "sales2", "admin1","admin2")
sma_female.cov <- lav_matrix_lower2full(
  c(67,
    72, 107,
    40, 55, 63,
    28, 38, 39, 35))
colnames(sma_female.cov) <- rownames(sma_female.cov) <- c("sales1", "sales2", "admin1","admin2")

model1 <- "
# measurement model
S_M =~ NA*sales1 + sales2
ADM =~ NA*admin1 + admin2
# factor Variances and covariance
S_M ~~ 1*S_M
ADM ~~ 1*ADM
S_M ~~ ADM
"
fit1 <-sem(model1, sample.cov=list(Group1=(sma_male.cov), Group2=(sma_female.cov)), sample.nobs=c(265, 300))
summary(fit1)

summary(sem(model1, sample.cov=(sma_male.cov), sample.nobs=265))
summary(sem(model1, sample.cov=(sma_female.cov), sample.nobs=300))
summary(sem(model1, sample.cov=list(Group1=cov2cor(sma_male.cov), Group2=cov2cor(sma_female.cov)), sample.nobs=c(265, 300)))

model1 <- "
# measurement model
S_M =~ sales1 + sales2
ADM =~ admin1 + admin2
# factor Variances and covariance
S_M ~~ S_M
ADM ~~ ADM
S_M ~~ ADM
"
model2 <- "
# measurement model
S_M =~ c(eq1,eq1)*sales1 + c(eq2,eq2)*sales2
ADM =~ c(eq3,eq3)*admin1 + c(eq4,eq4)*admin2
# factor Variances and covariance
S_M ~~ c(1,NA)*S_M
ADM ~~ c(1,NA)*ADM
S_M ~~ ADM
"
fit2 <-sem(model2, sample.cov=list(Group1=(sma_male.cov), Group2=(sma_female.cov)), sample.nobs=c(265, 300),std.lv=T)
summary(fit2,standardize=T)

model2 <- "
# measurement model
S_M =~ c(eq1,eq1)*sales1 + c(eq2,eq2)*sales2
ADM =~ c(eq3,eq3)*admin1 + c(eq4,eq4)*admin2
# factor Variances and covariance
S_M ~~ S_M
ADM ~~ ADM
S_M ~~ ADM
"
fit2 <-sem(model2, sample.cov=list(Group1=(sma_male.cov), Group2=(sma_female.cov)), sample.nobs=c(265, 300))
summary(fit2,standardize=T)


model3 <- "
# measurement model
S_M =~ c(eq1,eq1)*sales1 + c(eq2,eq2)*sales2
ADM =~ c(eq3,eq3)*admin1 + c(eq4,eq4)*admin2
# factor Variances and covariance
S_M ~~ c(1,1)*S_M
ADM ~~ c(1,1)*ADM
S_M ~~ ADM
"
fit3 <-sem(model3, sample.cov=list(Group1=(sma_male.cov), Group2=(sma_female.cov)), sample.nobs=c(265, 300),std.lv=T)
summary(fit3,standardize=T)

model3 <- "
# measurement model
S_M =~ c(eq1,eq1)*sales1 + c(eq2,eq2)*sales2
ADM =~ c(eq3,eq3)*admin1 + c(eq4,eq4)*admin2
# factor Variances and covariance
S_M ~~ c(eq11,eq11)*S_M
ADM ~~ c(eq22,eq22)*ADM
S_M ~~ ADM
"
fit3 <-sem(model3, sample.cov=list(Group1=(sma_male.cov), Group2=(sma_female.cov)), sample.nobs=c(265, 300))
summary(fit3,standardize=T)


model4 <- "
# measurement model
S_M =~ c(eq1,eq1)*sales1 + c(eq2,eq2)*sales2
ADM =~ c(eq3,eq3)*admin1 + c(eq4,eq4)*admin2
# factor Variances and covariance
S_M ~~ c(1,1)*S_M
ADM ~~ c(1,1)*ADM
S_M ~~ c(eq5,eq5)*ADM
"
fit4 <-sem(model4, sample.cov=list(Group1=(sma_male.cov), Group2=(sma_female.cov)), sample.nobs=c(265, 300),std.lv=T)
summary(fit4,standardize=T)

model4 <- "
# measurement model
S_M =~ c(eq1,eq1)*sales1 + c(eq2,eq2)*sales2
ADM =~ c(eq3,eq3)*admin1 + c(eq4,eq4)*admin2
# factor Variances and covariance
S_M ~~ c(eq11,eq11)*S_M
ADM ~~ c(eq22,eq22)*ADM
S_M ~~ c(eq5,eq5)*ADM
"
fit4 <-sem(model4, sample.cov=list(Group1=(sma_male.cov), Group2=(sma_female.cov)), sample.nobs=c(265, 300))
summary(fit4,standardize=T)


model5 <- "
# measurement model
S_M =~ c(eq1,eq1)*sales1 + c(eq2,eq2)*sales2
ADM =~ c(eq3,eq3)*admin1 + c(eq4,eq4)*admin2
# factor Variances and covariance
S_M ~~ c(1,1)*S_M
ADM ~~ c(1,1)*ADM
S_M ~~ c(eq5,eq5)*ADM
# error Variances
sales1 ~~ c(eq6,eq6)*sales1
sales2 ~~ c(eq7,eq7)*sales2
admin1 ~~ c(eq8,eq8)*admin1
admin2 ~~ c(eq9,eq9)*admin2
"
fit5 <-sem(model5, sample.cov=list(Group1=(sma_male.cov), Group2=(sma_female.cov)), sample.nobs=c(265, 300),std.lv=T)
summary(fit5,standardize=T)

#Model Comparisons
lavTestLRT(fit1, fit2, fit3, fit4, fit5)



?HolzingerSwineford1939
HSmodel <- "
# measurement model
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9 
# factor Variances and covariance
visual ~~ visual
textual ~~ textual
speed ~~ speed
visual ~~ textual
visual ~~ speed
textual ~~ speed
"
# configural invariance
fit1 <- sem(HSmodel, data = HolzingerSwineford1939, group = "school")
summary(fit1)


# weak invariance
fit2 <- sem(HSmodel, data = HolzingerSwineford1939, group = "school",
            group.equal = "loadings")
summary(fit2)


# strong invariance
fit3 <- sem(HSmodel, data = HolzingerSwineford1939, group = "school",
            group.equal = c("loadings","intercepts"))
summary(fit3)

fit3.1 <- sem(HSmodel, data = HolzingerSwineford1939, group = "school",
              group.equal = c("loadings","intercepts"),group.partial=c("visual=~x2","x7~1"))
summary(fit3.1)

# model comparison tests
lavTestLRT(fit1, fit2, fit3)

fit4 <- sem(HS.model, data = HolzingerSwineford1939, group = "school",
            group.equal = c("intercepts","means","loadings","lv.variances","lv.covariances","residuals"))
summary(fit4)
fit41 <-  cfa(HS.model, data = HolzingerSwineford1939,meanstructure = TRUE)
summary(fit41)



