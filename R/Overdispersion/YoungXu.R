# 13/4/19
# Reproduce example of table1 in Ýoung-Xu article
# to test overdispersion

workingDir<-"D:/Doctorat/Simplex/R/Overdispersion"
setwd(workingDir)
library(VGAM)
library(fitdistrplus)

dades <- read.table("Table1.YoungXu.csv", sep=",",stringsAsFactors = F,skip=6)
head(dades)
tail(dades)

#recodificar columnes
colnames(dades) <- c("study","ni","xi","perc")
dades$prop <- dades$xi/dades$ni
hist(dades$prop)
hist(dades$prop,breaks=20) #mmm expon...
hist(dades$xi)

fitdist(dades$prop, distr=dbeta, start=list(shape1=1,shape2=1),lower=c(0,0))
#         estimate Std. Error
# shape1  1.00000         NA
# shape2 30.72126         NA
fitdist(dades$xi, distr=dbetabinom.ab, start=list(size=13,shape1=1,shape2=1),lower=c(0,0))
fitdistr(dades$xi, densfun = dbetabinom.ab, start=list(shape1=1,shape2=1),size=13)

#       estimate Std. Error
# size   13.0000021         NA
# shape1  0.3529998         NA
# shape2  1.3022960         NA

#utilitzo el codi de l'exemple de la wiki betabinomial.R
M=dades$xi
F=dades$ni

m1=sum(M*F)/sum(F)
#4.013658
m2=sum((M^2*F))/sum(F)
#31.0563
n=max(M)

alpha <-  (n*m1-m2)/(n*(m2/m1-m1-1)+m1)
beta <-  ((n-m1)*(n-(m2/m1)))/(n*((m2/m1)-m1-1)+m1)
alpha
#0.5357241
beta
#1.199455

gamma <- 1/(1+alpha+beta)

gamma 0.3656068 #si és proper a 0 no hi ha overdispersion







a <- dades$xi
n <- max(a)

m1=sum(a)/length(a)
#6
m2=sum((a^2))/length(a)
#50
n=12

alpha <-  (n*m1-m2)/(n*((m2/m1)-m1-1)+m1)
beta <-  ((n-m1)*(n-(m2/m1)))/(n*((m2/m1)-m1-1)+m1)
alpha
#1
beta
#1
#NO!!
mu=alpha/(alpha+beta)
zita= 1/(alpha+beta) #+1???

mu
zita

#provar betabinomial!
source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params
#est params
est.all.params(dades$prop)
# s.mle.mu   s.mle.sig   s.zoip.mu  s.zoip.sig    b.mom.s1    b.mom.s2    b.mle.s1    b.mle.s2 
# NA          NA        0.03509619  5.67504707  1.07849900 29.51191132  1.36213966 37.17191925 
# binf.mle.s1 binf.mle.s2     n.mom.m    n.mom.sd     n.mle.m    n.mle.sd 
# NA          NA            0.03135368  0.03609890  0.03135368  0.03565595 
est.betabin.params(dades$prop)
# bb.mom.s1 bb.mom.s2     sim.n bb.mle.s1 bb.mle.s2    real.n 
# NaN       NaN        NA        NA        NA         0 
est.betabin.params(b=dades$prop,ni=dades$ni)

#MILLOR DISTRIBUCIÓ
#li passo xi, que és el que li passo a la fn
best.dist.betabin(dades$xi)
max(dades$xi)

#miro a veure quina dist s'ajusta més
best.dist(dades$prop)
# simplex.aic      beta.aic    normal.aic  simplex.ks.p     beta.ks.p   normal.ks.p 
#       NA            NA    -153.02185729            NA            NA    0.07369865 


#recordo com fer-ho
data("GasolineYield", package = "betareg")
hist(GasolineYield$yield,breaks=10)

## Table 1
gy <- betareg(yield ~ batch + temp, data = GasolineYield)
summary(gy)
