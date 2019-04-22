#proves amb els diferents paquets que permeten fer simplex (hi ha tb ZOIP.simplexregression.R) i altres mètodes
#provo: 
#1. dmp de minfi (f fisher)
#2. simplexreg de simplexreg
#3. beta reg inflated amb paquet gamlss
#10/9/18 sembla que la nova versió (5.0.6) del paquet gamlss.dist permet treballar amb simplex, la installo del zip
#error en intentar instal·lar-lo
#faig una quantile normalization
#he d'agafar la pglobal del model i no de la mu en simplexreg i ZOIP? check

projectDir <- "D:/Doctorat/1s year/Projecte defence/Presentation"
resultsDir <- "D:/Doctorat/Simplex/R"
rDir <- "D:/Doctorat/R"

###############################################
#agafo exemple minfiData
library(minfi)
library(minfiData)

data(MsetEx)
RSet.betas <- getBeta(MsetEx) 

hist(RSet.betas,freq=F) #?s el mateix hist que per les betas de MSet
range(RSet.betas,na.rm=TRUE) #[0,1] i tb hi ha nas

GSet <- preprocessQuantile(MsetEx)
GSet.b <- getBeta(GSet) 
hist(GSet.b,freq=F) #?s el mateix hist que per les betas de MSet
range(GSet.b,na.rm=TRUE) #(0,1) i ja no hi ha NAs, la manera és sempre fer un quantile per
#a no tenir dades inflades

########## 1. dmpFinder: per la presentació de la tesi 1st year #############
#dmpFinder, del paquet minfi, basat en reg lineal i F-test per a categ?riques
status <- pData(MsetEx)$status

# time1.0<-Sys.time()
# dmp <- dmpFinder(RSet.betas, pheno = status , type = "categorical")
# time2.0<-Sys.time()
# time.dmp<-time2.0-time1.0 #2.49mins!! 1,5 amb el pepino! si mirem nom?s un registre ?s autom?tic
# 
# res.dmp<-dmp[dmp$pval<0.05,]
# head(res.dmp)
# dim(res.dmp) #83406 #molt m?s laxe tot i que els 2 primers s?n exactament els mateixos
# 
# res.dmp.qval<-dmp[dmp$qval<0.05,] #14120 tamb? en s?n molts!
# dim(res.dmp.qval)

#quantile normalized
time1.0<-Sys.time()
  dmp <- dmpFinder(GSet.b, pheno = status , type = "categorical")
time2.0<-Sys.time()
time.dmp<-time2.0-time1.0 #2.49mins!! 1,5 amb el pepino! si mirem nom?s un registre ?s autom?tic

res.dmp<-dmp[dmp$pval<0.05,]
head(res.dmp)
dim(res.dmp) #92487

res.dmp.qval<-dmp[dmp$qval<0.05,] #17067 tamb? en s?n molts!
dim(res.dmp.qval)

save(res.dmp,res.dmp.qval,file=file.path(resultsDir,"dmp.res.RData"))


########## 1. simplexreg:  #############
library(simplexreg)
# In the simplexreg package, the function dsimplex gives the density function, psimplex provides
# the distribution function, qsimplex calculates the quantile function and rsimplex gives
# random numbers generated from the simplex distribution. They thus possess the same forms
# as other distribution functions in R:
  
#Pel que he entès de l'article del simplexreg i dels exemples que presenta, el fet de tenir els dos params
#mu i sigma, permet fer dos tipus d'assumpció sobre les dades: 
# model amb dispersió homogènia i model amb dispersió heterogènia
# si el model amb dispersió heterogènia no presenta coeficients significatius i AIC NO és molt
# més baixa que el model homogeni (sense modelar la dispersió) llavors podem dir que només amb el
# model amb dispersió hogeènia ja estem ajustant bé les dades però si no hauríem de tenir en compte 
# també el modelat de la dispersió. simplex reg només permet modelar la dispersió amb log
#tests del paquet
data("sdac", package = "simplexreg")
sim.glm1 <- simplexreg(rcd~ageadj+chemo, link = "logit", 
                       data = sdac)
summary(sim.glm1)
#     Coefficients (mean model with logit link):
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 1.100226   0.140683   7.821 5.26e-15 ***
#   ageadj      0.013575   0.006519   2.082   0.0373 *  
#   chemo       0.266092   0.124991   2.129   0.0333 *  

#per a estimar sigma sempre fa servir la link log segons l'ajut i ~1 default, es pot canviar fent
#Del help de simplexreg: Four types of function are available linking the regressors to the mean. However, for dispersion, the link function is restricted to logarithm function. 
# When modeling dispersion, the regressor modelling the dispersion parameter should be specified in a formula form of type y ~ x1 + x2 | z1 + z2 where z1 and z2 are linked to the dispersion parameter σ^2.

sim.glm11 <- simplexreg(rcd~ageadj+chemo|age, link = "logit", 
                       data = sdac)
summary(sim.glm11)
# Coefficients (mean model with logit link):
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 1.115550   0.141396   7.890 3.03e-15 ***
#   ageadj      0.013013   0.006452   2.017   0.0437 *  
#   chemo       0.251921   0.121807   2.068   0.0386 *  
#   
#   Coefficients (dispersion model with log link):
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  2.60750    0.36687   7.107 1.18e-12 ***
#   age         -0.01500    0.00688  -2.181   0.0292 *  

AIC(sim.glm1,sim.glm11) #millor el segon

hist(rcd) #una mica normal skewed cap a la dreta
table(sdac$chemo)
# 
# 0   1 
# 130 109 
hist(sdac$ageadj)
#adjusted age variable. age < 40 is set as the baseline age and other ages are adjusted by subtracting by 40

# cut<-10000
# betas2test<-RSet.betas[1:cut,]
cond<-as.factor(pData(MsetEx)$status)
# 
# #Simplex per totes les mostres
# #falla per valors 0 i 1 prenem 0.00001 com a epsilon i substituim
# table(RSet.betas==0) #10
# table(RSet.betas==1) #8, tampoc ?s un drama
# e<-0.00001
# 
# RSet.betas1<-RSet.betas
# for (i in 1:nrow(RSet.betas1)){
#   for (j in 1:ncol(RSet.betas1)){
#   if (!is.na(RSet.betas1[i,j])){
#   if (RSet.betas1[i,j]==0) RSet.betas1[i,j]=e else if (RSet.betas1[i,j]==1) RSet.betas1[i,j]=1-e
#    }
#   }
# }

# 
# time1<-Sys.time()
# lng<-nrow(RSet.betas1)
# p.s.all<-array(NA,dim=lng)
# for (i in 1:lng){
#   print(i)
#   y=RSet.betas1[i,]
#   m1 <- simplexreg(y ~ cond)
#   m1.sum<-summary(m1, save=TRUE)
#   p.s.all[i]<-m1.sum$coefficients$mean["condnormal","Pr(>|z|)"] 
# }
# time2<-Sys.time()
# (time.simplex<-time2-time1) #1.511832 hours, 
# simplex.res.all<-RSet.betas1[p.s.all<0.05,]
# dim(simplex.res.all) #148410

#amb dades normalitzades 
time1<-Sys.time()
lng<-nrow(GSet.b)
p.s.all<-array(NA,dim=lng)
for (i in 1:lng){
  print(i)
  y=GSet.b[i,]
  m1 <- simplexreg(y ~ cond)
  m1.sum<-summary(m1, save=TRUE)
  p.s.all[i]<-m1.sum$coefficients$mean["condnormal","Pr(>|z|)"] 
}
time2<-Sys.time()
(time.simplex<-time2-time1) #1,49 hours
simplex.res.all<-GSet.b[p.s.all<0.05,]
dim(simplex.res.all) #159927

save(simplex.res.all,file=file.path(resultsDir,"simplex.res.RData"))



#gamlss per totes les mostres
library(gamlss)
sink("NUL")

lng<-nrow(RSet.betas)
lng<-100000
p.all<-array(NA,dim=lng)
time1<-Sys.time()
for (i in 1:lng){
  y=RSet.betas[i,]
  m1 <- gamlss(y ~ cond,family=BEINF)
  m1.sum<-summary(m1, save=TRUE)
  p.all[i]<-m1.sum$pvalue["condnormal"] #check si ?s realment la p de la mu...
}
sink()
time2<-Sys.time()
(time.gamlss<-time2-time1) #Time difference of 13.42937 hours
gamlss.res.all<-RSet.betas[p.all<0.05,] # 3178 

#guardo objecte RData
#save(p.all,gamlss.res.all,file=file.path("F:/Doctorat/R","Gamlss.Minfidata.RData"))
load(file=file.path("D:/Doctorat/R","Gamlss.Minfidata.RData"))
dim(gamlss.res.all)

##### prova gamlss amb simplex (not inflated) ###########
# installar gamlss.dist 5.0.6 que dóna error! Ho he aconseguit: carregava gamlss.dist al namespace i 
# he aconseguit reinstal·lar-lo en la versió 5.0.6 dp gamlss 

time3<-Sys.time()
#lng<-nrow(RSet.betas)
p.s.all<-array(NA,dim=lng)
for (i in 1:lng){
  y=RSet.betas[i,]
  m1 <- gamlss(y ~ cond,family=SIMPLEX)
  m1.sum<-summary(m1, save=TRUE)
  p.s.all[i]<-m1.sum$pvalue["condnormal"] #check si ?s realment la p de la mu...
}
sink()
time4<-Sys.time()
(time.gamlss.simplex<-time4-time3) #Time difference of 13.42937 hours

pdf(file.path(resultsDir,"correl.gamlss.p.BEINF.SIMPLEX.pdf"))
  plot(p.all,p.s.all,main="all.p")
  plot(p.all,p.s.all,xlim=c(0,0.1),ylim=c(0,0.1),main="p.0.1")
dev.off()
#sembla que hi ha bona correl

############################# ZOIP #######################
#primer provo sobre les dades de simplexreg

data("sdac", package = "simplexreg")
sim.glm2 <- RM.ZOIP(formula.mu = rcd~ageadj+chemo, family="Simplex",
                    link = c("logit","identity","identity","identity"), data = sdac)
#sigma per defecte és ~1, i tb p0 i p1
summary(sim.glm2) 
#             Estimate Std. Error z value  Pr(>|z|)    
# X.Intercept. 1.1002248  0.1382182  7.9601 1.72e-15 ***
#   ageadj       0.0135749  0.0063455  2.1393  0.03241 *  
#   chemo        0.2660943  0.1235931  2.1530  0.03132 *  
  #molt semblant a simplexreg 

#ara modelo sigma també, igual que a simplexreg, que diu que utilitza com a link sempre log
#Del help de simplexreg: Four types of function are available linking the regressors to the mean. However, for dispersion, the link function is restricted to logarithm function. When modeling dispersion, the regressor modelling the dispersion parameter should be specified in a formula form of type y ~ x1 + x2 | z1 + z2 where z1 and z2 are linked to the dispersion parameter σ^2.

sim.glm22 <- RM.ZOIP(formula.mu = rcd~ageadj+chemo, formula.sigma = ~age, 
                     link = c("logit","log","identity","identity"), family="Simplex",data = sdac)
summary(sim.glm22)
# Fixed effects for logit(mu) 
# ---------------------------------------------------------------
#            Estimate Std. Error z value  Pr(>|z|)    
# X.Intercept. 1.1155503  0.2014173  5.5385 3.051e-08 ***
#   ageadj       0.0130135  0.0092742  1.4032    0.1606    
#   chemo        0.2519206  0.1807401  1.3938    0.1634    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# ---------------------------------------------------------------
#   Fixed effects for log(sigma) 
# ---------------------------------------------------------------
#                 Estimate Std. Error z value  Pr(>|z|)    
# X.Intercept..1  2.6075001  0.4520894  5.7677 8.038e-09 ***
#   age            -0.0150026  0.0089105 -1.6837   0.09224 .  

#upps, han deixat de ser significatius...mmm
AIC(sim.glm2,sim.glm22) #error

#ara provo amb les dades de les betes de minfi, normalitzades per quantile
time1<-Sys.time()
lng<-nrow(GSet.b)
p.z.all<-array(NA,dim=lng)
#primer muntar el data.frame, dp aniré creant un subdf amb la cond i la cpg a estudiar
df <- data.frame(cond,t(GSet.b))
for (i in 1:lng){
  print(i)
  cpg= colnames(df)[i+1]
  df.cpg <- df[,c("cond",cpg)]
  colnames(df.cpg) <- c("cond","cpg")
  m1 <- RM.ZOIP(formula.mu = cpg~cond,  
                link = c("logit","identity","identity","identity"), family="Simplex", data=df.cpg)
  #m1.sum<-summary(m1) #sembla que això no xuta
  #copio això de SummaryZOIP.R del paquet ZOIP, p de la mu
  estimate <- m1$par
  se       <- sqrt(diag(solve(m1$HM)))
  zvalue   <- estimate / se
  pvalue   <- 2 * stats::pnorm(abs(zvalue), lower.tail=F)
  p.z.all[i]<-pvalue["condnormal"]
}
time2<-Sys.time()
(time.ZOIP<-time2-time1) 
ZOIP.res.all<-GSet.b[p.z.all<0.05,]
dim(ZOIP.res.all) # 219270      6 la mitat!

save(ZOIP.res.all,file=file.path(resultsDir,"ZOIP.res.RData"))
