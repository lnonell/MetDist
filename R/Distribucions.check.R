#3/1/18 intento veure quina distribució tenen les dades
#https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
#passo la fn a SimulationFunctions per a poder utilitzar-la des de tots els data sets


library(simplexreg)
library(fitdistrplus) #conté la fn fitdist

b <-rbeta(1000,2,2)
plot(density(b))
mu <- mean(b)
var <- var(b)

descdist(b,discrete=F) #plot molt mono però realment no el puc aprofitar, retorna els params
# min, max, median, mean, sd, skewness and kurtosis però amb això no podem estimar

fit.b <- fitdistr(b, densfun=dbeta, start=list(shape1=mu,shape2=var))
plot(fit.b) #error

fit.b <- fitdist(b, distr=dbeta, start=list(shape1=mu,shape2=var)) #aqta és del pacque fitdistrplus
plot(fit.b) #ara sí 
fit.b$aic #-260.1639

fit.n <- fitdist(b, distr="norm")
plot(fit.n) #error
fit.n$aic #-172.9212

fit.s <-fitdist(b, distr=dsimplex, start=list(mu=0.5,sig=1))
fit.s$aic #-162.2187

########## Kolmogorov-Smirnov test
ks.b <- ks.test(r, "pbeta", fit.b$estimate[1],fit.b$estimate[2])
ks.b$p.value
#0.5371517
ks.s <- ks.test(r, psimplex, fit.s$estimate[1],fit.s$estimate[2])
ks.s$p.value
#0.0002323421

ks.n <- ks.test(r, pnorm, fit.n$estimate[1],fit.n$estimate[2])
ks.n$p.value
#0.1098184

#està clar que el millor és ks.b 
#HAURIA DE FER AIXÒ PER A TOTES LES DADES

#faig la fn
x <- runif(100)
best.dist(x)

best.dist <- function(x){
  library(fitdistrplus)
  library(simplexreg)
  mu <- mean(x)
  var <- var(x)
  
  fit.s <- try(fitdist(x, distr=dsimplex, start=list(mu=0.5,sig=2),optim.method="Nelder-Mead"),TRUE) #no entenc pq falla amb les dades guardades!!
  if(class(fit.s)!="try-error")  s.aic <- fit.s$aic else s.aic=NA
  
  fit.b <- try(fitdist(x, distr=dbeta, start=list(shape1=mu,shape2=var)),TRUE) #aqta és del paquet fitdistrplus
  if(class(fit.b)!="try-error")  b.aic <- fit.b$aic else b.aic=NA
  
  fit.n <- try(fitdist(x, distr="norm"),TRUE)
  if(class(fit.n)!="try-error")  n.aic <- fit.n$aic else n.aic=NA
  
  ########## Kolmogorov-Smirnov test
  ks.s <- try(ks.test(x, psimplex, fit.s$estimate[1],fit.s$estimate[2]),TRUE)
  if(class(ks.s)!="try-error") s.ks <- ks.s$p.value else s.ks=NA
  
  ks.b <- try(ks.test(x, "pbeta", fit.b$estimate[1],fit.b$estimate[2]),TRUE)
  if(class(ks.b)!="try-error") b.ks <- ks.b$p.value else b.ks=NA
  
  ks.n <- try(ks.test(x, pnorm, fit.n$estimate[1],fit.n$estimate[2]),TRUE)
  if(class(ks.n)!="try-error") n.ks <- ks.n$p.value else n.ks=NA
  
  v <- c(s.aic,b.aic,n.aic,s.ks,b.ks,n.ks)
  names(v) <- c("simplex.aic","beta.aic","normal.aic","simplex.ks.p","beta.ks.p","normal.ks.p")
  return(v)
}

#crear fn en paral·lel i anar creant vectors b, s, i n amb diferents paràmetres i anar guardant els resultats de la fn
b <-rbeta(200,2,2)
best.dist(b)


b <-rbeta(1000,0.5,0.5)
best.dist(b) #error

################
library(doParallel)
library(foreach)

#puc provar per a les dades de cada subset
#RRBS_216 
#load(file=file.path("D:/Doctorat/Simplex/Data/RRBS_216_Tissue_CL","simulated.objects.5percdml.2000.100.100.RData"))
load(file=file.path("D:/Doctorat/Simplex/Data/GSE50660_Smoking","simulated.objects.5percdml.2000.100.100.RData"))

N=nrow(cpgs.beta)
e=0.01
t1 <- Sys.time()
cl <- makeCluster(4,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
b.best.dist<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- round(cpgs.beta[i,],4)
  b.prime <- ifelse(b==1, 1-e, ifelse(b==0,0+e,b))
  best.dist(b.prime)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #18 secs super schnell!!

head(b.best.dist) #no estima bé simplex per a les RRBS, per GSE50660 sí!



N=nrow(cpgs.simplex)
e=0.01
t1 <- Sys.time()
cl <- makeCluster(4,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
s.best.dist<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- round(cpgs.simplex[i,],4)
  b.prime <- ifelse(b==1, 1-e, ifelse(b==0,0+e,b))
  best.dist(b.prime)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #2.88 mins

head(s.best.dist) #estima bé tot..!

N=nrow(cpgs.normal)
e=0.01
t1 <- Sys.time()
cl <- makeCluster(4,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
n.best.dist<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- round(cpgs.normal[i,],4)
  b.prime <- ifelse(b==1, 1-e, ifelse(b==0,0+e,b))
  best.dist(b.prime)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #

head(n.best.dist) #ha donat un error: sembla que hi ha una de les cpgs amb tot 0!!
#kk no estima bé ni simplex ni beta per a cap data set, pq surten de l'interval 0,1!

#fn per a trobar el millor dels aic i ks
best.dist.i <- b.best.dist[1,]

best.aic.ks <- function(best.dist.i){
  best.aic <- names(which.min(best.dist.i[1:3]))
  best.ks <- names(which.max(best.dist.i[4:6]))
  return(c(best.aic,best.ks))
}

best.aic.ks(best.dist.i)

#simplex
s.best <- apply(s.best.dist,1,best.aic.ks)
table(s.best[1,])
# beta.aic  normal.aic simplex.aic 
# 395         203        1402 
table(s.best[2,])
# beta.ks.p  normal.ks.p simplex.ks.p 
# 413          351         1236 

#betes
b.best <- apply(b.best.dist,1,best.aic.ks)

table(b.best[1,])
# beta.aic  normal.aic simplex.aic 
# 781         392         827 
table(b.best[2,])
# beta.ks.p  normal.ks.p simplex.ks.p 
# 685          550          765 

########################################
#estudio una mica els resultats
#simplex, miro alguns normals
which(s.best[1,]=="normal.aic")
#result.1    result.4   result.29   result.43...
hist(cpgs.simplex[1,])
hist(cpgs.simplex[4,])
#té sentit!
#simplex, miro algunes betes
which(s.best[1,]=="beta.aic")
#result.9   result.11   result.50   result.66   result.79..
hist(cpgs.simplex[9,])
hist(cpgs.simplex[11,])
#són com normals skewed a dreta o a esquerra
which(s.best[1,]=="simplex.aic")
#result.2    result.3    result.5    result.6    result.7    result.8   result.10
hist(cpgs.simplex[2,]) #té cues
hist(cpgs.simplex[3,]) #cua a la dreta
hist(cpgs.simplex[5,]) #bimodal
hist(cpgs.simplex[6,]) #bimodal
hist(cpgs.simplex[7,]) #bimodal
hist(cpgs.simplex[8,]) #bimodal

#les bimodals les estima millor amb la simplex!
#miro què passa amb els ks
table(s.best[,s.best[1,]=="simplex.aic"]) #??
table(s.best[,s.best[1,]=="beta.aic"]) #??
table(s.best[,s.best[1,]=="normal.aic"]) #??

#potser té més sentit aic que p del Kolmogorov-smirnov test

################################## proves ################################
#segons aconsellen s'hauria de generar una mostra amb les estimacions i agafar els estimates! INVIABLE
n.sims <- 10

stats <- replicate(n.sims, {      
  r <- rbeta(n = 1000, 2,2)
  estfit.b <- fitdist(r, "beta") # added to account for the estimated parameters
  as.numeric(ks.test(r, "pbeta", estfit.b$estimate[1],estfit.b$estimate[2])$statistic)
})

########### recomanent empirical cumulative distribution function (ecdf)


