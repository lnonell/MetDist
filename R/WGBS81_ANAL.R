#05/05/19 FINAL ANALYSIS SCRIPT WITH THE SAME STRUCTURE FOR ALL
#changed dir from Bueprint to Blueprint finally
# la resta de proves estan a WGBS_DEF.R i en els scripts previs

workingDir<-"D:/Doctorat/Simplex/MetDist/Data/WGBS_81_Blueprint"
setwd(workingDir)
codeDir <- "D:/Doctorat/Simplex/MetDist/R"

#############################################################################
########################### 0. Data load ####################################
#############################################################################

#ratios
load(file="meth.betas.poquesNAs.covg3.RData") 
dim(rnb.meth.f)
#204073     81

#coverage
load(file="rnb.covg.poquesNAs.covg3.RData") 
dim(rnb.covg.f)
#204073     81

#i ara all
rnb.all.f <- rnb.meth.f * rnb.covg.f

#pheno data
samples <- read.delim("samples.tsv", stringsAsFactors=FALSE)
dim(samples) # 81 52

table(samples$DONOR_SEX)
#     Female   Male 
# 4     47     30 

#######################################################
#functions load
source(file=file.path(codeDir,"call.functions.R")) #carrego directament la fn est.betabin.params
#to parallelize
library(doParallel)
library(foreach)

#############################################################################
####################### 1.Best distribution #################################
#############################################################################

N=nrow(rnb.meth.f)
cores=3
t1 <- Sys.time()
cl <- makeCluster(cores,type="PSOCK",outfile="output.txt") #poso només 1 pq està processant l'altre
registerDoParallel(cl)
best.dist.all<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- rnb.meth.f[i,]
  best.dist(b)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #2.617019 hours

head(best.dist.all) 

save(best.dist.all,file=file.path("best.dist.all.RData"))
load(file=file.path("best.dist.all.RData"))

dim(best.dist.all)
# 204073      6

wgbs.best.dist <- t(apply(best.dist.all,1,best.aic.ks))
table(wgbs.best.dist[,1])
# beta.aic  normal.aic simplex.aic 
# 108304       10163       85606 
table(wgbs.best.dist[,2])
# beta.ks.p  normal.ks.p simplex.ks.p 
# 79761        64951        59361 

#######################################################
#AKIbeta-binomial: AKI FER SUBSETTING

# N=nrow(rnb.meth.f)
# t1 <- Sys.time()
# cl <- makeCluster(6,type="PSOCK",outfile="output.txt") #poso només 1 pq està processant l'altre
# registerDoParallel(cl)
# best.dist.all.betabin<- foreach(i=1:N,.combine=rbind, .packages=c("VGAM","fitdistrplus")) %dopar% {
#   xi <- rnb.all.f[i,]
#   ni <- rnb.covg.f[i,]
#   best.dist.betabin(xi,ni)
# }
# stopCluster(cl)
# t2 <- Sys.time()
# t2-t1 #7.573625 hours 6 cores

#agafo els vells, els he anotat!!
# load(file=file.path("D:/Doctorat/Simplex/Data_OLD/WGBS_81_Bueprint","best.dist.betabin.filtered.RData"))
# head(best.dist.all.betabin) #ja hi ha els rownames
# best.dist.all.betabin <- best.dist.all.betabin[rownames(rnb.meth.f),]
# dim(best.dist.all.betabin)
# # 204073      2
# head(best.dist.all.betabin) #moltes NAs
# best.dist.all.betabin <- as.data.frame(best.dist.all.betabin)
# save(best.dist.all.betabin,file=file.path("best.dist.betabin.RData"))
load(file=file.path("best.dist.betabin.RData"))

#############################################################################
####################### 2.Estimate params for each distribution #############
#############################################################################
#aprofito els antics i faig subsetting 
# n <- nrow(rnb.meth.f)
# t1 <- Sys.time()
# cl <- makeCluster(4,type="PSOCK",outfile="output.txt")
# registerDoParallel(cl)
# rnb.est.params<- foreach(i=1:n,.combine=rbind, .packages=c("fitdistrplus","VGAM","ZOIP")) %dopar% {
#   #b <- rnb.meth[i,]
#   b <- rnb.meth.f[i,] #filtrant els que tenen menys de 30
#   est.all.params(b)
# }
# stopCluster(cl)
# t2 <- Sys.time()
# t2-t1 #1.24 hours 6 cores #2.69 hours 4cores
# 
# head(rnb.est.params)
# rnb.est.params <-as.data.frame(rnb.est.params)

# load(file=file.path("D:/Doctorat/Simplex/Data_OLD/WGBS_81_Bueprint","params.est.all.090319.RData"))
# dim(rnb.est.params)
# #213748     14
# rnb.est.params <- rnb.est.params[rownames(rnb.meth.f),]
# dim(rnb.est.params)
# # 204073      14
# head(rnb.est.params) 
#save(rnb.est.params,file=file.path("params.est.all.RData"))
load(file=file.path("params.est.all.RData"))

######################
#beta binomial

# N <- nrow(rnb.all.f)
# t1 <- Sys.time()
# cl <- makeCluster(6,type="PSOCK",outfile="output.txt")
# registerDoParallel(cl)
# rnb.betabin.est.params<- foreach(i=1:N,.combine=rbind, .packages=c("fitdistrplus","VGAM")) %dopar% {
#   xi <- rnb.all.f[i,]
#   ni <- rnb.covg.f[i,]
#   est.betabin.params(xi,ni)
# }
# stopCluster(cl)
# t2 <- Sys.time()
# t2-t1 #11.30993 hours
# head(rnb.betabin.est.params,25)

# aprofito les dades antigues i faig subsetting, doncs d'aquí no he canviat res
load(file=file.path("D:/Doctorat/Simplex/Data_OLD/WGBS_81_Bueprint","betabin.params.est.RData"))
dim(rnb.betabin.est.params)# 213748      6
head(rnb.betabin.est.params) #ja té els rownames!
rnb.betabin.est.params <- rnb.betabin.est.params[rownames(rnb.meth.f),]
dim(rnb.betabin.est.params)# 204073      6
#guardo
save(rnb.betabin.est.params,file="betabin.params.est.RData")
load(file="betabin.params.est.RData")
head(rnb.betabin.est.params)


#############################################################################
####################### 3.Simulations #######################################
#############################################################################

t1 <- Sys.time()
simulations.all <- fn.simulations(est.params=rnb.est.params,cond.n= c(3,5,10,30,100,500),cores=3)
t2 <- Sys.time()
t2-t1 #2.028116 hours (2 cores) 3.147326 hours (1 core) 56.85982 mins (3 cores)

load(file="simulated.cpgs.list.RData")

############ betabin
t1 <- Sys.time()
simulations.bb.all <- fn.simulations.betabin(est.params=rnb.betabin.est.params,cond.n= c(3,5,10,30,100,500),cores=6)
t2 <- Sys.time()
t2-t1 #3.713122 mins (6 cores)


#############################################################################
#################### 4. Specific data analyisis: SEX ########################
#############################################################################

samples2comp <- samples[samples$DONOR_SEX %in% c("Female","Male"),]
dim(samples2comp) #77
rnb.meth2comp <- rnb.meth.f[,samples2comp$sampleName]
dim(rnb.meth2comp) # 204073     77
table(colnames(rnb.meth2comp)==samples2comp$sampleName) #true!

cond <- as.factor(car::recode(samples2comp$DONOR_SEX,"'Female'=1;'Male'=2"))

t1 <- Sys.time()
rnb.meth2comp.models <- fn.models.parallel(rnb.meth2comp, cond1=cond, cores=3)
t2 <- Sys.time()
t2-t1 # 4.37942 hours (3 cores)

dim(rnb.meth2comp.models) # 204073     14
limma.res <- apply.limma(rnb.meth2comp,cond)

rnb.meth2comp.models.withlimma <- data.frame(rnb.meth2comp.models,p.limma=limma.res$P.Value)
# #save(rnb.meth2comp.models.withlimma,file="rnb.meth2comp.models.withlimma.RData")
# load(file="rnb.meth2comp.models.withlimma.RData")

# PER ARREGLAR L'AIC de simplex inflated he creat una fn
# #ho faig aquí, no havia guardat els models sense limma
# t1 <- Sys.time()
# zoip.models <- fn.models.parallel.onlyZOIP(rnb.meth2comp, cond1=cond, cores=2)
# t2 <- Sys.time()
# t2-t1 #  2.094571 hours (2 cores)
# zoip.models <- as.data.frame(zoip.models)
# 
# #matxambro i guardo
# rnb.meth2comp.models.withlimma$aic.sinf <- zoip.models$aic.sinf
#save(rnb.meth2comp.models.withlimma,file="rnb.meth2comp.models.withlimma.RData")
load(file="rnb.meth2comp.models.withlimma.RData")



#RENOMBRO PER APROFITAR EL QUE HI HAVIA
rnb.meth2comp.models <- rnb.meth2comp.models.withlimma

############# regressio betabinomial

all <-rnb.all.f[,samples2comp$sampleName]
covg <-rnb.covg.f[,samples2comp$sampleName]

t1 <- Sys.time()
rnb.overd<-fn.models.betabin.parallel(all=all,covg=covg, cond1=cond, cores=7)
t2 <- Sys.time()
t2-t1 # Time difference of 25.00339 mins
head(rnb.overd)

#RENOMBRO!!
rnb.meth2comp.betabin.models <- as.data.frame(rnb.overd)

#save(rnb.meth2comp.betabin.models,file="rnb.meth2comp.betabin.models.RData")
load(file="rnb.meth2comp.betabin.models.RData")

all.models.sex <- cbind(rnb.meth2comp.models.withlimma,rnb.meth2comp.betabin.models)
head(all.models.sex)

#ho guardo. IMPORTANT: ultima col es p.bb i aic.bb

#ajusto pvals
all.models.sex.adj.p <- as.data.frame(apply(all.models.sex,2,p.adjust))

all.models.sex.adj.p.s <- all.models.sex.adj.p[all.models.sex.adj.p[,1]<0.05 & 
                                                 !is.na(all.models.sex.adj.p[,1]),]
all.models.sex.adj.p.b <- all.models.sex.adj.p[all.models.sex.adj.p[,2]<0.05 & 
                                                 !is.na(all.models.sex.adj.p[,2]),]
all.models.sex.adj.p.sinf <- all.models.sex.adj.p[all.models.sex.adj.p[,3]<0.05 & 
                                                    !is.na(all.models.sex.adj.p[,3]),]
all.models.sex.adj.p.binf <- all.models.sex.adj.p[all.models.sex.adj.p[,4]<0.05 & 
                                                    !is.na(all.models.sex.adj.p[,4]),]
all.models.sex.adj.p.n <- all.models.sex.adj.p[all.models.sex.adj.p[,5]<0.05 & 
                                                 !is.na(all.models.sex.adj.p[,5]),]
all.models.sex.adj.p.l <- all.models.sex.adj.p[all.models.sex.adj.p[,6]<0.05 & 
                                                 !is.na(all.models.sex.adj.p[,6]),]
all.models.sex.adj.p.q <- all.models.sex.adj.p[all.models.sex.adj.p[,7]<0.05 & 
                                                 !is.na(all.models.sex.adj.p[,7]),]
all.models.sex.adj.p.limma <- all.models.sex.adj.p[all.models.sex.adj.p$p.limma<0.05 & 
                                                     !is.na(all.models.sex.adj.p$p.limma),]
all.models.sex.adj.p.bb <- all.models.sex.adj.p[all.models.sex.adj.p$p.bb<0.05 & 
                                                  !is.na(all.models.sex.adj.p$p.bb),]

#extrec els resultats (fn a SimulationFunctions.R)
res.chr.s <- res.chr(res=all.models.sex.adj.p.s)
res.chr.b <- res.chr(res=all.models.sex.adj.p.b)
res.chr.sinf <- res.chr(res=all.models.sex.adj.p.sinf)
res.chr.binf <- res.chr(res=all.models.sex.adj.p.binf)
res.chr.n <- res.chr(res=all.models.sex.adj.p.n)
res.chr.l <- res.chr(res=all.models.sex.adj.p.l)
res.chr.q <- res.chr(res=all.models.sex.adj.p.q)
res.chr.limma <- res.chr(res=all.models.sex.adj.p.limma)
res.chr.bb <- res.chr(res=all.models.sex.adj.p.bb)

names(res.chr.s)[2] <- "s"
names(res.chr.b)[2] <- "b"  
names(res.chr.sinf)[2] <- "sinf"  
names(res.chr.binf)[2] <- "binf"  
names(res.chr.n)[2] <- "n"  
names(res.chr.l)[2] <- "l"  
names(res.chr.q)[2] <- "q"  
names(res.chr.limma)[2] <- "limma"  
names(res.chr.bb)[2] <- "bb"  
require(plyr)

#hi ha algunes que no donen resultats
res.table <- join_all(list(res.chr.s,
                           res.chr.b,
                           res.chr.sinf,
                           res.chr.binf,
                           res.chr.n,
                           res.chr.l,
                           res.chr.q,
                           res.chr.limma,
                           res.chr.bb), by = 'chr', type = 'full')

res.table 

#     chr    s    b  sinf binf    n  l    q limma    bb
# 1   chr1  709   NA  2755   NA   NA NA   NA    NA   450
# 2  chr10  418    1  2133   NA   NA NA   NA    NA   267
# 3  chr11  559   NA  1845   NA   NA NA   NA    NA   512
# 4  chr12  343    2  1772    1    2 NA    1     2   245
# 5  chr13  304    1  1577   NA   NA NA   NA    NA   180
# 6  chr14  407    1  1860   NA   NA NA   NA    NA   241
# 7  chr15  474   NA  1563   NA   NA NA   NA    NA   382
# 8  chr16  386    2  1311   NA    1 NA   NA    NA   293
# 9  chr17  440    1  1553   NA   NA NA   NA    NA   348
# 10 chr18  306   NA  1443   NA   NA NA   NA    NA   153
# 11 chr19  545    1  2087   NA   NA NA    9    NA   369
# 12  chr2  764   NA  3315   NA   NA NA   NA    NA   530
# 13 chr20  331    2  1245   NA    1 NA   NA    NA   215
# 14 chr21   75   NA   340   NA   NA NA   NA    NA    49
# 15 chr22  156   NA   634   NA   NA NA    2    NA   111
# 16  chr3  337   NA  1777   NA   NA NA   NA    NA   214
# 17  chr4  431   NA  2440   NA   NA NA   NA    NA   248
# 18  chr5  838    1  3304   NA   NA NA   NA    NA   547
# 19  chr6  561    1  2174   NA   NA NA   NA    NA   423
# 20  chr7  526    2  2417   NA   NA NA   NA    NA   397
# 21  chr8  365    1  1800   NA   NA NA   NA    NA   232
# 22  chr9  339    1  1579   NA   NA NA    1    NA   198
# 23  chrX 8299 7822 12688 5680 5780  9 2025  6693 10508
# 24  chrY   10   15    25    5    7 NA    1     2    17

save(res.table,file="betadata.models.sex.res.table.RData") #li poso el mateix nom que als arrays

