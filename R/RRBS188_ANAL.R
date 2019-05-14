#05/05/19 FINAL ANALYSIS SCRIPT WITH THE SAME STRUCTURE FOR ALL
# la resta de proves estan a RRBS_DEF.R i en els scripts previs

workingDir<-"D:/Doctorat/Simplex/MetDist/Data/RRBS_188"
setwd(workingDir)
codeDir <- "D:/Doctorat/Simplex/MetDist/R"

#############################################################################
########################### 0. Data load ####################################
#############################################################################

#ratios
load(file="meth.betas.poquesNAs.covg3.RData") 
dim(rnb.meth.f)
#270569    188

#coverage
load(file="rnb.covg.poquesNAs.covg3.RData") 
dim(rnb.covg.f)
#270569    188

#i ara all
rnb.all.f <- rnb.meth.f * rnb.covg.f

#pheno data
samples <- read.table("samples.csv", sep=",", header=T, stringsAsFactors = F)
str(samples)

table(samples$patient_sex)
# f   m N/A 
# 63  96  29 
#according to exploratory analysis (cluster) there is a big difference in PCA in Tissue vs MSC
table(samples$patient_sex,samples$Tissue_vs_MSC)
#       Ewing_Tumor MSC
# f            56   4
# m            83   5
# N/A           1  23

#resulta que agafa patient_age com a factor...
samples$patient_age <- as.numeric(samples$patient_age)
range(samples$patient_age,na.rm=T)
# 0 62
hist(samples$patient_age,breaks=15)
hist(samples[samples$patient_sex=="f","patient_age"],breaks=15)
hist(samples[samples$patient_sex=="m","patient_age"],breaks=15)

t.test(samples[samples$patient_sex=="f","patient_age"],samples[samples$patient_sex=="m","patient_age"])
#sí que hi ha diferències...p-value = 0.0216

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
cores=4
t1 <- Sys.time()
cl <- makeCluster(cores,type="PSOCK",outfile="output.txt") #poso només 1 pq està processant l'altre
registerDoParallel(cl)
best.dist.all<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- rnb.meth.f[i,]
  best.dist(b)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #5.714317 hours 4 cores

head(best.dist.all) 
rownames(best.dist.all) <- rownames(rnb.meth.f)

#save(best.dist.all,file=file.path("best.dist.all.RData"))
load(file=file.path("best.dist.all.RData"))

dim(best.dist.all)
# 270569      6

rrbs.best.dist <- t(apply(best.dist.all,1,best.aic.ks))
table(rrbs.best.dist[,1])
# beta.aic  normal.aic simplex.aic 
# 43138        1949      225482 

#######################################################
#agafo els vells
load(file=file.path("D:/Doctorat/Simplex/Data_OLD/RRBS_188","best.dist.betabin.filtered.RData"))
head(best.dist.all.betabin) #si que hi ha els rownames!
dim(best.dist.all.betabin) #479057      2
best.dist.all.betabin <- best.dist.all.betabin[rownames(best.dist.all),]
dim(best.dist.all.betabin)
#  270569      2
head(best.dist.all.betabin) #moltes NAs
best.dist.all.betabin <- as.data.frame(best.dist.all.betabin)
#save(best.dist.all.betabin,file=file.path("best.dist.betabin.RData"))
load(file=file.path("best.dist.betabin.RData"))


#############################################################################
####################### 2.Estimate params for each distribution #############
#############################################################################
n <- nrow(rnb.meth.f)
t1 <- Sys.time()
cl <- makeCluster(7,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
rnb.est.params<- foreach(i=1:n,.combine=rbind, .packages=c("fitdistrplus","VGAM","ZOIP")) %dopar% {
  #b <- rnb.meth[i,]
  b <- rnb.meth.f[i,] #filtrant els que tenen menys de 30
  est.all.params(b)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #2.149716 hours (7 cores)
# 
head(rnb.est.params)
rnb.est.params <-as.data.frame(rnb.est.params)
rownames(rnb.est.params) <- rownames(rnb.meth.f)
save(rnb.est.params,file=file.path("params.est.all.RData"))

# #NO: HO HE FET REALMENT PQ limma sortia sospitosament alt i estava bé però he comparat i estava bé
# aprofito les dades antigues i faig subsetting, doncs d'aquí no he canviat res
#  load(file=file.path("D:/Doctorat/Simplex/Data_OLD/RRBS_188","params.est.all.RData"))
#  dim(rnb.est.params)# 479057      14
#  head(rnb.est.params) #ja té els rownames!
#  rnb.est.params <- rnb.est.params[rownames(rnb.meth.f),]
# # dim(rnb.est.params)
# # #  270569     14
# # head(rnb.est.params)
# # save(rnb.est.params,file=file.path("params.est.all.RData"))
# load(file=file.path("params.est.all.RData"))
# dim(rnb.est.params)
#  270569     14

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
load(file=file.path("D:/Doctorat/Simplex/Data_OLD/RRBS_188","betabin.params.est.RData"))
dim(rnb.betabin.est.params)# 479057      14
head(rnb.betabin.est.params) #ja té els rownames!
rnb.betabin.est.params <- rnb.betabin.est.params[rownames(rnb.meth.f),]
rnb.betabin.est.params <- as.data.frame(rnb.betabin.est.params)
#guardo
save(rnb.betabin.est.params,file="betabin.params.est.RData")
load(file="betabin.params.est.RData")


#############################################################################
####################### 3.Simulations #######################################
#############################################################################
t1 <- Sys.time()
simulations.all <- fn.simulations(est.params=rnb.est.params,cond.n= c(3,5,10,30,100,500),cores=7)
t2 <- Sys.time()
t2-t1 #2.028116 hours (2 cores) 3.054866 hours (1 core)  55.81003 mins (3 cores) 29.87096 mins (7 cores)


#############################################################################
#################### 4. Specific data analyisis: SEX ########################
#############################################################################

table(samples$patient_sex)
# f   m N/A 
# 63  96  29 
samples2comp <- samples[samples$patient_sex %in% c("f","m"),]
dim(samples2comp) #159
rnb.meth2comp <- rnb.meth.f[,samples2comp$sample_id]
dim(rnb.meth2comp) #270569    159
table(colnames(rnb.meth2comp)==samples2comp$sample_id) #true!

cond <- as.factor(car::recode(samples2comp$patient_sex,"'f'=1;'m'=2"))

t1 <- Sys.time()
rnb.meth2comp.models <- fn.models.parallel(rnb.meth2comp, cond1=cond, cores=4)
t2 <- Sys.time()
t2-t1 # 4.074444 hours (4 cores)
head(rnb.meth2comp.models)


limma.res <- apply.limma(rnb.meth2comp,cond)

# rnb.meth2comp.models.withlimma <- data.frame(rnb.meth2comp.models,p.limma=limma.res$P.Value)
# #save(rnb.meth2comp.models.withlimma,file="rnb.meth2comp.models.withlimma.RData")
# load(file="rnb.meth2comp.models.withlimma.RData")
# 
# # PER ARREGLAR L'AIC de simplex inflated he creat una fn
# # #ho faig aquí, no havia guardat els models sense limma
# t1 <- Sys.time()
# zoip.models <- fn.models.parallel.onlyZOIP(rnb.meth2comp, cond1=cond, cores=2)
# t2 <- Sys.time()
# t2-t1 # 2.486902 hours (2 cores)
# zoip.models <- as.data.frame(zoip.models)
# 
# # #matxambro i guardo
# rnb.meth2comp.models.withlimma$aic.sinf <- zoip.models$aic.sinf
# save(rnb.meth2comp.models.withlimma,file="rnb.meth2comp.models.withlimma.RData")
load(file="rnb.meth2comp.models.withlimma.RData")

############# regressio betabinomial
all <-rnb.all.f[,samples2comp$sample_id]
covg <-rnb.covg.f[,samples2comp$sample_id]
# # remove=c(64385,257324) nose pq aquestes es salten els try i retornen un error
# all[64385,] <- NA
# all[257324,] <- NA

t1 <- Sys.time()
rnb.overd<-fn.models.betabin.parallel(all=all,covg=covg, cond1=cond, cores=7)
t2 <- Sys.time()
t2-t1 # Time difference of 25.00339 mins
head(rnb.overd)

#RENOMBRO!!
rnb.meth2comp.betabin.models <- as.data.frame(rnb.overd)

#save(rnb.meth2comp.betabin.models,file="rnb.meth2comp.betabin.models.RData")
#load(file="rnb.meth2comp.betabin.models.RData")

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
#     chr   s   b sinf binf   n   l   q limma   bb
# 1   chr1   3   2   85   NA  NA  NA  NA    NA   NA
# 2  chr10   3   2   47   NA  NA  NA  NA    NA   NA
# 3  chr11   4   1   39   NA  NA  NA  NA    NA   NA
# 4  chr12   3   2   49   NA  NA  NA  NA    NA   NA
# 5  chr13   1   2   31   NA   1  NA  NA    NA   NA
# 6  chr15   3   1   35   NA  NA  NA  NA    NA    1
# 7  chr16   7   4   42   NA  NA  NA  NA    NA   NA
# 8  chr17   5   2   60   NA  NA  NA  NA    NA   NA
# 9  chr19   9   8   71   NA  NA  NA  NA    NA   NA
# 10  chr2   2   1   57   NA  NA  NA  NA    NA   NA
# 11 chr20   4  NA   51   NA  NA  NA  NA    NA    1
# 12 chr21   2   1   15   NA  NA  NA  NA    NA   NA
# 13 chr22   2   2   34   NA  NA  NA  NA    NA    1
# 14  chr3   3   2   47   NA  NA  NA  NA    NA    2
# 15  chr4   3   2   30   NA  NA  NA  NA    NA   NA
# 16  chr5   4   3   47   NA  NA  NA  NA    NA   NA
# 17  chr6   5   1   58   NA   2  NA  NA    NA   NA
# 18  chr7   1   1   56   NA   1  NA  NA    NA   NA
# 19  chr8   3   2   50   NA  NA  NA  NA    NA    1
# 20  chr9   2   2   36   NA  NA  NA  NA    NA    3
# 21  chrX 788 727 1434  535 601 234 311   683 7952
# 22  chrY   5   4    6    4   5   2  NA     4   21
# 23 chr14  NA   2   41   NA  NA  NA  NA    NA   NA
# 24 chr18  NA  NA   13   NA  NA  NA  NA    NA   NA
save(res.table,file="betadata.models.sex.res.table.RData") #li poso el mateix nom que als arrays

