#05/05/19 FINAL ANALYSIS SCRIPT WITH THE SAME STRUCTURE FOR ALL
# la resta de proves estan a GSE116339_DEF.R i en els scripts previs

workingDir<-"D:/Doctorat/Simplex/MetDist/Data/GSE116339_PBB"
setwd(workingDir)
codeDir <- "D:/Doctorat/Simplex/MetDist/R"

#############################################################################
########################### 0. Data load ####################################
#############################################################################

load(file="PM.filtered.RData") #atenció li he posat els noms més a baix pq no havia guardat els CGs
dim(PM.f) #198711    679 un cop eliminat els SNPs

load(file="pheno.RData")

#######################################################
#functions load
source(file=file.path(codeDir,"call.functions.R")) #carrego directament la fn est.betabin.params
#to parallelize
library(doParallel)
library(foreach)

#############################################################################
####################### 1.Best distribution #################################
#############################################################################

N=nrow(PM.f)
cores=4
t1 <- Sys.time()
cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
best.dist.all<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- PM.f[i,]
  best.dist(b)
}
stopCluster(cl)
return(best.dist.all)
t2 <- Sys.time()
t2-t1 #6.918498 hours

head(best.dist.all) 

save(best.dist.all,file=file.path("best.dist.all.RData"))
load(file=file.path("best.dist.all.RData"))
dim(best.dist.all)
#198711      6

gse11.best.dist <- t(apply(best.dist.all,1,best.aic.ks))
table(gse11.best.dist[,1])
# beta.aic  normal.aic simplex.aic 
# 31236       33927      133548 

table(gse11.best.dist[,2])
# beta.ks.p  normal.ks.p simplex.ks.p 
# 34779        43748       120184 

#############################################################################
####################### 2.Estimate params for each distribution #############
#############################################################################
#no ha canviat res, aprofito objecte antic

# N <- nrow(PM.f)
# cores=4
# t1 <- Sys.time()
# cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
# registerDoParallel(cl)
# est.params<- foreach(i=1:N,.combine=rbind, .packages=c("fitdistrplus","VGAM","ZOIP")) %dopar% {
#   b <- PM.f[i,]
#   est.all.params(b)
# }
# stopCluster(cl)
# t2 <- Sys.time()
# t2-t1 
# head(est.params) #les betainflated donen na pq no hi ha 0's ni 1's
# 
# est.params <-as.data.frame(est.params)

#save(est.params,file="params.est.all.210219.RData")
load(file="params.est.all.RData")
dim(est.params)
#198711     14

na.col <- function(x) sum(is.na(x))
apply(est.params, 2, na.col)

# s.mle.mu   s.mle.sig   s.zoip.mu  s.zoip.sig    b.mom.s1    b.mom.s2    b.mle.s1    b.mle.s2 binf.mle.s1 
# 168         168           0           0           0           0           1           1      198619 
# binf.mle.s2     n.mom.m    n.mom.sd     n.mle.m    n.mle.sd 
# 198619           0           0           0           0 

#ranges
apply(est.params, 2, range, na.rm=T)
#         s.mle.mu  s.mle.sig  s.zoip.mu s.zoip.sig    b.mom.s1    b.mom.s2    b.mle.s1    b.mle.s2 binf.mle.s1
# [1,] 0.03923843  0.2415039 0.03923187  0.2414835   0.2117655   0.1443674   0.5029186   0.3629927    0.413979
# [2,] 0.96207305 11.1077817 0.96207895 11.1122240 163.8193951 163.5935903 167.3989935 165.7837980   59.933499
#       binf.mle.s2    n.mom.m   n.mom.sd    n.mle.m   n.mle.sd
# [1,]   0.4155871 0.03858037 0.03000018 0.03858037 0.02997808
# [2,]  67.3928258 0.96306584 0.42159840 0.96306584 0.42128783

#############################################################################
####################### 3.Simulations #######################################
#############################################################################
t1 <- Sys.time()
simulations.all <- fn.simulations(est.params=est.params,cond.n= c(3,5,10,30,100,500),cores=2)
t2 <- Sys.time()
t2-t1 #1.341569 hours  1.647001 hours (2 cores)



#############################################################################
#################### 4. Specific data analyisis: SEX ########################
#############################################################################

cond <- as.factor(car::recode(pheno$`gender:ch1`,"'Female'=1;'Male'=2"))

t1 <- Sys.time()
beta.models <- fn.models.parallel(PM.f, cond1=cond, cores=4)
t2 <- Sys.time()
t2-t1 #8.224955 hours 4 cores
dim(beta.models)  #198711
head(beta.models)

# PER ARREGLAR L'AIC de simplex inflated he creat una fn
# t1 <- Sys.time()
# zoip.models <- fn.models.parallel.onlyZOIP(PM.f, cond1=cond, cores=3)
# t2 <- Sys.time()
# t2-t1 # 3.997382 hours (3 cores)
# zoip.models <- as.data.frame(zoip.models)
# 
# #save(beta.models,file="betadata.models.sex.RData") #he matxacat els betamodels de smoke!!!
# load(file="betadata.models.sex.RData")
# beta.models <- as.data.frame(beta.models)
# #matxambro i guardo
# beta.models$aic.sinf <- zoip.models$aic.sinf
# save(beta.models,file="betadata.models.sex.RData")
# load(file="betadata.models.sex.RData")
# 



limma.res <- apply.limma(PM.f,cond) #cannot allocate vector, ho parteixo en trossos
limma.res1 <- apply.limma(PM.f[1:50000,],cond) 
limma.res2 <- apply.limma(PM.f[50001:100000,],cond) 
limma.res3 <- apply.limma(PM.f[100001:150000,],cond) 
limma.res4 <- apply.limma(PM.f[150001:198711,],cond) 

limma.res <- rbind(limma.res1,
                   limma.res2,
                   limma.res3,
                   limma.res4)

beta.models.withlimma <- data.frame(beta.models,p.limma=limma.res$P.Value)
#per a no liar-la també carrego les dades i modifico l'objecte
#save(beta.models.withlimma,file="betadata.models.sex.withlimma.RData")
# beta.models.withlimma <- data.frame(beta.models,p.limma=limma.res$P.Value)
# save(beta.models.withlimma,file="betadata.models.sex.withlimma.RData")
# load(file="betadata.models.sex.withlimma.RData")
# beta.models.withlimma$aic.sinf <- zoip.models$aic.sinf
# save(beta.models.withlimma,file="betadata.models.sex.withlimma.RData")

load(file="betadata.models.sex.withlimma.RData")

beta.models.adj.p <- apply(beta.models.withlimma,2,p.adjust)

# load(file="PM.f.annot.RData")
# dim(PM.f.annot) 198711    684
# # 
# #beta.models.anot <- cbind(beta.models,PM.f.annot[680:684])
# beta.models.adj.p.anot <- cbind(beta.models.adj.p,PM.f.annot[680:684])
save(beta.models.adj.p.anot,file="betadata.models.adj.p.anot.sex.RData") 

beta.models.adj.p.anot.p.s <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.s<0.05,]
beta.models.adj.p.anot.p.b <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.b<0.05,]
beta.models.adj.p.anot.p.sinf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.sinf<0.05,]
beta.models.adj.p.anot.p.binf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.binf<0.05,]
beta.models.adj.p.anot.p.n <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.n<0.05,]
beta.models.adj.p.anot.p.l <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.l<0.05,]
beta.models.adj.p.anot.p.q <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.q<0.05,]
beta.models.adj.p.anot.p.limma <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.limma<0.05,]

# 
# #sembla que hi ha moltes al chr X
res.chr.s <- as.data.frame(table(beta.models.adj.p.anot.p.s$chr))
res.chr.b <- as.data.frame(table(beta.models.adj.p.anot.p.b$chr))
res.chr.sinf <- as.data.frame(table(beta.models.adj.p.anot.p.sinf$chr))
res.chr.binf <- as.data.frame(table(beta.models.adj.p.anot.p.binf$chr))
res.chr.n <- as.data.frame(table(beta.models.adj.p.anot.p.n$chr))
res.chr.l <- as.data.frame(table(beta.models.adj.p.anot.p.l$chr))
res.chr.q <- as.data.frame(table(beta.models.adj.p.anot.p.q$chr))
res.chr.limma <- as.data.frame(table(beta.models.adj.p.anot.p.limma$chr))

names(res.chr.s)[2] <- "s"
names(res.chr.b)[2] <- "b"  
names(res.chr.sinf)[2] <- "sinf"  
names(res.chr.binf)[2] <- "binf"  
names(res.chr.n)[2] <- "n"  
names(res.chr.l)[2] <- "l"  
names(res.chr.q)[2] <- "q" 
names(res.chr.limma)[2] <- "limma" 

require(plyr)
res.table <- join_all(list(res.chr.s,
                           res.chr.b,
                           res.chr.sinf,
                           res.chr.binf,
                           res.chr.n,
                           res.chr.l,
                           res.chr.q,
                           res.chr.limma), by = 'Var1', type = 'full')

res.table 

#     Var1     s     b  sinf  binf     n    l     q limma
# 1   chr1  4493  4541  4998  4471  4407 4202  1198  4486
# 2  chr10  1978  1997  2239  1966  1945 1843   467  1973
# 3  chr11  2582  2601  2871  2562  2537 2401   661  2573
# 4  chr12  2329  2333  2587  2301  2276 2158   561  2321
# 5  chr13   897   907  1010   888   881  845   178   895
# 6  chr14  1577  1595  1734  1577  1549 1468   426  1587
# 7  chr15  1386  1399  1540  1380  1363 1293   363  1382
# 8  chr16  1771  1778  1968  1746  1727 1637   441  1752
# 9  chr17  2577  2585  2849  2543  2527 2401   717  2552
# 10 chr18   674   685   772   665   656  622   167   676
# 11 chr19  1861  1864  2139  1834  1823 1732   497  1846
# 12  chr2  3379  3422  3722  3360  3313 3123   797  3375
# 13 chr20  1383  1382  1500  1364  1347 1286   327  1375
# 14 chr21   510   507   563   500   494  459   101   502
# 15 chr22  1080  1081  1191  1064  1048 1004   321  1067
# 16  chr3  2543  2574  2853  2533  2495 2359   591  2547
# 17  chr4  1609  1627  1813  1601  1564 1473   353  1610
# 18  chr5  2105  2125  2417  2091  2066 1938   473  2104
# 19  chr6  2559  2568  2840  2529  2471 2357   596  2548
# 20  chr7  1967  1991  2226  1948  1917 1807   450  1958
# 21  chr8  1816  1842  2048  1805  1782 1683   415  1814
# 22  chr9  1369  1378  1527  1353  1343 1269   322  1362
# 23  chrX 12382 12392 12416 12384 12373 7841 12021 12394
# 24  chrY    86    85    86    85    85   28    80    85

save(res.table,file="betadata.models.sex.res.table.RData")


