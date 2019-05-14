#05/05/19 FINAL ANALYSIS SCRIPT WITH THE SAME STRUCTURE FOR ALL
# la resta de proves estan a GSE50660_DEF.R i en els scripts previs

workingDir <- "D:/Doctorat/Simplex/MetDist/Data/GSE50660_Smoking"
setwd(workingDir)
codeDir <- "D:/Doctorat/Simplex/MetDist/R"

#############################################################################
########################### 0. Data load ####################################
#############################################################################

load(file="pD.all.RData") #phenoData
load(file="beta.filtered.RData") #carrega les betes preprocessades i filtrades que són les que farem servir aquí
dim(beta) #125950    464

#Condition to analyze
sex <- pD.all$`gender:ch1` 

#######################################################
#functions load
source(file=file.path(codeDir,"call.functions.R")) #carrego directament la fn est.betabin.params
#to parallelize
library(doParallel)
library(foreach)

#############################################################################
####################### 1.Best distribution #################################
#############################################################################

N=nrow(beta)
cores=3
t1 <- Sys.time()
cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
best.dist.all<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- beta[i,]
  best.dist(b)
}
stopCluster(cl)
return(best.dist.all)
t2 <- Sys.time()
t2-t1 #3.952837 hours 3 cores

head(best.dist.all) 

#save(best.dist.all,file=file.path("best.dist.all.RData"))
load(file=file.path("best.dist.all.RData"))
dim(best.dist.all)
#125950      6

gse50660.best.dist <- t(apply(best.dist.all,1,best.aic.ks))
table(gse50660.best.dist[,1])
# beta.aic  normal.aic simplex.aic 
# 19737       32179       74034 
#genial!!!

#############################################################################
####################### 2.Estimate params for each distribution #############
#############################################################################
#no ha canviat la fn, aprofito els antics

# N <- nrow(beta)
# cores=4
# t1 <- Sys.time()
# cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
# registerDoParallel(cl)
# est.params<- foreach(i=1:N,.combine=rbind, .packages=c("fitdistrplus","VGAM","ZOIP")) %dopar% {
#   b <- beta[i,]
#   est.all.params(b)
# }
# stopCluster(cl)
# t2 <- Sys.time()
# t2-t1 #4.878 hours
# head(est.params) #les betainflated donen na pq no hi ha 0's ni 1's

# est.params <-as.data.frame(est.params)

#save(est.params,file="params.est.all.120119.RData")
load(file="params.est.all.RData")
dim(est.params)
# 125950     14

na.col <- function(x) sum(is.na(x))
apply(est.params, 2, na.col)
# s.mle.mu   s.mle.sig   s.zoip.mu  s.zoip.sig    b.mom.s1    b.mom.s2    b.mle.s1    b.mle.s2 binf.mle.s1 
# 31          31           0           0           0           0           0           0      125903 
# binf.mle.s2     n.mom.m    n.mom.sd     n.mle.m    n.mle.sd 
# 125903           0           0           0           0 

#ranges
apply(est.params, 2, range, na.rm=T)
#         s.mle.mu s.mle.sig  s.zoip.mu s.zoip.sig    b.mom.s1   b.mom.s2    b.mle.s1    b.mle.s2 binf.mle.s1
# [1,] 0.03675831 0.2413557 0.03675583  0.2413354   0.2108814   0.188386   0.5140319   0.5199418    1.133821
# [2,] 0.96840259 6.0728849 0.96839279  6.0717148 163.6719375 163.913653 188.3485588 189.6722137   57.708526
#       binf.mle.s2    n.mom.m   n.mom.sd    n.mle.m   n.mle.sd
# [1,]   0.8942635 0.03117263 0.03000025 0.03117263 0.02996791
# [2,]  33.7265609 0.97070543 0.32087425 0.97070543 0.32052829

#############################################################################
####################### 3.Simulations #######################################
#############################################################################
t1 <- Sys.time()
simulations.all <- fn.simulations(est.params=est.params,cond.n= c(3,5,10,30,100,500),cores=2)
t2 <- Sys.time()
t2-t1 #1.36002 hours (2 cores) 1.559195 hours (2 cores)

load(file="simulated.cpgs.list.RData")
length(cpgs.list)
load(file="simulated.models.list.RData")
length(res.list)
res.list[[1]][[1]]

#############################################################################
#################### 4. Specific data analyisis: SEX ########################
#############################################################################

sex.r <- car::recode(sex,"'Female'=1;'Male'=2")
table(sex.r)
#   1   2 
# 137 327 

table(pD.all$geo_accession==colnames(beta)) #yes
cond <- as.factor(sex.r) #la fn està preparada per a que la cond sigui 1 i 2
t1 <- Sys.time()
beta.models <- fn.models.parallel(beta, cond1=cond, cores=6)
t2 <- Sys.time()
t2-t1 #5.393927 hours (3cores)
dim(beta.models) #125950
# 

# PER ARREGLAR L'AIC de simplex inflated he creat una fn
# cond <- as.factor(sex.r) #la fn està preparada per a que la cond sigui 1 i 2
# t1 <- Sys.time()
# zoip.models <- fn.models.parallel.onlyZOIP(beta, cond1=cond, cores=2)
# t2 <- Sys.time()
# t2-t1 # 45.32981 mins (6 cores)
# zoip.models <- as.data.frame(zoip.models)
# 
# #save(beta.models,file="betadata.models.sex.RData") #he matxacat els betamodels de smoke!!!
# load(file="betadata.models.sex.RData")
# beta.models <- as.data.frame(beta.models)
# #matxambro i guardo
# beta.models$aic.sinf <- zoip.models$aic.sinf
# save(beta.models,file="betadata.models.sex.RData")
load(file="betadata.models.sex.RData")

limma.res <- apply.limma(beta,cond) 

beta.models.withlimma <- data.frame(beta.models,p.limma=limma.res$P.Value)
save(beta.models.withlimma,file="betadata.models.sex.withlimma.RData")

beta.models.adj.p <- apply(beta.models.withlimma,2,p.adjust)

#anotacions
# load(file="beta.anot.RData")
# dim(beta.anot)

# beta.models.anot <- cbind(beta.models,beta.anot[465:469])
# beta.models.adj.p.anot <- cbind(beta.models.adj.p,beta.anot[465:469])
save(beta.models.adj.p.anot,file="betadata.models.adj.p.anot.sex.RData") 

beta.models.adj.p.anot.p.s <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.s<0.05,]
beta.models.adj.p.anot.p.b <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.b<0.05,]
beta.models.adj.p.anot.p.sinf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.sinf<0.05,]
beta.models.adj.p.anot.p.binf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.binf<0.05,]
beta.models.adj.p.anot.p.n <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.n<0.05,]
beta.models.adj.p.anot.p.l <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.l<0.05,]
beta.models.adj.p.anot.p.q <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.q<0.05,]
beta.models.adj.p.anot.p.limma <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.limma<0.05,]

#puc intentar fer un plot per chr amb KaryoploteR

# # #sembla que hi ha moltes al chr X
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

res.table #molts resultats a X!!
#     Var1    s    b sinf binf    n    l    q limma
# 1   chr1  298  287  416  273  267  226  106   271
# 2  chr10  146  136  198  128  125  103   42   128
# 3  chr11  124  119  177  107  107   85   42   111
# 4  chr12  139  130  190  125  123  103   41   124
# 5  chr13   59   53   96   50   49   44   25    51
# 6  chr14   97   87  130   85   85   67   26    89
# 7  chr15   86   84  117   83   83   73   41    83
# 8  chr16  124  111  163  104  105   97   49   111
# 9  chr17  184  175  248  167  163  142   57   167
# 10 chr18   42   36   61   34   34   32   16    34
# 11 chr19  186  168  253  159  160  135   59   165
# 12  chr2  220  198  280  182  183  154   70   184
# 13 chr20   64   58   85   56   55   41   20    58
# 14 chr21   31   23   40   22   22   18    8    24
# 15 chr22   42   40   59   40   39   35   16    40
# 16  chr3  148  136  198  127  127  103   40   129
# 17  chr4  119  115  177  107  108   90   40   107
# 18  chr5  152  153  215  147  141  125   41   146
# 19  chr6  257  243  337  231  222  188   68   238
# 20  chr7  149  127  217  120  122  104   48   126
# 21  chr8   97   84  140   77   80   59   28    81
# 22  chr9   61   59   91   57   58   53   25    58
# 23  chrX 8097 8120 8113 8108 8095 4970 7927  8113
# 24  chrY  375  375  375  375  375  238  359   375

save(res.table,file="betadata.models.sex.res.table.RData")




