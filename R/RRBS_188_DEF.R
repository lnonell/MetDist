#21/4/19
#Les dades les analitzo amb RnBeads i les filtro al cluster samscratch:\Simplex\RRBS_188
#faig l'analisi completa al cluster (samscratch\Simplex\RRBS_188\reports_anal), triga unes 24h en fer import,
#QC, exploratory i anàlisi. resultats al cluster, tot i que dona un error a l'hora de fer l'analisi i no dona
#resultats de Differential methylation
#les simulacions i altres objectes els faig al cluster, especialment els de beta bin triguen més d'un dia!
# Aquests generats al cluster amb diferents scripts pq trigaven molt
# betabin.params.est.RData
# best.dist.betabin.filtered.RData
# best.dist.all.filtered.RData
# simulated.cpgs.list.RData
# simulated.models.list.RData
# 1/5/19 
# les simulated, com que a simplex inflated donaven tot NAs, 
# les he corregut aquí al portàtil eliminant el p.sinf.gamlss i canvio la fn que sumaritza en taules
# 5/5/19 filtro les que tenen poc coverage (seguint consells de l0article del mètode MACAU)
# aixo es l'utlim que faig abans de crear RRBS188_ANAL
# 12/5/19 poso els rownames dels objectes rnb.est.params, rnb.betabin.est.params i best.dist.all.betabin
# per a poder fer el subsetting des de RRBS_188_ANAL.R

workingDir<-"D:/Doctorat/Simplex/Data_OLD/RRBS_188"
setwd(workingDir)

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

#de moment ho deixo així...

#resulta que agafa patient_age com a factor...
samples$patient_age <- as.numeric(samples$patient_age)
range(samples$patient_age,na.rm=T)
# 0 62
hist(samples$patient_age,breaks=15)
hist(samples[samples$patient_sex=="f","patient_age"],breaks=15)
hist(samples[samples$patient_sex=="m","patient_age"],breaks=15)

t.test(samples[samples$patient_sex=="f","patient_age"],samples[samples$patient_sex=="m","patient_age"])
#sí que hi ha diferències...p-value = 0.0216

#########################################################################
#########0. carrego les dades #############################
#########################################################################
library(RnBeads)

# load(file="rnb.meth.f.annot.RData") #filtrades i anotades al cluster
# dim(rnb.meth.f) #3.059.804     188 inviable
# #
# annot <- rownames(rnb.meth.f)
# cpg_type <- sapply(strsplit(annot,split="_"),function(x) x[3])
# table(cpg_type)
# # Island Open Sea    Shelf    Shore
# # 470589  2478548    49220    61447
#
# #en comparacio amb WGBS hi ha moltes més islands aquí...normal pq son regions seleccionades
# #millor cobertes
#
# sel.ssi<- grep("Island|Shelf|Shore",rownames(rnb.meth.f))
# length(sel.ssi)
# rnb.meth.ff <- rnb.meth.f[sel.ssi,]
# head(rnb.meth.ff)
# #
# dim(rnb.meth.ff) #581.256 188
# #
# # #elimino les que tenen moltes NAs
# na.col <- function(x) sum(is.na(x))
# rnb.meth.na.sample <- apply(rnb.meth.ff, 2, na.col)
# range(rnb.meth.na.sample) #122647 372008
# rnb.meth.na.cpg <- apply(rnb.meth.ff[,samples$patient_sex %in% c("f","m")], 1, na.col) #hauria de fer-ho pels que tenen sex
# range(rnb.meth.na.cpg) #0 186
# rnb.meth.f <-rnb.meth.ff[rnb.meth.na.cpg<129,] #159-30!! (159 son les que tenen sex)
# dim(rnb.meth.f) # 479057

# #guardo
# save(rnb.meth.f,file="meth.betas.poquesNAs.RData") #les anoto al punt 3 per a fer l'analisi!
load(file="meth.betas.poquesNAs.RData") #AQUESTES SON LA BASE DE L'ANALISI
dim(rnb.meth.f) #479057


# #filtro les mateixes pel coverage tb
# load(file="rnb.covg.f.annot.RData") 
# dim(rnb.covg.f)
# head(rnb.covg.f)
# #matxaco objecte directament
# rnb.covg.f <- rnb.covg.f[rownames(rnb.meth.f),]
# dim(rnb.covg.f) #479057    188
# # comprovo
# identical(rownames(rnb.meth.f),rownames(rnb.covg.f)) #TRUE
# 
# save(rnb.covg.f,file="rnb.covg.poquesNAs.RData") 
load(file="rnb.covg.poquesNAs.RData") #AQUESTES SON LA BASE DE L'ANALISI beta binomial
dim(rnb.covg.f)
# 479057    188

#miro coverage
range(rnb.covg.f,na.rm=T) #0 5183
covg.m <- apply(rnb.covg.f,1,mean)
sum(is.na(covg.m)) #0
(covg.m.q <-quantile(covg.m,probs = seq(0, 1, 0.05)))
# 0%          5%         10%         15%         20%         25%         30%         35%         40% 
# 0.1808511   0.5797872   0.8351064   1.1010638   1.3829787   1.6914894   2.0212766   2.3670213   2.7340426 
# 45%         50%         55%         60%         65%         70%         75%         80%         85% 
# 3.1223404   3.5319149   3.9946809   4.5106383   5.1010638   5.8031915   6.6436170   7.7180851   9.1755319 
# 90%         95%        100% 
# 11.3297872  15.4585106 717.7500000
#és una kk!!
#filtro les que tenen menys de 3
rnb.covg.f <- rnb.covg.f[covg.m>3,]
dim(rnb.covg.f)
#270569    188
save(rnb.covg.f,file="rnb.covg.poquesNAs.covg3.RData") 

#i ara les altres
rnb.meth.f <- rnb.meth.f[covg.m>3,] 
dim(rnb.meth.f)
#270569    188
save(rnb.meth.f,file="meth.betas.poquesNAs.covg3.RData") 

#i ara all
rnb.all.f <- rnb.meth.f * rnb.covg.f

#############################################################################
####################### 1.Estimate params for each distribution #############
#############################################################################
#necessary functions
source(file=file.path("D:/Doctorat/Simplex/MetDist/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params
library(doParallel)
library(foreach)

##################### all params (but betabin, més a baix for beta-binomial)
#Estimacions amb les betes (rnb.meth) de totes les distribucions considerades
#amb els filtrats per NAs!!
n <- nrow(rnb.meth.f)
t1 <- Sys.time()
cl <- makeCluster(7,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
rnb.est.params<- foreach(i=1:n,.combine=rbind, .packages=c("fitdistrplus","VGAM","ZOIP")) %dopar% {
  #b <- rnb.meth[i,] #sobre les dades amb com a mínim 30 valors no NA, pq simulacions sortien un xurro
  b <- rnb.meth.f[i,]
  est.all.params(b)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #3.827254 hours with 7 cores

head(rnb.est.params)

rnb.est.params <-as.data.frame(rnb.est.params)
#save(rnb.est.params,file="params.est.all.RData")
load(file="params.est.all.RData")
# rownames(rnb.est.params) <- rownames(rnb.meth.f)
# save(rnb.est.params,file="params.est.all.RData")

na.col <- function(x) sum(is.na(x))
apply(rnb.est.params, 2, na.col)

apply(rnb.est.params, 2, range, na.rm=T)
#       s.mle.mu s.mle.sig s.zoip.mu s.zoip.sig    b.mom.s1    b.mom.s2  b.mle.s1  b.mle.s2  binf.mle.s1  binf.mle.s2
# [1,] 0.1103813  2.130462 0.1105003   2.130461 0.004347147 0.004347147 0.2082989 0.2084554 3.134143e-01 3.530879e-01
# [2,] 0.8886737 19.721473 0.8885252  19.698741 3.120151619 3.094272926 3.5121126 3.0024909 1.589884e+08 1.589705e+08
#       n.mom.m  n.mom.sd    n.mle.m n.mle.sd
# [1,] 0.04166667 0.2000000 0.04166667 0.197478
# [2,] 0.95833333 0.5080005 0.95833333 0.500000

#AKI pero ja he fet anàlisi sex
######################
#beta binomial

#estimació de params: de les filtrades!
N <- nrow(rnb.all.f)
t1 <- Sys.time()
cl <- makeCluster(3,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
rnb.betabin.est.params<- foreach(i=1:N,.combine=rbind, .packages=c("fitdistrplus","VGAM")) %dopar% {
  xi <- rnb.all.f[i,]
  ni <- rnb.covg.f[i,]
  est.betabin.params(xi,ni)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #5.036301 hours
head(rnb.betabin.est.params,25)

dim(rnb.betabin.est.params)
#guardo
save(rnb.betabin.est.params,file="betabin.params.est.RData")
load(file="betabin.params.est.RData")
# head(rnb.betabin.est.params)
# rownames(rnb.betabin.est.params) <- rownames(rnb.meth.f)
# save(rnb.betabin.est.params,file="betabin.params.est.RData")


##########################################################
#2. Select vector of parameters and generate simulated data (in this case we also study)
##########################################################

#enlloc de fer les simulacions per a 100 i dp per 10, creo una fn a SimulationFunctions i ho faig per varies Ns

t1 <- Sys.time()
simulations.all <- fn.simulations(est.params=rnb.est.params,cond.n= c(3,5,10,30,100,500),cores=7)
t2 <- Sys.time() #Time difference of 35.89495 mins CHECK sinf!!!
t2-t1 #

load(file="simulated.cpgs.list.RData")
length(cpgs.list) #6 un per cada cond.n
length(cpgs.list[[1]]) #3: simplex, beta, normal
dim(cpgs.list[[1]][[1]]) #2000    6
head(cpgs.list[[1]][[1]])

load(file="simulated.models.list.RData")
length(res.list) #6 un per cada cond.n
length(res.list[[1]]) #3: simplex, beta, normal
dim(res.list[[1]][[1]]) #2000    6
head(res.list[[1]][[1]]) #sembla ok

#mirem si hi ha molts NAs
table(is.na(res.list[[1]][[1]]))
table(is.na(res.list[[1]][[2]]))
table(is.na(res.list[[1]][[3]])) #6704
table(is.na(res.list[[2]][[1]]))
table(is.na(res.list[[2]][[2]]))
table(is.na(res.list[[2]][[3]])) #7777
#fins aquí hi ha uns 2000-2500 NAs però normal mogollon
table(is.na(res.list[[3]][[1]])) #675
table(is.na(res.list[[3]][[2]])) #133
table(is.na(res.list[[3]][[3]])) #6932
table(is.na(res.list[[4]][[1]])) #
table(is.na(res.list[[4]][[2]])) #
table(is.na(res.list[[4]][[3]])) #
table(is.na(res.list[[5]][[1]])) #
table(is.na(res.list[[5]][[2]])) #
table(is.na(res.list[[5]][[3]])) #
table(is.na(res.list[[6]][[1]])) #
table(is.na(res.list[[6]][[2]])) #
table(is.na(res.list[[6]][[3]])) #

#fer el mateix per a la beta binomial
rnb.betabin.est.params <- as.data.frame(rnb.betabin.est.params)
t1 <- Sys.time()
simulations.betabin.all <- fn.simu.betabin(est.params=rnb.betabin.est.params[1:10,],cond.n= c(3,5,10,30,100,500),cores=3)
t2 <- Sys.time()
t2-t1 #3.1237159 hours 3 cores (afegint gamlss, abans era poc més de la mitat)

head(simulations.betabin.all)
 
#############################################################################
########### 3. Specific data analyisis #####################################
#############################################################################
#Aplicar tots els models a les dades i veure què surt (amb quina comparació???)
# ############################
# 3.1 sexe entre les no tractades
# ############################

table(samples$sex)

#selecciono sex in "m","f" 
samples2comp <- samples[samples$patient_sex %in% c("f","m"),]
dim(samples2comp) #159

#load(file="meth.betas.poquesNAs.RData") 
#colnames(rnb.meth.f)

rnb.meth2comp <- rnb.meth.f[,samples2comp$sample_id]
dim(rnb.meth2comp)

table(colnames(rnb.meth2comp)==samples2comp$sample_id) #true!

cond <- as.factor(car::recode(samples2comp$patient_sex,"'f'=1;'m'=2"))
t1 <- Sys.time()
rnb.meth2comp.models <- fn.models.parallel(rnb.meth2comp, cond1=cond, cores=6)
t2 <- Sys.time()
t2-t1 # 8.649863 hours

head(rnb.meth2comp.models)
#save(rnb.meth2comp.models,file="rnb.meth2comp.models.RData")
load(file="rnb.meth2comp.models.RData")

limma.res <- apply.limma(rnb.meth2comp,cond)
#comprovo dims i ordres i afegeixo
head(rnb.meth2comp.models)
head(limma.res)
dim(rnb.meth2comp.models)
#[1]  479057      8
dim(limma.res)
# 479057      6

rnb.meth2comp.models.withlimma <- data.frame(rnb.meth2comp.models,p.limma=limma.res$P.Value)
save(rnb.meth2comp.models.withlimma,file="rnb.meth2comp.models.withlimma.RData")
load(file="rnb.meth2comp.models.withlimma.RData")

#RENOMBRO PER APROFITAR EL QUE HI HAVIA
rnb.meth2comp.models <- rnb.meth2comp.models.withlimma

############# regressio betabinomial
all <-rnb.all.f[,samples2comp$sample_id]
covg <-rnb.covg.f[,samples2comp$sample_id]
# remove=c(64385,257324) nose pq aquestes es salten els try i retornen un error
all[64385,] <- NA
all[257324,] <- NA



#primer miro overdispersion, he d'extraure coefs
t1 <- Sys.time()
rnb.overd<-fn.betabin.od.parallel(all=all,covg=covg, cond1=cond, cores=7)
t2 <- Sys.time()
t2-t1 # Time difference of 49.51469 mins
head(rnb.overd)



t1 <- Sys.time()
rnb.meth2comp.betabin.model <- fn.models.betabin.parallel(all=all,covg=covg, cond1=cond, cores=7)
t2 <- Sys.time()
t2-t1 # Time difference of 48.83685 mins 7 cores

head(rnb.meth2comp.betabin.models) #moltes NAs!
rnb.meth2comp.betabin.models <- as.data.frame(rnb.meth2comp.betabin.models)
table(rnb.meth2comp.betabin.models[rnb.meth2comp.betabin.models$p.bb<0.05,])

save(rnb.meth2comp.betabin.models,file="rnb.meth2comp.betabin.models.RData")
load(file="rnb.meth2comp.betabin.models.RData")

all.models.sex <- cbind(rnb.meth2comp.models.withlimma,rnb.meth2comp.betabin.models)
head(all.models.sex)

#ho guardo. IMPORTANT: ultima col és el param phi de betabinomial de la fn betabin de aod

#ajusto pvals
all.models.sex.adj.p <- apply(all.models.sex,2,p.adjust)

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
all.models.sex.adj.p.limma <- all.models.sex.adj.p[all.models.sex.adj.p[,9]<0.05 & 
                                                     !is.na(all.models.sex.adj.p[,9]),]
all.models.sex.adj.p.bb <- all.models.sex.adj.p[all.models.sex.adj.p[,10]<0.05 & 
                                                  !is.na(all.models.sex.adj.p[,10]),]

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
# 1   chr1 198 132  940   NA 158  NA  NA    NA   91
# 2  chr10  88  72  483   NA  64  NA  NA    NA   37
# 3  chr11  96  67  493   NA  89  NA  NA    NA   60
# 4  chr12  92  63  467   NA  70  NA  NA    NA   47
# 5  chr13  45  37  223   NA  32  NA  NA    NA   25
# 6  chr14  55  46  338   NA  63  NA  NA    NA   22
# 7  chr15  52  39  268   NA  35  NA  NA    NA   16
# 8  chr16 141 106  639   NA 115  NA  NA    NA  106
# 9  chr17 143 113  706   NA 115  NA  NA    NA   80
# 10 chr18  38  21  190   NA  30  NA  NA    NA   25
# 11 chr19 192 146  911   NA 152  NA  NA    NA  167
# 12  chr2 110  95  638   NA  88  NA  NA    NA   52
# 13 chr20  69  54  347   NA  57  NA  NA    NA   39
# 14 chr21  34  21  167   NA  29  NA  NA    NA   31
# 15 chr22  82  65  367   NA  69  NA  NA    NA   53
# 16  chr3  81  52  407   NA  65  NA  NA    NA   32
# 17  chr4  77  55  389   NA  58  NA  NA    NA   51
# 18  chr5  98  67  494   NA  81  NA  NA    NA   44
# 19  chr6  82  59  450   NA  56  NA  NA    NA   39
# 20  chr7 117  89  610   NA 107  NA  NA    NA   81
# 21  chr8  79  56  390   NA  54  NA  NA    NA   59
# 22  chr9 114  78  514   NA  95  NA  NA    NA   54
# 23  chrX 865 782 1852  518 658 214 309   644 9512
# 24  chrY  10   6   19    4   7   2  NA     4   32
save(res.table,file="betadata.models.sex.res.table.RData") #li poso el mateix nom que als arrays

#DSS i la comparació amb els mateixos resultats però complet està més a baix

############## CHECK a partir d'aqui

rnb.meth2comp.nas <-apply(rnb.meth2comp,1,function(x) sum(is.na(x))) 
rnb.meth2comp.models.withnas <- cbind(rnb.meth2comp.models,rnb.meth2comp.nas)
tail(rnb.meth2comp.models.withnas)

rnb.meth2comp.models.adj.p <- apply(rnb.meth2comp.models,2,p.adjust)
(rnb.meth2comp.models.adj.p.sign <- apply(rnb.meth2comp.models.adj.p,2,find.dml.n))
# p.s        p.b     p.sinf     p.binf        p.n        p.l        p.q p.s.gamlss 
# 469        171       3781          0        110          0          0        549 
head(rnb.meth2comp.models.adj.p)

genes.adj.p <- annot.rnbeads(rnb.meth2comp.models.adj.p)
length(genes.adj.p) #38003
# a veure quant son no NA
sum(!is.na(genes.adj.p)) #24895

# provo ara d incloure els promotors
genes.pr.adj.p <- annot.rnbeads(rnb.meth2comp.models.adj.p,promoters = T)
length(genes.pr.adj.p) #38003
#a veure quant son no NA
sum(!is.na(genes.pr.adj.p)) #6673, els hem d'afegir

#ho ajunto en una llista, si la pos està buida afegeixo el que hagi trobat el promotor

genes.all.adj.p <- genes.adj.p
for (i in 1:length(genes.adj.p)){
  if (is.na(genes.all.adj.p[i])) genes.all.adj.p[i] <- genes.pr.adj.p[i]
}

sum(!is.na(genes.all.adj.p)) #28392 
res.annot <- as.data.frame(rnb.meth2comp.models.adj.p[!is.na(genes.all.adj.p),])
dim(res.annot) ##28392 
res.annot$gene <- as.vector(genes.all.adj.p[!is.na(genes.all.adj.p)])
res.annot$chr <- unlist(lapply(strsplit(rownames(res.annot),split="_"), function(l) l[[1]]))
table(res.annot$chr)
# chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5 
# 1239   384  5336   537   207   379   294  1227  1010   723  6799   941   377   218   499   439   592  1035 
# chr6  chr7  chr8  chr9  chrX  chrY 
# 1301  1145   458  3030   219     3 

#nomes 219 annots al chrX!
#carrego les dades de Singmann i les comparo
load(file=file.path("D:/Doctorat/Simplex/Data/Paper_Singmann","sing.genes.chr.RData"))
sing.genes <- unique(unlist(strsplit(sing.genes.chr$UCSC_RefGene_Name,split=";")))
length(sing.genes) #679

# interseco amb el que teniem 
length(intersect(res.annot$gene,sing.genes)) # 216 nomes!!!

length(intersect(res.annot.s$gene,sing.genes)) #14
length(intersect(res.annot.b$gene,sing.genes)) #5
length(intersect(res.annot.sinf$gene,sing.genes)) #62
length(intersect(res.annot.b$gene,sing.genes)) #5
length(intersect(res.annot.binf$gene,sing.genes)) #0
length(intersect(res.annot.n$gene,sing.genes)) #5
length(intersect(res.annot.l$gene,sing.genes)) #0
length(intersect(res.annot.q$gene,sing.genes)) #0
# molt pobres les coincidències


##############################################################
#############anoto, per això he de carregar l'rnb.set
#creo la fn a R/Rnbeads.annot.fns.R

# load("RnBiseqSet.RData")
# library(RnBeads)
# library(RnBeads.hg38) #això conté les anotacions
# 
# sum(rnb.annotation.size(assembly="hg38")) #58802720
# 
# rnb.region.types(assembly="hg38")
# #"tiling"     "genes"      "promoters"  "cpgislands"
# 
# rnb.get.annotation(type = "CpG", assembly = "hg38")
# rnb.get.annotation(type = "promoters", assembly = "hg38")
# rnb.get.annotation(type = "genes", assembly = "hg38")
# rnb.get.annotation(type = "tiling", assembly = "hg38")
# rnb.get.annotation(type = "cpgislands", assembly = "hg38")
# rnb.get.annotation(type = "Open Sea", assembly = "hg38") #error
# rnb.get.annotation(type = "shore", assembly = "hg38") #error
# 
# annot.genes <- rnb.get.annotation(type = "genes", assembly = "hg38")
# #GRangeslist amb 24 GRanges, un per cada chr
# annot.genes.gr <- unlist(annot.genes)
# annot.genes.gr #60070 ranges

# library(SummarizedExperiment)
# annot.cpg <- rnb.get.annotation(type = "CpG", assembly = "hg38") 
# annot.cpg.1 <- annot.genes[[1]] #hi ha la mateixa info que a genes!!!
# The total number of dinucleotides annotated in HG19 is 28,217,009 represented 
# both on the forward and reverse DNA strands. per aixo dona tot error!!

# annot.cpg.com <- setClass(annot.cpg,"CompressedGRangesList")
# annot.cpg.com <- updateObject(annot.cpg, verbose=TRUE) #rror
# annot.cpg.gr <- unlist(annot.cpg) #error


annot <- annotation(rnb.set, type="sites")
dim(annot) #62387     8 be!
head(annot)

table(annot$`CGI Relation`)

# Open Sea    Shelf    Shore   Island 
# 38542     1735     2142    19968 

#potser ens hauríem de centrar només en Islands i shores en algunes regions com els promotors??
#com a mínim eliminar els Open Sea...

table(annot$End-annot$Start)
# 1 
# 62387 

annot.col <- paste(annot$Chromosome, annot$Start,annot$"CGI Relation",sep="_" )

#anoto models
rownames(rnb.meth2comp.models) <- annot.col
#aprofito per anotar les dades
rownames(rnb.meth) <- annot.col

rnb.meth2comp.models.adj.p <- apply(rnb.meth2comp.models,2,p.adjust)

rnb.meth2comp.models.adj.p[rnb.meth2comp.models.adj.p[,1]<0.05 & !is.na(rnb.meth2comp.models.adj.p[,1]),]
rnb.meth2comp.models.adj.p[rnb.meth2comp.models.adj.p[,2]<0.05 & !is.na(rnb.meth2comp.models.adj.p[,2]),]
rnb.meth2comp.models.adj.p[rnb.meth2comp.models.adj.p[,3]<0.05 & !is.na(rnb.meth2comp.models.adj.p[,3]),]
rnb.meth2comp.models.adj.p[rnb.meth2comp.models.adj.p[,4]<0.05 & !is.na(rnb.meth2comp.models.adj.p[,4]),]
rnb.meth2comp.models.adj.p[rnb.meth2comp.models.adj.p[,5]<0.05 & !is.na(rnb.meth2comp.models.adj.p[,5]),]
rnb.meth2comp.models.adj.p[rnb.meth2comp.models.adj.p[,6]<0.05 & !is.na(rnb.meth2comp.models.adj.p[,6]),]
rnb.meth2comp.models.adj.p[rnb.meth2comp.models.adj.p[,7]<0.05 & !is.na(rnb.meth2comp.models.adj.p[,7]),]


######################################################################
############ 4. Estimació de la millor dist per a cada cpg ###########
######################################################################
#CHECK: NO ESTÀ FET SOBRE LES FILTRADES PER NA's, tornar a fer??: sí, sobre les filtrades
source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params
#dim(rnb.meth) #62387   216

library(doParallel)
library(foreach)

N=nrow(rnb.meth.f) #amb tots els valors, sense eliminar NAs!!!
range(rnb.meth.f,na.rm=T)
#  0 1 ara sí que hi ha 0's i 1's!

#provo
rnb.meth.f <- rnb.meth.f[covg.m>3,]
e=0.001
t1 <- Sys.time()
cl <- makeCluster(6,type="PSOCK",outfile="output.txt") #poso només 1 pq està processant l'altre
registerDoParallel(cl)
best.dist.all<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- round(unlist(rnb.meth.f[i,]),4)
  #in this case we need to remove NAs
  b.prime <-b[!is.na(b)]
  b.prime <- ifelse(b.prime==1, 1-e, ifelse(b.prime==0,0+e,b.prime))
  best.dist(b.prime)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #7.55 6cores

head(best.dist.all, 10) 

#save(best.dist.all,file=file.path("best.dist.all.RData"))
#load(file=file.path("best.dist.all.RData"))
#save(best.dist.all,file=file.path("best.dist.all.filtered.RData"))
load(file=file.path("best.dist.all.filtered.RData"))


### beta-binomial: especific de NGS!!!
#en aquest cas li hem de passar els reads rnb.all FILTRATS!!!
N=nrow(rnb.meth.f)
t1 <- Sys.time()
cl <- makeCluster(3,type="PSOCK",outfile="output.txt") #poso només 1 pq està processant l'altre
registerDoParallel(cl)
best.dist.all.betabin<- foreach(i=1:N,.combine=rbind, .packages=c("VGAM","fitdistrplus")) %dopar% {
  xi <- rnb.all.f[i,]
  ni <- rnb.covg.f[i,]
  best.dist.betabin(xi,ni)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #5.076667 hours

#save(best.dist.all.betabin,file=file.path("best.dist.betabin.filtered.RData"))
load(file=file.path("best.dist.betabin.filtered.RData"))
rownames(best.dist.all.betabin) <- rownames(rnb.meth.f)
#save(best.dist.all.betabin,file=file.path("best.dist.betabin.filtered.RData"))

################################# DSS #############################
# aplico DSS, basat en beta-binomial, a veure que passa
# library(DSS)
# 
# #seleccionem les mostres
# samples2comp <- samples.new[samples.new$sex %in% c("F","M") & samples.new$treatment=="None",]
# dim(samples2comp) #162
# 
# covg.test <- rnb.covg.f[,samples2comp$sampleName]
# all.test <- rnb.all.f[,samples2comp$sampleName]
# meth.test <- rnb.meth.f[,samples2comp$sampleName]
# 
# #no hi pot haver NAs a M!!
# 
# all.test.nona <- na.omit(all.test)
# covg.test.nona <- na.omit(covg.test)
# meth.test.nona <- meth.test[rownames(all.test.nona),]
# dim(all.test.nona) # 11125
# dim(meth.test.nona) # 11125
# dim(covg.test.nona) # 38003
# 
# all.equal(rownames(all.test.nona),rownames(covg.test.nona))
# #i ara treiem els que són NA de all a covg
# 
# covg.test.nona <- covg.test.nona[rownames(all.test.nona),]
# dim(covg.test.nona) # 11125
# 
# annot <- strsplit(rownames(covg.test.nona),split="_")
# chr <-unlist(lapply(annot, function(l) l[[1]]))
# pos <-unlist(lapply(annot, function(l) l[[2]]))
# fn <-unlist(lapply(annot, function(l) l[[3]]))
# 
# table(chr)
# # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5 
# # 541   203  1873   238    82   174   136   491   425   264  2137   514   219    90   214   186   282   521 
# # chr6  chr7  chr8  chr9  chrX 
# # 556   525   243  1103   108
# #n'hi ha poquets a X
# 
# 
# BSseq.obj <- BSseq(chr = chr, pos = as.numeric(pos),
#                    M = all.test.nona,
#                    Cov = covg.test.nona,
#                    sampleNames = colnames(all.test.nona))
# 
# dmlTest <- DMLtest(BSseq.obj, group1=colnames(all.test.nona)[samples2comp$sex=="F"], 
#                    group2=colnames(all.test.nona)[samples2comp$sex=="M"])
# Sys.time() 
# save(dmlTest,file="DSS.dmlTest.RData")
# #load(file="DSS.dmlTest.RData")
# 
# 
# #seleccionem les significatives
# dmlTest.sig <- dmlTest[dmlTest$fdr<0.05,]
# dim(dmlTest.sig) #67 11
# #poquets i no n'hi ha cap a X!
# 
# #################################
# #torno a fer les anàlisis d'abans de sexe, només amb els complete
# rnb.meth2comp <- rnb.meth.f[,samples2comp$sampleName]
# dim(rnb.meth2comp) # 38003   162
# 
# table(colnames(meth.test.nona)==samples2comp$sampleName) #true!
# 
# cond <- as.factor(car::recode(samples2comp$sex,"'F'=1;'M'=2"))
# 
# t1 <- Sys.time()
# meth.test.nona.models <- fn.models.parallel(meth.test.nona, cond1=cond, cores=6)
# t2 <- Sys.time()
# t2-t1 # 11.02239 mins (6 cores)
# 
# limma.nona.res <- apply.limma(meth.test.nona,cond)
# 
# meth.test.nona.models.withlimma <- data.frame(meth.test.nona.models,p.limma=limma.nona.res$P.Value)
# 
# save(meth.test.nona.models.withlimma,file="meth.test.nona.models.sex.RData")
# #load(file="meth.test.nona.models.sex.RData")
# 
# ############ regressio betabinomial
# all <-all.test.nona
# covg <-covg.test.nona
# 
# t1 <- Sys.time()
# meth.test.nona.betabin.models <- fn.models.betabin.parallel(all=all.test.nona,covg=covg.test.nona, cond1=cond, cores=3)
# t2 <- Sys.time()
# t2-t1 # 3.222671 hours
# 
# #no les matxambro, les guardo directament,
# save(meth.test.nona.betabin.models,file="meth.test.nona.complete.betabin.models.sex.RData")
# load(file="meth.test.nona.complete.betabin.models.sex.RData")
# 
# 
# limma.nona.res <- apply.limma(meth.test.nona,cond)
# meth.test.nona.models.withlimma <- data.frame(meth.test.nona.models,p.limma=limma.nona.res$P.Value)
# 
# save(meth.test.nona.models.withlimma,file="meth.test.nona.complete.models.sex.RData")
# load(file="meth.test.nona.complete.models.sex.RData")
# 
# 
# ########### comparo dss amb la resta
# meth.test.nona.models.withlimma.adj.p <- apply(meth.test.nona.models.withlimma,2,p.adjust)
# 
# #comprovo ordre
# all.equal(rownames(meth.test.nona.models.withlimma.adj.p),rownames(dmlTest))
# #no
# dmlTest.s <- dmlTest[rownames(meth.test.nona.models.withlimma.adj.p),]
# all.equal(rownames(meth.test.nona.models.withlimma.adj.p),rownames(dmlTest.s)) #TRUE
# 
# 
# plot(meth.test.nona.models.withlimma$p.s,dmlTest.s$fdr) #buff
# cor.test(meth.test.nona.models.withlimma$p.s,dmlTest.s$fdr) #0.43
# plot(meth.test.nona.models.withlimma$p.s,dmlTest.s$fdr,xlim=c(0,0.1),ylim=c(0,0.1)) 
# 
# cor.test(meth.test.nona.models.withlimma$p.b,dmlTest.s$fdr) #0.44l...
# plot(meth.test.nona.models.withlimma$p.b,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.sinf,dmlTest.s$fdr) #0.68
# plot(meth.test.nona.models.withlimma$p.sinf,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.binf,dmlTest.s$fdr) #0.32
# plot(meth.test.nona.models.withlimma$p.binf,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.n,dmlTest.s$fdr) #0.68
# plot(meth.test.nona.models.withlimma$p.n,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.l,dmlTest.s$fdr) #0.55
# plot(meth.test.nona.models.withlimma$p.l,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.q,dmlTest.s$fdr) #0.67
# plot(meth.test.nona.models.withlimma$p.q,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.limma,dmlTest.s$fdr) #0.23
# plot(meth.test.nona.models.withlimma$p.limma,dmlTest.s$fdr) #bu
# #FATAL LES CORRELS!!
# 
# ############# regressio betabinomial
# all <-all.test.nona
# covg <-covg.test.nona
# 
# #AKI
# t1 <- Sys.time()
# meth.test.nona.betabin.models <- fn.models.betabin.parallel(all=all.test.nona,covg=covg.test.nona, cond1=cond, cores=6)
# t2 <- Sys.time()
# t2-t1 # 10.08667 mins
# 
# meth.test.nona.betabin.models<-as.data.frame(meth.test.nona.betabin.models)
# #save(meth.test.nona.betabin.models,file="meth.test.nona.betabin.models.sex.RData")
# load(file="meth.test.nona.betabin.models.sex.RData")
# 
# meth.test.nona.betabin.models <- as.data.frame(meth.test.nona.betabin.models)
# #no tinc molt clar això de l'overdispersion
# dim(meth.test.nona.betabin.models[meth.test.nona.betabin.models$p.bb.od<0.05,])
# #totes signif, vol dir que hi ha overdispersion a tot arreu o que no n'hi ha o què??
# 
# identical(rownames(meth.test.nona.models.withlimma),rownames(meth.test.nona.betabin.models))
# 
# all.models.sex <- cbind(meth.test.nona.models.withlimma,meth.test.nona.betabin.models)
# head(all.models.sex)
# 
# #ho guardo. IMPORTANT: ultima col és el param phi de betabinomial de la fn betabin de aod
# # save(all.models.sex,file="all.models.sex.RData")
# # plot(all.models.sex$p.s,all.models.sex$p.b)
# # cor.test(all.models.sex$p.s,all.models.sex$p.b) 
# # cor.test(all.models.sex$p.s,all.models.sex$p.bb) 
# #spearman surt una miqueta més...
# 
# #ajusto pvals
# all.models.sex.adj.p <- apply(all.models.sex,2,p.adjust)
# 
# all.models.sex.adj.p.s <- all.models.sex.adj.p[all.models.sex.adj.p[,1]<0.05 & 
#                                                  !is.na(all.models.sex.adj.p[,1]),]
# all.models.sex.adj.p.b <- all.models.sex.adj.p[all.models.sex.adj.p[,2]<0.05 & 
#                                                  !is.na(all.models.sex.adj.p[,2]),]
# all.models.sex.adj.p.sinf <- all.models.sex.adj.p[all.models.sex.adj.p[,3]<0.05 & 
#                                                     !is.na(all.models.sex.adj.p[,3]),]
# all.models.sex.adj.p.binf <- all.models.sex.adj.p[all.models.sex.adj.p[,4]<0.05 & 
#                                                     !is.na(all.models.sex.adj.p[,4]),]
# all.models.sex.adj.p.n <- all.models.sex.adj.p[all.models.sex.adj.p[,5]<0.05 & 
#                                                  !is.na(all.models.sex.adj.p[,5]),]
# all.models.sex.adj.p.l <- all.models.sex.adj.p[all.models.sex.adj.p[,6]<0.05 & 
#                                                  !is.na(all.models.sex.adj.p[,6]),]
# all.models.sex.adj.p.q <- all.models.sex.adj.p[all.models.sex.adj.p[,7]<0.05 & 
#                                                  !is.na(all.models.sex.adj.p[,7]),]
# all.models.sex.adj.p.limma <- all.models.sex.adj.p[all.models.sex.adj.p[,9]<0.05 & 
#                                                      !is.na(all.models.sex.adj.p[,9]),]
# all.models.sex.adj.p.bb <- all.models.sex.adj.p[all.models.sex.adj.p[,10]<0.05 & 
#                                                   !is.na(all.models.sex.adj.p[,10]),]
# 
# 
# #extrec els resultats (fn més a munt)
# res.chr.s <- res.chr(res=all.models.sex.adj.p)
# res.chr.b <- res.chr(res=all.models.sex.adj.p.b)
# res.chr.sinf <- res.chr(res=all.models.sex.adj.p.sinf)
# res.chr.binf <- res.chr(res=all.models.sex.adj.p.binf)
# res.chr.n <- res.chr(res=all.models.sex.adj.p.n)
# res.chr.l <- res.chr(res=all.models.sex.adj.p.l)
# res.chr.q <- res.chr(res=all.models.sex.adj.p.q)
# res.chr.limma <- res.chr(res=all.models.sex.adj.p.limma)
# res.chr.bb <- res.chr(res=all.models.sex.adj.p.bb)
# 
# names(res.chr.s)[2] <- "s"
# names(res.chr.b)[2] <- "b"  
# names(res.chr.sinf)[2] <- "sinf"  
# names(res.chr.binf)[2] <- "binf"  
# names(res.chr.n)[2] <- "n"  
# names(res.chr.l)[2] <- "l"  
# names(res.chr.q)[2] <- "q"  
# names(res.chr.limma)[2] <- "limma"  
# names(res.chr.bb)[2] <- "bb"  
# require(plyr)
# res.table <- join_all(list(res.chr.s,
#                            res.chr.b,
#                            res.chr.sinf,
#                            res.chr.binf,
#                            res.chr.n,
#                            res.chr.l,
#                            res.chr.q,
#                            res.chr.limma,
#                            res.chr.bb), by = 'chr', type = 'full')
# 
# res.table #NO HI HA QUASI RESULTATS: KK
