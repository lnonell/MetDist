#18/1/19
# a partir de les dades processades al cluster (samscratch/WGBS_81_Blueprint) amb RnBeads i posteriorment 
# partides en 4
# al cluster hi ha els scripts en què es processen les dades i després es particiona, els .R i .sh
# anotades i filtrades (cadascuna de les 4 parts) i 
# ajuntades: D:\Doctorat\Simplex\Data\WGBS_81_Bueprint\rnb.meth.f.annot
# basat en l'script de RRBS_216
# faig el filtrar del coverage en base al filtrat del rnb.meth.f!

#revisem els samples.csv
workingDir<-"D:/Doctorat/Simplex/MetDist/Data/WGBS_81_Bueprint"
setwd(workingDir)

#dades originals, dp de filtrar al cluster
# load(file="rnb.meth.f.annot.RData")
# dim(rnb.meth.f)
#2.905.277      81 moltes més files que la resta, a veure si puc fer l'anàlisi!

#load(file="rnb.covg.RData")
samples <- read.delim("D:/Doctorat/Simplex/MetDist/Data/WGBS_81_Bueprint/samples.tsv", stringsAsFactors=FALSE)
dim(samples) # 81 52

# str(samples)
# table(samples$sampleGroup)
# # Bcell     cancer         DC       eryt       gran       megK         Mf       mono         NK      other 
# # 10         10          2          2          7          1         19          8          2          4 
# # plasma progenitor      Tcell 
# # 2          1         13 
# 
# table(samples$cellTypeShort)
# # AML_BM    AML_PBMC    Bcell_gc   Bcell_mem Bcell_naive   Bcell_pre          DC   Endo_prol   Endo_rest 
# # 3           3           1           3           5           1           2           2           2 
# # Eosi_mat       Eryth   Leuk_myel      Mac_f0      Mac_f1      Mac_f2          MK        Mono         MPP 
# # 1           2           4           8           6           5           1           8           1 
# # Neut_mat          NK        Plas        TCD4   TCD4_cenM   TCD4_effM        TCD8   TCD8_cenM   TCD8_effM 
# # 6           2           2           3           1           1           4           1           1 
# # TCD8_term        Treg 
# # 1           1 
# 
# table(samples$DISEASE)
# table(samples$BIOMATERIAL_TYPE)
# # Primary Cell Primary Cell Culture 
# # 52                   29 
# table(samples$CELL_LINE)
# table(samples$CELL_TYPE) #molts
# table(samples$TISSUE_TYPE)
# Bone marrow   Cord blood       Tonsil Venous blood 
#       11           22            1           47 

table(samples$DONOR_SEX)
#     Female   Male 
# 4     47     30 
#############################################################################
####################### 0.Filtro, elimino els open sea ######################
#############################################################################
# annot <- rownames(rnb.meth.f)
# cpg_type <- sapply(strsplit(annot,split="_"),function(x) x[3])
# table(cpg_type)
# # Island Open Sea    Shelf    Shore 
# # 179519  2649162    31513    45083 

#select shelf shore or island
# sel.ssi<- grep("Island|Shelf|Shore",annot)
# length(sel.ssi) #256115
# rnb.meth.ff <- rnb.meth.f[sel.ssi,]
# head(rnb.meth.ff)
# 
# #renombo l'objecte per a poder aprofitar el codi
# rnb.meth <- rnb.meth.ff
# dim(rnb.meth) #256115     81
# #save(rnb.meth,file="meth.betas.filtered.RData")
# 
# #DATA
# load(file="meth.betas.filtered.RData")
# dim(rnb.meth)  #256115     81 si!
# 
# #hi ha molts NAs, això s'hauria de reportar
# na.col <- function(x) sum(is.na(x))
# rnb.meth.na.sample <- apply(rnb.meth, 2, na.col)
# range(rnb.meth.na.sample) #29283 82867
# 
# rnb.meth.na.cpg <- apply(rnb.meth, 1, na.col)
# range(rnb.meth.na.cpg) #0 79
# 
# hist(rnb.meth.na.cpg,breaks=100)
# #hauria de filtrar aquelles CpGs que tenen pocs valors
# #miro aquells que tenen més de 30
# #rnb.meth.f <-rnb.meth[rnb.meth.na.cpg<30,] #MALAMENT!!
# rnb.meth.f <-rnb.meth[rnb.meth.na.cpg<52,] #81-30!!
# dim(rnb.meth.f) #213748
# #comprovo
# range(apply(rnb.meth.f, 1, na.col))
# table(rnb.meth.f[1,])

#save(rnb.meth.f,file="meth.betas.poquesNAs.RData") #aquestes estan anotades!
load(file="meth.betas.poquesNAs.RData")
dim(rnb.meth.f) #213748     81

# sel.rows <- rownames(rnb.meth)
# #a partir d'aquí seleccionaré les files del covg
# 
# load(file="rnb.covg1.annot.RData")
# rnb.covg1.f <- rnb.covg1[intersect(rownames(rnb.covg1),sel.rows),]
# head(rnb.covg1.f)
# rm(rnb.covg1)
# load(file="rnb.covg2.annot.RData")
# rnb.covg2.f <- rnb.covg2[intersect(rownames(rnb.covg2),sel.rows),]
# head(rnb.covg2.f)
# rm(rnb.covg2)
# load(file="rnb.covg3.annot.RData")
# rnb.covg3.f <- rnb.covg3[intersect(rownames(rnb.covg3),sel.rows),]
# head(rnb.covg3.f)
# rm(rnb.covg3)
# load(file="rnb.covg4.annot.RData")
# rnb.covg4.f <- rnb.covg4[intersect(rownames(rnb.covg4),sel.rows),]
# 
# head(rnb.covg4.f)
# rm(rnb.covg4)

# rnb.covg.f <- rbind(rnb.covg1.f,
#                   rnb.covg2.f,
#                   rnb.covg3.f,
#                   rnb.covg4.f)
# 
# table(rownames(rnb.covg)==rownames(rnb.meth)) #si!
# 
# #save(rnb.covg.f, file="rnb.covg.filtered.RData")
# load(file="rnb.covg.filtered.RData")
# 
# #renombro tb per a que sigui tot igual
# rnb.covg <- rnb.covg.f
# rm(rnb.covg.f)
# dim(rnb.covg) #256115     81

#i ara filtro de nou per a eliminar les NAs que ja he eliminat a rnb.meth
#filtro les mateixes pel coverage tb
# rnb.covg.f <- rnb.covg[rnb.meth.na.cpg<52,]
# dim(rnb.covg.f) #213748     81 yes!
#save(rnb.covg.f,file="rnb.covg.poquesNAs.RData") 
load(file="rnb.covg.poquesNAs.RData")
dim(rnb.covg.f) #213748     81
#i ara all
rnb.all.f <- rnb.meth.f * rnb.covg.f

table(rownames(rnb.meth.f)==rownames(rnb.covg.f)) #just in case:TRUE
#############################################################################
####################### 1.Estimate params for each distribution #############
#############################################################################
#necessary functions
source(file=file.path("D:/Doctorat/Simplex/MetDist/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params
library(doParallel)
library(foreach)

#ho faig sobre tots els CpGs, que són 62387
#rnb.all <- rnb.meth * rnb.covg
# #ccomprovo
# dim(rnb.all)
# rnb.all[1:5,1:5]
# range(rnb.all[1,],na.rm=T) #coincideix amb rnb, deu estar bé ;-)

##################### all params (but betabin, més a baix for beta-binomial)
#Estimacions amb les betes (rnb.meth) de totes les distribucions considerades
#fn def amb tot
n <- nrow(rnb.meth.f)
t1 <- Sys.time()
cl <- makeCluster(4,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
rnb.est.params<- foreach(i=1:n,.combine=rbind, .packages=c("fitdistrplus","VGAM","ZOIP")) %dopar% {
  #b <- rnb.meth[i,]
  b <- rnb.meth.f[i,] #filtrant els que tenen menys de 30
  est.all.params(b)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #1.24 hours 6 cores #2.69 hours 4cores

head(rnb.est.params)

rnb.est.params <-as.data.frame(rnb.est.params)
# save(rnb.est.params,file="params.est.all.RData")
# load(file="params.est.all.RData") ##renombro withNAs!
#save(rnb.est.params,file="params.est.all.020319.RData") #sense NAs!
#save(rnb.est.params,file="params.est.all.050319.RData") #canviant e de 0.001 a 0.01
save(rnb.est.params,file="params.est.all.090319.RData") #estaven malament agafades les NAs!
load(file="params.est.all.RData")

na.col <- function(x) sum(is.na(x))
apply(rnb.est.params, 2, na.col)
# s.mle.mu   s.mle.sig   s.zoip.mu  s.zoip.sig    b.mom.s1    b.mom.s2    b.mle.s1    b.mle.s2 binf.mle.s1 
# 0           0           8           8           0           0           1           1      168700 
# binf.mle.s2     n.mom.m    n.mom.sd     n.mle.m    n.mle.sd 
# 168700           0           0           0           0 

apply(rnb.est.params, 2, range, na.rm=T)
# s.mle.mu s.mle.sig  s.zoip.mu s.zoip.sig  b.mom.s1  b.mom.s2  b.mle.s1  b.mle.s2  binf.mle.s1  binf.mle.s2
# [1,] 0.04960212  1.218596 0.04956221   1.218596 -0.248997 -0.248997 0.1289177 0.1289263 2.845620e-01 3.114245e-01
# [2,] 0.94967629 63.273728 0.95043659  63.150645  3.052946  3.054790 6.8885565 7.2263283 8.239343e+06 4.122168e+06
# n.mom.m  n.mom.sd    n.mle.m  n.mle.sd
# [1,] 0.04166667 0.2000001 0.04166667 0.1416667
# [2,] 0.95833333 0.7071068 0.95833333 0.5000000

length(which(rnb.est.params$s.zoip.sig==0.0000)) # 0
length(which(rnb.est.params$s.zoip.sig==0)) # 0

##########################################################
#2. Select vector of parameters and generate simulated data (in this case we also study)
##########################################################
#enlloc de fer les simulacions per a 100 i dp per 10, creo una fn a SimulationFunctions i ho faig per varies Ns
cond.n= c(3,5,10,30,100,500)
t1 <- Sys.time()
simulations.all <- fn.simulations(est.params=rnb.est.params,cond.n= cond.n,cores=4)
t2 <- Sys.time()
t2-t1 #1.81804 hours 4 cores amb gamlss.simplex

load(file="simulated.cpgs.list.RData")
length(cpgs.list) #6 un per cada cond.n
length(cpgs.list[[1]]) #3: simplex, beta, normal
dim(cpgs.list[[1]][[1]]) #2000    6
head(cpgs.list[[1]][[1]])
head(cpgs.list[[3]][[1]])
hist(cpgs.list[[3]][[1]][1,1:10],breaks=10)
hist(cpgs.list[[3]][[1]][1,11:20],breaks=10)
hist(cpgs.list[[4]][[1]][1,11:20],breaks=10)
hist(cpgs.list[[4]][[1]][1,1:10],breaks=10)
#si que sembla que hi ha diferències...
head(sort(apply(cpgs.list[[5]][[2]],1,min,na.rm=T))) #totes les dades normals tenen valors negatius i més grans d'u
tail(sort(apply(cpgs.list[[5]][[2]],1,max,na.rm=T)))

load(file="simulated.models.list.RData")
length(res.list) #6 un per cada cond.n
length(res.list[[1]]) #3: simplex, beta, normal
dim(res.list[[1]][[1]]) #2000    6
head(res.list[[1]][[1]]) #sembla ok
head(res.list[[3]][[1]])
head(res.list[[3]][[3]]) #sembla que si és normal simplex i beta no són capaços d'ajustar

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
table(is.na(res.list[[4]][[1]])) #756
table(is.na(res.list[[4]][[2]])) #99
table(is.na(res.list[[4]][[3]])) #7808
table(is.na(res.list[[5]][[1]])) #705
table(is.na(res.list[[5]][[2]])) #164
table(is.na(res.list[[5]][[3]])) #7992
table(is.na(res.list[[6]][[1]])) #864
table(is.na(res.list[[6]][[2]])) #382
table(is.na(res.list[[6]][[3]])) #8000 la mitat!

######################
#beta binomial

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
t2-t1 #11.30993 hours
head(rnb.betabin.est.params,25)
head(rnb.betabin.est.params,25) #molts negatius amb mom i les altres no se si me les puc creure..
range(rnb.betabin.est.params[,3],na.rm=T) #2.002180e+00 5.878612e+11
hist(rnb.betabin.est.params[,3],breaks=20)
range(rnb.betabin.est.params[,4],na.rm=T) #3.343218e-02 1.278933e+08
hist(rnb.betabin.est.params[,4],breaks=20)
range(rnb.betabin.est.params[,5],na.rm=T) #3.421720e-01 3.633867e+11
hist(rnb.betabin.est.params[,5],breaks=20)
#tot i que hi ha valors extrems, la majoria son petits
#la veritat és que són molt estranys els valors, no sé si ens en podem refiar...

dim(rnb.betabin.est.params)
#guardo
save(rnb.betabin.est.params,file="betabin.params.est.RData")

 
#############################################################################
########### 3. Specific data analyisis: sexe  #########
#############################################################################

#selecciono sex 
samples2comp <- samples[samples$DONOR_SEX %in% c("Female","Male"),]
dim(samples2comp) #77

load(file="meth.betas.poquesNAs.RData")
dim(rnb.meth.f) #213748     81

rnb.meth2comp <- rnb.meth.f[,samples2comp$sampleName]
dim(rnb.meth2comp) # 213748     77 ok!

table(colnames(rnb.meth2comp)==samples2comp$sampleName) #true!

cond <- as.factor(car::recode(samples2comp$DONOR_SEX,"'Female'=1;'Male'=2"))
# t1 <- Sys.time()
# rnb.meth2comp.models <- fn.models(rnb.meth2comp, cond1=cond)
# t2 <- Sys.time()
# 
# cond <- as.factor(car::recode(samples2comp$sex,"'F'=1;'M'=2"))
t1 <- Sys.time()
rnb.meth2comp.models <- fn.models.parallel(rnb.meth2comp, cond1=cond, cores=6)
t2 <- Sys.time()
t2-t1 # 3.016305 hours (6 cores)

# save(rnb.meth2comp.models,file="rnb.meth2comp.models.RData")
# load(file="rnb.meth2comp.models.RData")

limma.res <- apply.limma(rnb.meth2comp,cond)
#comprovo dims i ordres i afegeixo
head(rnb.meth2comp.models)
head(limma.res)
dim(rnb.meth2comp.models)
#[1]  213748     8
dim(limma.res)
#[1]  213748    6

rnb.meth2comp.models.withlimma <- data.frame(rnb.meth2comp.models,p.limma=limma.res$P.Value)
#save(rnb.meth2comp.models.withlimma,file="rnb.meth2comp.models.withlimma.RData")
load(file="rnb.meth2comp.models.withlimma.RData")

#RENOMBRO PER APROFITAR EL QUE HI HAVIA
rnb.meth2comp.models <- rnb.meth2comp.models.withlimma

############# regressio betabinomial
all <-rnb.all.f[,samples2comp$sampleName]
covg <-rnb.covg.f[,samples2comp$sampleName]

t1 <- Sys.time()
rnb.meth2comp.betabin.models <- fn.models.betabin.parallel(all=all,covg=covg, cond1=cond, cores=6)
t2 <- Sys.time()
t2-t1 # Time difference of 25.00339 mins
head(rnb.meth2comp.betabin.models)
rnb.meth2comp.betabin.models <- as.data.frame(rnb.meth2comp.betabin.models)
table(rnb.meth2comp.betabin.models[rnb.meth2comp.betabin.models$p.bb<0.05,])

save(rnb.meth2comp.betabin.models,file="rnb.meth2comp.betabin.models.RData")

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

#extrec els resultats (fn més a munt)
res.chr.s <- res.chr(res=all.models.sex.adj.p)
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
res.table <- join_all(list(res.chr.s,
                           res.chr.b,
                           res.chr.sinf,
                           res.chr.binf,
                           res.chr.n,
                           res.chr.l,
                           res.chr.q,
                           res.chr.limma,
                           res.chr.bb), by = 'chr', type = 'full')

res.table #molt bé la taula, bb mogollon i només a X
#     chr     s    b  sinf binf    n  l    q limma   bb
# 1   chr1 15227    2  2771   NA    7 NA   NA    NA   NA
# 2  chr10  8778    1  2128   NA   NA NA   NA    NA   NA
# 3  chr11  9627    2  1848   NA    1 NA    1    NA   NA
# 4  chr12  7734    2  1786    1    2 NA    1     2   NA
# 5  chr13  4856    1  1576   NA   NA NA   NA    NA   NA
# 6  chr14  6365    2  1848   NA   NA NA   NA    NA   NA
# 7  chr15  5519    2  1576   NA    2 NA   NA    NA   NA
# 8  chr16  9807    5  1315   NA    4 NA   NA    NA   NA
# 9  chr17 12501    1  1569   NA    1 NA   NA    NA   NA
# 10 chr18  4751   NA  1445   NA   NA NA   NA    NA   NA
# 11 chr19 16692    3  2076   NA    4 NA    4    NA   NA
# 12  chr2 13298    1  3326   NA    1 NA   NA    NA   NA
# 13 chr20  6203    2  1245   NA    1 NA   NA    NA    1
# 14 chr21  4026   NA   345   NA   NA NA   NA    NA   NA
# 15 chr22  5814    1   642   NA   NA NA    1    NA   NA
# 16  chr3  6685   NA  1790   NA    2 NA   NA    NA   NA
# 17  chr4  7414    1  2436   NA    1 NA   NA    NA   NA
# 18  chr5 12476    2  3317   NA    1 NA   NA    NA   NA
# 19  chr6  7782    1  2184   NA   NA NA   NA    NA   NA
# 20  chr7 10616    4  2411   NA   NA NA   NA    NA   NA
# 21  chr8  7729    1  1817   NA    1 NA   NA    NA   NA
# 22  chr9  9383    3  1601   NA    2 NA    1    NA   NA
# 23  chrX 20035 7798 12683 5661 5749  8 2025  6964 7091
# 24  chrY   430   15    25    5    7 NA    1     3   NA

save(res.table,file="betadata.models.sex.res.table.RData") #li poso el mateix nom que als arrays

#DSS i la comparació amb els mateixos resultats però complet està més a baix
######################################################################
############ 4. Estimació de la millor dist per a cada cpg ###########
######################################################################
#CHECK: NO ESTÀ FET SOBRE LES FILTRADES PER NA's, tornar a fer: si sobre els filtrats!!! 
source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params
dim(rnb.meth) #62387   216

library(doParallel)
library(foreach)

N=nrow(rnb.meth.f)
range(rnb.meth.f,na.rm=T)
#  0 1 ara sí que hi ha 0's i 1's!
e=0.001
t1 <- Sys.time()
cl <- makeCluster(3,type="PSOCK",outfile="output.txt") #poso només 1 pq està processant l'altre
registerDoParallel(cl)
best.dist.all<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- round(unlist(rnb.meth.f[i,]),6)
  #in this case we need to remove NAs
  b.prime <-b[!is.na(b)]
  b.prime <- ifelse(b.prime==1, 1-e, ifelse(b.prime==0,0+e,b.prime))
  best.dist(b.prime)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #2.494726 hours 3cores

head(best.dist.all, 10) 

# save(best.dist.all,file=file.path("best.dist.all.RData"))
# load(file=file.path("best.dist.all.RData"))
#save(best.dist.all,file=file.path("best.dist.all.filtered.RData"))
load(file=file.path("best.dist.all.filtered.RData"))


# WGBS81.best.dist <- t(apply(best.dist.all,1,best.aic.ks))
# 
# table(WGBS81.best.dist[,1]) #mirar què passa, sembla que no fa be best.aic.ks!!! Que hi havia resultats que eren NAs i 
# #eren NAs tota la fila i 
# # beta.aic  normal.aic simplex.aic 
# # 198014       14938       43163 
# table(WGBS81.best.dist[,2])
# beta.ks.p  normal.ks.p simplex.ks.p 
# 126144       105801        24170 

### beta-binomial: especific de NGS!!!
#en aquest cas li hem de passar els reads rnb.all
N=nrow(rnb.meth.f)
t1 <- Sys.time()
cl <- makeCluster(6,type="PSOCK",outfile="output.txt") #poso només 1 pq està processant l'altre
registerDoParallel(cl)
best.dist.all.betabin<- foreach(i=1:N,.combine=rbind, .packages=c("VGAM","fitdistrplus")) %dopar% {
  xi <- rnb.all.f[i,]
  ni <- rnb.covg.f[i,]
  best.dist.betabin(xi,ni)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #7.573625 hours 6 cores

head(best.dist.all.betabin) #moltes NAs
save(best.dist.all.betabin,file=file.path("best.dist.betabin.filtered.RData"))
#load(file=file.path("best.dist.betabin.filtered.RData"))


# #les afegeixo a les best.dist.all, serà un data.frame: VELL
# best.dist.with.betabin <- data.frame(best.dist.all[,1:3],betabin.aic=best.dist.all.betabin[,1],
#                                      best.dist.all[,4:6],betabin.ks.p=best.dist.all.betabin[,2])
# save(best.dist.with.betabin,file=file.path("best.dist.with.betabin.RData"))
# # load(file=file.path("best.dist.with.betabin.RData"))

###################################################################
################################# DSS #############################
# aplico DSS, basat en beta-binomial, a veure que passa
# library(DSS)
# 
# #seleccionem les mostres
# samples2comp <- samples[samples$DONOR_SEX %in% c("Female","Male"),]
# dim(samples2comp) #77
# 
# covg.test <- rnb.covg.f[,samples2comp$sampleName]
# all.test <- rnb.all.f[,samples2comp$sampleName]
# meth.test <- rnb.meth.f[,samples2comp$sampleName]
# 
# #no hi pot haver NAs a M!!
# all.test.nona <- na.omit(all.test)
# covg.test.nona <- na.omit(covg.test)
# meth.test.nona <- meth.test[rownames(all.test.nona),]
# dim(all.test.nona) #102426
# dim(covg.test.nona) #213748
# dim(meth.test.nona) #102426
# 
# all.equal(rownames(all.test.nona),rownames(covg.test.nona))
# #i ara treiem els que són NA de all a covg
# 
# covg.test.nona <- covg.test.nona[rownames(all.test.nona),]
# dim(covg.test.nona) #102426 
# 
# #annot
# annot <- strsplit(rownames(covg.test.nona),split="_")
# chr <-unlist(lapply(annot, function(l) l[[1]]))
# pos <-unlist(lapply(annot, function(l) l[[2]]))
# fn <-unlist(lapply(annot, function(l) l[[3]]))
# 
# BSseq.obj <- BSseq(chr = chr, pos = as.numeric(pos),
#                   M = all.test.nona,
#                   Cov = covg.test.nona,
#                   sampleNames = colnames(all.test.nona))
# 
# dmlTest <- DMLtest(BSseq.obj, group1=colnames(all.test.nona)[samples2comp$DONOR_SEX=="Female"], 
#                    group2=colnames(all.test.nona)[samples2comp$DONOR_SEX=="Male"])
# Sys.time() #1-2 minuts!!
# 
# dim(dmlTest) #102426
# save(dmlTest,file="DSS.dmlTest.RData")
# load(file="DSS.dmlTest.RData")
# 
# 
# #seleccionem les significatives
# dmlTest.sig <- dmlTest[dmlTest$fdr<0.05,]
# dim(dmlTest.sig) #43263 
# 
# table(dmlTest.sig$chr)
# # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3 
# # 1250  2012  1589  1416  1256  1634  1345  1246  1682  1144  2119  2889  1143   347   755  1546 
# # chr4  chr5  chr6  chr7  chr8  chr9  chrX 
# # 2075  3070  1732  1990  1657  1450  7916 
# 
# table(dmlTest$chr)
# # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3 
# # 4614  4594  4927  3699  2429  3366  2799  4558  6384  2433  9023  6606  3104  1190  2968  3333 
# # chr4  chr5  chr6  chr7  chr8  chr9  chrX 
# # 3606  6712  3831  5243  4193  4472  8342 
# 
# #n'hi ha moltes significatives a X en comparació amb els totals, 
# table(dmlTest.sig$chr)/table(dmlTest$chr)
# 
# # chr1     chr10     chr11     chr12     chr13     chr14     chr15     chr16     chr17 
# # 0.2709146 0.4379626 0.3225086 0.3828062 0.5170852 0.4854427 0.4805288 0.2733655 0.2634712 
# # chr18     chr19      chr2     chr20     chr21     chr22      chr3      chr4      chr5 
# # 0.4702014 0.2348443 0.4373297 0.3682345 0.2915966 0.2543801 0.4638464 0.5754298 0.4573897 
# # chr6      chr7      chr8      chr9      chrX 
# # 0.4521013 0.3795537 0.3951824 0.3242397 0.9489331 
# 
# #################################
# #torno a fer les anàlisis d'abans de sexe, només amb els complete
# rnb.meth2comp <- rnb.meth.f[,samples2comp$sampleName]
# dim(rnb.meth2comp) # 213748     77 ok!
# 
# table(colnames(meth.test.nona)==samples2comp$sampleName) #true!
# 
# cond <- as.factor(car::recode(samples2comp$DONOR_SEX,"'Female'=1;'Male'=2"))
# 
# t1 <- Sys.time()
# meth.test.nona.models <- fn.models.parallel(meth.test.nona, cond1=cond, cores=6)
# t2 <- Sys.time()
# t2-t1 # 1.296551 hours (6 cores)
# 
# limma.nona.res <- apply.limma(meth.test.nona,cond)
# 
# meth.test.nona.models.withlimma <- data.frame(meth.test.nona.models,p.limma=limma.nona.res$P.Value)
# dim(meth.test.nona.models.withlimma) #102426
# 
# #save(meth.test.nona.models.withlimma,file="meth.test.nona.models.sex.RData")
# load(file="meth.test.nona.models.sex.RData")
# 
# ############# regressio betabinomial
# all <-all.test.nona
# covg <-covg.test.nona
# 
# t1 <- Sys.time()
# meth.test.nona.betabin.models <- fn.models.betabin.parallel(all=all.test.nona,covg=covg.test.nona, cond1=cond, cores=3)
# t2 <- Sys.time()
# t2-t1 # Time difference of 15.60176 mins
# 
# #no les matxambro, les guardo directament,
# save(meth.test.nona.betabin.models,file="meth.test.nona.betabin.models.sex.RData")
# load(file="meth.test.nona.betabin.models.sex.RData")
# 
# ########### comparo dss amb la resta
# meth.test.nona.models.withlimma.adj.p <- apply(meth.test.nona.models.withlimma,2,p.adjust)
# 
# #comprovo ordre
# all.equal(rownames(meth.test.nona.models.withlimma.adj.p),rownames(dmlTest))
# #no
# dmlTest.s <- dmlTest[rownames(meth.test.nona.models.withlimma.adj.p),]
# 
# plot(meth.test.nona.models.withlimma$p.s,dmlTest.s$fdr) #buff
# cor.test(meth.test.nona.models.withlimma$p.s,dmlTest.s$fdr) #0.7 i signif pero al plot es veu fatal...
# plot(meth.test.nona.models.withlimma$p.s,dmlTest.s$fdr,xlim=c(0,0.1),ylim=c(0,0.1)) 
# 
# cor.test(meth.test.nona.models.withlimma$p.b,dmlTest.s$fdr) #0.83 i signif pero al plot es veu fatal...
# plot(meth.test.nona.models.withlimma$p.b,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.sinf,dmlTest.s$fdr) #0.68
# plot(meth.test.nona.models.withlimma$p.sinf,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.binf,dmlTest.s$fdr) #0.83
# plot(meth.test.nona.models.withlimma$p.binf,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.n,dmlTest.s$fdr) #0.68
# plot(meth.test.nona.models.withlimma$p.n,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.l,dmlTest.s$fdr) #0.89
# plot(meth.test.nona.models.withlimma$p.l,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.q,dmlTest.s$fdr) #0.67
# plot(meth.test.nona.models.withlimma$p.q,dmlTest.s$fdr) #bu
# 
# cor.test(meth.test.nona.models.withlimma$p.limma,dmlTest.s$fdr) #0.85
# plot(meth.test.nona.models.withlimma$p.limma,dmlTest.s$fdr) #bu
# 
# 
# ############# regressio betabinomial
# all <-all.test.nona
# covg <-covg.test.nona
# 
# cond <- as.factor(car::recode(samples2comp$DONOR_SEX,"'Female'=1;'Male'=2"))
# 
# t1 <- Sys.time()
# meth.test.nona.betabin.models <- fn.models.betabin.parallel(all=all.test.nona,covg=covg.test.nona, cond1=cond, cores=6)
# t2 <- Sys.time()
# t2-t1 # 10.08667 mins
# 
# #save(meth.test.nona.betabin.models,file="meth.test.nona.betabin.models.sex.RData")
# load(file="meth.test.nona.betabin.models.sex.RData")
# meth.test.nona.betabin.models<-as.data.frame(meth.test.nona.betabin.models)
# 
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
# #save(all.models.sex,file="all.models.sex.RData")
# plot(all.models.sex$p.s,all.models.sex$p.b)
# cor.test(all.models.sex$p.s,all.models.sex$p.b) #0.7177047 
# cor.test(all.models.sex$p.s,all.models.sex$p.bb) #0.5754344 
# cor.test(all.models.sex$p.b,all.models.sex$p.bb) #0.7853257
# #spearman surt una miqueta més...
# 
# #ajusto pvals
# all.models.sex.adj.p <- apply(all.models.sex,2,p.adjust)
# 
# all.models.sex.adj.p.s <- all.models.sex.adj.p[all.models.sex.adj.p[,1]<0.05 & 
#                                                              !is.na(all.models.sex.adj.p[,1]),]
# all.models.sex.adj.p.b <- all.models.sex.adj.p[all.models.sex.adj.p[,2]<0.05 & 
#                                                              !is.na(all.models.sex.adj.p[,2]),]
# all.models.sex.adj.p.sinf <- all.models.sex.adj.p[all.models.sex.adj.p[,3]<0.05 & 
#                                                                 !is.na(all.models.sex.adj.p[,3]),]
# all.models.sex.adj.p.binf <- all.models.sex.adj.p[all.models.sex.adj.p[,4]<0.05 & 
#                                                                 !is.na(all.models.sex.adj.p[,4]),]
# all.models.sex.adj.p.n <- all.models.sex.adj.p[all.models.sex.adj.p[,5]<0.05 & 
#                                                              !is.na(all.models.sex.adj.p[,5]),]
# all.models.sex.adj.p.l <- all.models.sex.adj.p[all.models.sex.adj.p[,6]<0.05 & 
#                                                              !is.na(all.models.sex.adj.p[,6]),]
# all.models.sex.adj.p.q <- all.models.sex.adj.p[all.models.sex.adj.p[,7]<0.05 & 
#                                                              !is.na(all.models.sex.adj.p[,7]),]
# all.models.sex.adj.p.limma <- all.models.sex.adj.p[all.models.sex.adj.p[,9]<0.05 & 
#                                                                  !is.na(all.models.sex.adj.p[,9]),]
# all.models.sex.adj.p.bb <- all.models.sex.adj.p[all.models.sex.adj.p[,10]<0.05 & 
#                                                  !is.na(all.models.sex.adj.p[,10]),]
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
# res.table #molt bé la taula, bb mogollon i només a X
# 
# #     chr    s    b sinf binf    n  l    q limma   bb
# # 1   chr1 4614   NA 1006   NA   NA NA   NA    NA   NA
# # 2  chr10 4594   NA 1685   NA   NA NA   NA    NA   NA
# # 3  chr11 4927   NA 1268   NA   NA NA    1    NA   NA
# # 4  chr12 3699   NA 1202   NA   NA NA   NA    NA   NA
# # 5  chr13 2429   NA 1110   NA   NA NA   NA    NA   NA
# # 6  chr14 3366   NA 1344   NA   NA NA   NA    NA   NA
# # 7  chr15 2799   NA 1180   NA   NA NA   NA    NA   NA
# # 8  chr16 4558   NA  939   NA   NA NA   NA    NA   NA
# # 9  chr17 6384   NA 1188   NA   NA NA   NA    NA   NA
# # 10 chr18 2433   NA  986   NA   NA NA   NA    NA   NA
# # 11 chr19 9023   NA 1638   NA   NA NA    4    NA   NA
# # 12  chr2 6606   NA 2390   NA   NA NA   NA    NA   NA
# # 13 chr20 3104    1  888   NA    1 NA   NA    NA    1
# # 14 chr21 1190   NA  246   NA   NA NA   NA    NA   NA
# # 15 chr22 2968   NA  539   NA   NA NA    1    NA   NA
# # 16  chr3 3333   NA 1220   NA   NA NA   NA    NA   NA
# # 17  chr4 3606   NA 1831   NA   NA NA   NA    NA   NA
# # 18  chr5 6712   NA 2484   NA   NA NA   NA    NA   NA
# # 19  chr6 3831   NA 1448   NA   NA NA   NA    NA   NA
# # 20  chr7 5243   NA 1629   NA   NA NA   NA    NA   NA
# # 21  chr8 4193   NA 1338   NA   NA NA   NA    NA   NA
# # 22  chr9 4472   NA 1193   NA   NA NA    1    NA   NA
# # 23  chrX 8342 5166 7017 4066 4110 41 1539  4472 7091
# 
# #save(res.table,file="betadata.nona.models.sex.res.table.RData") #li poso el mateix nom que als arrays
# load(file="betadata.nona.models.sex.res.table.RData") #li poso el mateix nom que als arrays
