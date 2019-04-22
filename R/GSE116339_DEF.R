#12/1/19
#script def, la resta de proves són a l'altre script
#go to line 216!!

workingDir<-"D:/Doctorat/Simplex/MetDist/Data/GSE116339_PBB"
setwd(workingDir)

load(file="PM.filtered.RData") #atenció li he posat els noms més a baix pq no havia guardat els CGs
dim(PM.f) #198711    679 un cop eliminat els SNPs
load(file="pheno.RData")

#############################################################################
####################### 1.Estimate params for each distribution #############
#############################################################################

#necessary functions
source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params
library(doParallel)
library(foreach)

#1. generate estimations for all distributions from PM.f 
######################
n <- nrow(PM.f)
t1 <- Sys.time()
cl <- makeCluster(6,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
rnb.est.params<- foreach(i=1:n,.combine=rbind, .packages=c("fitdistrplus","VGAM","ZOIP")) %dopar% {
  b <- PM.f[i,]
  est.all.params(b)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #4.41h
head(rnb.est.params)

colnames(rnb.est.params)

#treure rnb de l'objecte!!
est.params <-as.data.frame(rnb.est.params)
#save(est.params,file="params.est.all.120119.RData")
save(est.params,file="params.est.all.210219.RData")
load(file="params.est.all.210219.RData")

#study NAs
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

##########################################################
#2. Select vector of parameters and generate simulated data
##########################################################
#enlloc de fer les simulacions per a 100 i dp per 10, creo una fn a SimulationFunctions i ho faig per varies Ns
#la fn guarda dos objectes amb llistes
# simulated.cpgs.list.RData: conté cpgs.list
# simulated.models.list.RData: conte res.list

t1 <- Sys.time()
simulations.all <- fn.simulations(est.params=est.params,cond.n= c(3,5,10,30,100,500),cores=4)
t2 <- Sys.time()
t2-t1 #1.184927 hours (4 cores)

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


############################################################################
########### 3.1 Specific data analyisis: sexe  ##############################
#############################################################################
str(pheno)

table(pheno$characteristics_ch1)
table(pheno$`gender:ch1`)
# # Female   Male 
# # 399    280
# table(pheno$`gender:ch1`,pheno$characteristics_ch1.2) #si!
# table(pheno$`tissue:ch1`)
# #miro totes les vars de les characteristics
# table(pheno$characteristics_ch1)
# table(pheno$characteristics_ch1.1) #tots peripheral blood
# table(pheno$characteristics_ch1.3) #age cont
# hist(as.numeric(pheno$characteristics_ch1.3)) #en mesos?
# range(as.numeric(pheno$characteristics_ch1.3)) #1 570
# table(pheno$characteristics_ch1.4) #ln(totalpbb) cont
# table(pheno$characteristics_ch1.5) #pbb-153
# pbb.153.d <- as.numeric(pheno$characteristics_ch1.5)>1
# table(pbb.153.d) 
# # FALSE  TRUE 
# # 8   671 
# table(pheno$characteristics_ch1.6) #pbb-77
# pbb.77.d <- as.numeric(pheno$characteristics_ch1.6)>1
# table(pbb.77.d) 
# # FALSE  TRUE 
# # 616    63 
# table(pheno$characteristics_ch1.7) #pbb-101
# pbb.101.d <- as.numeric(pheno$characteristics_ch1.7)>1
# table(pbb.101.d) 
# # pbb.101.d
# # FALSE  TRUE 
# # 608    71 
# table(pheno$characteristics_ch1.8) #pbb-180
# pbb.180.d <- as.numeric(pheno$characteristics_ch1.8)>1
# table(pbb.180.d) 
# # FALSE  TRUE 
# # 676     3 
# 
# #res dicotòmic, mmm, hauré d'agafar el gender
# 
cond <- as.factor(car::recode(pheno$`gender:ch1`,"'Female'=1;'Male'=2"))
t1 <- Sys.time()
beta.models <- fn.models.parallel(PM.f, cond1=cond, cores=3)
t2 <- Sys.time()
t2-t1 #10.75 hours
dim(beta.models)  #257029      7
# 
# #no hi ha els noms, faig una xapusa per no haver-ho de fer tot denou!!
# cg.f <- PM$V1[as.numeric(rownames(PM.f))]
# rownames(PM.f) <-cg.f
# rownames(beta.models) <- cg.f
# 
# #rownames(beta.models) <-rownames(beta)
# #save(beta.models,file="betadata.models.RData")
# load(file="betadata.models.RData")
# #com que estava fet abans de filtrar els snps, faig una altra xapusa per a seleccionar els cpgs
# beta.models <- beta.models[rownames(PM.f),]
# dim(beta.models) #198711

save(beta.models,file="betadata.models.sex.RData") #ho havia guardat sense sex
load(file="betadata.models.sex.RData")

#AKI Falta afegir limma: cannot allocate vector
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
save(beta.models.withlimma,file="betadata.models.sex.withlimma.RData")
load(file="betadata.models.sex.withlimma.RData")

beta.models.adj.p <- apply(beta.models.withlimma,2,p.adjust)

###################################################################
################# anotacions i estudi dels resultats
###################################################################
############## fn per anotar ###################################
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #hi ha les b2 b3 i b4. Sembla que la b4 és l'última
# 
# annotFile=  getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #getAnnotation es de minfi
# colnames(annotFile)

#no aconsegueixo que funcioni paral·lelitzant!!!
# get.cpg.annot <- function(df, annotFile,cores=3){
#   library(doParallel)
#   library(foreach)
#   
#   #rownames are the cpg id
#    cpgs <- rownames(df)
#    lng <- length(cpgs)
#   # genes <- array(NA,dim=lng)
#   # chr <- array(NA,dim=lng) 
#   # pos <- array(NA,dim=lng) 
#   # group <- array(NA,dim=lng) 
#   # island_rel <- array(NA,dim=lng)
#   #ho he de paral·lelitzar o serà la locura
#   cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
#   registerDoParallel(cl)
#   
#   get.annot <- function(cpgi, annotFile){
#     genes <- annotFile[cpgi,"UCSC_RefGene_Name"]
#     chr <- annotFile[cpgi,"chr"]
#     pos <- annotFile[cpgi,"pos"]
#     group <- annotFile[cpgi,"UCSC_RefGene_Group"]
#     island_rel <- annotFile[cpgi,"Relation_to_Island"]
#     c(cpgi,genes,chr,pos,group,island_rel)
#   }  
#   
#   cpg.annot<- foreach(i=1:lng,.combine=rbind) %dopar% {
#     cpgi <- cpgs[i]
#     get.annot(cpgi,annotFile)
#   }
#   stopCluster(cl)
#   colnames(cpg.annot) <- c("cpgi","genes","chr","pos","group","island_rel") 
#   df1 <- data.frame(df,cpg.annot)
#   return(df1)
# }   

# get.cpg.annot <- function(df){
#   #rownames are the cpg id
#   cpgs <- rownames(df)
#   lng <- length(cpgs)
#   genes <- array(NA,dim=lng)
#   chr <- array(NA,dim=lng) 
#   pos <- array(NA,dim=lng) 
#   group <- array(NA,dim=lng) 
#   island_rel <- array(NA,dim=lng)
#   #ho posarem a cada pas del bucle, per a que sigui m?s optimitzat
#   for (i in 1:lng){
#     print(i)
#     genes[i] <- annotFile[cpgs[i],"UCSC_RefGene_Name"]
#     chr[i] <- annotFile[cpgs[i],"chr"]
#     pos[i] <- annotFile[cpgs[i],"pos"]
#     group[i] <- annotFile[cpgs[i],"UCSC_RefGene_Group"]
#     island_rel[i] <- annotFile[cpgs[i],"Relation_to_Island"]
#   }
#   df1 <- data.frame(df,genes,chr,pos,group,island_rel)
#   return(df1)
# }

#prova
#PM.f.annot <- get.cpg.annot(PM.f[1:5,])

########################## resum de tot ###########################
#anoto: fetes sobre les  257029 abans d'eliminar els SNPs i renombro l'objecte RData amb els SNPs
# t1 <- Sys.time()
#   PM.f.annot <- get.cpg.annot(PM.f)
# t2 <- Sys.time()
# t2-t1 #18h!

# table(colnames(PM.f)==colnames(PM.f.annot)[1:679]) #si
# annots <- PM.f.annot[,680:684]
# head(annots)
# annots.f <- annots[rownames(PM.f),]
# dim(annots.f) #198711
# 
# #miro que els rownames siguin els mateixos abans d'ajuntar
# table(rownames(annots.f)==rownames(PM.f)) #TRUE
# 
# PM.f.annot <- data.frame(PM.f,annots.f)
# #save(PM.f.annot,file="PM.f.annot.RData")
load(file="PM.f.annot.RData")
# 
#beta.models.anot <- cbind(beta.models,PM.f.annot[680:684])
beta.models.adj.p.anot <- cbind(beta.models.adj.p,PM.f.annot[680:684])
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

res.table #molts resultats a X!!
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

find.gene <- function(res, gene=gene){
  res.genes <- res[grep(gene,res$genes),]
  n1 <- nrow(res.genes)
  print(n1)
  return(n1)
}
# 
# #l'article parla d'alguns gens, anem a veure
beta.models.adj.p.anot[grep("EZH2",beta.models.adj.p.anot $genes),] #2 s,b limma pero 3 sinf
beta.models.adj.p.anot [grep("H3K9",beta.models.adj.p.anot $genes),] #none
beta.models.adj.p.anot [grep("H19",beta.models.adj.p.anot $genes),] #molta cosa

#############################################################################
########### 3.2 Specific data analyisis: PBB  ##############################
#############################################################################
str(pheno)

#PBB, que és la var d'exposició són valors contínuus
table(pheno$`pbb-101:ch1`) #608 0s
hist(as.numeric(pheno$`pbb-101:ch1`),breaks=100)

########## PBB
table(pheno$`ln(totalpbb):ch1`) 
hist(as.numeric(pheno$`ln(totalpbb):ch1`),breaks=100)
#1890 CpG sites associated with total PBB levels
dim(pheno[pheno$`ln(totalpbb):ch1`>0,]) #166 individuus amb nivell és gran que 0
cond.pbb <- as.factor(ifelse(pheno$`ln(totalpbb):ch1`>0,2,1))

t1 <- Sys.time()
beta.models.pbb <- fn.models.parallel (PM.f, cond1=cond.pbb, cores=6)
t2 <- Sys.time()
t2-t1 #7,515 hours

#falta aplicar limma
limma.res <- apply.limma(PM.f,cond.pbb) #retorna la top table amb l'ordre original

#Proves de models diferents
#intento mirar si es pot fer model cont
# cpgi <- PM.f[1,]
# #hist(as.numeric(cpgi))
# 
# cond <- cond.pbb[!is.na(cpgi)]
# age <- as.numeric(pheno$"age:ch1")[!is.na(cpgi)]
# sex <- pheno$"gender:ch1"[!is.na(cpgi)]
# cpgi <-cpgi[!is.na(cpgi)]
# 
# m.s <- simplexreg(cpgi ~ cond)
# summary(m.s, save=FALSE)
# 
# m.s1 <- simplexreg(cpgi ~ cond + age + sex)
# summary(m.s1, save=FALSE) #ja no és significatiu!
# 
# 
# #condició contínua
# cpgi <- PM.f[1,]
# cond <- as.numeric(pheno$`ln(totalpbb):ch1`)[!is.na(cpgi)]
# age <- as.numeric(pheno$"age:ch1")[!is.na(cpgi)]
# sex <- pheno$"gender:ch1"[!is.na(cpgi)]
# cpgi <-cpgi[!is.na(cpgi)]
# 
# m.s <- simplexreg(cpgi ~ cond)
# summary(m.s, save=FALSE) #significatiu però una kk d'estimate
# 
# m.s1 <- simplexreg(cpgi ~ cond + age + sex)
# summary(m.s1, save=FALSE) #ja no és significatiu!

#provo els models
beta.models <- beta.models.pbb
#save(beta.models,file="betadata.models.pbb.RData")

load(file="betadata.models.pbb.RData")
#falta aplicar limma
limma.res <- apply.limma(PM.f,cond.pbb) #retorna la top table amb l'ordre original

#comprovo dims i ordres i afegeixo
dim(beta.models)
#[1] 198711      7
dim(limma.res)
#[1] 198711      6

beta.models.withlimma <- data.frame(beta.models,p.limma=limma.res$P.Value)
save(beta.models.withlimma,file="betadata.models.withlimma.pbb.RData")

#RENOMBRO per a que ensenyi tb els resultats de limma
beta.models <- beta.models.withlimma

beta.models.adj.p <- apply(beta.models,2,p.adjust)
(beta.models.adj.p.sign <- apply(beta.models.adj.p,2,find.dml.n))
# p.s     p.b  p.sinf  p.binf     p.n     p.l     p.q p.limma 
# 148     195    4532     161     343     130       4     344 
#millor, serà més fàcil de veure!
load(file="PM.f.annot.RData")

beta.models.anot <- cbind(beta.models,PM.f.annot[680:684])
beta.models.adj.p.anot <- cbind(beta.models.adj.p,PM.f.annot[680:684])

beta.models.adj.p.anot.p.s <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.s<0.05,]
beta.models.adj.p.anot.p.b <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.b<0.05,]
beta.models.adj.p.anot.p.sinf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.sinf<0.05,]
beta.models.adj.p.anot.p.binf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.binf<0.05 & 
                                                          !is.na(beta.models.adj.p.anot$p.binf),] #no se que passa amb binf
dim(beta.models.adj.p.anot.p.binf) #161
beta.models.adj.p.anot.p.n <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.n<0.05,]
beta.models.adj.p.anot.p.l <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.l<0.05,]
beta.models.adj.p.anot.p.q <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.q<0.05,]
beta.models.adj.p.anot.p.limma <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.limma<0.05,]

#sembla que hi ha moltes X
table(beta.models.adj.p.anot.p.s$chr) #105
table(beta.models.adj.p.anot.p.b$chr) #156
table(beta.models.adj.p.anot.p.sinf$chr) #3902
table(beta.models.adj.p.anot.p.binf$chr) #126
table(beta.models.adj.p.anot.p.n$chr) #306
table(beta.models.adj.p.anot.p.l$chr) #114
table(beta.models.adj.p.anot.p.q$chr) #1
table(beta.models.adj.p.anot.p.limma$chr) #312
#SI!

#mirar gens de l'article!
#ARNT and ESR2
beta.models.adj.p.anot[grep("ARNT",beta.models.adj.p.anot$genes),]
beta.models.adj.p.anot[grep("ESR2",beta.models.adj.p.anot$genes),]
#res significatiu!

#top gen de CTD amb interactions amb Polybrominated Biphenyls
beta.models.adj.p.anot[grep("PPARG",beta.models.adj.p.anot$genes),]
beta.models.adj.p.anot[grep("THRB",beta.models.adj.p.anot$genes),]
beta.models.adj.p.anot[grep("THRA",beta.models.adj.p.anot$genes),]
#RES!
#sembla que surten moltes coses al chr X
table(cond.pbb,pheno$`gender:ch1`)
# cond.pbb Female Male
# 1    326  187
# 2     73   93
beta.models.adj.p.anot[grep("TBX3",beta.models.adj.p.anot$genes),]


#altres gens de Disgenet associats amb Thyroid diseases (232 gens)
beta.models.adj.p.anot[grep("DIO2",beta.models.adj.p.anot$genes),]
beta.models.adj.p.anot[grep("ID3",beta.models.adj.p.anot$genes),]

sort(beta.models.adj.p.anot.p.s$genes[beta.models.adj.p.anot.p.s$genes!=""])
#ABCD1, MECP2 només de p.s i de Disgenet: ho hauria de tornar a mirar

plot.hist.dens.2cond <- function(models,pdf.name){
  pdf(pdf.name)
  for (i in 1:nrow(models)){
    ki <- rownames(models)[i]
    h.b <- hist(as.numeric(PM.f[ki,]),breaks=seq(0,1,0.05),freq=FALSE,main=ki,xlab="", ylim=c(0,16)) 
    lines(density(as.numeric(PM.f[ki,cond.pbb==1]),na.rm=T),col="green")
    lines(density(as.numeric(PM.f[ki,cond.pbb==2]),na.rm=T),col="red") 
    #lines(density(beta[ki,smok==2]),col="red")
  }
  dev.off()
}  

plot.hist.dens.2cond(models=beta.models.adj.p.anot.p.s,pdf.name="DE.adj.p.simplex.hist.pbb.pdf")
plot.hist.dens.2cond(beta.models.adj.p.anot.p.b,"DE.adj.p.beta.hist.pbb.pdf")
plot.hist.dens.2cond(beta.models.adj.p.anot.p.sinf,"DE.adj.p.simplexinflated.hist.pbb.pdf")
plot.hist.dens.2cond(models=beta.models.adj.p.anot.p.binf,"DE.adj.p.betainflated.hist.pbb.pdf")
plot.hist.dens.2cond(beta.models.adj.p.anot.p.n,"DE.adj.p.normal.hist.pbb.pdf")
plot.hist.dens.2cond(beta.models.adj.p.anot.p.l,"DE.adj.p.logistic.hist.pbb.pdf")
plot.hist.dens.2cond(beta.models.adj.p.anot.p.q,"DE.adj.p.quartil.hist.pbb.pdf")


######################################################################
############ 4. Estimació de la millor dist per a cada cpg ###########
######################################################################

source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params
dim(PM.f)

library(doParallel)
library(foreach)

N=nrow(PM.f)
range(PM.f,na.rm=T)
# 0.001394743 0.999601009 tampoc elimino!

#e=0.01
t1 <- Sys.time()
cl <- makeCluster(4,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
best.dist.all<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- round(unlist(PM.f[i,]),4)
  #in this case we need to remove NAs
  b.prime <-b[!is.na(b)]
  #b.prime <- ifelse(b==1, 1-e, ifelse(b==0,0+e,b))
  #best.dist(b.prime) #aqui no cal treure la inflation pq les dades estan norm i no cotnenen 0s ni 1s
  best.dist(b.prime)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #7.12 Hours

head(best.dist.all) #

save(best.dist.all,file=file.path("best.dist.all.RData"))

GSE116339.best.dist <- t(apply(best.dist.all,1,best.aic.ks))

table(GSE116339.best.dist[,1])
# beta.aic  normal.aic simplex.aic 
# 31236       33927      133548 
#si!

table(GSE116339.best.dist[,2])
# beta.ks.p  normal.ks.p simplex.ks.p 
# 34779        43748       120184 