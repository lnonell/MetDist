#12/1/19
#script def, la resta de proves són als altres script
#DMR dóna els mateixos resultats que lm, makes sense!


workingDir<-"D:/Doctorat/Simplex/MetDist/Data/GSE50660_Smoking"
setwd(workingDir)

#load(file="GMSet.RData") #carrega dades i phenoData pD.all

#load(file="beta.RData") #carrega les betes preprocessades i sense filtrar
load(file="pD.all.RData") #phenoData
load(file="beta.filtered.RData") #carrega les betes preprocessades i filtrades que són les que farem servir aquí
#atenció que és el mateix nom d'objecte: beta
dim(beta) #125950    464

#condició que utilitzarem en tots els mètodes
#smok <- pD.all$"smoking (0, 1 and 2, which represent never, former and current smokers):ch1"
sex <- pD.all$`gender:ch1` 
#############################################################################
####################### 1.Estimate params for each distribution #############
#############################################################################

#necessary functions
source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params
library(doParallel)
library(foreach)

##########################################################
#1. generate estimations for all distributions from beta (filtered!)
##########################################################
n <- nrow(beta)
t1 <- Sys.time()
cl <- makeCluster(4,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
est.params<- foreach(i=1:n,.combine=rbind, .packages=c("fitdistrplus","VGAM","ZOIP")) %dopar% {
  b <- beta[i,]
  est.all.params(b)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #4.878 hours
head(est.params) #les betainflated donen na pq no hi ha 0's ni 1's

est.params <-as.data.frame(est.params)
#save(est.params,file="params.est.all.120119.RData")
load(file="params.est.all.120119.RData")

#NAs
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

##########################################################
#2. Select vector of parameters and generate simulated data
##########################################################

######################## FN per a fer-ho tot!!!#####################
#enlloc de fer les simulacions per a 100 i dp per 10, creo una fn a SimulationFunctions i ho faig per varies Ns

t1 <- Sys.time()
simulations.all <- fn.simulations(est.params=est.params,cond.n= c(3,5,10,30,100,500),cores=2)
t2 <- Sys.time()
t2-t1 #2.028116 hours (2 cores)

load(file="simulated.cpgs.list.RData")
length(cpgs.list) #6 un per cada cond.n
length(cpgs.list[[1]]) #3: simplex, beta, normal
dim(cpgs.list[[1]][[1]]) #2000    6
head(cpgs.list[[1]][[1]])
hist(cpgs.list[[3]][[1]][1,1:10],breaks=10)
hist(cpgs.list[[3]][[1]][1,11:20],breaks=10)
hist(cpgs.list[[4]][[1]][1,11:20],breaks=10)
hist(cpgs.list[[4]][[1]][1,1:10],breaks=10)
hist(cpgs.list[[6]][[3]][1,1:500],breaks=20) #aquest no el pot ajustar simplex ni beta
hist(cpgs.list[[6]][[3]][1,501:1000],breaks=20) #aquest no el pot ajustar simplex ni beta, té negatius!!
range(cpgs.list[[5]][[3]][7,1:100],breaks=20) #aquest no el pot ajustar simplex ni beta, té negatius!!
range(cpgs.list[[5]][[3]][7,101:200],breaks=20) #aquest no el pot ajustar simplex ni beta
#miro els anteriors
range(cpgs.list[[5]][[3]][1,1:200])

load(file="simulated.models.list.RData")
length(res.list) #6 un per cada cond.n
length(res.list[[1]]) #3: simplex, beta, normal
dim(res.list[[1]][[1]]) #2000    6
head(res.list[[1]][[1]]) #sembla ok
head(res.list[[1]][[1]]) #sembla ok
head(res.list[[3]][[1]]) #sembla que no hi ha NAs
head(res.list[[3]][[3]],15) #algunes ditribucions normals no són capaces de ser ajustades per una simplex o beta
head(res.list[[6]][[3]],15)
table(is.na(res.list[[1]][[1]]))
table(is.na(res.list[[1]][[2]]))
table(is.na(res.list[[1]][[3]]))
table(is.na(res.list[[2]][[1]]))
table(is.na(res.list[[2]][[2]]))
table(is.na(res.list[[2]][[3]]))
#fins aquí hi ha uns 2000-2220 NAs, són els quantiles, no deu ser capaç d'estimar res
table(is.na(res.list[[3]][[1]])) #no
table(is.na(res.list[[3]][[2]])) #no
table(is.na(res.list[[3]][[3]])) #432
table(is.na(res.list[[4]][[1]])) #no
table(is.na(res.list[[4]][[2]])) #no
table(is.na(res.list[[4]][[3]])) #536
table(is.na(res.list[[5]][[1]])) #no
table(is.na(res.list[[5]][[2]])) #no
table(is.na(res.list[[5]][[3]])) #720
table(is.na(res.list[[6]][[1]])) #no
table(is.na(res.list[[6]][[2]])) #no
table(is.na(res.list[[6]][[3]])) #1152

#miro si les normals estan a (0,1)
head(sort(apply(cpgs.list[[6]][[3]],1,min,na.rm=T))) #totes les dades normals tenen valors negatius i més grans d'u
tail(sort(apply(cpgs.list[[6]][[3]],1,max,na.rm=T)))
#però això no semblaafectar als resultas

head(cpgs.list[[3]][[3]]) #el 6 té un neg
head(res.list[[3]][[3]]) #efectivament el 6 té NAs a les simplex i betes...


#comento tota la resta fins a les avaluacions dels models
models.eval(res.list[[1]][[1]],dml.r=100, alpha=0.05, adjust=TRUE) #sembla que sí, pero no son pitjors amb n petita

#passo les simulacions antigues fetes per 100 i 10 a OLDCodeSimulations.100.10.R del datadir


#################################################################################
####################### 3. Specific data analysis ###############################
#################################################################################

####3.1 comparacio fumar amb fns de tots els models #############

#amb la fn que prova tots els models, utilitzo var dicotòmica: aquesta en ppi no la faig servir
#agrupo 1 i 2 per a tenir una variable dicotòmica
# smok.b <- car::recode(smok,"2=1") #ajunto els que han fumat amb els que fumen. tot i que a nivell de metilació
# #això potser no t'e sentit pq volem veure què passa quan fumes i deixes de fumar, no?
# table(smok.b)
# #   0   1 
# 
# source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params

# table(pD.all$geo_accession==colnames(beta)) #yes
# 
# cond <- as.factor(smok.b+1) #la fn està preparada per a que la cond sigui 1 i 2
# t1 <- Sys.time()
# beta.models <- fn.models.parallel(beta, cond1=cond, cores=3)
# t2 <- Sys.time()
# t2-t1 #4.368636 hours
# dim(beta.models) #125950 
# 
# rownames(beta.models) <-rownames(beta)
# #save(beta.models,file="betadata.models.RData")
# #falta aplicar limma
# limma.res <- apply.limma(beta,cond)
# 
# #comprovo dims i ordres i afegeixo
# head(beta.models)
# head(limma.res)
# 
# dim(beta.models)
# #[1] 125950      7
# dim(limma.res)
# #[1]  125950      6
# 
# beta.models.withlimma <- data.frame(beta.models,p.limma=limma.res$P.Value)
# 
# #RENOMBRO per a que ensenyi tb els resultats de limma
# beta.models <- beta.models.withlimma
# 
# apply(beta.models,2,function(x) sum(is.na(x)))
# # p.s     p.b  p.sinf  p.binf     p.n     p.l     p.q p.limma 
# # 0       0       0       1       0       0       0       0 
# #molt bé, no hi ha NA's!
# 
# (beta.models.sign <- apply(beta.models,2,find.dml.n))
# # p.s     p.b  p.sinf  p.binf     p.n     p.l     p.q p.limma 
# # 40267   37092   46339   36951   37775   37394   26376   37588 
# beta.models.adj.p <- apply(beta.models,2,p.adjust)
# (beta.models.adj.p.sign <- apply(beta.models.adj.p,2,find.dml.n))
# # p.s     p.b  p.sinf  p.binf     p.n     p.l     p.q p.limma 
# # 776     408    3720     346     331     203     183     370 
# 
# #anoto, carrego les betes que ja havia anotat (GSE50660analysis.R)
# load(file="beta.anot.RData")
# 
# beta.models.anot <- cbind(beta.models,beta.anot[465:469])
# beta.models.adj.p.anot <- cbind(beta.models.adj.p,beta.anot[465:469])
# 
# beta.models.adj.p.anot.p.s <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.s<0.05,] #776
# beta.models.adj.p.anot.p.b <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.b<0.05,] #408
# beta.models.adj.p.anot.p.sinf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.sinf<0.05,] #3720
# beta.models.adj.p.anot.p.binf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.binf<0.05,] #347
# #hi ha una fila amb NAs, no sé d'on ha sortit. De moment l'esborro
# beta.models.adj.p.anot.p.binf <- beta.models.adj.p.anot.p.binf[-171,]
# beta.models.adj.p.anot.p.n <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.n<0.05,] #331
# beta.models.adj.p.anot.p.l <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.l<0.05,] #203
# beta.models.adj.p.anot.p.q <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.q<0.05,] #183
# beta.models.adj.p.anot.p.limma <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.limma<0.05,] #370
# 
# #ara ordenats
# beta.models.adj.p.anot.p.s.s<-beta.models.adj.p.anot.p.s[order(beta.models.adj.p.anot.p.s$p.s),]
# beta.models.adj.p.anot.p.b.s<-beta.models.adj.p.anot.p.b[order(beta.models.adj.p.anot.p.b$p.b),]
# beta.models.adj.p.anot.p.sinf.s<-beta.models.adj.p.anot.p.sinf[order(beta.models.adj.p.anot.p.sinf$p.sinf),]
# beta.models.adj.p.anot.p.binf.s<-beta.models.adj.p.anot.p.binf[order(beta.models.adj.p.anot.p.binf$p.binf),]
# beta.models.adj.p.anot.p.n.s<-beta.models.adj.p.anot.p.n[order(beta.models.adj.p.anot.p.n$p.n),]
# beta.models.adj.p.anot.p.l.s<-beta.models.adj.p.anot.p.l[order(beta.models.adj.p.anot.p.l$p.l),]
# beta.models.adj.p.anot.p.q.s<-beta.models.adj.p.anot.p.q[order(beta.models.adj.p.anot.p.q$p.q),]
# beta.models.adj.p.anot.p.limma.s<-beta.models.adj.p.anot.p.limma[order(beta.models.adj.p.anot.p.limma$p.limma),]
# 
# find.ahrr.gpr15 <- function(res){
#   res.ahrr <- res[grep("AHRR",res$genes),]
#   n1 <- nrow(res.ahrr)
#   print(n1)
#   res.gpr15 <- res[grep("GPR15",res$genes),]
#   n2 <- nrow(res.gpr15)
#   print(n2)
#   return(list(res.ahrr,res.gpr15))
# }
# 
# g.p.s <- find.ahrr.gpr15(beta.models.adj.p.anot.p.s) #7 0
# g.p.b <- find.ahrr.gpr15(beta.models.adj.p.anot.p.b) #5 0
# g.p.sinf <- find.ahrr.gpr15(beta.models.adj.p.anot.p.sinf) #7 0
# g.p.binf <- find.ahrr.gpr15(beta.models.adj.p.anot.p.binf) #5 0
# g.p.n <- find.ahrr.gpr15(beta.models.adj.p.anot.p.n) #5 0
# g.p.l <- find.ahrr.gpr15(beta.models.adj.p.anot.p.l) #5 0
# g.p.q <- find.ahrr.gpr15(beta.models.adj.p.anot.p.q) #2 0
# g.p.limma <- find.ahrr.gpr15(beta.models.adj.p.anot.p.limma) 
# 
# intersect(rownames(g.p.s[[1]]),rownames(g.p.b[[1]])) #les 5 del beta
# setdiff(rownames(g.p.s[[1]]),rownames(g.p.b[[1]])) # "cg09854184" "cg26703534"
# intersect(rownames(g.p.sinf[[1]]),rownames(g.p.s[[1]])) #6
# setdiff(rownames(g.p.s[[1]]),rownames(g.p.sinf[[1]])) #"cg26703534"
# intersect(rownames(g.p.binf[[1]]),rownames(g.p.b[[1]])) #les mateixes
# setdiff(rownames(g.p.s[[1]]),rownames(g.p.n[[1]])) #"cg09854184" "cg26703534"
# setdiff(rownames(g.p.s[[1]]),rownames(g.p.l[[1]]))
# setdiff(rownames(g.p.s[[1]]),rownames(g.p.q[[1]])) #"cg09854184" "cg03991871" "cg05575921" "cg26703534" "cg25648203"
# 
# #són els piquets baixos els que canvien
# sort(beta["cg09854184",smok==0]) #8 mostres <0.1 1 mostra 0.5 i la resta 0.8-09
# sort(beta["cg09854184",smok==1]) #3<0.1 2 0.4 i 2 0.6, la resta 0.8-0.9
# sort(beta["cg09854184",smok==2]) #tots 0.8-0.9

#cg26703534 la veritat és que presenta poques diferències
#"cg03991871" força diferent pero no al quantil 0.8?
#"cg05575921" molt diferent, hauria de ser significatiu, és la top de simplex!!
#"cg25648203" com cg03991871

#estrany que no hi hagi cap en GPR15

# #provo sense ajustar
# beta.models.anot.p.s <- beta.models.anot[beta.models.anot$p.s<0.05,]
# beta.models.anot.p.b <- beta.models.anot[beta.models.anot$p.b<0.05,]
# beta.models.anot.p.sinf <- beta.models.anot[beta.models.anot$p.sinf<0.05,]
# beta.models.anot.p.binf <- beta.models.anot[beta.models.anot$p.binf<0.05,]
# beta.models.anot.p.n <- beta.models.anot[beta.models.anot$p.n<0.05,]
# beta.models.anot.p.l <- beta.models.anot[beta.models.anot$p.l<0.05,]
# beta.models.anot.p.q <- beta.models.anot[beta.models.anot$p.q<0.05,]
# 
# g.p.s <- find.ahrr.gpr15(beta.models.anot.p.s) #23 15
# g.p.b <- find.ahrr.gpr15(beta.models.anot.p.b) #22 13
# g.p.sinf <- find.ahrr.gpr15(beta.models.anot.p.sinf) #25 16
# g.p.binf <- find.ahrr.gpr15(beta.models.anot.p.binf) #22 13
# g.p.n <- find.ahrr.gpr15(beta.models.anot.p.n) #23 14
# g.p.l <- find.ahrr.gpr15(beta.models.anot.p.l) #23 13
# g.p.q <- find.ahrr.gpr15(beta.models.anot.p.q) #14 9
# 
# g.all <- find.ahrr.gpr15(beta.models.anot) #34 40 aquest és el total de CpGs per AHRR i GPR15

#clar però surt una tercera part significatiu...
# beta.models.anot["cg05575921",] #super signif en tots els models
# 
# #de la taula1 de l'article, que guardo
# known.cpg <- c("cg06644428","cg05951221","cg21566642","cg01940273","cg19859270","cg05575921","cg21161138",
#                "cg06126421","cg03636183","cg27241845","cg03991871","cg25648203","cg19572487","cg25189904",
#                "cg23480021","cg21611682")
# 
# beta.models.adj.p.anot[known.cpg,1:8] #no troba "cg19859270" (GPR15) "cg21161138" "cg21611682", la resta supersignificatius
# 
# new.cpg <- c("cg03329539","cg13193840","cg02657160","cg14817490","cg24090911","cg24859433","cg15342087")
# beta.models.adj.p.anot[new.cpg,1:8] #no troba "cg02657160" 
# 
# conf.cell.cpg <- c("cg17024919")
# beta.models.adj.p.anot[conf.cell.cpg,1:8] #no signif
# beta.models.anot[conf.cell.cpg,1:8] #aqui si...
# 
# new.cpg.2 <- c("cg11660018","cg20295214","cg03547355","cg23079012","cg22717080","cg02451831")
# beta.models.adj.p.anot[new.cpg.2,1:8] #només troba "cg11660018","cg23079012" i signif

# ################ 3.2 el mateix pero comparant former smokers i non smokers per si de cas ####
# #ho deixo a l'script GSE50660analysis pq no millora, 
#  
# new.cpg.2 <- c("cg11660018","cg20295214","cg03547355","cg23079012","cg22717080","cg02451831")
# beta.models.0.2.adj.p.anot[new.cpg.2,1:8] #només troba "cg11660018","cg23079012" i signif

############################ anàlisi dels resultats i plots ######################
################### distribució del cg05575921 super signif en tots els models

# pdf("global.dist.pdf")
#   hist(beta,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="") 
#   lines(density(beta[,smok==0]),col="green")
#   lines(density(beta[,smok==1]),col="orange") 
#   lines(density(beta[,smok==2]),col="red")
# dev.off()
# #segons això no hi ha gaires diferències
# 
# hist(beta["cg05575921",],breaks=100,ylim=c(0,50))
# plot(density(beta["cg05575921",]),ylim=c(0,50))
# 
# pdf("cg05575921.AHRR.pdf")
# plot(density(beta["cg05575921",smok==0]),main="cg05575921", col="green",xlim=c(0,1))
# lines(density(beta["cg05575921",smok==1]),col="orange") 
# lines(density(beta["cg05575921",smok==2]),col="red") 
# dev.off()
# 
# pdf("cg07599136.AHRR.pdf")
# plot(density(beta["cg07599136",smok==0]),main="cg07599136", col="green",xlim=c(0,1))
# lines(density(beta["cg07599136",smok==1]),col="orange") 
# lines(density(beta["cg07599136",smok==2]),col="red") 
# dev.off()
# 
# #miro cpgs per mostrar
# known.cpg.in <- intersect(dimnames(beta)[[1]],known.cpg)
# length(known.cpg.in)
# 
# pdf("known.cpg.hist.pdf")
# for (i in 1:length(known.cpg.in)){
#   ki <- known.cpg.in[i]
#   h.b <- hist(beta[ki,],breaks=seq(0,1,0.05),freq=FALSE,main=ki,xlab="") 
#   lines(h.b$mids,h.b$density,col="black",lty=2,lwd=1)
# }
# dev.off()
# 
# #curiosament sembla força normal
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# 
# annotFile=  getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) #getAnnotation es de minfi
# #els gens que ens interessen són AHRR i GPR15, miro quines CpGs contenen
# annot.AHRR <- annotFile[grep("AHRR",annotFile$UCSC_RefGene_Name),]
# dim(annot.AHRR) #149
# annot.GPR15 <- annotFile[grep("GPR15",annotFile$UCSC_RefGene_Name),]
# dim(annot.GPR15) #145
# 
# GPR15.AHR.in <- intersect(dimnames(beta)[[1]],unique(rownames(annot.AHRR),rownames(annot.GPR15)))
# length(GPR15.AHR.in) #34 només

# pdf("GPR15.AHR.in.hist.pdf")
# for (i in 1:length(GPR15.AHR.in)){
#   ki <- GPR15.AHR.in[i]
#   h.b <- hist(beta[ki,],breaks=seq(0,1,0.05),freq=FALSE,main=ki,xlab="") 
#   lines(h.b$mids,h.b$density,col="black",lty=2,lwd=1)
# }
# dev.off()
# 
# pdf("GPR15.AHR.in.hist.xungo.pdf")
# for (i in 1:length(GPR15.AHR.in)){
#   ki <- GPR15.AHR.in[i]
#   hist(beta[ki,],breaks=100,freq=FALSE,main=ki,xlab="") 
# }
# dev.off()

#agafo els 200 primers...
# N=200
# pdf("first.200.hist.pdf")
# for (i in 1:N){
#   h.b <- hist(beta[i,],breaks=seq(0,1,0.05),freq=FALSE,main=dimnames(beta)[[1]][i],xlab="") 
#   lines(h.b$mids,h.b$density,col="black",lty=2,lwd=1)
# }
# dev.off()

#agafo els DE, p.adj i ordenats, faig una fn per a fer-los tots: 3 condicions de variable
# plot.hist.dens <- function(models,pdf.name){
#   pdf(pdf.name)
#   for (i in 1:nrow(models)){
#     ki <- rownames(models)[i]
#     h.b <- hist(beta[ki,],breaks=seq(0,1,0.05),freq=FALSE,main=ki,xlab="", ylim=c(0,12)) 
#     lines(density(beta[ki,smok==0]),col="green")
#     lines(density(beta[ki,smok==1]),col="orange") 
#     lines(density(beta[ki,smok==2]),col="red")
#   }
#   dev.off()
# }  
#   
# plot.hist.dens(beta.models.adj.p.anot.p.s.s,"DE.adj.p.simplex.hist.pdf")
# plot.hist.dens(beta.models.adj.p.anot.p.b.s,"DE.adj.p.beta.hist.pdf")
# plot.hist.dens(beta.models.adj.p.anot.p.sinf.s,"DE.adj.p.simplexinflated.hist.pdf")
# plot.hist.dens(beta.models.adj.p.anot.p.binf.s,"DE.adj.p.betainflated.hist.pdf")
# plot.hist.dens(beta.models.adj.p.anot.p.n.s,"DE.adj.p.normal.hist.pdf")
# plot.hist.dens(beta.models.adj.p.anot.p.l.s,"DE.adj.p.logistic.hist.pdf")
# plot.hist.dens(beta.models.adj.p.anot.p.q.s,"DE.adj.p.quartil.hist.pdf")
# 
# beta["cg01533021",smok==0]
# beta["cg01533021",smok==1]
# beta["cg01533021",smok==2]
# 
# #simplex inflated n'hi ha molts, miro com es distribueixen
# sinf.spec <- setdiff(rownames(beta.models.adj.p.anot.p.sinf.s),rownames(beta.models.adj.p.anot.p.s.s))
# plot.hist.dens(beta.models.adj.p.anot.p.sinf.s[sinf.spec,],"DE.adj.p.simplexinflated.specific.hist.pdf")
# #miro els plots i semblen absolutament falsos positius!!!
# 
# #ara nomes dos condicions
# plot.hist.dens.2cond <- function(models,pdf.name){
#   pdf(pdf.name)
#   for (i in 1:nrow(models)){
#     ki <- rownames(models)[i]
#     h.b <- hist(beta[ki,],breaks=seq(0,1,0.05),freq=FALSE,main=ki,xlab="", ylim=c(0,16)) 
#     lines(density(beta[ki,smok==0]),col="green")
#     lines(density(beta[ki,smok %in% c(1,2)]),col="red") 
#     #lines(density(beta[ki,smok==2]),col="red")
#   }
#   dev.off()
# }  
# 
# plot.hist.dens.2cond(beta.models.adj.p.anot.p.s.s,"DE.adj.p.simplex.hist.2cond.pdf")
# plot.hist.dens.2cond(beta.models.adj.p.anot.p.b.s,"DE.adj.p.beta.hist.2cond.pdf")
# plot.hist.dens.2cond(beta.models.adj.p.anot.p.sinf.s,"DE.adj.p.simplexinflated.hist.2cond.pdf")
# plot.hist.dens.2cond(beta.models.adj.p.anot.p.binf.s,"DE.adj.p.betainflated.hist.2cond.pdf")
# plot.hist.dens.2cond(beta.models.adj.p.anot.p.n.s,"DE.adj.p.normal.hist.2cond.pdf")
# plot.hist.dens.2cond(beta.models.adj.p.anot.p.l.s,"DE.adj.p.logistic.hist.2cond.pdf")
# plot.hist.dens.2cond(beta.models.adj.p.anot.p.q.s,"DE.adj.p.quartil.hist.2cond.pdf")
# 
# 
# pdf("GPR15.AHR.in.dens.groups.pdf")
# for (i in 1:length(GPR15.AHR.in)){
#   ki <- GPR15.AHR.in[i]
#   plot(density(beta[ki,smok==0]),main=ki, col="green",xlim=c(0,1))
#   lines(density(beta[ki,smok==1]),col="orange") 
#   lines(density(beta[ki,smok==2]),col="red")
# }
# dev.off()

#no acabo de veure quins plotar, no veig alguns que tinguin dos pics als extrems potser justament
#els que presenten més diferències
# cg05575921 #super signif en tots els models 
# cg07599136 #interessant, fumadors perfils diferents
# cg16219322 #tot i tenir el mateix perfil en els tres grups menys extrem com més es fuma
# en altres data sets ho havia vist...potser efecte de la normalització??

#############################################################################################
################ 3.4 dmpfinder de minfi (F-test of linear regression) ####
#############################################################################################

library(minfi)

time1.0<-Sys.time()
dmp <- dmpFinder(beta, pheno = cond , type = "categorical") #
#dmp és un data.frame
time2.0<-Sys.time()
time.dmp<-time2.0-time1.0 #9.64 secs pepino!


#### results ###
res.dmp<-dmp[dmp$pval<0.05,]
head(res.dmp)
dim(res.dmp) #37775! 

res.dmp.qval<-dmp[dmp$qval<0.05,] #20873 són molts!
dim(res.dmp.qval)  #28733 són molts igualment!

#provo d'ajustar jo el pval
head(dmp)
dmp$adj.p <- p.adjust(dmp$pval)
head(dmp)
res.dmp.adj.p<-dmp[dmp$adj.p<0.05,] #331 ara si!

length(intersect(rownames(res.dmp.adj.p),rownames(beta.models.adj.p.anot.p.s))) #327
length(intersect(rownames(res.dmp.adj.p),rownames(beta.models.adj.p.anot.p.b))) #330
length(intersect(rownames(res.dmp.adj.p),rownames(beta.models.adj.p.anot.p.sinf))) #322
length(intersect(rownames(res.dmp.adj.p),rownames(beta.models.adj.p.anot.p.binf))) #321
length(intersect(rownames(res.dmp.adj.p),rownames(beta.models.adj.p.anot.p.n))) #331 #is in fact the same model!!!
length(intersect(rownames(res.dmp.adj.p),rownames(beta.models.adj.p.anot.p.l))) #203
length(intersect(rownames(res.dmp.adj.p),rownames(beta.models.adj.p.anot.p.q))) #40

#############################################################################################
################ 3.5 SEX ####
#############################################################################################

sex.r <- car::recode(sex,"'Female'=1;'Male'=2") #ajunto els que han fumat amb els que fumen. tot i que a nivell de metilació
#això potser no t'e sentit pq volem veure què passa quan fumes i deixes de fumar, no?
table(sex.r)
#   1   2 
# 137 327 

source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params

table(pD.all$geo_accession==colnames(beta)) #yes

cond <- as.factor(sex.r) #la fn està preparada per a que la cond sigui 1 i 2
t1 <- Sys.time()
beta.models <- fn.models.parallel(beta, cond1=cond, cores=2)
t2 <- Sys.time()
t2-t1 #7.764886 hours 2 cores
dim(beta.models) #125950 

rownames(beta.models) <-rownames(beta)
#save(beta.models,file="betadata.models.sex.RData") #he matxacat els betamodels de smoke!!!
load(file="betadata.models.sex.RData")

#AKI falta aplicar limma: cannot allocate vector: FET
limma.res <- apply.limma(beta,cond) 

beta.models.withlimma <- data.frame(beta.models,p.limma=limma.res$P.Value)
save(beta.models.withlimma,file="betadata.models.sex.withlimma.RData")
beta.models.adj.p <- apply(beta.models.withlimma,2,p.adjust)

#anotacions
load(file="beta.anot.RData")

beta.models.anot <- cbind(beta.models,beta.anot[465:469])
beta.models.adj.p.anot <- cbind(beta.models.adj.p,beta.anot[465:469])
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
#     Var1    s    b sinf binf    n    l    q limma
# 1   chr1  298  287  416  273  267  226  106   271
# 2  chr10  146  136  198  128  125  103   42   128
# 3  chr11  123  118  176  106  106   84   41   111
# 4  chr12  139  130  190  125  123  103   41   124
# 5  chr13   59   53   96   50   49   44   25    51
# 6  chr14   97   87  130   85   85   67   26    89
# 7  chr15   86   84  117   83   83   73   41    83
# 8  chr16  124  111  163  104  105   97   49   111
# 9  chr17  184  175  248  167  163  142   57   167
# 10 chr18   42   36   61   34   34   32   16    34
# 11 chr19  186  168  253  159  160  135   59   165
# 12  chr2  220  198  280  182  183  154   70   184
# 13 chr20   65   58   86   56   55   41   20    58
# 14 chr21   31   23   40   22   22   18    8    24
# 15 chr22   42   40   59   40   39   35   16    40
# 16  chr3  148  136  198  127  127  103   40   129
# 17  chr4  119  115  177  107  108   90   40   107
# 18  chr5  152  153  215  147  141  125   41   146
# 19  chr6  257  243  338  231  222  188   68   238
# 20  chr7  149  127  217  120  122  104   48   126
# 21  chr8   97   84  140   77   80   59   28    81
# 22  chr9   61   59   91   57   58   53   25    58
# 23  chrX 8097 8120 8113 8108 8095 4970 7927  8113
# 24  chrY  375  375  375  375  375  238  359   375
save(res.table,file="betadata.models.sex.res.table.RData")

######################################################################
############ 4. Estimació de la millor dist per a cada cpg ###########
######################################################################
source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params
dim(beta)

library(doParallel)
library(foreach)

N=nrow(beta)
#e=0.01
t1 <- Sys.time()
cl <- makeCluster(4,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
best.dist.all<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","fitdistrplus")) %dopar% {
  b <- round(beta[i,],4)
  #b.prime <- ifelse(b==1, 1-e, ifelse(b==0,0+e,b))
  #best.dist(b.prime) #aqui no cal treure la inflation pq les dades estan norm i no cotnenen 0s ni 1s
  best.dist(b)
}
stopCluster(cl)
t2 <- Sys.time()
t2-t1 #2.23 hours

head(best.dist.all) 
#ho guardo
save(best.dist.all,file=file.path("best.dist.all.RData"))
load(file=file.path("best.dist.all.RData"))

gse50660.best.dist <- t(apply(best.dist.all,1,best.aic.ks))

table(gse50660.best.dist[,1])
# beta.aic  normal.aic simplex.aic 
# 19737       32179       74034 
#genial!!!

#provo amb VGAM, doncs en ppi és la que farem servir: NO TE psimplex!! HE DE FER SERVIR simplexreg per a fer ks test
# N=nrow(beta)
# #e=0.01
# t1 <- Sys.time()
# cl <- makeCluster(4,type="PSOCK",outfile="output.txt")
# registerDoParallel(cl)
# best.dist.all.VGAM<- foreach(i=1:N,.combine=rbind, .packages=c("VGAM","fitdistrplus")) %dopar% {
#   b <- round(beta[i,],4)
#   #b.prime <- ifelse(b==1, 1-e, ifelse(b==0,0+e,b))
#   #best.dist(b.prime) #aqui no cal treure la inflation pq les dades estan norm i no cotnenen 0s ni 1s
#   best.dist.VGAM(b)
# }
# stopCluster(cl)
# t2 <- Sys.time()
# t2-t1 #22mins hours
# 
# head(best.dist.all.VGAM) #exactament igual excepte simplex.ks.p
# best.dist.all.VGAM <- t(apply(best.dist.all.VGAM,1,best.aic.ks))
# table(best.dist.all.VGAM[,1])
# 

