#27/4/19 Estudi de l'efecte de la normalització de les dades en les distribucions
#les dades les carrego al cluster des dels fitxers idat que he baixat de GEO de nou
#al cluster he 
# 1. Normalitzat el RGSet que carrega amb els 6 metodes:
#  ?preprocessRaw #RGChannelSet
#  ?preprocessIllumina #RGChannelSet
#  ?preprocessSWAN #RGChannelSet
#  ?preprocessQuantile #An object of class RGChannelSet or [Genomic]MethylSet, aquest sí que el podem aplicar!
#  ?preprocessNoob #RGChannelSet
#  ?preprocessFunnorm #RGCChannelSet
# 2. Eliminat SNPs
# 3. Extret i guardat les betes resultants de cada normalització
# ara les he de filtrar aquí igual que les originals i mirar les distribucions

library(GEOquery)
library(minfi)

workingDir <- "D:/Doctorat/Simplex/MetDist/Data/GSE116339_PBB"
setwd(workingDir)
resultsDir <- "D:/Doctorat/Simplex/MetDist/Data/Summary"

##################################################################
# pheno data
load(file="pheno.RData")

################# load data #####################################
load(file="PM.filtered.RData") #atenció li he posat els noms més a baix pq no havia guardat els CGs
dim(PM.f) #198711
load(file="RSet.r.b.f.RData")
dim(RSet.r.b.f)
length(intersect(rownames(PM.f),rownames(RSet.r.b.f))) #198711 si!
#i ara la resta
load(file="RSet.fn.b.f.RData")
#RSet.fn.b.f <- as.data.frame(RSet.fn.b.f)
dim(RSet.fn.b.f) #198711
load(file="RSet.q.b.f.RData")
#RSet.q.b.f <- as.data.frame(RSet.q.b.f)
dim(RSet.q.b.f) #198711
load(file="RSet.i.b.f.RData")
#RSet.i.b.f <- as.data.frame(RSet.i.b.f)
dim(RSet.i.b.f) #198711
load(file="RSet.s.b.f.RData")
#RSet.s.b.f <- as.data.frame(RSet.s.b.f)
dim(RSet.s.b.f) #198711
load(file="RSet.n.b.f.RData")
#RSet.n.b.f <- as.data.frame(RSet.n.b.f)
dim(RSet.n.b.f) #198711

table(pheno$`gender:ch1`)
# Female   Male 
# 399    280 
cond.b <- as.factor(car::recode(pheno$`gender:ch1`,"'Female'=1;'Male'=2"))

colors <- c("green4", "blue", "red2", "orange")
col.lines<-c("brown","mediumorchid2")
# col.lines<-c("green","red")
#FIGURA: SI NO LA VULL COLOURED TREURE ELS BORDERS I TORNAR A POSAR EL RED GREEN COM A col.lines
pdf(file=file.path(resultsDir,"SuppFig5.GSE116339.Normaliz.coloured.270419.pdf"))
par(mfrow=c(3,2))
  ylim=c(0,3)
  hist(RSet.r.b.f,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="Raw",ylab="",ylim= ylim, border=colors[2])
  hist(RSet.i.b.f,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="Illumina",ylab="",ylim= ylim, border=colors[2])
  hist(RSet.q.b.f,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="Quantile",ylab="",ylim= ylim, border=colors[2])
  hist(RSet.fn.b.f,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="Funnorm",ylab="",ylim= ylim, border=colors[2])
  hist(RSet.s.b.f,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="SWAN",ylab="",ylim= ylim, border=colors[2])
  hist(RSet.n.b.f,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="noob",ylab="",ylim= ylim, border=colors[2])
  # lines(density(RSet.r.b.f[,cond.b==1],na.rm=T),col=col.lines[1])
  # lines(density(RSet.r.b.f[,cond.b==2],na.rm=T),col=col.lines[2])
dev.off()

#sospito que la que tinc és quantile 
cor.test(RSet.q.b.f,as.matrix(PM.f))

#heatmap amb les correls??
cor.test(RSet.r.b.f,RSet.i.b.f) #0.9719538 
cor.test(RSet.r.b.f,RSet.q.b.f) #0.9706348 
cor.test(RSet.r.b.f,RSet.fn.b.f) #0.9706348 
cor.test(RSet.r.b.f,RSet.s.b.f) #0.9945702 
cor.test(RSet.r.b.f,RSet.n.b.f) #0.9789752 

cor.test(RSet.i.b.f,RSet.q.b.f) #0.9699878 
cor.test(RSet.i.b.f,RSet.fn.b.f) #0.9461845
cor.test(RSet.i.b.f,RSet.s.b.f) #0.986376 
cor.test(RSet.i.b.f,RSet.n.b.f) #0.9618818 

cor.test(RSet.q.b.f,RSet.fn.b.f) #0.9523339 
cor.test(RSet.q.b.f,RSet.s.b.f) #0.9783623 
cor.test(RSet.q.b.f,RSet.n.b.f) #0.9554284

cor.test(RSet.fn.b.f,RSet.s.b.f) #0.9667648 
cor.test(RSet.fn.b.f,RSet.n.b.f) #0.9910817 

cor.test(RSet.s.b.f,RSet.n.b.f) #0.9767127 



