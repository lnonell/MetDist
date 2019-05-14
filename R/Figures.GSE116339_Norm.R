#27/4/19 Estudi de l'efecte de la normalització de les dades en les distribucions
# supp figure 5 i supp figure 6
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
# BESTDIST generades al cluster
#12/05/19: no canvia res, ho deixo com estava, D:\Doctorat\Simplex\MetDist\Data\Summary_OLD
# IMPORTANT: si cal modfiicar res caldrà canviar els directoris

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

########################################################
################ BEST DIST Supp figure 6 
#la mateixa fn que a Figures.BestDist.R
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

colors <- c("#E69F00", "#56B4E9", "#999999")

#
aic.plot <- function(best.dist.all,tit="A"){
  ll <- nrow(best.dist.all)
  aic.plot <- data.frame(x=1:ll,aic=best.dist.all[,1], model="simplex")
  aic.plot <- rbind(aic.plot,
                    data.frame(x=1:ll,aic=best.dist.all[,2], model="beta"),
                    data.frame(x=1:ll,aic=best.dist.all[,3], model="normal"))
  aic.plot$model <- factor(aic.plot$model,levels=c("simplex","beta","normal")) #aixo no estava a la fn original!
  
  
  p <-ggplot(aic.plot) + geom_smooth(aes(x=x, y=aic, color=model)) + xlab("CpG")+ ylab("AIC")+
    theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=7,angle=45), axis.text.y = element_text(size=8),
          axis.title = element_text(size=8), legend.key.size = unit(0.8,"line"))+
    scale_color_manual(values=colors) + ggtitle(tit)
  p
}

#ara les dades, grades al cluster
load(file=file.path(GSE116339_data,"best.dist.r.RData"))
r <- aic.plot(best.dist.r,tit="A")
r

load(file=file.path(GSE116339_data,"best.dist.i.RData"))
i <- aic.plot(best.dist.i,tit="B")
i

load(file=file.path(GSE116339_data,"best.dist.q.RData"))
q <- aic.plot(best.dist.q,tit="C")
q

load(file=file.path(GSE116339_data,"best.dist.fn.RData"))
fn <- aic.plot(best.dist.fn,tit="D")
fn

load(file=file.path(GSE116339_data,"best.dist.s.RData"))
s <- aic.plot(best.dist.s,tit="E")
s

load(file=file.path(GSE116339_data,"best.dist.n.RData"))
n <- aic.plot(best.dist.n,tit="F")
n

#per a tenir una única llegenda
mylegend<-g_legend(r)

#png(file=file.path(resultsDir,"Fig2.modelest.datasets.png"), width=480, height=480,res=100) #param res és l'important!!
ga <- grid.arrange(arrangeGrob(r + theme(legend.position="none"), #potser treure això de la legend??
                               i + theme(legend.position="none"),
                               q + theme(legend.position="none"),
                               fn + theme(legend.position="none"), 
                               s + theme(legend.position="none"), 
                               n + theme(legend.position="none"),
                               nrow=3),
                   mylegend, nrow=2,heights=c(10, 1))
#dev.off()
ggsave(file=file.path(resultsDir,"SuppFig6.GSE116339.norm.bestdist.pdf"), ga, width = 14, height = 14, units = "cm")


