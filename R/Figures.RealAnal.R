#Figure 3: Global distribution for all data sets and the variable sex
#23/4/19 change variable, all to sex and canvio RRBS216 a RRBS188
#Supp Figure 7: Venn diagram, en aquest cas no sé si té gaire sentit
#Supp table 6 amb els resultats...res.table per a cada data set

#11/05/19 ho actualitzo tot
#Figura 5: AIC dels models

#header amb tot
source("D:/Doctorat/Simplex/MetDist/R/Figures.HeaderScript.R")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


######################################################################
########### Figure 3: Distribució de les betes per cada data set amb les condicions######
######################################################################

#fer fn i per a cada data set grar un plot, que després ajuntarem
dens.plot <- function(betes,cond,cond.names=c("cond1","cond2"),tit="A"){
    #just in case 
    betes <- as.matrix(betes)
    betes.cond <- data.frame(x=as.vector(betes[,cond==1]),cond=cond.names[1])
    betes.cond <- rbind(betes.cond,
                      data.frame(x=as.vector(betes[,cond==2]),cond=cond.names[2]))
    ggplot(betes.cond, aes(x=x)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
    geom_density(alpha=.2,aes(color=cond, fill=cond))  + ggtitle(tit)
}

# a <- dens.plot(betes=beta,cond,cond.names=c("Non Smokers","Smokers"),tit="A")
# a
# tit <- c("A","B")
# windowsFonts(A = windowsFont("Arial"))
#textplot(tit,cex=4,family="A")


#########################
#GSE50660
load(file=file.path(GSE50660_data,"beta.filtered.RData")) #carrega les betes preprocessades i filtrades que són les que farem servir aquí
load(file=file.path(GSE50660_data,"pD.all.RData")) #phenoData
betes.a <- beta
rm(beta)
sex <- pD.all$`gender:ch1` 
sex.r <- car::recode(sex,"'Female'=1;'Male'=2") #ajunto els que han fumat amb els que fumen. tot i que a nivell de metilació
#això potser no t'e sentit pq volem veure què passa quan fumes i deixes de fumar, no?
table(sex.r)
#   1   2 
# 137 327 
cond.a <-as.factor(sex.r) #la fn està preparada per a que la cond sigui 1 i 2

#########################
#GSE116339
load(file=file.path(GSE116339_data,"PM.filtered.RData")) 
load(file=file.path(GSE116339_data,"pheno.RData"))

dim(PM.f) #198711  
betes.b <- as.matrix(PM.f)
rm(PM.f)
table(pheno$`gender:ch1`)
# Female   Male 
# 399    280 
cond.b <- as.factor(car::recode(pheno$`gender:ch1`,"'Female'=1;'Male'=2"))

hist(as.matrix(PM.f),breaks=seq(0,1,0.02),freq=FALSE,main="",xlab="")
lines(density(as.matrix(PM.f[,cond.pbb==1]),na.rm=T),col="green")
lines(density(as.matrix(PM.f[,cond.pbb==2]),na.rm=T),col="red")

#########################
#RRBS188 enlloc de RRBS216
load(file=file.path(RRBS188_data,"meth.betas.poquesNAs.covg3.RData"))
dim(rnb.meth.f) # 270569    188

samples <- read.table(file=file.path(RRBS188_data,"samples.csv"), sep=",", header=T, stringsAsFactors = F)
samples2comp <- samples[samples$patient_sex %in% c("f","m"),]
dim(samples2comp) #159

betes.c <- rnb.meth.f[,samples2comp$sample_id]
rm(rnb.meth.f)
cond.c <- as.factor(car::recode(samples2comp$patient_sex,"'f'=1;'m'=2"))
table(cond.c)
# 1  2 
# 63 96 
#########################
#WGBS81
load(file=file.path(WGBS81_data,"meth.betas.poquesNAs.covg3.RData"))
dim(rnb.meth.f) #204073     81
samples <- read.delim(file=file.path(WGBS81_data,"samples.tsv"), stringsAsFactors=FALSE)

samples2comp <- samples[samples$DONOR_SEX %in% c("Female","Male"),]
dim(samples2comp) #77
betes.d <- rnb.meth.f[,samples2comp$sampleName]
rm(rnb.meth.f)
cond.d <- as.factor(car::recode(samples2comp$DONOR_SEX,"'Female'=1;'Male'=2"))
table(cond.d)

colors <- c("green4", "blue", "red2", "orange")
col.lines<-c("brown","mediumorchid2")
# col.lines<-c("green","red")
#FIGURA: SI NO LA VULL COLOURED TREURE ELS BORDERS I TORNAR A POSAR EL RED GREEN COM A col.lines
pdf(file=file.path(resultsDir,"Fig3.datasetsdist.coloured.110519.pdf"))
  par(mfrow=c(2,2))
  hist(betes.a,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="",ylab="", border=colors[1])
    lines(density(betes.a[,cond.a==1],na.rm=T),col=col.lines[1])
    lines(density(betes.a[,cond.a==2],na.rm=T),col=col.lines[2])
  # title("A",outer=T,line=1)
  
  hist(betes.b,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="",ylab="", border=colors[2])
    lines(density(betes.b[,cond.b==1],na.rm=T),col=col.lines[1])
    lines(density(betes.b[,cond.b==2],na.rm=T),col=col.lines[2])
  # title("B",outer=T)
  
  hist(betes.c,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="",ylab="", border=colors[3])
    lines(density(betes.c[,cond.c==1],na.rm=T),col=col.lines[1])
    lines(density(betes.c[,cond.c==2],na.rm=T),col=col.lines[2])
  # title("C",outer=T)
  
  hist(betes.d,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="",ylab="", ylim=c(0,4), border=colors[4])
    lines(density(betes.d[,cond.d==1],na.rm=T),col=col.lines[1])
    lines(density(betes.d[,cond.d==2],na.rm=T),col=col.lines[2])
  # title("D",outer=T)
dev.off()


########################################################################################
######### Supp Figure 7: Comparacions de cada dataset: Venn diagrams dels significatius #####
########################################################################################
#posarés els mateixos colors que a les simulacions, que són aquests
#faig unes proves amb la inflated... dóna moooolts resultats, no les faig servir

scales::show_col(scales::hue_pal()(6))
# i corresponen a: beta, simplex, normal,logistic,quantile, limma
# color <- c("#F8766D" "#B79F00" "#00BA38" "#00BFC4" "#619CFF" "#F564E3")
# elimino la logística que ja veig que és un xurro: he de canviar els dos primers d'orde que aquí tinc simplex
color <- c("#B79F00","#F8766D","#00BA38","#619CFF","#F564E3")

#fn per fer el Venn
Venn5D <- function(list.1=l1,list.2=l2,list.3=l3,list.4=l4,list.5=l5,listNames=l.names,
                   filename=file.path(resultsDir,"Venn.GSE50.100.100"),CatCex=0.8, CatDist=rep(0.1, 5)){
  
  # require(Vennerable)
  require(VennDiagram)
 # require(colorfulVennPlot) #Per generar plots amb colors diferents
 # require(RColorBrewer)
  
  # cols <- c(brewer.pal(8,"Pastel1"), brewer.pal(8,"Pastel2"))  #Fixar els colors per als venn diagrams amb 4 condicions
  # cols <- cols[c(8,2,3,15,5,6,7,1,9,10,11,12,13,14,4,16)]
  #remove na's just in case
  list.1 <- list.1[!is.na(list.1)]
  list.2 <- list.2[!is.na(list.2)]
  list.3 <- list.3[!is.na(list.3)]
  list.4 <- list.4[!is.na(list.4)]
  list.5 <- list.5[!is.na(list.5)]
  
  list.venn<-list(list.1,list.2,list.3,list.4,list.5)
  names(list.venn)<-listNames
  
  pdf(filename)
  venn.plot <- draw.quintuple.venn(
    area1 = length(list.1),
    area2 = length(list.2),
    area3 = length(list.3),
    area4 = length(list.4),
    area5 = length(list.5),
    n12 = length(intersect(list.1,list.2)),
    n13 = length(intersect(list.1,list.3)),
    n14 = length(intersect(list.1,list.4)),
    n15 = length(intersect(list.1,list.5)),
    n23 = length(intersect(list.2,list.3)),
    n24 = length(intersect(list.2,list.4)),
    n25 = length(intersect(list.2,list.5)),
    n34 = length(intersect(list.3,list.4)),
    n35 = length(intersect(list.3,list.5)),
    n45 = length(intersect(list.4,list.5)),
    n123 = length(intersect(list.1,intersect(list.2,list.3))),
    n124 = length(intersect(list.1,intersect(list.2,list.4))),
    n125 = length(intersect(list.1,intersect(list.2,list.5))),
    n134 = length(intersect(list.1,intersect(list.3,list.4))),
    n135 = length(intersect(list.1,intersect(list.3,list.5))),
    n145 = length(intersect(list.1,intersect(list.4,list.5))),
    n234 = length(intersect(list.2,intersect(list.3,list.4))),
    n235 = length(intersect(list.2,intersect(list.3,list.5))),
    n245 = length(intersect(list.2,intersect(list.4,list.5))),
    n345 = length(intersect(list.3,intersect(list.4,list.5))),
    n1234 = length(intersect(list.1,intersect(list.2,intersect(list.3,list.4)))),
    n1235 = length(intersect(list.1,intersect(list.2,intersect(list.3,list.5)))),
    n1245 = length(intersect(list.1,intersect(list.2,intersect(list.4,list.5)))),
    n1345 = length(intersect(list.1,intersect(list.3,intersect(list.4,list.5)))),
    n2345 = length(intersect(list.2,intersect(list.3,intersect(list.4,list.5)))),
    n12345 = length(intersect(list.1,intersect(list.2,intersect(list.3,intersect(list.4,list.5))))),
    category = listNames,
    fill = color,
    cat.cex = CatCex,
    cat.dist = CatDist,
    margin = 0.05,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
            1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    ind = TRUE)
  
  dev.off()          
  
}

venn5D_betamodels <-function(beta.models.adj.p,file=file.path(resultsDir,"Venn.GSE50.smoking.pdf")){
  l1 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p[,1]<0.05) & beta.models.adj.p[,1]<0.05,])
  l2 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p[,2]<0.05) & beta.models.adj.p[,2]<0.05,])
  l3 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p[,3]<0.05) & beta.models.adj.p[,3]<0.05,])
  l4 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p[,4]<0.05) & beta.models.adj.p[,4]<0.05,])
  l5 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p[,5]<0.05) & beta.models.adj.p[,5]<0.05,])
  l.names=c("simplex","beta","normal","quantile","limma")
  #fn definida més a munt
  Venn5D(list.1=l1,list.2=l2,list.3=l3,list.4=l4,list.5=l5,listNames=l.names,
         filename=file, CatCex=0.8, CatDist=rep(0.1, 5))
}

#GSE50660
load(file=file.path(GSE50660_data,file="betadata.models.sex.withlimma.RData"))
head(beta.models.withlimma)
beta.models.adj.p <- as.data.frame(apply(beta.models.withlimma,2,p.adjust))
#i ara selecciono els params que vull plotar
to.plot <-beta.models.adj.p[,c("p.s","p.b","p.n","p.q","p.limma")]
venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.GSE50.sex.080719.pdf"))
# to.plot <-beta.models.adj.p[,c("p.sinf","p.b","p.n","p.q","p.limma")]
# venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.GSE50.sex.sinf.050719.pdf"))

# #miro de quin chr son els 281 per si estan majoritariament a X però NO!
# l1 <- rownames(to.plot[!is.na(to.plot[,1]<0.05) & to.plot[,1]<0.05,])
# l2 <- rownames(to.plot[!is.na(to.plot[,2]<0.05) & to.plot[,2]<0.05,])
# l3 <- rownames(to.plot[!is.na(to.plot[,3]<0.05) & to.plot[,3]<0.05,])
# l4 <- rownames(to.plot[!is.na(to.plot[,4]<0.05) & to.plot[,4]<0.05,])
# l5 <- rownames(to.plot[!is.na(to.plot[,5]<0.05) & to.plot[,5]<0.05,])
# 
# spec.s<-setdiff(l1,c(l2,l3,l4,l5))
# 
# load(file=file.path(GSE50660_data,"beta.anot.RData"))
# dim(beta.anot)
# 
# table(beta.anot[spec.s,"chr"])
# # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3 
# # 29    15    13    11     7     7     3    14    13     7    17    25     5     7     3    20 
# # chr4  chr5  chr6  chr7  chr8  chr9  chrX 
# # 8     7    20    23    12     2    12 

#GSE116339
load(file=file.path(GSE116339_data,file="betadata.models.sex.withlimma.RData"))
head(beta.models.withlimma)
beta.models.adj.p <- as.data.frame(apply(beta.models.withlimma,2,p.adjust))
to.plot <-beta.models.adj.p[,c("p.s","p.b","p.n","p.q","p.limma")]
venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.GSE116339.sex.110519.pdf"))
# to.plot <-beta.models.adj.p[,c("p.sinf","p.b","p.n","p.q","p.limma")]
# venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.GSE116339.sex.sinf.050719.pdf"))


#RRBS188 enlloc de RRBS216
load(file=file.path(RRBS188_data,file="rnb.meth2comp.models.withlimma.RData")) 
head(rnb.meth2comp.models.withlimma)
#li poso el mateix nom i així ho aprofito tot!
beta.models.adj.p <- as.data.frame(apply(rnb.meth2comp.models.withlimma,2,p.adjust))
to.plot <-beta.models.adj.p[,c("p.s","p.b","p.n","p.q","p.limma")]
venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.RRBS188.sex.110519.pdf"))
# to.plot <-beta.models.adj.p[,c("p.sinf","p.b","p.n","p.q","p.limma")]
# venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.RRBS188.sex.sinf.050719.pdf"))


#WGBS81
load(file=file.path(WGBS81_data,file="rnb.meth2comp.models.withlimma.RData")) 
head(rnb.meth2comp.models.withlimma)
#ccarrega objecte rnb.meth2comp.lineage.models
#li poso el mateix nom i així ho aprofito tot!
beta.models.adj.p <- as.data.frame(apply(rnb.meth2comp.models.withlimma,2,p.adjust))
to.plot <-beta.models.adj.p[,c("p.s","p.b","p.n","p.q","p.limma")]
venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.WGBS81.sex.110519.pdf"))
# to.plot <-beta.models.adj.p[,c("p.sinf","p.b","p.n","p.q","p.limma")]
# venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.WGBS81.sex.sinf.050719.pdf"))

########################################################################################
######################## Figura 5: AIC dels models   ######################################
########################################################################################
#5/7/19 elimino logistic
library(scales)
colors <- hue_pal()(6)
colors <- colors[c(1:3,5)]


aic.plot <- function(beta.models.withlimma,tit="A"){
  ll <- nrow(beta.models.withlimma)
  aic.plot <- data.frame(x=1:ll,aic=beta.models.withlimma$aic.sinf, model="simplex") #inflated!!
  aic.plot <- rbind(aic.plot,
                    data.frame(x=1:ll,aic=beta.models.withlimma$aic.b, model="beta"),
                    data.frame(x=1:ll,aic=beta.models.withlimma$aic.n, model="normal"),
                   # data.frame(x=1:ll,aic=beta.models.withlimma$aic.l, model="logistic"),
                    data.frame(x=1:ll,aic=beta.models.withlimma$aic.q, model="quantile"))
  
  aic.plot$model <- factor(aic.plot$model, levels=c("beta", "simplex","normal","quantile"))
  
  p <-ggplot(aic.plot) + geom_smooth(aes(x=x, y=aic, color=model)) + xlab("CpG")+ ylab("AIC")+
    theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=7,angle=45), axis.text.y = element_text(size=8),
          axis.title = element_text(size=8), legend.key.size = unit(0.8,"line"))+
     scale_color_manual(values=colors) + 
    ggtitle(tit)
  p
}

aic.bb.plot <- function(beta.models.withlimma,betabin.models,tit="A"){
  ll <- nrow(beta.models.withlimma)
  aic.plot <- data.frame(x=1:ll,aic=beta.models.withlimma$aic.sinf, model="simplex") #inflated!!
  aic.plot <- rbind(aic.plot,
                    data.frame(x=1:ll,aic=beta.models.withlimma$aic.b, model="beta"),
                    data.frame(x=1:ll,aic=beta.models.withlimma$aic.n, model="normal"),
                    data.frame(x=1:ll,aic=beta.models.withlimma$aic.l, model="logistic"),
                    data.frame(x=1:ll,aic=beta.models.withlimma$aic.q, model="quantile"),
                    data.frame(x=1:ll,aic=betabin.models$aic.bb, model="beta-binomial"))
  
  aic.plot$model <- factor(aic.plot$model, levels=c("beta", "simplex","normal","logistic","quantile","beta-binomial"))
  
  p <-ggplot(aic.plot) +  geom_smooth(aes(x=x, y=aic, color=model)) + xlab("CpG")+ylab("AIC")+
    theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=7,angle=45), axis.text.y = element_text(size=8),
          axis.title = element_text(size=8), legend.key.size = unit(0.8,"line"))+
   #scale_color_manual(values=colors2) + 
    ggtitle(tit)
  p
}


load(file=file.path(GSE50660_data,"betadata.models.sex.withlimma.RData"))
head(beta.models.withlimma)
a <- aic.plot(beta.models.withlimma,tit="A")
a

load(file=file.path(GSE116339_data,"betadata.models.sex.withlimma.RData"))
head(beta.models.withlimma)
b <- aic.plot(beta.models.withlimma,tit="B")
b

load(file=file.path(RRBS188_data,"rnb.meth2comp.models.withlimma.RData")) 
head(rnb.meth2comp.models.withlimma)
c <- aic.plot(rnb.meth2comp.models.withlimma,tit="C")
c
load(file=file.path(RRBS188_data,"rnb.meth2comp.betabin.models.RData"))
head(rnb.meth2comp.betabin.models)
c.bb <- aic.bb.plot(rnb.meth2comp.models.withlimma,rnb.meth2comp.betabin.models,tit="C")
c.bb

load(file=file.path(WGBS81_data,"rnb.meth2comp.models.withlimma.RData")) 
head(rnb.meth2comp.models.withlimma)
d <- aic.plot(rnb.meth2comp.models.withlimma,tit="D")
d

load(file=file.path(WGBS81_data,"rnb.meth2comp.betabin.models.RData"))
head(rnb.meth2comp.betabin.models)
d.bb <- aic.bb.plot(rnb.meth2comp.models.withlimma,rnb.meth2comp.betabin.models,tit="D")
d.bb

#per a tenir una única llegenda
mylegend<-g_legend(a)

ga <- grid.arrange(arrangeGrob(a + theme(legend.position="none"), 
                               b + theme(legend.position="none"),
                               c + theme(legend.position="none"),
                               d + theme(legend.position="none"), 
                               nrow=2),
                   mylegend, nrow=2,heights=c(10, 1))

ggsave(file=file.path(resultsDir,"Fig5.AIC.sex.datasets.050719.pdf"), ga, width = 14, height = 14, units = "cm")

mylegend<-g_legend(d)
ga <- grid.arrange(arrangeGrob(a + theme(legend.position="none"), 
                               b + theme(legend.position="none"),
                               c.bb + theme(legend.position="none"),
                               d.bb + theme(legend.position="none"), 
                               nrow=2),
                   mylegend, nrow=2,heights=c(10, 1))
ggsave(file=file.path(resultsDir,"Fig5.AIC.sex.datasets.betabin.120519.pdf"), ga, width = 14, height = 14, units = "cm")


########################################################################################
##### SuppTable 6: Comparacions de cada dataset: res.table at FDR 0.05 per chr     #####
########################################################################################
load(file=file.path(GSE50660_data,file="betadata.models.sex.res.table.RData"))
rt.gse50 <-res.table
load(file=file.path(GSE116339_data,file="betadata.models.sex.res.table.RData"))
rt.gse11 <-res.table
load(file=file.path(RRBS188_data,file="betadata.models.sex.res.table.RData"))
rt.rrbs <-res.table
load(file=file.path(WGBS81_data,file="betadata.models.sex.res.table.RData"))
rt.wgbs <-res.table

#ho ordeno, ajunto tot i ho poso en forma de taula, transposada
chr.l <- paste0("chr",c(1:22,"X","Y"))
#per a poder ordenar poso els rownames
rownames(rt.gse50) <- rt.gse50$Var1
rt.gse50.s <- rt.gse50[chr.l,]
rownames(rt.gse11) <- rt.gse11$Var1
rt.gse11.s <- rt.gse11[chr.l,]
rownames(rt.rrbs) <- rt.rrbs$chr
rt.rrbs.s <- rt.rrbs[chr.l,]
rownames(rt.wgbs) <- rt.wgbs$chr
rt.wgbs.s <- rt.wgbs[chr.l,]

#trasposo i ajunto
rt.all <- rbind(data.frame(dataset="A450k-smoking",t(rt.gse50.s)[-1,]),
          data.frame(dataset="EPIC-PBB",t(rt.gse11.s)[-1,]),
          data.frame(dataset="RRBS-ES",t(rt.rrbs.s)[-1,]),
          data.frame(dataset="WGBS-BLUE",t(rt.wgbs.s)[-1,]))

rt.all
write.csv2(rt.all,file=file.path(resultsDir,"SuppTable6.Results.sex.110519.csv"), row.names = T, quote = F)
