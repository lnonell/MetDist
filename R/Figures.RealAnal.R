#Figure 4: Global distribution for all data sets

#header amb tot
source("D:/Doctorat/Simplex/R/Figures.HeaderScript.R")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


######################################################################
########### Figure 4: Distribució de les betes per cada data set amb les condicions######
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
smok <- pD.all$"smoking (0, 1 and 2, which represent never, former and current smokers):ch1"
smok.b <- car::recode(smok,"2=1") #ajunto els que han fumat amb els que fumen. tot i que a nivell de metilació
table(smok.b)
#   0   1 
cond.a <- as.factor(smok.b+1) #la fn està preparada per a que la cond sigui 1 i 2

#########################
#GSE116339
load(file=file.path(GSE116339_data,"PM.filtered.RData")) 
load(file=file.path(GSE116339_data,"pheno.RData"))

dim(PM.f) #198711  
betes.b <- as.matrix(PM.f)
rm(PM.f)
cond.b <- ifelse(pheno$`ln(totalpbb):ch1`>0,2,1)


hist(as.matrix(PM.f),breaks=seq(0,1,0.02),freq=FALSE,main="",xlab="")
lines(density(as.matrix(PM.f[,cond.pbb==1]),na.rm=T),col="green")
lines(density(as.matrix(PM.f[,cond.pbb==2]),na.rm=T),col="red")

#########################
#RRBS216
load(file=file.path(RRBS216_data,"meth.betas.poquesNAs.RData"))
dim(rnb.meth.f) #38003   216
betes.c<-rnb.meth
rm(rnb.meth)

samples.new <- read.table(file=file.path(RRBS216_data,"samples.new.csv"), sep="\t", header=T, stringsAsFactors = F)
samples2comp <- samples.new[samples.new$lineage %in% c("endoderm","mesoderm"),]
dim(samples2comp) #136

betes.c <- rnb.meth.f[,samples2comp$sampleName]
rm(rnb.meth.f)
cond.c <- as.factor(car::recode(samples2comp$lineage,"'endoderm'=1;'mesoderm'=2"))

#########################
#WGBS81
load(file=file.path(WGBS81_data,"meth.betas.poquesNAs.RData"))
dim(rnb.meth.f) #213748     81
rm(rnb.meth)
samples <- read.delim(file=file.path(WGBS81_data,"samples.tsv"), stringsAsFactors=FALSE)

samples2comp <- samples[samples$DONOR_SEX %in% c("Female","Male"),]
dim(samples2comp) #77
betes.d <- rnb.meth.f[,samples2comp$sampleName]
rm(rnb.meth.f)
cond.d <- as.factor(car::recode(samples2comp$DONOR_SEX,"'Female'=1;'Male'=2"))

colors <- c("green4", "blue", "red2", "orange")
col.lines<-c("brown","mediumorchid2")
# col.lines<-c("green","red")
#FIGURA: SI NO LA VULL COLOURED TREURE ELS BORDERS I TORNAR A POSAR EL RED GREEN COM A col.lines
pdf(file=file.path(resultsDir,"Fig4.datasetsdist.coloured.pdf"))
  par(mfrow=c(2,2))
  hist(betes.a,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="", border=colors[1])
    lines(density(betes.a[,cond.a==1],na.rm=T),col=col.lines[1])
    lines(density(betes.a[,cond.a==2],na.rm=T),col=col.lines[2])
  # title("A",outer=T,line=1)
  
  hist(betes.b,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="", border=colors[2])
    lines(density(betes.b[,cond.b==1],na.rm=T),col=col.lines[1])
    lines(density(betes.b[,cond.b==2],na.rm=T),col=col.lines[2])
  # title("B",outer=T)
  
  hist(betes.c,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="", border=colors[3])
    lines(density(betes.c[,cond.c==1],na.rm=T),col=col.lines[1])
    lines(density(betes.c[,cond.c==2],na.rm=T),col=col.lines[2])
  # title("C",outer=T)
  
  hist(betes.d,breaks=seq(0,1,0.05),freq=FALSE,main="",xlab="", ylim=c(0,4), border=colors[4])
    lines(density(betes.d[,cond.d==1],na.rm=T),col=col.lines[1])
    lines(density(betes.d[,cond.d==2],na.rm=T),col=col.lines[2])
  # title("D",outer=T)
dev.off()


########################################################################################
######### Figure 5: Comparacions de cada dataset: Venn diagrams dels significatius #####
########################################################################################
#posarés els mateixos colors que a les simulacions, que són aquests
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
load(file=file.path(GSE50660_data,file="betadata.models.withlimma.RData"))
head(beta.models.withlimma)
beta.models.adj.p <- as.data.frame(apply(beta.models.withlimma,2,p.adjust))
#i ara selecciono els params que vull plotar
to.plot <-beta.models.adj.p[,c("p.s","p.b","p.n","p.q","p.limma")]
venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.GSE50.smoking.pdf"))

#GSE116339
load(file=file.path(GSE116339_data,file="betadata.models.withlimma.pbb.RData"))
head(beta.models.withlimma)
beta.models.adj.p <- as.data.frame(apply(beta.models.withlimma,2,p.adjust))
to.plot <-beta.models.adj.p[,c("p.s","p.b","p.n","p.q","p.limma")]
venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.GSE116339.pbb.pdf"))

#RRBS216
load(file=file.path(RRBS216_data,file="rnb.meth2comp.lineage.models.withlimma.RData")) 
#ccarrega objecte rnb.meth2comp.lineage.models
#li poso el mateix nom i així ho aprofito tot!
beta.models.adj.p <- as.data.frame(apply(rnb.meth2comp.lineage.models.withlimma,2,p.adjust))
to.plot <-beta.models.adj.p[,c("p.s","p.b","p.n","p.q","p.limma")]

venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.RRBS216.lineage.pdf"))

#RRBS216
load(file=file.path(RRBS216_data,file="rnb.meth2comp.lineage.models.withlimma.RData")) 
#ccarrega objecte rnb.meth2comp.lineage.models
#li poso el mateix nom i així ho aprofito tot!
beta.models.adj.p <- as.data.frame(apply(rnb.meth2comp.lineage.models.withlimma,2,p.adjust))
to.plot <-beta.models.adj.p[,c("p.s","p.b","p.n","p.q","p.limma")]
venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.RRBS216.lineage.pdf"))

#WGBS81
load(file=file.path(WGBS81_data,file="rnb.meth2comp.models.withlimma.RData")) 
#ccarrega objecte rnb.meth2comp.lineage.models
#li poso el mateix nom i així ho aprofito tot!
beta.models.adj.p <- as.data.frame(apply(rnb.meth2comp.models.withlimma,2,p.adjust))
to.plot <-beta.models.adj.p[,c("p.s","p.b","p.n","p.q","p.limma")]
venn5D_betamodels(to.plot,file=file.path(resultsDir,"Venn.WGBS81.lineage.pdf"))

