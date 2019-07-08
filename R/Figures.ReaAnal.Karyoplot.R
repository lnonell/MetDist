#Figure 6 EF CIRCOS per a tot el genoma i zoom a X amb karyoplotR, density (abans era la 5)
#primer faig els chrX de cada data set amb karyoplotR i al final el CIRCOS

#La munto al pptx fent instantànies amb zoom 410% al fox i activant regles per a fer-ho tot igual: 
# 1 circos i giro 90º per a tenir la X a la dreta
# 2 Cada chrX de cada data set fet amb KaryoploteR faig zoom 410 amb el fox i selecciono de x=(1,5,17) y=(4,5,14,5)
# 3 canvio el tamany a pptx per a 3,5 x 13 cada track (elimino nom chrX)
# Deixo les altres proves comentades per si de cas

#05/07/19 elimino la logística!!!

#header amb tot
source("D:/Doctorat/Simplex/MetDist/R/Figures.HeaderScript.R")

#load(file=file.path(GSE50660_data,"betadata.models.adj.p.anot.sex.RData")) #per si de cas

# library(qqman) #provo amb un manhattan..mmmm
# beta.models.adj.p.anot$chr1 <- gsub("chr","",beta.models.adj.p.anot$chr)
# beta.models.adj.p.anot$chr1 <- as.numeric(ifelse(beta.models.adj.p.anot$chr1=="X",23,
#                                                  ifelse(beta.models.adj.p.anot$chr1=="Y",24,
#                                                         beta.models.adj.p.anot$chr1)))
# beta.models.adj.p.anot$p <- ifelse(beta.models.adj.p.anot$p.s<1e-10,1e-10,beta.models.adj.p.anot$p.s)
# 
# df4qqman <- beta.models.adj.p.anot[,c("genes","chr1","pos","pos","p")]
# manhattan(df4qqman, chr="chr1", bp="pos", snp="genes", p="p" )


# fn.karyo <- function(beta.models.adj.p.anot,p.th=0.05,cex=0.6){
#   #plot amb KaryoploteR
#   library(karyoploteR)
#   #i els colors que he fet servir a la fig 3
#   library(scales)
#   colors <- hue_pal()(6)
#   
#   # beta.models.adj.p.anot.p.s <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.s<p.th,]
#   # beta.models.adj.p.anot.p.b <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.b<p.th,]
#   # beta.models.adj.p.anot.p.sinf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.sinf<p.th,]
#   # beta.models.adj.p.anot.p.binf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.binf<p.th,]
#   # beta.models.adj.p.anot.p.n <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.n<p.th,]
#   # beta.models.adj.p.anot.p.l <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.l<p.th,]
#   # beta.models.adj.p.anot.p.q <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.q<p.th,]
#   # beta.models.adj.p.anot.p.limma <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.limma<p.th,]
#   
#   #preparo dades
#   df4kar <- beta.models.adj.p.anot[,c("chr","pos","p.b","p.s","p.n","p.l","p.q","p.limma")]
#   df4kar.1 <- df4kar[df4kar$p.b<p.th,c("chr","pos","p.b")]
#   df4kar.2 <- df4kar[df4kar$p.s<p.th,c("chr","pos","p.s")]
#   df4kar.3 <- df4kar[df4kar$p.n<p.th,c("chr","pos","p.n")]
#   df4kar.4 <- df4kar[df4kar$p.l<p.th,c("chr","pos","p.l")]
#   df4kar.5 <- df4kar[df4kar$p.q<p.th,c("chr","pos","p.q")]
#   df4kar.6 <- df4kar[df4kar$p.limma<p.th,c("chr","pos","p.limma")]
#   #plot
#   kp <- plotKaryotype(genome="hg19", plot.type=1)
#   kpDataBackground(kp, data.panel = 1, r0=0, r1=0.16)
#   kpDataBackground(kp, data.panel = 1, r0=0.33, r1=0.50)
#   kpDataBackground(kp, data.panel = 1, r0=.67, r1=0.84)
#   kpPoints(kp, chr=df4kar.1$chr, x=df4kar.1$pos, y=df4kar.1$p.b,data.panel = 1, r0=0,  r1=0.16,pch=".",col=colors[1])  
#   kpPoints(kp, chr=df4kar.2$chr, x=df4kar.2$pos, y=df4kar.2$p.s,data.panel = 1, r0=0.17,  r1=0.33,pch=".",col=colors[2])  
#   kpPoints(kp, chr=df4kar.3$chr, x=df4kar.3$pos, y=df4kar.3$p.n,data.panel = 1, r0=0.34,  r1=0.50,pch=".",col=colors[3])  
#   kpPoints(kp, chr=df4kar.4$chr, x=df4kar.4$pos, y=df4kar.4$p.l,data.panel = 1, r0=0.51,  r1=0.66,pch=".",col=colors[4])  
#   kpPoints(kp, chr=df4kar.5$chr, x=df4kar.5$pos, y=df4kar.5$p.q,data.panel = 1, r0=0.67,  r1=0.84,pch=".",col=colors[5])  
#   kpPoints(kp, chr=df4kar.6$chr, x=df4kar.6$pos, y=df4kar.6$p.limma,data.panel = 1, r0=0.85,  r1=1,pch=".",col=colors[6])  
# }  

#i ara només chrX: OLD amb els punts
# fn.chrx<- function(beta.models.adj.p.anot,p.th=0.05,cex=0.6){
#   library(karyoploteR)
#   library(scales)
#   colors <- hue_pal()(6)
#   e<-1e-10 #queden enganxats aqí pq hi ha molts 0's
#   #preparo dades
#   df4kar <- beta.models.adj.p.anot[,c("chr","pos","p.b","p.s","p.n","p.l","p.q","p.limma")]
#   df4kar.1 <- df4kar[df4kar$p.b<p.th,c("chr","pos","p.b")]
#   df4kar.2 <- df4kar[df4kar$p.s<p.th,c("chr","pos","p.s")]
#   df4kar.3 <- df4kar[df4kar$p.n<p.th,c("chr","pos","p.n")]
#   df4kar.4 <- df4kar[df4kar$p.l<p.th,c("chr","pos","p.l")]
#   df4kar.5 <- df4kar[df4kar$p.q<p.th,c("chr","pos","p.q")]
#   df4kar.6 <- df4kar[df4kar$p.limma<p.th,c("chr","pos","p.limma")]
#   
#   kpPoints(kp, chr=df4kar.1$chr, x=df4kar.1$pos, y=-log(df4kar.1$p.b+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.1$p.b+e))),data.panel = 1,r0=0,r1=0.16,pch=".",col=colors[1])  
#   kpPoints(kp, chr=df4kar.2$chr, x=df4kar.2$pos, y=-log(df4kar.2$p.s+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.2$p.s+e))),data.panel = 1,r0=0.17,r1=0.33,pch=".",col=colors[2])  
#   kpPoints(kp, chr=df4kar.3$chr, x=df4kar.3$pos, y=-log(df4kar.3$p.n+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.3$p.n+e))),data.panel = 1, r0=0.33,  r1=0.50,pch=".",col=colors[3])  
#   kpPoints(kp, chr=df4kar.4$chr, x=df4kar.4$pos, y=-log(df4kar.4$p.l+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.4$p.l+e))),data.panel = 1, r0=0.51,  r1=0.66,pch=".",col=colors[4])  
#   kpPoints(kp, chr=df4kar.5$chr, x=df4kar.5$pos, y=-log(df4kar.5$p.q+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.5$p.q+e))),data.panel = 1, r0=0.67,  r1=0.84,pch=".",col=colors[5])  
#   kpPoints(kp, chr=df4kar.6$chr, x=df4kar.6$pos, y=-log(df4kar.6$p.limma+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.6$p.limma+e))),data.panel = 1, r0=0.85,  r1=1,pch=".",col=colors[6])  
#   kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.1)), label.margin=0.01, r0=0,  r1=0.16,col=colors[1],cex=cex)
#   kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.2)), label.margin=0.01, r0=0.17,  r1=0.33,col=colors[2],cex=cex)
#   kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.3)), label.margin=0.01, r0=0.33,  r1=0.50,col=colors[3],cex=cex)
#   kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.4)), label.margin=0.01, r0=0.51,  r1=0.66,col=colors[4],cex=cex)
#   kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.5)), label.margin=0.01, r0=0.67,  r1=0.84,col=colors[5],cex=cex)
#   kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.6)), label.margin=0.01, r0=0.85,  r1=1,col=colors[6],cex=cex)
# 
# }

#1/5/19 ZOOM A X amb KaryoploteR, amb els density
# fn.chrx<- function(beta.models.adj.p.anot,p.th=0.05,cex=0.6){
#   library(karyoploteR)
#   library(scales)
#   colors <- hue_pal()(6)
#   #deixo el cex tot i que no el faig servir per res
#   e<-1e-10 #queden enganxats aqí pq hi ha molts 0's
#   #preparo dades, necessitaré un granges amb start i end
#   #ho converteixo directament
#   df4kar <- beta.models.adj.p.anot[!is.na(beta.models.adj.p.anot$chr),c("chr","pos","p.b","p.s","p.n","p.l","p.q","p.limma")]
#   #he d'eliminar els que tenen chr NA!
# 
#   df4kar.1 <- toGRanges(df4kar[!is.na(df4kar$p.b) & df4kar$p.b<p.th,c("chr","pos","pos","p.b")])
#   df4kar.2 <- toGRanges(df4kar[!is.na(df4kar$p.s) & df4kar$p.s<p.th,c("chr","pos","pos","p.s")])
#   df4kar.3 <- toGRanges(df4kar[!is.na(df4kar$p.n) & df4kar$p.n<p.th,c("chr","pos","pos","p.n")])
#   df4kar.4 <- toGRanges(df4kar[!is.na(df4kar$p.l) & df4kar$p.l<p.th,c("chr","pos","pos","p.l")])
#   df4kar.5 <- toGRanges(df4kar[!is.na(df4kar$p.q) & df4kar$p.q<p.th,c("chr","pos","pos","p.q")])
#   df4kar.6 <- toGRanges(df4kar[!is.na(df4kar$p.limma) & df4kar$p.limma<p.th,c("chr","pos","pos","p.limma")])
# 
#   kp <- plotKaryotype(genome="hg19", plot.type=1,chromosomes=c("chrX"))
#   kpPlotDensity(kp,data=df4kar.1, data.panel = 1,r0=0,r1=0.16,col=colors[1])
#   kpPlotDensity(kp,data=df4kar.2, data.panel = 1,r0=0.17,r1=0.33,col=colors[2])
#   kpPlotDensity(kp,data=df4kar.3, data.panel = 1,r0=0.34,r1=0.50,col=colors[3])
#   kpPlotDensity(kp,data=df4kar.4, data.panel = 1,r0=0.51,r1=0.66,col=colors[4])
#   kpPlotDensity(kp,data=df4kar.5, data.panel = 1,r0=0.67,r1=0.84,col=colors[5])
#   kpPlotDensity(kp,data=df4kar.6, data.panel = 1,r0=0.85,r1=1,col=colors[6])
# 
# }

#elimino logistica, he d'ajustar els tracks
fn.chrx<- function(beta.models.adj.p.anot,p.th=0.05,cex=0.6){
  library(karyoploteR)
  library(scales)
  colors <- hue_pal()(6)
  #deixo el cex tot i que no el faig servir per res
  e<-1e-10 #queden enganxats aqí pq hi ha molts 0's
  #preparo dades, necessitaré un granges amb start i end
  #ho converteixo directament
#  df4kar <- beta.models.adj.p.anot[!is.na(beta.models.adj.p.anot$chr),c("chr","pos","p.b","p.s","p.n","p.l","p.q","p.limma")]
  df4kar <- beta.models.adj.p.anot[!is.na(beta.models.adj.p.anot$chr),c("chr","pos","p.b","p.s","p.n","p.q","p.limma")]
  #he d'eliminar els que tenen chr NA!
  
  df4kar.1 <- toGRanges(df4kar[!is.na(df4kar$p.b) & df4kar$p.b<p.th,c("chr","pos","pos","p.b")])
  df4kar.2 <- toGRanges(df4kar[!is.na(df4kar$p.s) & df4kar$p.s<p.th,c("chr","pos","pos","p.s")])
  df4kar.3 <- toGRanges(df4kar[!is.na(df4kar$p.n) & df4kar$p.n<p.th,c("chr","pos","pos","p.n")])
#  df4kar.4 <- toGRanges(df4kar[!is.na(df4kar$p.l) & df4kar$p.l<p.th,c("chr","pos","pos","p.l")])
  df4kar.5 <- toGRanges(df4kar[!is.na(df4kar$p.q) & df4kar$p.q<p.th,c("chr","pos","pos","p.q")])
  df4kar.6 <- toGRanges(df4kar[!is.na(df4kar$p.limma) & df4kar$p.limma<p.th,c("chr","pos","pos","p.limma")])
  
  kp <- plotKaryotype(genome="hg19", plot.type=1,chromosomes=c("chrX"))
  kpPlotDensity(kp,data=df4kar.1, data.panel = 1,r0=0,r1=0.20,col=colors[1])
  kpPlotDensity(kp,data=df4kar.2, data.panel = 1,r0=0.21,r1=0.40,col=colors[2])
  kpPlotDensity(kp,data=df4kar.3, data.panel = 1,r0=0.41,r1=0.60,col=colors[3])
#  kpPlotDensity(kp,data=df4kar.4, data.panel = 1,r0=0.51,r1=0.66,col=colors[4])
  kpPlotDensity(kp,data=df4kar.5, data.panel = 1,r0=0.61,r1=0.80,col=colors[5])
  kpPlotDensity(kp,data=df4kar.6, data.panel = 1,r0=0.81,r1=1,col=colors[6])
  
}


load(file=file.path(GSE50660_data,"betadata.models.adj.p.anot.sex.RData")) 
# pdf(file=file.path(resultsDir,"Fig5.Sex.Karyo.10.min8.GSE50.pdf"))
# #png(file=file.path(resultsDir,"Fig5.Sex.Karyo.GSE50.png"))
#  fn.karyo(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
# dev.off()

pdf(file=file.path(resultsDir,"Fig6.Sex.ChrX.10.min8.dens.GSE50.050719.pdf"))
#png(file=file.path(resultsDir,"Fig5.Sex.ChrX.GSE50.png")) #no funciona, reportat a KaryoplotR
 fn.chrx(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()

#
load(file=file.path(GSE116339_data,"betadata.models.adj.p.anot.sex.RData")) 
# pdf(file=file.path(resultsDir,"Fig5.Sex.Karyo.10.min8.GSE11.pdf"))
#  fn.karyo(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
# dev.off()
pdf(file=file.path(resultsDir,"Fig6.Sex.ChrX.10.min8.dens.GSE11.050719.pdf"))
 fn.chrx(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()

load(file=file.path(WGBS81_data,"rnb.meth2comp.models.withlimma.RData"))
#head(rnb.meth2comp.models.withlimma)
#ho converteixo a format similar
annot <- rownames(rnb.meth2comp.models.withlimma)
cpg_chr <- sapply(strsplit(annot,split="_"),function(x) x[1])
cpg_pos <- sapply(strsplit(annot,split="_"),function(x) x[2])
cpg_type <- sapply(strsplit(annot,split="_"),function(x) x[3])

#li dono el mateix nom: beta.models.adj.p.anot
beta.models.adj.p.anot <- data.frame(rnb.meth2comp.models.withlimma,chr=cpg_chr,pos=as.numeric(cpg_pos))
# pdf(file=file.path(resultsDir,"Fig5.Sex.Karyo.10.min8.WGBS.pdf"))
#  fn.karyo(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
# dev.off()
pdf(file=file.path(resultsDir,"Fig6.Sex.ChrX.10.min8.dens.WGBS.050719.pdf"))
 fn.chrx(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()

load(file=file.path(RRBS188_data,"rnb.meth2comp.models.withlimma.RData"))
#head(rnb.meth2comp.models.withlimma)
#ho converteixo a format similar
annot <- rownames(rnb.meth2comp.models.withlimma)
cpg_chr <- sapply(strsplit(annot,split="_"),function(x) x[1])
cpg_pos <- sapply(strsplit(annot,split="_"),function(x) x[2])
cpg_type <- sapply(strsplit(annot,split="_"),function(x) x[3])

#li dono el mateix nom: beta.models.adj.p.anot
beta.models.adj.p.anot <- data.frame(rnb.meth2comp.models.withlimma,chr=cpg_chr,pos=as.numeric(cpg_pos))
# pdf(file=file.path(resultsDir,"Fig5.Sex.Karyo.10.min8.RRBS188.pdf"))
#  fn.karyo(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
# dev.off()
pdf(file=file.path(resultsDir,"Fig6.Sex.ChrX.10.min8.dens.RRBS188.050719.pdf"))
 fn.chrx(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()


#CIRCOS FIG PPAL
library(OmicCircos)
#he de fer 6 tracks per data set

library(scales)
colors <- hue_pal()(6)

#5/7/19: elimino logistica
fn.generate.tracks <- function(beta.models.adj.p.anot,p.th=0.05){
  #preparo dades
  colnms <- c("chr","Start","End","Value")
#  df4kar <- beta.models.adj.p.anot[,c("chr","pos","p.b","p.s","p.n","p.l","p.q","p.limma")]
  df4kar <- beta.models.adj.p.anot[,c("chr","pos","p.b","p.s","p.n","p.q","p.limma")]
  track<-list(NA)
  #no se pq a gse11 hi ha alguns que no tenen chr
  track[[1]] <- df4kar[!is.na(df4kar$p.b) & !is.na(df4kar$chr) & df4kar$p.b<p.th,c("chr","pos","pos","p.b")]
  track[[2]] <- df4kar[!is.na(df4kar$p.s) & !is.na(df4kar$chr) & df4kar$p.s<p.th,c("chr","pos","pos","p.s")]
  track[[3]] <- df4kar[!is.na(df4kar$p.n) & !is.na(df4kar$chr) & df4kar$p.n<p.th,c("chr","pos","pos","p.n")]
 # track[[4]] <- df4kar[!is.na(df4kar$p.l) & !is.na(df4kar$chr) & df4kar$p.l<p.th,c("chr","pos","pos","p.l")]
  track[[4]] <- df4kar[!is.na(df4kar$p.q) & !is.na(df4kar$chr) & df4kar$p.q<p.th,c("chr","pos","pos","p.q")]
  track[[5]] <- df4kar[!is.na(df4kar$p.limma) & !is.na(df4kar$chr) & df4kar$p.limma<p.th,c("chr","pos","pos","p.limma")]
  for (i in 1:5) colnames(track[[i]]) <- colnms
  return(track)
}
# fn.generate.tracks.X <- function(beta.models.adj.p.anot,p.th=0.05){
#   #preparo dades
#   colnms <- c("chr","Start","End","Value")
#   df4kar <- beta.models.adj.p.anot[,c("chr","pos","p.b","p.s","p.n","p.l","p.q","p.limma")]
#   track<-list(NA)
#   #no se pq a gse11 hi ha alguns que no tenen chr
#   track[[1]] <- df4kar[!is.na(df4kar$p.b) & !is.na(df4kar$chr) & df4kar$p.b<p.th & df4kar$chr=="chrX",c("chr","pos","pos","p.b")]
#   track[[2]] <- df4kar[!is.na(df4kar$p.s) & !is.na(df4kar$chr) & df4kar$p.s<p.th & df4kar$chr=="chrX",c("chr","pos","pos","p.b")]
#   track[[3]] <- df4kar[!is.na(df4kar$p.n) & !is.na(df4kar$chr) & df4kar$p.n<p.th & df4kar$chr=="chrX",c("chr","pos","pos","p.b")]
#   track[[4]] <- df4kar[!is.na(df4kar$p.l) & !is.na(df4kar$chr) & df4kar$p.l<p.th & df4kar$chr=="chrX",c("chr","pos","pos","p.b")]
#   track[[5]] <- df4kar[!is.na(df4kar$p.q) & !is.na(df4kar$chr) & df4kar$p.q<p.th & df4kar$chr=="chrX",c("chr","pos","pos","p.b")]
#   track[[6]] <- df4kar[!is.na(df4kar$p.limma) & !is.na(df4kar$chr) & df4kar$p.limma<p.th & df4kar$chr=="chrX",c("chr","pos","pos","p.b")]
#   for (i in 1:6) colnames(track[[i]]) <- colnms
#   return(track)
# }

p.th=10^-8
load(file=file.path(GSE50660_data,"betadata.models.adj.p.anot.sex.RData"))   
track.gse50 <- fn.generate.tracks(beta.models.adj.p.anot,p.th=p.th)
#track.gse50.X <- fn.generate.tracks.X(beta.models.adj.p.anot,p.th=p.th)

load(file=file.path(GSE116339_data,"betadata.models.adj.p.anot.sex.RData")) 
track.gse11 <- fn.generate.tracks(beta.models.adj.p.anot,p.th=p.th)
#track.gse11.X <- fn.generate.tracks.X(beta.models.adj.p.anot,p.th=p.th)

load(file=file.path(RRBS188_data,"rnb.meth2comp.models.withlimma.RData"))
annot <- rownames(rnb.meth2comp.models.withlimma)
cpg_chr <- sapply(strsplit(annot,split="_"),function(x) x[1])
cpg_pos <- sapply(strsplit(annot,split="_"),function(x) x[2])
cpg_type <- sapply(strsplit(annot,split="_"),function(x) x[3])
beta.models.adj.p.anot <- data.frame(rnb.meth2comp.models.withlimma,chr=cpg_chr,pos=as.numeric(cpg_pos))
track.rrbs <- fn.generate.tracks(beta.models.adj.p.anot,p.th=p.th)
#track.rrbs.X <- fn.generate.tracks.X(beta.models.adj.p.anot,p.th=p.th)

load(file=file.path(WGBS81_data,"rnb.meth2comp.models.withlimma.RData"))
annot <- rownames(rnb.meth2comp.models.withlimma)
cpg_chr <- sapply(strsplit(annot,split="_"),function(x) x[1])
cpg_pos <- sapply(strsplit(annot,split="_"),function(x) x[2])
cpg_type <- sapply(strsplit(annot,split="_"),function(x) x[3])
beta.models.adj.p.anot <- data.frame(rnb.meth2comp.models.withlimma,chr=cpg_chr,pos=as.numeric(cpg_pos))
track.wgbs <- fn.generate.tracks(beta.models.adj.p.anot,p.th=p.th)
#track.wgbs.X <- fn.generate.tracks.X(beta.models.adj.p.anot,p.th=p.th)

#revisar: 
# 1.hg38/19: hg38 no hi és...
# 2. Zoom chrx en una altra figura per composar
# 3. Fons de cada color del data set
# 4. provar de posar lletres (per cada dataset) o nums amb els totals
# 5. el mateix a 10^-8

#per a fer el bg de cada track
data("UCSC.hg19.chr")
ref <- UCSC.hg19.chr;
ref.d <- c();
for (chr in c(1:22, "X", "Y")){
  chr.s <- paste0("chr", chr);
  ref.i <- which(ref[,1]==chr.s);
  ref.s <- ref[ref.i,];
  ref.d <- rbind(ref.d, c(chr.s, 1, ref.s[length(ref.i),3]))
}

#elimino logística, com que era el 4rt element elimino el 4rt color
library(scales)
colors <- hue_pal()(6)[c(1:3,5:6)]

pdf(file=file.path(resultsDir,"Fig6.Circos.bg.10menys8.050719.pdf")) 
    col.bg <- c("darkseagreen1","lightskyblue1","rosybrown1","lightgoldenrod")
    #  col.bg <- rep("white",4)
    par(mar=c(2, 1, 2, 2));
    plot(c(0,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="");
    circos(R=400, cir="hg19", W=10,   type="chr", print.chr.lab=T, scale=F);
    circos(R=355, cir="hg19", W=30, mapping=ref.d, type="arc2",  B=F, col= col.bg[1], lwd=25, scale=F);
    for (i in 1:5){
      circos(R=390-7*i, cir="hg19", W=4,  mapping=track.gse50[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
    }
   
    circos(R=315, cir="hg19", W=30, mapping=ref.d, type="arc2",  B=F, col=col.bg[2], lwd=25, scale=F);
    for (i in 1:5){
      circos(R=350-7*i, cir="hg19", W=4,  mapping=track.gse11[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
    }
   
     circos(R=275, cir="hg19", W=30, mapping=ref.d, type="arc2",  B=F, col=col.bg[3], lwd=25, scale=F);
    for (i in 1:5){
      circos(R=310-7*i, cir="hg19", W=4,  mapping=track.rrbs[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
    }
    
     circos(R=235, cir="hg19", W=30, mapping=ref.d, type="arc2",  B=F, col=col.bg[4], lwd=25, scale=F);
    for (i in c(1:5)){
      circos(R=270-7*i, cir="hg19", W=4,  mapping=track.wgbs[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
    }
dev.off()  


# pdf(file=file.path(resultsDir,"Fig6.Circos.bg.10menys8.120519.pdf")) #l'ultim dona error track.wgbs[[4]] no te elements
# #pdf(file=file.path(resultsDir,"Fig5.Circos.bgwhite.10menys8.pdf")) #l'ultim dona error track.wgbs[[4]] no te elements
# # color del background de cada pista es pot fer amb h enlloc de s però tampoc es un density 
#   col.bg <- c("darkseagreen1","lightskyblue1","rosybrown1","lightgoldenrod")
# #  col.bg <- rep("white",4)
#   par(mar=c(2, 1, 2, 2));
#   plot(c(0,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="");
#   circos(R=400, cir="hg19", W=10,   type="chr", print.chr.lab=T, scale=F);
#   circos(R=360, cir="hg19", W=25, mapping=ref.d, type="arc2",  B=F, col= col.bg[1], lwd=25, scale=F);
#   for (i in 1:6){
#     circos(R=390-5*i, cir="hg19", W=4,  mapping=track.gse50[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
#   }
#   circos(R=320, cir="hg19", W=25, mapping=ref.d, type="arc2",  B=F, col=col.bg[2], lwd=25, scale=F);
#   for (i in 1:6){
#     circos(R=350-5*i, cir="hg19", W=4,  mapping=track.gse11[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
#   }
#   circos(R=280, cir="hg19", W=25, mapping=ref.d, type="arc2",  B=F, col=col.bg[3], lwd=25, scale=F);
#   for (i in 1:6){
#     circos(R=310-5*i, cir="hg19", W=4,  mapping=track.rrbs[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
#   }
#   circos(R=240, cir="hg19", W=25, mapping=ref.d, type="arc2",  B=F, col=col.bg[4], lwd=25, scale=F);
#   for (i in c(1:3,5:6)){
#     circos(R=270-5*i, cir="hg19", W=4,  mapping=track.wgbs[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
#   }
# dev.off()  
  
#per al powerpoint: NO HO NECESSITO PQ HO FAIG AMB EL CUENTAGOTAS!
col2rgb(col.bg)
#       [,1] [,2] [,3] [,4]
# red    193  176  255  238
# green  255  226  193  221
# blue   193  255  193  130
#provo de fer el mateix a X, sembla que no es pot fer amb un zoom, està reportat...https://www.biostars.org/p/129864/
# par(mar=c(2, 2, 2, 2));
# plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="");
# zoom <- c("chr1", "chr1", 1,155270000, 0, 360);
# circos(R=400, cir="hg19", W=4,   type="chr", print.chr.lab=T, scale=T, zoom=zoom);
#data(UCSC.hg19) #aquesta ?s la bona, d'aqu? treurem X i Y i recalculem angle, per? sembla que no l'enten per tant no tinc clar que sigui aquest el bo
# data(UCSC.hg19.chr)
# chrX.i <- which(UCSC.hg19.chr[,1]=="chrX");
# chrX   <- UCSC.hg19.chr[chrX.i,];
# ## segment data of chromosome X
# seg.c   <- segAnglePo(chrX, seg="chrX");
# 
# pdf(file=file.path(resultsDir,"Fig5.Circos.X.bgwhite.10menys8.pdf")) #l'ultim dona error track.wgbs[[4]] no te elements
# pdf(file=file.path(resultsDir,"Fig5.Circos.X.bg.10menys8.pdf")) #l'ultim dona error track.wgbs[[4]] no te elements
#   col.bg <- c("darkseagreen1","lightskyblue1","rosybrown1","lightgoldenrod")
#   par(mar=c(2, 2, 2, 2))
#   plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="");
#   circos(R=400, type="chr2", cir=seg.c, mapping=chrX, print.chr.lab=TRUE, W=4, scale=T);
#   
#   circos(R=360, cir=seg.c, W=25, mapping=chrX, type="arc2",  B=F, col= col.bg[1], lwd=25, scale=F);
#   for (i in 1:6){
#     circos(R=390-5*i, cir=seg.c, W=4,  mapping=track.gse50.X[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
#   }
#   circos(R=320, cir=seg.c, W=25, mapping=chrX, type="arc2",  B=F, col=col.bg[2], lwd=25, scale=F);
#   for (i in 1:6){
#     circos(R=350-5*i, cir=seg.c, W=4,  mapping=track.gse11.X[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
#   }
#   circos(R=280, cir=seg.c, W=25, mapping=chrX, type="arc2",  B=F, col=col.bg[3], lwd=25, scale=F);
#   for (i in 1:6){
#     circos(R=310-5*i, cir=seg.c, W=4,  mapping=track.rrbs.X[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
#   }
#   circos(R=240, cir=seg.c, W=25, mapping=chrX, type="arc2",  B=F, col=col.bg[4], lwd=25, scale=F);
#   for (i in c(1:3,5:6)){
#     circos(R=270-5*i, cir=seg.c, W=4,  mapping=track.wgbs.X[[i]],   col.v=4,    type="s",B=FALSE, lwd=0.1, col=colors[i],cex=0.1,cutoff=0);
#   }
# dev.off()  
