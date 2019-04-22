#Figure 5 per substituir els Venn
#KaryoplotR with results, all genome and zoom in chrX

#header amb tot
source("D:/Doctorat/Simplex/R/Figures.HeaderScript.R")

load(file=file.path(GSE50660_data,"betadata.models.adj.p.anot.sex.RData")) #per si de cas

library(qqman) #provo amb un manhattan..mmmm
beta.models.adj.p.anot$chr1 <- gsub("chr","",beta.models.adj.p.anot$chr)
beta.models.adj.p.anot$chr1 <- as.numeric(ifelse(beta.models.adj.p.anot$chr1=="X",23,
                                                 ifelse(beta.models.adj.p.anot$chr1=="Y",24,
                                                        beta.models.adj.p.anot$chr1)))
beta.models.adj.p.anot$p <- ifelse(beta.models.adj.p.anot$p.s<1e-10,1e-10,beta.models.adj.p.anot$p.s)

df4qqman <- beta.models.adj.p.anot[,c("genes","chr1","pos","pos","p")]
manhattan(df4qqman, chr="chr1", bp="pos", snp="genes", p="p" )


fn.karyo <- function(beta.models.adj.p.anot,p.th=0.05,cex=0.6){
  #plot amb KaryoploteR
  library(karyoploteR)
  #i els colors que he fet servir a la fig 3
  library(scales)
  colors <- hue_pal()(6)
  
  # beta.models.adj.p.anot.p.s <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.s<p.th,]
  # beta.models.adj.p.anot.p.b <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.b<p.th,]
  # beta.models.adj.p.anot.p.sinf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.sinf<p.th,]
  # beta.models.adj.p.anot.p.binf <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.binf<p.th,]
  # beta.models.adj.p.anot.p.n <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.n<p.th,]
  # beta.models.adj.p.anot.p.l <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.l<p.th,]
  # beta.models.adj.p.anot.p.q <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.q<p.th,]
  # beta.models.adj.p.anot.p.limma <- beta.models.adj.p.anot[beta.models.adj.p.anot$p.limma<p.th,]
  
  #preparo dades
  df4kar <- beta.models.adj.p.anot[,c("chr","pos","p.b","p.s","p.n","p.l","p.q","p.limma")]
  df4kar.1 <- df4kar[df4kar$p.b<p.th,c("chr","pos","p.b")]
  df4kar.2 <- df4kar[df4kar$p.s<p.th,c("chr","pos","p.s")]
  df4kar.3 <- df4kar[df4kar$p.n<p.th,c("chr","pos","p.n")]
  df4kar.4 <- df4kar[df4kar$p.l<p.th,c("chr","pos","p.l")]
  df4kar.5 <- df4kar[df4kar$p.q<p.th,c("chr","pos","p.q")]
  df4kar.6 <- df4kar[df4kar$p.limma<p.th,c("chr","pos","p.limma")]
  #plot
  kp <- plotKaryotype(genome="hg19", plot.type=1)
  kpDataBackground(kp, data.panel = 1, r0=0, r1=0.16)
  kpDataBackground(kp, data.panel = 1, r0=0.33, r1=0.50)
  kpDataBackground(kp, data.panel = 1, r0=.67, r1=0.84)
  kpPoints(kp, chr=df4kar.1$chr, x=df4kar.1$pos, y=df4kar.1$p.b,data.panel = 1, r0=0,  r1=0.16,pch=".",col=colors[1])  
  kpPoints(kp, chr=df4kar.2$chr, x=df4kar.2$pos, y=df4kar.2$p.s,data.panel = 1, r0=0.17,  r1=0.33,pch=".",col=colors[2])  
  kpPoints(kp, chr=df4kar.3$chr, x=df4kar.3$pos, y=df4kar.3$p.n,data.panel = 1, r0=0.34,  r1=0.50,pch=".",col=colors[3])  
  kpPoints(kp, chr=df4kar.4$chr, x=df4kar.4$pos, y=df4kar.4$p.l,data.panel = 1, r0=0.51,  r1=0.66,pch=".",col=colors[4])  
  kpPoints(kp, chr=df4kar.5$chr, x=df4kar.5$pos, y=df4kar.5$p.q,data.panel = 1, r0=0.67,  r1=0.84,pch=".",col=colors[5])  
  kpPoints(kp, chr=df4kar.6$chr, x=df4kar.6$pos, y=df4kar.6$p.limma,data.panel = 1, r0=0.85,  r1=1,pch=".",col=colors[6])  
}  

#i ara només chrX
fn.chrx<- function(beta.models.adj.p.anot,p.th=0.05,cex=0.6){
  library(karyoploteR)
  library(scales)
  colors <- hue_pal()(6)
  e<-1e-10 #queden enganxats aqí pq hi ha molts 0's
  #preparo dades
  df4kar <- beta.models.adj.p.anot[,c("chr","pos","p.b","p.s","p.n","p.l","p.q","p.limma")]
  df4kar.1 <- df4kar[df4kar$p.b<p.th,c("chr","pos","p.b")]
  df4kar.2 <- df4kar[df4kar$p.s<p.th,c("chr","pos","p.s")]
  df4kar.3 <- df4kar[df4kar$p.n<p.th,c("chr","pos","p.n")]
  df4kar.4 <- df4kar[df4kar$p.l<p.th,c("chr","pos","p.l")]
  df4kar.5 <- df4kar[df4kar$p.q<p.th,c("chr","pos","p.q")]
  df4kar.6 <- df4kar[df4kar$p.limma<p.th,c("chr","pos","p.limma")]
  
  kp <- plotKaryotype(genome="hg19", plot.type=1,chromosomes=c("chrX"))
  kpPoints(kp, chr=df4kar.1$chr, x=df4kar.1$pos, y=-log(df4kar.1$p.b+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.1$p.b+e))),data.panel = 1,r0=0,r1=0.16,pch=".",col=colors[1])  
  kpPoints(kp, chr=df4kar.2$chr, x=df4kar.2$pos, y=-log(df4kar.2$p.s+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.2$p.s+e))),data.panel = 1,r0=0.17,r1=0.33,pch=".",col=colors[2])  
  kpPoints(kp, chr=df4kar.3$chr, x=df4kar.3$pos, y=-log(df4kar.3$p.n+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.3$p.n+e))),data.panel = 1, r0=0.33,  r1=0.50,pch=".",col=colors[3])  
  kpPoints(kp, chr=df4kar.4$chr, x=df4kar.4$pos, y=-log(df4kar.4$p.l+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.4$p.l+e))),data.panel = 1, r0=0.51,  r1=0.66,pch=".",col=colors[4])  
  kpPoints(kp, chr=df4kar.5$chr, x=df4kar.5$pos, y=-log(df4kar.5$p.q+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.5$p.q+e))),data.panel = 1, r0=0.67,  r1=0.84,pch=".",col=colors[5])  
  kpPoints(kp, chr=df4kar.6$chr, x=df4kar.6$pos, y=-log(df4kar.6$p.limma+e), ymin=-log(p.th), ymax=ceiling(max(-log(df4kar.6$p.limma+e))),data.panel = 1, r0=0.85,  r1=1,pch=".",col=colors[6])  
  kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.1)), label.margin=0.01, r0=0,  r1=0.16,col=colors[1],cex=cex)
  kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.2)), label.margin=0.01, r0=0.17,  r1=0.33,col=colors[2],cex=cex)
  kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.3)), label.margin=0.01, r0=0.33,  r1=0.50,col=colors[3],cex=cex)
  kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.4)), label.margin=0.01, r0=0.51,  r1=0.66,col=colors[4],cex=cex)
  kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.5)), label.margin=0.01, r0=0.67,  r1=0.84,col=colors[5],cex=cex)
  kpAddLabels(kp, labels=paste0("N=",nrow(df4kar.6)), label.margin=0.01, r0=0.85,  r1=1,col=colors[6],cex=cex)

}

load(file=file.path(GSE50660_data,"betadata.models.adj.p.anot.sex.RData")) 
res.chrx.GSE50 <- beta.models.adj.p.anot[beta.models.adj.p.anot$chr=="chrX",]
dim(res.chrx.GSE50) #potser hauria de buscar la proporcio de cada model signif
load(file=file.path(GSE50660_data,"betadata.models.sex.res.table.RData")) 
res.table.GSE50 <-res.table

pdf(file=file.path(resultsDir,"Fig5.Sex.Karyo.10.min8.GSE50.pdf"))
#png(file=file.path(resultsDir,"Fig5.Sex.Karyo.GSE50.png"))
fn.karyo(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()

pdf(file=file.path(resultsDir,"Fig5.Sex.ChrX.10.min8.GSE50.pdf"))
#png(file=file.path(resultsDir,"Fig5.Sex.ChrX.GSE50.png")) #no funciona, reportat a KaryoplotR
fn.chrx(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()

#aquestes no són!!
load(file=file.path(GSE116339_data,"betadata.models.adj.p.anot.sex.RData")) 
pdf(file=file.path(resultsDir,"Fig5.Sex.Karyo.10.min8.GSE11.pdf"))
fn.karyo(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()
pdf(file=file.path(resultsDir,"Fig5.Sex.ChrX.10.min8.GSE11.pdf"))
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
pdf(file=file.path(resultsDir,"Fig5.Sex.Karyo.10.min8.WGBS.pdf"))
fn.karyo(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()
pdf(file=file.path(resultsDir,"Fig5.Sex.ChrX.10.min8.WGBS.pdf"))
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
pdf(file=file.path(resultsDir,"Fig5.Sex.Karyo.10.min8.RRBS188.pdf"))
fn.karyo(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()
pdf(file=file.path(resultsDir,"Fig5.Sex.ChrX.10.min8.RRBS188.pdf"))
fn.chrx(beta.models.adj.p.anot,p.th=1e-8,cex=0.6)
dev.off()



# png(file=file.path(resultsDir,"Fig5.Sex.ChrX.png")) #no funciona, reportat a KaryoplotR ho he
# de fer per a cada data set
# layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
#   fn.chrx(beta.models.adj.p.anot)
#   fn.chrx(beta.models.adj.p.anot)
#   fn.chrx(beta.models.adj.p.anot)
#   fn.chrx(beta.models.adj.p.anot)
# dev.off()

