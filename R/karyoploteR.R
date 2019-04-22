# R 3.4.3
# 24/1/19 Lara
# Proves amb el paquet
# vull fer uns gràfics pel projecte de la Blanca del panell, per a generar imatges de regions "a demanda"
# Sembla que el plot.type=5 no va bé i que quan fas l'axis, no funciona si fas un zoom
# l'nstal·lo des de bioc però l'última versio, doncs sembla que el problema de l'axis està resolt
library(BiocManager)
BiocManager::install("karyoploteR")

wd<-"Z:/R Code/karyoploteR"
setwd(wd)

library(karyoploteR)

#el primer és pintar el cariotip, escolir regions, pden ser vàries i tb fer zoom

#plot.type
kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes=c("chr1")) #Horizontal ideograms with a single data panel above them
kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chr1")) #Horizontal ideograms with two data panels, one above and one below them
kp <- plotKaryotype(genome="hg19", plot.type=3, chromosomes=c("chr1")) #Horizontal ideograms with all chromosomes in a single line with two data panels, one above and one below them
kp <- plotKaryotype(genome="hg19", plot.type=4, chromosomes=c("chr1")) #Horizontal ideograms with all chromosomes in a single line with one data panel above
kp <- plotKaryotype(genome="hg19", plot.type=5, chromosomes=c("chr1")) #Horizontal ideograms with all chromosomes in a single line with one data panel below them

kpAddBaseNumbers(kp) #add posicio bases
kpAddCytobandLabels(kp) #add cytobands!!

#puc fer zoom i tot el plot quedarà en aquesta regió, definir una regió en un GRanges
(gr <-GRanges(seqnames ="chr1", ranges =  IRanges(start=1,end=100000000)))
kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chr1"), zoom=gr)

#i ara es poden afegir tracks, com sempre
kpPoints(kp, chr = data.points$chr, x=data.points$pos, y=data.points$value, col=rainbow(240))


#exemple tutorial

rand.data <- createRandomRegions(genome="hg19", nregions=1000, length.mean=1, length.sd=0,
                                 mask=NA, non.overlapping=TRUE) 

## Attaching package: 'Biostrings'
## The following object is masked from 'package:base':
## 
##     strsplit
rand.data <- toDataframe(sort(rand.data))
rand.data <- cbind(rand.data, y=runif(n=1000, min=-1, max=1))

#Select some data points as "special ones"
sel.data <- rand.data[c(7, 30, 38, 52),] 
head(rand.data)

rand.data.chr1 <- rand.data[rand.data$chr=="chr1",]
isSorted(rand.data.chr1$start) #yes! 

#proves doncs tinc la sostipta que algun plot.type no xuta, efectivament plot.type=5 no va 
pdf("karyoploteR.DifferentTypes.pdf")
#1
kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes=c("chr1", "chr2", "chr3"))
kpDataBackground(kp, data.panel = 1, r0=0, r1=0.45) #important especificar r0 i r1 doncs és l'alçada del track
kpAxis(kp, ymin=-1, ymax=1, r0=0.05, r1=0.4, col="gray50", cex=0.5) #en aquest exemple no especifica el panell
kpPoints(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="black", pch=".", cex=1.5)
kpPoints(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="red")
kpText(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
       ymin=-1, ymax=1, r0=0.05, r1=0.4, labels=c("A", "B", "C", "D"), col="red",
       pos=4, cex=0.8)
kpAbline(kp, h=c(-1,0,1), col="#DDAAAA", ymin=-1, ymax=1, r0=0.05, r1=0.4, lty=2)
#2
kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chr1", "chr2", "chr3"))
kpDataBackground(kp, data.panel = 1, r0=0, r1=0.45) #
kpAxis(kp, ymin=-1, ymax=1, r0=0.05, r1=0.4, col="gray50", cex=0.5) #en aquest exemple no especifica el panell
kpPoints(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="black", pch=".", cex=2)
kpPoints(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="red")
kpText(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
       ymin=-1, ymax=1, r0=0.05, r1=0.4, labels=c("A", "B", "C", "D"), col="red",
       pos=4, cex=0.8)
kpAbline(kp, h=0, col="#DDAAAA", ymin=-1, ymax=1, r0=0.05, r1=0.4, lty=2)
#3
kp <- plotKaryotype(genome="hg19", plot.type=3, chromosomes=c("chr1", "chr2", "chr3"))
kpDataBackground(kp, data.panel = 1, r0=0, r1=0.45) #
kpAxis(kp, ymin=-1, ymax=1, r0=0.05, r1=0.4, col="gray50", cex=0.5) #en aquest exemple no especifica el panell
kpPoints(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="black", pch=".", cex=2)
kpPoints(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="red")
kpText(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
       ymin=-1, ymax=1, r0=0.05, r1=0.4, labels=c("A", "B", "C", "D"), col="red",
       pos=4, cex=0.8)
kpAbline(kp, h=0, col="#DDAAAA", ymin=-1, ymax=1, r0=0.05, r1=0.4, lty=2)
#4
kp <- plotKaryotype(genome="hg19", plot.type=4, chromosomes=c("chr1", "chr2", "chr3"))
kpDataBackground(kp, data.panel = 1, r0=0, r1=0.45) #
kpAxis(kp, ymin=-1, ymax=1, r0=0.05, r1=0.4, col="gray50", cex=0.5) #en aquest exemple no especifica el panell
kpPoints(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="black", pch=".", cex=2)
kpPoints(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="red")
kpText(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
       ymin=-1, ymax=1, r0=0.05, r1=0.4, labels=c("A", "B", "C", "D"), col="red",
       pos=4, cex=0.8)
kpAbline(kp, h=0, col="#DDAAAA", ymin=-1, ymax=1, r0=0.05, r1=0.4, lty=2)
#5: és el que sembla que està malament
kp <- plotKaryotype(genome="hg19", plot.type=5, chromosomes=c("chr1", "chr2", "chr3"))
kpDataBackground(kp, data.panel = 1, r0=0, r1=0.45) #
kpAxis(kp, ymin=-1, ymax=1, r0=0.05, r1=0.4, col="gray50", cex=0.5) #en aquest exemple no especifica el panell
kpPoints(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="black", pch=".", cex=2)
kpPoints(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="red")
kpText(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
       ymin=-1, ymax=1, r0=0.05, r1=0.4, labels=c("A", "B", "C", "D"), col="red",
       pos=4, cex=0.8)
kpAbline(kp, h=0, col="#DDAAAA", ymin=-1, ymax=1, r0=0.05, r1=0.4, lty=2)
dev.off()


#aquí sembla que l'axis no el mostra bé
(gr <-GRanges(seqnames ="chr1", ranges =  IRanges(start=1,end=400000000)))
kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chr1"), zoom=gr)
kpDataBackground(kp, data.panel = 1, r0=0, r1=0.45) #
kpAxis(kp, ymin=-1, ymax=1, r0=0.05, r1=0.4, col="gray50", cex=0.5) #en aquest exemple no especifica el panell
kpPoints(kp, chr=rand.data.chr1$chr, x=rand.data.chr1$start, y=rand.data.chr1$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="black", pch=".", cex=2)
kpPoints(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, col="red")
kpText(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
       ymin=-1, ymax=1, r0=0.05, r1=0.4, labels=c("A", "B", "C", "D"), col="red",
       pos=4, cex=0.8)
kpAbline(kp, h=0, col="#DDAAAA", ymin=-1, ymax=1, r0=0.05, r1=0.4, lty=2)


############## això és de la vignette!
#Upper part: data.panel=1
kpDataBackground(kp, data.panel = 1, r0=0.5, r1=1)
kpAxis(kp, ymin=-1, ymax=1, r0=0.5, r1=1, col="gray50", cex=0.5, numticks = 5)
kpAbline(kp, h=c(-0.5, 0, 0.5), col="gray50", ymin=-1, ymax=1, r0=0.5, r1=1)
kpLines(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
        col="#AA88FF", ymin=-1, ymax=1, r0=0.5, r1=1)
#Use kpSegments to add small tic to the line
kpSegments(kp, chr=rand.data$chr, x0=rand.data$start, x1=rand.data$start,
           y0=rand.data$y-0.1, y1=rand.data$y+0.1,
           col="#8866DD", ymin=-1, ymax=1, r0=0.5, r1=1)
#Plot the same line but inverting the data by pssing a r0>r1
kpLines(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
        col="#FF88AA", ymin=-1, ymax=1, r0=1, r1=0.5)


#Lower part: data.panel=2
kpDataBackground(kp, r0=0, r1=0.29, color = "#EEFFEE", data.panel = 2)
kpAxis(kp, col="#AADDAA", ymin=-1, ymax=1, r0=0, r1=0.29, data.panel = 2,
       numticks = 2, cex=0.5, tick.len = 0)
kpAbline(kp, h=0, col="#AADDAA", ymin=-1, ymax=1, r0=0, r1=0.29, data.panel = 2)
kpBars(kp, chr=rand.data$chr, x0=rand.data$start, x1=rand.data$end, y1 = rand.data$y,
       col="#AADDAA", ymin=-1, ymax=1, r0=0, r1=0.29, data.panel = 2, border="#AADDAA" )

kpDataBackground(kp, r0=0.34, r1=0.63, color = "#EEEEFF", data.panel = 2)
kpAxis(kp, col="#AAAADD", ymin=-1, ymax=1, r0=0.34, r1=0.63, data.panel = 2, 
       numticks = 2, cex=0.5, tick.len = 0)
kpAbline(kp, h=0, col="#AAAADD", ymin=-1, ymax=1, r0=0.34, r1=0.63, data.panel = 2)
kpSegments(kp, chr=rand.data$chr, x0=rand.data$start, x1=rand.data$end, 
           y0=rand.data$y-0.2, y1=rand.data$y, 
           col="#AAAADD", ymin=-1, ymax=1, r0=0.34, r1=0.63, data.panel = 2, lwd=2)

kpDataBackground(kp, r0=0.68, r1=0.97, color = "#FFEEEE", data.panel = 2)
kpAxis(kp, col="#DDAAAA", ymin=-1, ymax=1, r0=0.68, r1=0.97, data.panel = 2,
       numticks = 2, cex=0.5, tick.len = 0)
kpPoints(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
         col="#DDAAAA", ymin=-1, ymax=1, r0=0.68, r1=0.97, data.panel = 2, pch=".", cex=3)


###################### BAM FILES
library(pasillaBamSubset) #A package with 2 example bam files
un1.bam.file <- untreated1_chr4() # get the name of the first bam
window.size <- 1e4 #compute the density with 10kb windows
kp <- plotKaryotype(genome="dm6", chromosomes="chr4") #The pasilla data comes from drosophila
kp <- kpPlotBAMDensity(kp, data = un1.bam.file, window.size = window.size, r0=0.5, r1=1, ymax=50000, col="darkorange")


#Or normalizing by the number of mapped reads
kp <- plotKaryotype(genome="dm6", chromosomes="chr4") #The pasilla data comes from drosophila
kp <- kpAddBaseNumbers(kp, tick.dist = 1e5)
kp <- kpPlotBAMDensity(kp, data = un1.bam.file, window.size = window.size, normalize=TRUE, r0=0.5, r1=1, ymax=0.2, col="darkorange")

kp <- plotKaryotype(genome="hg19",  chromosomes="chr1")
kp <- kpPlotBAMDensity(kp, data = file.path(bamDir,bam.file), window.size = window.size, normalize=FALSE, r0=0, r1=1, ymax=300, col="darkorange")

#realment el problema és el chr, el nostre no té chr
gr <-GRanges(seqnames = "4", ranges =  IRanges(start=c(1,1000001), end = c(1000000,1348131)))

windows.1 <-GRanges(seqnames = "1", ranges =  ranges(windows))

######################### Miro com són els bam files
library(Rsamtools)

test <- scanBam(un1.bam.file)
p = ScanBamParam(what=c("rname", "pos"))
table(as.data.frame(scanBam(un1.bam.file, param=p))$rname)
# chr2L   chr2R   chr3L   chr3R    chr4    chrM    chrX chrYHet 
# 0       0       0       0  204355       0       0       0 

test.1 <- BamFile(un1.bam.file)
seqlevels(test.1)

test2 <- scanBam(file=file.path(bamDir,bam.file))
table(as.data.frame(scanBam(file=file.path(bamDir,bam.file), param=p))$rname) #no tenen chr els meus!!
test2.1 <- BamFile(file=file.path(bamDir,bam.file))
seqlevels(test2.1)


dens <- countBam(data, param = ScanBamParam(which = gr))$records


######################## tot genoma 
kp <- plotKaryotype("hg19", plot.type=1, chromosomes=c("chr1", "chr2"))

all.regs <- GRanges()

nreps <- 20
for(i in 1:nreps) {
  regs <- createRandomRegions(nregions = 100, length.mean = 10000000, length.sd = 1000000,
                              non.overlapping = TRUE, genome = "hg19", mask=NA)
  all.regs <- c(all.regs, regs)
  kpPlotRegions(kp, regs, r0 = (i-1)*(0.8/nreps), r1 = (i)*(0.8/nreps), col="#AAAAAA")
}

kpPlotCoverage(kp, all.regs, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF")
kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)


#Example 2: Do the same with a single bigger set of possibly overlapping regions

kp <- plotKaryotype("hg19", plot.type=1)

regs <- createRandomRegions(nregions = 1000, length.mean = 10000000, length.sd = 1000000,
                            non.overlapping = FALSE, genome = "hg19", mask=NA)

kpPlotRegions(kp, regs, r0 = 0, r1 = 0.8, col="#AAAAAA")

kpPlotCoverage(kp, regs, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF", border=NA)
kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)


