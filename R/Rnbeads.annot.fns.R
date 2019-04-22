#annot RnBeads

library(RnBeads)
library(RnBeads.hg38) #això conté les anotacions

sum(rnb.annotation.size(assembly="hg38")) #58802720

rnb.region.types(assembly="hg38")
#"tiling"     "genes"      "promoters"  "cpgislands"

rnb.get.annotation(type = "CpG", assembly = "hg38")
rnb.get.annotation(type = "promoters", assembly = "hg38")
rnb.get.annotation(type = "genes", assembly = "hg38")
rnb.get.annotation(type = "tiling", assembly = "hg38")
rnb.get.annotation(type = "cpgislands", assembly = "hg38")
rnb.get.annotation(type = "Open Sea", assembly = "hg38") #error
rnb.get.annotation(type = "shore", assembly = "hg38") #error

annot.genes <- rnb.get.annotation(type = "genes", assembly = "hg38")
#GRangeslist amb 24 GRanges, un per cada chr
annot.genes.gr <- unlist(annot.genes)
annot.genes.gr #60070 ranges

# library(SummarizedExperiment)
# annot.cpg <- rnb.get.annotation(type = "CpG", assembly = "hg38") 
# annot.cpg.1 <- annot.genes[[1]] #hi ha la mateixa info que a genes!!!
# The total number of dinucleotides annotated in HG19 is 28,217,009 represented 
# both on the forward and reverse DNA strands. per aixo dona tot error!!

# annot.cpg.com <- setClass(annot.cpg,"CompressedGRangesList")
# annot.cpg.com <- updateObject(annot.cpg, verbose=TRUE) #rror
# annot.cpg.gr <- unlist(annot.cpg) #error

annot.genes.gr.df <- as.data.frame(annot.genes.gr)
dim(annot.genes.gr.df) #60070    11

#fn to return gene symbol given a chr and position
annot.rnbeads <- function(dat, promoters=F, flank.l=NULL, flank.r=NULL){
  #dat data containing  rownames that I are in format chrN_pos_fn
  #promoters uses promoters function of GRanges to obtain promoters, complementary to the gene info
  #flank.l flanking nt to the left of the annotated genes to check
  #flank.r flanking nt to the left of the annotated genes to check

  aux.fn <- function(ri){ #given a rowname in format chrN_pos_fn it returns the corresponding gene or NA
    annot <- strsplit(ri,split="_")
    chr <-unlist(lapply(annot, function(l) l[[1]]))
    pos <-unlist(lapply(annot, function(l) l[[2]]))
    cpg.gr <- GRanges(seqnames=chr,
                  ranges=IRanges(start=as.numeric(pos),width=1))
    if (promoters) genes <- subsetByOverlaps(promoters(annot.genes.gr),cpg.gr)$symbol[1]
    else genes <- subsetByOverlaps(annot.genes.gr,cpg.gr)$symbol[1] #agafo el primer per si de cas...
    if(length(genes)==0) genes <- NA
    return(genes)
  }
  
  gene.list <- lapply(rownames(dat),aux.fn)
  return(unlist(gene.list))
}


