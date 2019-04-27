#header for package
apply.limma <- function(cpgs, cond){
  #first transform to M's
  cpgs.M <- log2(cpgs/(1-cpgs))
  
  #agafo els noms per a poder ordenar dp
  nams <- rownames(cpgs)
  library(limma)
  design<-model.matrix(~0+cond)
  fit<-lmFit(cpgs.M,design) #model could also be fitted to the batch corrected data
  contrast.matrix<-makeContrasts(cond2-cond1,levels=design)
  fit2<-contrasts.fit(fit,contrast.matrix)
  fite<-eBayes(fit2)
  top.table<-topTable(fite,coef=1,number=Inf,adjust="BH")
  top.table.s <- top.table[nams,] #ordenat com al principi!
  return(top.table.s)
}
