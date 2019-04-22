models.eval <- function(models,dml.r=100,alpha=0.05, adjust=TRUE){
  #for simulated data
  #models: matrix with each column a model tested
  #dml.r: number of rows dml (always at the beginning), the rest should be non dml
  #alpha: signif
  
  #adjust
  if (adjust)  models <- apply(models,2,p.adjust)
  
  #for each column
  models.dml <- models[1:dml.r,]
  models.nondml <- models[(dml.r+1):nrow(models),]
  
  nas <- apply(models,2,function(x) sum(is.na(x)))
  
  TP <- apply(models.dml,2,find.dml.n,alpha=0.05)
  FP <- apply(models.nondml,2,find.dml.n,alpha=0.05)
  TN <- 1900-FP
  FN <- dml.r-TP
  
  sens <- TP/(TP+FN)
  spec <- TN/(TN+FP)
  
  #find.dml <- function(x) {y <- x[!is.na(x) & x<0.05]; length(intersect(y,))
  
  jaccard <-TP/(TP+FP+FN)
  
  eval.meas <- rbind(nas, sens, spec, TP, TN, FP, FN,jaccard)  
  return(eval.meas)
  
}
