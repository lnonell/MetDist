best.aic.ks.with.betabin <- function(best.dist.i){
  best.aic <- names(which.min(round(best.dist.i[1:4],4)))
  best.ks <- names(which.max(round(best.dist.i[5:8],4)))
  if (is.null(best.aic))  best.aic <- NA
  if (is.null(best.ks))  best.ks <- NA
  return(c(best.aic,best.ks))
}