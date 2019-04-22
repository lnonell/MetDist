best.aic.ks <- function(best.dist.i){
  best.aic <- names(which.min(round(best.dist.i[1:3],4)))
  best.ks <- names(which.max(round(best.dist.i[4:6],4)))
  if (is.null(best.aic))  best.aic <- NA
  if (is.null(best.ks))  best.ks <- NA
  return(c(best.aic,best.ks))
}