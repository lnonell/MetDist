best.dist <- function(x){
  library(fitdistrplus)
  library(simplexreg) 
  mu <- mean(x,na.rm=T)
  var <- var(x,na.rm=T)
  
  fit.s <- try(fitdist(x, distr=dsimplex, start=list(mu=0.5,sig=2),optim.method="Nelder-Mead"),TRUE) #no entenc pq falla amb les dades guardades!!
  if(class(fit.s)!="try-error")  s.aic <- fit.s$aic else s.aic=NA
  
  fit.b <- try(fitdist(x, distr=dbeta, start=list(shape1=mu,shape2=var)),TRUE) #aqta és del paquet fitdistrplus
  if(class(fit.b)!="try-error")  b.aic <- fit.b$aic else b.aic=NA
  
  fit.n <- try(fitdist(x, distr="norm"),TRUE)
  if(class(fit.n)!="try-error")  n.aic <- fit.n$aic else n.aic=NA
  
  ########## Kolmogorov-Smirnov test
  ks.s <- try(ks.test(x, "psimplex", fit.s$estimate[1],fit.s$estimate[2]),TRUE)
  if(class(ks.s)!="try-error") s.ks <- ks.s$p.value else s.ks=NA
  
  ks.b <- try(ks.test(x, "pbeta", fit.b$estimate[1],fit.b$estimate[2]),TRUE)
  if(class(ks.b)!="try-error") b.ks <- ks.b$p.value else b.ks=NA
  
  ks.n <- try(ks.test(x, pnorm, fit.n$estimate[1],fit.n$estimate[2]),TRUE)
  if(class(ks.n)!="try-error") n.ks <- ks.n$p.value else n.ks=NA
  
  v <- c(s.aic,b.aic,n.aic,s.ks,b.ks,n.ks)
  names(v) <- c("simplex.aic","beta.aic","normal.aic","simplex.ks.p","beta.ks.p","normal.ks.p")
  return(v)
}
