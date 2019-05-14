best.dist.betabin <- function(xi,ni){
  require(fitdistrplus)
  require(VGAM)

  xi.prime <- round(xi[!is.na(xi)],0)
  ni.prime <- ni[!is.na(xi)]
  
  n<- max(xi.prime)
  
  fit.b <- try(fitdist(xi.prime, distr=dbetabinom.ab, start=list(size=n, shape1=1,shape2=1),lower=c(0,0)),TRUE) #aqta és del paquet fitdistrplus
  if(class(fit.b)!="try-error")  b.aic <- fit.b$aic else b.aic=NA
  
  ########## Kolmogorov-Smirnov test
  ks.b <- try(ks.test(b.prime, "pbetabinom.ab", fit.b$estimate[1],fit.b$estimate[2],fit.b$estimate[3]),TRUE)
  if(class(ks.b)!="try-error") b.ks <- ks.b$p.value else b.ks=NA
  
  v <- c(b.aic,b.ks)
  names(v) <- c("betabin.aic","betabin.ks.p")
  return(v)
}

