est.all.params <- function(b){
  #mirar si afegir aquí binf de gamlss
  require(fitdistrplus)
  require(VGAM)
  require(ZOIP)
  
  #NAs removal
  b.1 <-b[!is.na(b)]
  e<-0.01 #molt important per l'ajust!! quan e=0.001 dona unes estimacions de les cues brutals per la simplex 
  b.prime <- ifelse(b.1==1, 1-e, ifelse(b.1==0,0+e,b.1))
  
  mu=mean(b.prime)
  var=var(b.prime)
  
  params <-rep(NA,14)
  names(params) <- c("s.mle.mu","s.mle.sig","s.zoip.mu","s.zoip.sig","b.mom.s1","b.mom.s2",
                     "b.mle.s1","b.mle.s2","binf.mle.s1","binf.mle.s2","n.mom.m","n.mom.sd","n.mle.m","n.mle.sd")
  
  #simplex: MLE
  
  params[1:2]<-tryCatch(fitdist(b.prime, distr=dsimplex, start=list(mu=0.5,dispersion=1), optim.method="Nelder-Mead")$estimate,error=function(err) NA)
  
  #simplex: ZOIP
  betai.df <- data.frame(b.prime)
  mod <- try(RM.ZOIP(formula.mu = b.prime ~ 1, data = betai.df, family = "Simplex"),silent=T)
  if(class(mod)!="try-error"){
    params[3] <- coef(mod)$`Parameters.mu`
    params[4] <- sqrt(coef(mod)$Parameters.sigma) ##ho he comprovat manualment, generant random dist
  }
  
  #beta: moments
  params[5] <-  ((1 - mu) / var - 1 / mu) * mu ^ 2
  params[6] <- params[5] * (1 / mu - 1)
  #beta: MLE
  params[7:8]<-tryCatch(fitdist(b.prime, distr=dbeta, start=list(shape1=mu,shape2=var))$estimate,error=function(err) NA)
  
  #now dist that can be inflated
  #beta inflated: recently added :-)
  mu=mean(b.1)
  var=var(b.1)
  
  f.1 <- sum(b.1 == 1)/length(b.1)
  f.0 <- sum(b.1 == 0)/length(b.1)
  
  params[9:10]<-tryCatch(fitdist(b.1, distr=dzoabeta,  
                                 start=list(shape1=mu,shape2=var,pobs0=f.0,pobs1=f.1))$estimate,error=function(err) NA)
  #normal: direct moments
  params[11] = mu
  params[12] = sqrt(var)
  
  #normal: MLE
  params[13:14]<-tryCatch(fitdist(b.1, distr=dnorm)$estimate,error=function(err) NA)
  
  return(params)
  
}
