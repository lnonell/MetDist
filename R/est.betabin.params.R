est.betabin.params <- function(xi,ni){
  #b should be a vector with the number of methylated reads by sample, no beta values!
  require(fitdistrplus)
  require(VGAM)
  
  #NAs removal 
  xi.prime <- round(xi[!is.na(xi)],0)
  ni.prime <- ni[!is.na(xi)]
  
  n<- max(xi.prime) #aqui estava xi.prime però hauria de ser ni.prime? sembla que no canvia
  betabin.params <-rep(NA,6)
  betabin.params[6] <- n
  names(betabin.params) <- c("bb.mom.s1","bb.mom.s2","sim.n","bb.mle.s1","bb.mle.s2","real.n")
  
  #moments 
  m1=sum(xi.prime*ni.prime)/sum(ni.prime)
  m2=sum(xi.prime^2*ni.prime)/sum(ni.prime)
  
  s1 <-  (n*m1-m2)/(n*(m2/m1-m1-1)+m1)
  s2 <-  ((n-m1)*(n-(m2/m1)))/(n*((m2/m1)-m1-1)+m1)
  #those are alpha and beta
  
  betabin.params[1] <- s1
  betabin.params[2] <- s2
  
  # if(s1<0 | is.na(s1)) s1<-1
  # if(s2<0 | is.na(s2)) s2<-1
  #fitdistr, using estimated params through moments 
  betabin.params[3:5]<-tryCatch(fitdist(xi.prime, distr=dbetabinom.ab,  
                                        start=list(size=n,shape1=1,shape2=1),lower=c(0,0))$estimate,error=function(err) NA)
  return(betabin.params)
}
