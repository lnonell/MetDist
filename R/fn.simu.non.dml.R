fn.simu.non.dml <- function(mod="beta",p1.v, p2.v, N=1000, N.cond1=100, N.cond2=100){
  library(VGAM) #abans simplexreg
  #mod= simplex, beta, normal
  #p1.v = vector of param1 of distribution to sample from
  #p2.v = vector of param2 of distribution to sample from 
  #N number of probes (CpG sites) to generate
  #N.cond1 = number of samples belonging to condition1  
  #N.cond2 = number of samples belonging to condition2
  #falta opció per a que siguin iguals i per tant N.cond2 agafi de la mateixa dist (amb igual params)
  
  #STUDY TO GENERATE P1.V & P2.V HERE INSTEAD OF PASSING AS PARAMETER
  
  #in this case 
  N.cond = N.cond1+N.cond2
  sample.names = c(paste("SC1",1:N.cond1,sep="_"),paste("CS2",1:N.cond2,sep="_"))
  row.names = paste0("l",1:N)
  cpgs <- matrix(NA, nrow=N, ncol= (N.cond), dimnames=list(row.names,sample.names))
  
  if (mod=="simplex"){
    
    for (i in 1:N) {
      #sample in the numbers an the take corresponding mu1 and sig1
      j <-sample(1:length(p1.v),1)
      mu1 <- p1.v[j]
      sig1 <- p2.v[j]
      cpgs[i,]<- rsimplex(N.cond,mu1,sig1)
    }
    
  } else if (mod=="beta"){
    
    for (i in 1:N) {
      j <-sample(1:length(p1.v),1)
      s1 <- p1.v[j]
      s2 <- p2.v[j]
      cpgs[i,] <- rbeta(N.cond,s1,s2)
      
    }
    
  } else if (mod=="normal"){
    
    for (i in 1:N) {
      j <-sample(1:length(p1.v),1)
      m1 <- p1.v[j]
      s1 <- p2.v[j]
      cpgs[i,] <- rnorm(N.cond,m1,s1)
    }
    
  }  
  return(cpgs)
  
}