fn.simulations <- function(est.params=est.params,cond.n= c(3,5,10,30,100,500),cores=6){
  #he de guardar els models, si no dp no els podré estudiar
  res.list<-NULL
  cpgs.list <-NULL
  i=1
  for (n in cond.n){
    # 1 generate simulations
    #simplex
    print(n)
    p1.v <- est.params$s.zoip.mu[which(est.params$s.zoip.sig>0)] #tot i que aquí no cal filtrar
    p2.v <- est.params$s.zoip.sig[which(est.params$s.zoip.sig>0)] #sembla que estigui malament però no!
    cpgs.simplex.dml <- fn.simu.dml(mod="simplex", p1.v=p1.v, p2.v=p2.v, N=100, N.cond1=n, N.cond2=n)
    cpgs.simplex.non.dml <- fn.simu.non.dml(mod="simplex", p1.v=p1.v, p2.v=p2.v, N=1900, N.cond1=n, N.cond2=n)
    cpgs.simplex <-rbind(cpgs.simplex.dml,cpgs.simplex.non.dml) 
    
    #beta
    p1.v <- est.params$b.mle.s1[!is.na(est.params$b.mle.s1)]
    p2.v <- est.params$b.mle.s2[!is.na(est.params$b.mle.s1)]
    cpgs.beta.dml <- fn.simu.dml(mod="beta", p1.v=p1.v, p2.v=p2.v, N=100, N.cond1=n, N.cond2=n)
    cpgs.beta.non.dml <- fn.simu.non.dml(mod="beta", p1.v=p1.v, p2.v=p2.v, N=1900, N.cond1=n, N.cond2=n)
    cpgs.beta <-rbind(cpgs.beta.dml,cpgs.beta.non.dml)
    
    #normal
    p1.v <- est.params$n.mle.m[!is.na(est.params$n.mle.m)]
    p2.v <- est.params$n.mle.sd[!is.na(est.params$n.mle.m)]
    cpgs.normal.dml <- fn.simu.dml(mod="normal", p1.v=p1.v, p2.v=p2.v, N=100, N.cond1=n, N.cond2=n)
    cpgs.normal.non.dml <- fn.simu.non.dml(mod="normal", p1.v=p1.v, p2.v=p2.v, N=1900, N.cond1=n, N.cond2=n)
    cpgs.normal <-rbind(cpgs.normal.dml,cpgs.normal.non.dml)
    
    # 2 run models in parallel!!
    cond <- as.factor(c(rep(1,n),rep(2,n))) #sempre 1 i 2, les fns estan preparades per això
    simplex.models <- fn.models.parallel(cpgs.simplex, cond1=cond, cores=cores)
    beta.models <- fn.models.parallel(cpgs.beta, cond1=cond, cores=cores)
    normal.models <- fn.models.parallel(cpgs.normal, cond1=cond, cores=cores)
    # 3 limma
    cpgs.simplex.limma <- apply.limma(cpgs.simplex,cond) #retorna la toptable ordenada
    cpgs.beta.limma <- apply.limma(cpgs.beta,cond) #retorna la toptable ordenada
    cpgs.normal.limma <- apply.limma(cpgs.normal,cond) #retorna la toptable ordenada
    
    #ho afegeixo a cada model 
    cpgs.simplex.models.withlimma <- data.frame(simplex.models,p.limma=cpgs.simplex.limma$P.Value)
    cpgs.beta.models.withlimma <- data.frame(beta.models,p.limma=cpgs.beta.limma$P.Value)
    cpgs.normal.models.withlimma <- data.frame(normal.models,p.limma=cpgs.normal.limma$P.Value)
    
    #guardo els models
    cpgs <- list(cpgs.simplex=cpgs.simplex, cpgs.beta=cpgs.beta, cpgs.normal=cpgs.normal)
    
    models.list <- list(cpgs.simplex.models=cpgs.simplex.models.withlimma, 
                        cpgs.beta.models=cpgs.beta.models.withlimma, 
                        cpgs.normal.models=cpgs.normal.models.withlimma)
    
    cpgs.list[[i]] <- cpgs
    res.list[[i]] <- models.list
    
    i=i+1
  }
  
  names(cpgs.list) <- cond.n 
  save(cpgs.list,file="simulated.cpgs.list.RData")
  names(res.list) <- cond.n  
  save(res.list,file="simulated.models.list.RData")
  #ho guardo tot i retorno els models
  return(res.list)  
}
