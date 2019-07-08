fn.simulations.betabin <- function(est.params=est.params,cond.n= c(3,5,10,30,100,500),cores=6){
  #en aquest cas faig moments i no mle
  require(VGAM)
  res.list<-NULL
  cpgs.list <-NULL
  #en aquest cas també cal guardar les sizes
  ns.list <-NULL
  
  i=1
  for (n in cond.n){
    # 1 generate simulations
    #betabin
    print(i)
    print(n)
    N.cond1=n
    N.cond2=n
    
    p1.v <- est.params$bb.mom.s1[!is.na(est.params$bb.mom.s1) & est.params$bb.mom.s1>0 ] #n'hi ha 58 negatius
    p2.v <- est.params$bb.mom.s2[!is.na(est.params$bb.mom.s1) & est.params$bb.mom.s1>0 ]
    n.v <- est.params$real.n[!is.na(est.params$bb.mom.s1) & est.params$bb.mom.s1>0 ]
    
    N.cond = N.cond1+N.cond2
   
    N=100
    sample.names = c(paste("SC1",1:N.cond1,sep="_"),paste("SC2",1:N.cond2,sep="_"))
    row.names = paste0("dml",1:N)
    cpgs.dml <- matrix(NA, nrow=N, ncol= (N.cond), dimnames=list(row.names,sample.names))
    ns.dml <- matrix(NA, nrow=N, ncol= (N.cond), dimnames=list(row.names,sample.names))
    for (j in 1:N) {
      k <-sample(1:length(p1.v),1)
      s1 <- p1.v[k]
      s2 <- p2.v[k]
      n.1 <- n.v[k]
      cpgi.1 <- rbetabinom.ab(n=N.cond1,size=n.1, s1,s2)
      #he de generar les n's entre el valor de cpgi.1 i ni.1, que és el maxim
      ni.1 <- unlist(lapply(cpgi.1,function(x) if (x==n.1) x else sample(x:n.1,1))) #aqui si poso ni no compto
      # que de vegades el coverage és molt més gran, hauria de fer servir el coverage?
      #si x=n.1 fa el que li dona la gana i retorna un valor inferior!
      
      k <-sample(1:length(p1.v),1)
      s1 <- p1.v[k]
      s2 <- p2.v[k]
      n.2 <- n.v[k] #potser la n hauria de ser la mateixa??
      cpgi.2 <- rbetabinom.ab(n=N.cond1,size=n.2, s1,s2) 
      ni.2 <- unlist(lapply(cpgi.2,function(x) if (x==n.2) x else sample(x:n.2,1)))
      
      cpgs.dml[j,] <- c(cpgi.1,cpgi.2)
      ns.dml[j,] <- c(ni.1,ni.2) #poso el màxim pq si no podria donar un error però no se si te gaire sentit
      if(any((ns.dml[j,]-cpgs.dml[j,])<0)) print (paste(j,"KK",n.1,n.2,sep="_"))
    }
    
    #ara nondml
    N=1900
    sample.names = c(paste("SC1",1:N.cond1,sep="_"),paste("SC2",1:N.cond2,sep="_"))
    row.names = paste0("l",1:N)
    cpgs.nondml <- matrix(NA, nrow=N, ncol= (N.cond), dimnames=list(row.names,sample.names))
    ns.nondml <- matrix(NA, nrow=N, ncol= (N.cond), dimnames=list(row.names,sample.names))
    for (j in 1:N) {
      k <-sample(1:length(p1.v),1)
      s1 <- p1.v[k]
      s2 <- p2.v[k]
      n.1 <- n.v[k]
      cpgs.nondml[j,] <- rbetabinom.ab(n=N.cond,size=n.1, s1,s2)
      ns.nondml[j,] <- unlist(lapply(cpgs.nondml[j,],function(x) if (x==n.1) x else sample(x:n.1,1)))
      if(any((ns.nondml[j,]-cpgs.nondml[j,])<0)) print (paste(j,"KK",n.1,sep="_"))
    }
    # ajuntar i return
    cpgs <-rbind(cpgs.dml,cpgs.nondml)
    ns.all <- rbind(ns.dml,ns.nondml)
    
    # 2 run models in parallel!!
    cond <- as.factor(c(rep(1,n),rep(2,n))) #sempre 1 i 2, les fns estan preparades per això
    #aquí es on tinc dues opcions, fer el covg tal i com l'he calculat o posant el real
    #si poso el real hauria de carregar les dades reals de coverage i agafar la k que correspongui...
    #complicat
    bb.models <- fn.models.betabin.parallel(all=cpgs,covg=ns.all, cond1=cond, cores=cores)
 
    #guardo les cpgs, els models i les ns
    cpgs <- list(cpgs.bb=cpgs)
    ns <- list(ns.bb=ns.all)
    models.list <- list(cpgs.bb.models=bb.models)
    
    cpgs.list[[i]] <- cpgs
    ns.list[[i]] <- ns
    res.list[[i]] <- models.list
    
    i=i+1
  }
  
  names(cpgs.list) <- cond.n 
  save(cpgs.list,file="simulated.cpgs.bb.list.RData")
  names(ns.list) <- cond.n 
  save(ns.list,file="simulated.ns.bb.list.RData")
  names(res.list) <- cond.n  
  save(res.list,file="simulated.models.bb.list.RData")
  #ho guardo tot i retorno els models
  return(res.list)  
}
