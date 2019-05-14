##################################################################################
############################### functions to reuse ###############################
##################################################################################
#l'unica fn que te incorporat el parallelisme es fn.model.parallel
##################################################################################
############################### best distribution  ###############################
##################################################################################

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


#to get the best
best.aic.ks <- function(best.dist.i){
  best.aic <- names(which.min(round(best.dist.i[1:3],4)))
  best.ks <- names(which.max(round(best.dist.i[4:6],4)))
  if (is.null(best.aic))  best.aic <- NA
  if (is.null(best.ks))  best.ks <- NA
  return(c(best.aic,best.ks))
}


best.dist.betabin <- function(xi,ni){
  #en ppi ni no ho necessito
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

best.aic.ks.with.betabin <- function(best.dist.i){
  best.aic <- names(which.min(round(best.dist.i[1:4],4)))
  best.ks <- names(which.max(round(best.dist.i[5:8],4)))
  if (is.null(best.aic))  best.aic <- NA
  if (is.null(best.ks))  best.ks <- NA
  return(c(best.aic,best.ks))
}

##################################################################################
######################### estimation of dist parameters ##########################
##################################################################################

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


#DEF betabinomial, only for sequencing: ho canvio a fitdistrplus
est.betabin.params <- function(xi,ni){
  #b should be a vector with the number of methylated reads by sample, no beta values!
  require(fitdistrplus)
  require(VGAM)
  
  #NAs removal 
  xi.prime <- round(xi[!is.na(xi)],0)
  ni.prime <- ni[!is.na(xi)]
  
  n<- max(xi.prime)
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

##################################################################################
############################# apply assessed models ##############################
##################################################################################

fn.models.parallel <- function(cpgs,cond1,cores=4){
  #should we obtain regression coefs too?
  #definitivament elimino gamlss
  
  library(doParallel)
  library(foreach)
  
  require(simplexreg)
  require(betareg)
  require(ZOIP) #simplex inflated
  require(gamlss) #beta inflated
  require(quantreg)
  
  N = nrow(cpgs)
  
  cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
  registerDoParallel(cl)
  cpgs.models<- foreach(i=1:N,.combine=rbind, .packages=c("simplexreg","betareg","gamlss","ZOIP","quantreg")) %dopar% {
    cpgi <- cpgs[i,]
    
    #remove NAs for all #we should also remove values out of [0,1]...this affects to the normal simulated data
    e<-0.00001
    cond <- cond1[!is.na(cpgi) & cpgi>e & cpgi<1-e]
    cpgi <-cpgi[!is.na(cpgi) & cpgi>e & cpgi<1-e]
    
    p.s <- NA
    p.b <- NA
    p.sinf <- NA
    p.binf <- NA
    p.n <- NA
    p.l <- NA
    p.q <- NA
    # p.s.gamlss <- NA
    
    #NON INFLATED models
    #simplex regression with simplexreg
    m.s <- try(simplexreg(cpgi ~ cond),silent=T)
    if(class(m.s)!="try-error"){
      m.s.sum<-summary(m.s)
      p.s <- m.s.sum$coefficients$mean["cond2","Pr(>|z|)"] 
    }
    
    # beta regression with betareg
    m.b <-  try(betareg(cpgi ~ cond),silent=T)
    if(class(m.b)!="try-error"){
      m.b.sum<-summary(m.b)
      p.b <- m.b.sum$coefficients$mean["cond2","Pr(>|z|)"] 
    }
    
    #simplex inflated with ZOIP package, need tot transform into a df
    df <- data.frame(cond,cpgi)
    m.sinf <- try(RM.ZOIP(formula.mu = cpgi~cond,  
                          link = c("logit","identity","identity","identity"), family="Simplex", data=df),silent=T)
    if(class(m.sinf)!="try-error"){
      estimate <- m.sinf$par #from ZOIP package, to extract p from mu param
      se       <- try(sqrt(diag(solve(m.sinf$HM))),silent=T)
      if(class(se)!="try-error"){
        zvalue   <- estimate / se
        pvalue   <- 2 * stats::pnorm(abs(zvalue), lower.tail=F)
        p.sinf<-pvalue["cond2"]
      }
    }  
    
    
    #beta inflated with gamlss package
    #sink("out.txt") #no tinc ous que no surti output!!!
    m.binf <- try(gamlss(cpgi ~ cond, family = BEINF),silent=T) #so suppress gamlss messages, que és molt pesat!
    #sink(NULL)
    if(class(m.binf)!="try-error"){
      m.binf.sum<-tryCatch(summary(m.binf, save=TRUE), error=function(err) NA)
      p.prov <- m.binf.sum$pvalue["cond2"]
      if(!is.null(p.prov)) p.binf<- p.prov
    }
    
    #Normal regression
    m.n <-  try(lm(cpgi ~ cond),silent=T)
    if(class(m.n)!="try-error"){
      m.n.sum<-summary(m.n)
      print(m.n.sum)
      p.n <- tryCatch(m.n.sum$coefficients["cond2","Pr(>|t|)"] , error=function(err) NA)
    }
    # #beta-binomial?
    
    # 
    #logistic, hem de canviar les condicions, logística vol 0s i 1s
    cond2 <- as.factor(as.numeric(cond)-1)
    m.l <-  try(glm(cond2 ~ cpgi, family=binomial),silent=T)
    if(class(m.l)!="try-error"){
      m.l.sum<-summary(m.l)
      p.l <-tryCatch( m.l.sum$coefficients["cpgi","Pr(>|z|)"], error=function(err) NA)
      #no sé pq en alguns casos només hi ha intercept!
    }
    
    #quantile regression with rq, 75%
    m.q <-  try(rq(cpgi ~ cond, tau=.75),silent=T)
    if(class(m.q)!="try-error"){
      m.q.sum<-tryCatch(summary(m.q, se="ker")[[3]], error=function(err) NA)
      p.q <- tryCatch(m.q.sum["cond2","Pr(>|t|)"] , error=function(err) NA)
    }
    
    # m.sgamlss <- try(gamlss(cpgi ~ cond, family = SIMPLEX),silent=T) #so suppress gamlss messages, que és molt pesat!
    # #sink(NULL)
    # if(class(m.sgamlss)!="try-error"){
    #   m.sgamlss.sum<-tryCatch(summary(m.sgamlss, save=TRUE), error=function(err) NA)
    #   p.prov <- m.sgamlss.sum$pvalue["cond2"]
    #   if(!is.null(p.prov)) p.s.gamlss<- p.prov
    # }
    
    # c(p.s,p.b,p.sinf,p.binf,p.n,p.l,p.q,p.s.gamlss)
    c(p.s,p.b,p.sinf,p.binf,p.n,p.l,p.q)
  }
  stopCluster(cl)
  #colnames(cpgs.models) <- c("p.s","p.b","p.sinf","p.binf","p.n","p.l","p.q","p.s.gamlss") #potser hauríem de retornar el num de NAs
  colnames(cpgs.models) <- c("p.s","p.b","p.sinf","p.binf","p.n","p.l","p.q") #potser hauríem de retornar el num de NAs
  rownames(cpgs.models) <- rownames(cpgs)
  return(cpgs.models)
}  


apply.limma <- function(cpgs, cond){
  library(limma)
  
  #first transform to M's
  cpgs.M <- log2(cpgs/(1-cpgs))
  
  #agafo els noms per a poder ordenar dp
  nams <- rownames(cpgs)
  
  design<-model.matrix(~0+cond)
  fit<-lmFit(cpgs.M,design) #model could also be fitted to the batch corrected data
  contrast.matrix<-makeContrasts(cond2-cond1,levels=design)
  fit2<-contrasts.fit(fit,contrast.matrix)
  fite<-eBayes(fit2)
  top.table<-topTable(fite,coef=1,number=Inf,adjust="BH")
  top.table.s <- top.table[nams,] #ordenat com al principi!
  return(top.table.s)
}

#proves beta-binomial
fn.betabin.od.parallel <- function(all,covg,cond1,cores=4){
  #test the overdispersion param, phi
  #https://www.r-bloggers.com/binary-beta-beta-binomial/
  library(aod)
  library(doParallel)
  library(foreach)
  
  N = nrow(all)
  
  cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
  registerDoParallel(cl)
  cpgs.models.bb<- foreach(i=1:N,.combine=rbind, .packages=c("aod")) %dopar% {
    y <- all[i,]
    n <- covg[i,]
    #y1 <- ifelse(n==y,y-1,y)
    p.bb <- NA
    phi.bb <- NA
    phi.p.bb <- NA
    df <-data.frame(y,n,cond1)
    df <- df[complete.cases(df),] #only complete cases!
    m.s <- try(betabin(cbind(y, n - y) ~ cond1, ~ 1,data=df),silent=T)
    if(class(m.s)!="try-error"){
      #if(class(m.s) %in% "glimML"){
      m.s.sum<-summary(m.s)
      p.bb <- m.s.sum@Coef["cond12","Pr(> |z|)"] 
      phi.bb <- m.s.sum@Phi["phi.(Intercept)","Estimate"]
      phi.p.bb <- m.s.sum@Phi["phi.(Intercept)","Pr(> z)"]
    } 
    c(p.bb,phi.bb,phi.p.bb)
  }
  stopCluster(cl)
  colnames(cpgs.models.bb) <- c("p.bb","phi.bb","phi.p.bb") 
  rownames(cpgs.models.bb) <- rownames(all)
  return(cpgs.models.bb)
}  

#aquestes dues funcions són per provar els params d'overdispersion)
fn.betabin.od.cond.parallel <- function(all,covg,cond1,cores=4){
  #test the overdispersion param, assumeixo que hi ha overdisp diferent a cada grup (p de phi <0.05)
  #https://www.r-bloggers.com/binary-beta-beta-binomial/
  library(aod)
  library(doParallel)
  library(foreach)
  
  N = nrow(all)
  
  cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
  registerDoParallel(cl)
  cpgs.models.bb<- foreach(i=1:N,.combine=rbind, .packages=c("aod")) %dopar% {
    y <- all[i,]
    n <- covg[i,]
    #y1 <- ifelse(n==y,y-1,y)
    p.bb <- NA
    df <-data.frame(y,n,cond1)
    df <- df[complete.cases(df),] #only complete cases!
    m.s <- try(betabin(cbind(y, n - y) ~ cond1, ~ cond1,data=df),silent=T)
    if(class(m.s)!="try-error"){
      #if(class(m.s) %in% "glimML"){
      m.s.sum<-summary(m.s)
      p.bb <- m.s.sum@Coef["cond12","Pr(> |z|)"] 
    } 
    c(p.bb)
  }
  stopCluster(cl)
  colnames(cpgs.models.bb) <- c("p.bb") 
  rownames(cpgs.models.bb) <- rownames(all)
  return(cpgs.models.bb)
}  

fn.betabin.od.cond.nood.parallel <- function(all,covg,cond1,cores=4){
  #test the overdispersion param, assumeixo que hi ha overdisp diferent a cada grup (p de phi <0.05)
  #https://www.r-bloggers.com/binary-beta-beta-binomial/
  library(aod)
  library(doParallel)
  library(foreach)
  
  N = nrow(all)
  
  cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
  registerDoParallel(cl)
  cpgs.models.bb<- foreach(i=1:N,.combine=rbind, .packages=c("aod")) %dopar% {
    y <- all[i,]
    n <- covg[i,]
    #y1 <- ifelse(n==y,y-1,y)
    p.bb <- NA
    df <-data.frame(y,n,cond1)
    df <- df[complete.cases(df),] #only complete cases!
    m.s <- try(betabin(cbind(y, n - y) ~ cond1, ~ cond1,data=df,fixpar = list(c(3, 4), c(0, 0))),silent=T)
    if(class(m.s)!="try-error"){
      #if(class(m.s) %in% "glimML"){
      m.s.sum<-summary(m.s)
      p.bb <- m.s.sum@Coef["cond12","Pr(> |z|)"] 
    } 
    c(p.bb)
  }
  stopCluster(cl)
  colnames(cpgs.models.bb) <- c("p.bb") 
  rownames(cpgs.models.bb) <- rownames(all)
  return(cpgs.models.bb)
}  

fn.models.betabin.parallel <- function(all,covg,cond1,cores=4){
  #in this case I also get the overdispersion param, phi
  #https://www.r-bloggers.com/binary-beta-beta-binomial/
  library(aod)
  library(doParallel)
  library(foreach)
  
  N = nrow(all)
  
  cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
  registerDoParallel(cl)
  cpgs.models.bb<- foreach(i=1:N,.combine=rbind, .packages=c("aod")) %dopar% {
    y <- all[i,]
    n <- covg[i,]
    #y1 <- ifelse(n==y,y-1,y)
    p.bb <- NA
    p.bb.od <- NA
    df <-data.frame(y,n,cond1)
    df <- df[complete.cases(df),] #only complete cases!
    m.s <- try(betabin(cbind(y, n - y) ~ cond1, ~ 1,data=df),silent=T)
    if(class(m.s)!="try-error"){
      #if(class(m.s) %in% "glimML"){
      m.s.sum<-summary(m.s)
      p.bb <- m.s.sum@Coef["cond12","Pr(> |z|)"] 
      p.bb.od <- m.s.sum@Phi["phi.(Intercept)","Pr(> z)"]
    } 
    c(p.bb,p.bb.od)
  }
  stopCluster(cl)
  colnames(cpgs.models.bb) <- c("p.bb","p.bb.od") 
  rownames(cpgs.models.bb) <- rownames(all)
  return(cpgs.models.bb)
}  

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

##################################################################################
############################### functions for simulations ########################
##################################################################################
#it also calls fn models parallel

#function to obtain simulations, for dml
fn.simu.dml <- function(mod="beta", p1.v, p2.v, N=1000, N.cond1=100, N.cond2=100){
  library(VGAM) #abans feia servir betareg
  #mod= simplex, beta, normal
  #p1.v = vector of param1 of distribution to sample from
  #p2.v = vector of param2 of distribution to sample from
  #N number of probes (CpG sites) to generate
  #N.cond1 = number of samples belonging to condition1
  #N.cond2 = number of samples belonging to condition2
  #falta opció per a que siguin iguals i per tant N.cond2 agafi de la mateixa dist (amb igual params)

  #STUDY TO GENERATE P1.V & P2.V HERE INSTEAD OF PASSING AS PARAMETER
  N.cond = N.cond1+N.cond2
  sample.names = c(paste("SC1",1:N.cond1,sep="_"),paste("CS2",1:N.cond2,sep="_"))
  row.names = paste0("dml",1:N)
  cpgs <- matrix(NA, nrow=N, ncol= (N.cond), dimnames=list(row.names,sample.names))

    if (mod=="simplex"){

    for (i in 1:N) {
      #sample in the numbers an the take corresponding mu1 and sig1
      j <-sample(1:length(p1.v),1)
      mu1 <- p1.v[j]
      sig1 <- p2.v[j]
      cpgi.1 <- rsimplex(N.cond1,mu1,sig1)

      k <-sample(1:length(p1.v),1)
      mu2 <- p1.v[k]
      sig2 <- p2.v[k]
      cpgi.2 <- rsimplex(N.cond2,mu2,sig2)

      cpgs[i,] <- c(cpgi.1,cpgi.2)
    }

  } else if (mod=="beta"){

    for (i in 1:N) {
      j <-sample(1:length(p1.v),1)
      s1 <- p1.v[j]
      s2 <- p2.v[j]
      cpgi.1 <- rbeta(N.cond1,s1,s2)

      k <-sample(1:length(p1.v),1)
      s1 <- p1.v[k]
      s2 <- p2.v[k]
      cpgi.2 <- rbeta(N.cond2,s1,s2)

      cpgs[i,] <- c(cpgi.1,cpgi.2)
    }

  } else if (mod=="normal"){

    for (i in 1:N) {
      j <-sample(1:length(p1.v),1)
      m1 <- p1.v[j]
      s1 <- p2.v[j]
      cpgi.1 <- rnorm(N.cond1,m1,s1)

      k <-sample(1:length(p1.v),1)
      m2 <- p1.v[k]
      s2 <- p2.v[k]
      cpgi.2 <- rnorm(N.cond2,m2,s2)

      cpgs[i,] <- c(cpgi.1,cpgi.2)
    }

  }
    return(cpgs)
}

#now for non.dml, from the same dist
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


#fn per a fer simulacions amb diferents n
fn.simulations <- function(est.params=est.params,cond.n= c(3,5,10,30,100,500),cores=6){
  #he de guardar els models, si no dp no els podré estudiar
  res.list<-NULL
  cpgs.list <-NULL
  i=1
  for (n in cond.n){
    # 1 generate simulations
    #simplex
    print(n)
    #aqui agafava zoip pq era el mateix en les proves a GSE50660paramest.R i tenia menys NAs
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

##################################################################################
############################### altres fns auxiliars #############################
##################################################################################

find.dml.n <- function(x,alpha=0.05) length(x[!is.na(x) & x<alpha])


#per si la poso: https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient
sorensen <- function(TP,FP,FN) {Sor <- 2*TP/(2*TP+FP+FN); return(Sor)}

#function to estimate the best distribution that fits the data (non inflated!!): simplex, beta or normal
#per estimar beta-binomial, hauria de tenir el coverage també

#per als resultats de les comparacions
res.chr <-function(res=rnb.meth2comp.models.adj.p.s){
  annot <- strsplit(rownames(res),split="_")
  chr <-unlist(lapply(annot, function(l) l[[1]]))
  return(as.data.frame(table(chr)))
}

##################################################################################
############################### VELLES O NO USADES ##############################
##################################################################################
#no faig les simulacions amb la beta-binomial al final
# fn.simulations.betabinomial <- function(est.params=est.params.betabin,cond.n= c(3,5,10,30,100,500),cores=6){
#   
#   res.list<-NULL
#   cpgs.list <-NULL
#   i=1
#   for (n in cond.n){
#     # 1 generate simulations
#     #simplex
#     print(n)
#     #beta bin, en aquest cas agafo mom pq mle surt molt raro
#     n <- est
#     p1.v <- est.params$bb.mle.s1[!is.na(est.params$bb.mle.s1)]
#     p2.v <- est.params$bb.mle.s2[!is.na(est.params$bb.mle.s1)]
#     cpgs.betabin <- fn.simu.betabin(mod="betabin",n=n, p1.v=p1.v, p2.v=p2.v, N=100, N.cond1=n, N.cond2=n)
#     
#     # 2 run models in parallel!!
#     cond <- as.factor(c(rep(1,n),rep(2,n))) #sempre 1 i 2, les fns estan preparades per això
#     cpgs.betabin.models <- fn.models.betabin.parallel( cpgs.betabin, cond1=cond, cores=cores)
#   
#   #guardo els models
#     cpgs <- list(cpgs.betabin=cpgs.betabin)
#     
#     models.list <- list(cpgs.simplex.models=cpgs.simplex.models.withlimma, 
#                         cpgs.beta.models=cpgs.beta.models.withlimma, 
#                         cpgs.normal.models=cpgs.normal.models.withlimma)
#     
#     cpgs.list[[i]] <- cpgs
#     res.list[[i]] <- models.list
#     
#     i=i+1
#   }
#   
#   names(cpgs.list) <- cond.n 
#   save(cpgs.list,"simulated.cpgs.list.RData")
#   names(res.list) <- cond.n  
#   save(res.list,"simulated.models.list.RData")
#   #ho guardo tot i retorno els models
#   return(res.list)  
# }

# best.dist.betabin.OLD <- function(b){
#   #b should be a vector with the number of methylated reads by sample, no beta values!
#   #malament: N should be the total number of measures not the max!
#   require(fitdistrplus)
#   require(VGAM)
#   
#   #NAs removal and round to integer
#   b.prime <- round(b[!is.na(b)],0)
#   n <- max(b.prime)
#   betabin.params <-rep(NA,4)
#   
#   mu <- mean(b.prime)
#   var <- var(b.prime) 
#   
#   #from wikipedia!
#   s1 <-  (n*mu-var)/(n*((var/mu)-mu-1)+mu)
#   s2 <-  ((n-mu)*(n-(var/mu)))/(n*((var/mu)-mu-1)+mu)
#   
#   fit.b <- try(fitdist(b.prime, distr=dbetabinom.ab, start=list(size=n, shape1=s1,shape2=s2)),TRUE) #aqta és del paquet fitdistrplus
#   if(class(fit.b)!="try-error")  b.aic <- fit.b$aic else b.aic=NA
# 
#     ########## Kolmogorov-Smirnov test
#   ks.b <- try(ks.test(b.prime, "pbetabinom.ab", fit.b$estimate[1],fit.b$estimate[2],fit.b$estimate[3]),TRUE)
#   if(class(ks.b)!="try-error") b.ks <- ks.b$p.value else b.ks=NA
#   
#   v <- c(b.aic,b.ks)
#   names(v) <- c("betabin.aic","betabin.ks.p")
#   return(v)
# }
