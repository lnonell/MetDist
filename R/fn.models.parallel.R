fn.models.parallel <- function(cpgs,cond1,cores=4){
  #should we obtain regression coefs too?
  
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
    p.s.gamlss <- NA
    
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
    
    m.sgamlss <- try(gamlss(cpgi ~ cond, family = SIMPLEX),silent=T) #so suppress gamlss messages, que és molt pesat!
    #sink(NULL)
    if(class(m.sgamlss)!="try-error"){
      m.sgamlss.sum<-tryCatch(summary(m.sgamlss, save=TRUE), error=function(err) NA)
      p.prov <- m.sgamlss.sum$pvalue["cond2"]
      if(!is.null(p.prov)) p.s.gamlss<- p.prov
    }
    
    c(p.s,p.b,p.sinf,p.binf,p.n,p.l,p.q,p.s.gamlss)
  }
  stopCluster(cl)
  colnames(cpgs.models) <- c("p.s","p.b","p.sinf","p.binf","p.n","p.l","p.q","p.s.gamlss") #potser hauríem de retornar el num de NAs
  rownames(cpgs.models) <- rownames(cpgs)
  return(cpgs.models)
}  