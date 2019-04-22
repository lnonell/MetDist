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
    y1 <- ifelse(n==y,y-1,y)
    p.bb <- NA
    p.bb.od <- NA
    df <-data.frame(y1,n,cond1)
    df <- df[complete.cases(df),] #only complete cases!
    m.s <- try(betabin(cbind(y1, n - y1) ~ cond1, ~ 1,data=df),silent=T)
    if(class(m.s)!="try-error"){
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
