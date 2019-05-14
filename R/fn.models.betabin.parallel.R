fn.models.betabin.parallel <- function(all,covg,cond1,cores=4){
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
    aic.bb <- NA
    df <-data.frame(y,n,cond1)
    df <- df[complete.cases(df),] #only complete cases!
    m.s <- try(betabin(cbind(y, n - y) ~ cond1, ~ cond1,data=df),silent=T)
    if(class(m.s)!="try-error"){
      #if(class(m.s) %in% "glimML"){
      m.s.sum<-summary(m.s)
      p.bb <- m.s.sum@Coef["cond12","Pr(> |z|)"] 
      aic.bb <- AIC(m.s)@istats[1,"AIC"]
    } 
    c(p.bb,aic.bb)
  }
  stopCluster(cl)
  colnames(cpgs.models.bb) <- c("p.bb","aic.bb") 
  rownames(cpgs.models.bb) <- rownames(all)
  return(cpgs.models.bb)
}  
