#Figures and tables from results of simulations
#Figure4: jaccard of all data sets for non inflated models N=3-500
#Supp figure3: jaccard of all data sets for non inflated models N=3,5,10 and 30
#Table 2: Summary of performance measures for simulations N=100
#Supplementary tables 3-7:  Summary of performance measures for simulations N=3,5,10,30,500
#11/05/19 ho actualitzo tot
#5/7/19 elimino logística

#header amb tot
source("D:/Doctorat/Simplex/MetDist/R/Figures.HeaderScript.R")


########################################################################################
############### Figure 3: Resultats de les simulacions: TOTES #######################
########################################################################################

#Fn per plotar les corbes
#la fn sense les inflated està abaix!!!
#per a fer nomes de 3 a 30 descomentar línia scale_x.... dins la fn  fn
#fn sense les inflated (i comptant que als evals hi ha gamlss per simplex), es diu igual!!!
# plot.eval.curves <- function(evals,measure=jaccard,ylim_0=0,ylim_1=1){
#   measure <- enquo(measure)
#   
#   eval.all <- cbind (do.call("rbind",evals),N=as.numeric(rep(names(res.list),each=9)))
#   eval.all<- data.frame(eval.all,model.p=rownames(eval.all),stringsAsFactors = F)
#   
#   #aquí elimino les inflated 
#   eval.all <- eval.all[!(eval.all$model.p %in% c("p.binf","p.sinf","p.s.gamlss")),]
#   
#   eval.all$model.p <- factor(recode(eval.all$model.p,  p.b="beta", p.s="simplex",
#                                     p.n="normal", p.l="logistic", 
#                                     p.q="quantile",p.limma="limma"),
#                              levels=c("beta", "simplex","normal","logistic","quantile","limma"))
#   
#   p <- ggplot(data=eval.all,aes(x=N, y=!!measure, colour=model.p)) +
#     geom_point(size = 1)+ geom_line(size=0.6) +ylim(ylim_0,ylim_1)
#   p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
#              axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
#              axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
#    scale_x_continuous(breaks=as.numeric(names(res.list)))
#  # scale_x_continuous(breaks=c(3,5,10,30),limits=c(0,30)) #la poso i la trec segons necessitat
# }    
# plot.eval.curves <- function(evals,measure=jaccard,ylim_0=0,ylim_1=1){
#   measure <- enquo(measure)
# 
#   eval.all <- cbind (do.call("rbind",evals),N=as.numeric(rep(names(res.list),each=8)))
#   eval.all<- data.frame(eval.all,model.p=rownames(eval.all))
#   eval.all$model.p <- factor(recode(eval.all$model.p,  p.b="beta", p.binf="beta.inf",p.s="simplex",
#                                     p.sinf="simplex.inf", p.n="normal", p.l="logistic",
#                                     p.q="quantile",p.limma="limma"),
#                              levels=c("beta", "beta.inf","simplex","simplex.inf","normal","logistic","quantile","limma"))
# 
#   p <- ggplot(data=eval.all,aes(x=N, y=!!measure, colour=model.p)) +
#     geom_point(size = 1)+ geom_line(size=0.6) +ylim(ylim_0,ylim_1)
#   p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
#              axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
#              axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
#   scale_x_continuous(breaks=as.numeric(names(res.list)))
#  # scale_x_continuous(breaks=c(3,5,10,30),limits=c(0,30)) #la poso i la trec segons necessitat
# }

#proves amb gamlss
# plot.eval.curves <- function(evals,measure=jaccard,ylim_0=0,ylim_1=1){
#   measure <- enquo(measure)
#   
#   eval.all <- cbind (do.call("rbind",evals),N=as.numeric(rep(names(res.list),each=9))) #each 9
#   eval.all<- data.frame(eval.all,model.p=rownames(eval.all),stringsAsFactors = F)
#   
#   #aquí elimino les inflated 
#  # eval.all <- eval.all[!(eval.all$model.p %in% c("p.binf","p.sinf")),]
#   
#   eval.all$model.p <- factor(recode(eval.all$model.p,  p.b="beta", p.binf="beta.inf",p.s="simplex",
#                                     p.sinf="simplex.inf", p.n="normal", p.l="logistic",
#                                     p.q="quantile",p.limma="limma",p.s.gamlss="sim.gamlss"),
#                              levels=c("beta", "beta.inf","simplex","simplex.inf","normal",
#                                       "logistic","quantile","limma","sim.gamlss"))
#   
#   
#   # eval.all$model.p <- factor(recode(eval.all$model.p,  p.b="beta", p.s="simplex",
#   #                                   p.n="normal", p.l="logistic", 
#   #                                   p.q="quantile",p.limma="limma",p.s.gamlss="sim.gamlss"),
#   #                            levels=c("beta", "simplex","normal","logistic","quantile","limma", "sim.gamlss"))
#   
#   p <- ggplot(data=eval.all,aes(x=N, y=!!measure, colour=model.p)) +
#     geom_point(size = 1)+ geom_line(size=0.6) +ylim(ylim_0,ylim_1)
#   p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
#              axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
#              axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
#     scale_x_continuous(breaks=as.numeric(names(res.list)))
#   #scale_x_continuous(breaks=c(3,5,10,30),limits=c(0,30)) #la poso i la trec segons necessitat
# }    

#load(file=file.path(GSE50660_data,"simulated.models.list.RData"))

#fn to eval and flatten data sets
eval_flatten <- function(data){
  is <- c(3,5,10,30,100,500)
  js <- c("simplex","beta","normal")
  evals <- NULL
  for (i in 1:6){ #N de cada simulació
    for (j in 1:3){
      mod <- t(models.eval(models=res.list[[i]][[j]], dml.r=100,alpha=0.05, adjust=TRUE))
      evals <- rbind(evals,
                     data.frame(mod, model.p=rownames(mod),N=is[i],dist=js[j],stringsAsFactors = F))
    }
  }
  return(evals)
}  

load(file=file.path(GSE50660_data,"simulated.models.list.RData"))
gse50.ef <- eval_flatten(res.list)
load(file=file.path(GSE116339_data,"simulated.models.list.RData"))
gse11.ef <- eval_flatten(res.list)
load(file=file.path(RRBS188_data,"simulated.models.list.RData"))
rrbs.ef <- eval_flatten(res.list)
load(file=file.path(WGBS81_data,"simulated.models.list.RData"))
wgbs.ef <- eval_flatten(res.list)



eval.all <- rbind(data.frame(gse50.ef,dataset="a450k-smoking"),
                   data.frame(gse11.ef,dataset="aEPIC-PBB"),
                   data.frame(rrbs.ef,dataset="RRBS-ES"),
                   data.frame(wgbs.ef,dataset="WGBS-BLUE"))
str(eval.all)

eval.all <- eval.all[(eval.all$model.p %in% c("p.b","p.s","p.n","p.q","p.limma")),]

eval.all$model <- factor(recode(eval.all$model.p,  p.b="beta", p.s="simplex",
                                  p.n="normal", 
                                  p.q="quantile",p.limma="limma"),
                           levels=c("beta", "simplex","normal","quantile","limma"))

eval.all$dist <- factor(eval.all$dist, levels=c("simplex","beta","normal"))

p <- ggplot(data=eval.all,aes(x=N, y=jaccard, colour=model)) +
  geom_point(size = 1)+ geom_line(size=0.6) +ylim(0,1)
p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
           axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
           axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
  scale_x_continuous(breaks=as.numeric(names(res.list)))

all <- p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
           axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
           axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
  scale_x_continuous(breaks=as.numeric(names(res.list))) + facet_grid (dist ~ dataset)

all

ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.facets.050719.png"), all, width = 20, height = 16, units = "cm")
ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.facets.050719.pdf"), all, width = 20, height = 16, units = "cm")

#faig subset de 0 a 30
eval.all.0.30 <- eval.all[eval.all$N %in% c(3,5,10,30),]

#i faig el mateix plot
p <- ggplot(data=eval.all.0.30,aes(x=N, y=jaccard, colour=model)) +
  geom_point(size = 1)+ geom_line(size=0.6) +ylim(0,1)
p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
           axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
           axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
  scale_x_continuous(breaks=as.numeric(names(res.list)))

all.0.30 <- p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
                  axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
                  axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
  scale_x_continuous(breaks=as.numeric(names(res.list))) + facet_grid (dist ~ dataset)

all.0.30

ggsave(file=file.path(resultsDir,"SuppFig3.jaccard.simulations.0.30.facets.050719.png"), all.0.30, width = 20, height = 16, units = "cm")
ggsave(file=file.path(resultsDir,"SuppFig3.jaccard.simulations.0.30.facets.050719.pdf"), all.0.30, width = 20, height = 16, units = "cm")


###############################################################################
################## provo d'afegir la beta-binomial per a NGS ##################
###############################################################################

#provo els beta.bin
eval_flatten_bb <- function(data){
  is <- c(3,5,10,30,100,500)
  js <- c("beta-binomial")
  evals <- NULL
  for (i in 1:6){ #N de cada simulació
    for (j in 1){
      mod <- t(models.eval(models=res.list[[i]][[j]], dml.r=100,alpha=0.05, adjust=TRUE))
      evals <- rbind(evals,
                     data.frame(mod, model.p=rownames(mod),N=is[i],dist=js[j],stringsAsFactors = F))
    }
  }
  return(evals)
}  

load(file=file.path(RRBS188_data,"simulated.models.bb.list.RData"))
rrbs.ef.bb <- eval_flatten_bb(res.list)
load(file=file.path(WGBS81_data,"simulated.models.bb.list.RData"))
wgbs.ef.bb <- eval_flatten_bb(res.list)

eval.all.t <- rbind(data.frame(gse50.ef,dataset="a450k-smoking"),
                    data.frame(gse11.ef,dataset="aEPIC-PBB"),
                    data.frame(rrbs.ef,dataset="RRBS-ES"),
                    data.frame(rrbs.ef.bb,dataset="RRBS-ES"),
                    data.frame(wgbs.ef,dataset="WGBS-BLUE"),
                    data.frame(wgbs.ef.bb,dataset="WGBS-BLUE"))

eval.all <- eval.all.t[(eval.all.t$model.p %in% c("p.b","p.s","p.n","p.q","p.limma","p.bb")),]

eval.all$model <- factor(recode(eval.all$model.p,  p.b="beta", p.s="simplex",
                                p.n="normal", p.q="quantile",p.limma="limma",p.bb="beta-binomial"),
                         levels=c("beta", "simplex","normal","quantile","limma","beta-binomial"))

eval.all$dist <- factor(eval.all$dist, levels=c("simplex","beta","normal","beta-binomial"))

p <- ggplot(data=eval.all,aes(x=N, y=jaccard, colour=model)) +
  geom_point(size = 1)+ geom_line(size=0.6) +ylim(0,1)
p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
           axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
           axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
  scale_x_continuous(breaks=as.numeric(names(res.list)))

all <- p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
                  axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
                  axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
  scale_x_continuous(breaks=as.numeric(names(res.list))) + facet_grid (dist ~ dataset)

all


ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.withbetabin.facets.050719.png"), all, width = 20, height = 16, units = "cm")
#ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.facets.120519.pdf"), all, width = 20, height = 16, units = "cm")

#faig subset de 0 a 30
eval.all.0.30 <- eval.all[eval.all$N %in% c(3,5,10,30),]

#i faig el mateix plot
p <- ggplot(data=eval.all.0.30,aes(x=N, y=jaccard, colour=model)) +
  geom_point(size = 1)+ geom_line(size=0.6) +ylim(0,1)
p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
           axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
           axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
  scale_x_continuous(breaks=as.numeric(names(res.list)))

all.0.30 <- p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
                       axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
                       axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
  scale_x_continuous(breaks=as.numeric(names(res.list))) + facet_grid (dist ~ dataset)

all.0.30

ggsave(file=file.path(resultsDir,"SuppFig3.jaccard.simulations.withbetabin.0.30.facets.050719.png"), all.0.30, width = 20, height = 16, units = "cm")

# #provo de fer només la part dels conjunts grans, a veure si es veu més clar: no queda gaire bé...
# eval.all.100.500 <- eval.all[eval.all$N %in% c(100,500),]
# 
# #i faig el mateix plot
# p <- ggplot(data=eval.all.100.500,aes(x=N, y=jaccard, colour=model)) +
#   geom_point(size = 1)+ geom_line(size=0.6) +ylim(0,1)
# p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
#            axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
#            axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
#   scale_x_continuous(breaks=as.numeric(names(res.list)))
# 
# all.100.500 <- p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
#                        axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
#                        axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
#   scale_x_continuous(breaks=as.numeric(names(res.list))) + facet_grid (dist ~ dataset)
# 
# all.100.500


#####################################################################################
########################### TAble 2 ###############################################
#####################################################################################
#primer extrec tota la info 
#sselecciono les columnes..per eliminar p.s.gamlss
selec <- c("p.s","p.b","p.sinf","p.binf","p.n","p.l","p.q","p.limma" )


load(file=file.path(GSE50660_data,"simulated.models.list.RData"))
evals.GSE50.s <- NULL
evals.GSE50.b <- NULL
evals.GSE50.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.GSE50.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
  evals.GSE50.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
  evals.GSE50.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
}

load(file=file.path(GSE116339_data,"simulated.models.list.RData"))
evals.GSE11.s <- NULL
evals.GSE11.b <- NULL
evals.GSE11.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.GSE11.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
  evals.GSE11.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
  evals.GSE11.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
}

#aqui de fet no cal fer el selec pq no hi ha p.s.gamlss
load(file=file.path(RRBS188_data,"simulated.models.list.RData"))
evals.RRBS.s <- NULL
evals.RRBS.b <- NULL
evals.RRBS.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.RRBS.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
  evals.RRBS.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
  evals.RRBS.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
}

load(file=file.path(WGBS81_data,"simulated.models.list.RData"))
evals.WGBS81.s <- NULL
evals.WGBS81.b <- NULL
evals.WGBS81.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.WGBS81.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
  evals.WGBS81.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
  evals.WGBS81.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))[selec,]
}

create.table <- function (n){
  tab.s <-cbind(evals.GSE50.s[[n]],evals.GSE11.s[[n]],
          evals.RRBS.s[[n]],evals.WGBS81.s[[n]])
  tab.b <-cbind(evals.GSE50.b[[n]],evals.GSE11.b[[n]],
                evals.RRBS.b[[n]],evals.WGBS81.b[[n]])
  tab.n <-cbind(evals.GSE50.n[[n]],evals.GSE11.n[[n]],
                evals.RRBS.n[[n]],evals.WGBS81.n[[n]])
  tab <- rbind(tab.s,tab.b,tab.n)
  rownames(tab) <- recode(rownames(tab),  p.b="beta", p.binf="beta.inf",p.s="simplex",
                                              p.sinf="simplex.inf", p.n="normal", p.l="logistic",
                                              p.q="quantile",p.limma="limma")
  return(tab)
}

eval.3 <-create.table(1)
write.csv2(eval.3,file=file.path(resultsDir,"Performances.N3.csv"), row.names = T, quote = F)
eval.5 <-create.table(2)
write.csv2(eval.5,file=file.path(resultsDir,"Performances.N5.csv"), row.names = T, quote = F)
eval.10 <-create.table(3)
write.csv2(eval.10,file=file.path(resultsDir,"Performances.N10.csv"), row.names = T, quote = F)
eval.30 <-create.table(4)
write.csv2(eval.30,file=file.path(resultsDir,"Performances.N30.csv"), row.names = T, quote = F)
eval.100 <-create.table(5)
write.csv2(eval.100,file=file.path(resultsDir,"Performances.N100.csv"), row.names = T, quote = F)
eval.500 <-create.table(6)
write.csv2(eval.500,file=file.path(resultsDir,"Performances.N500.csv"), row.names = T, quote = F)


