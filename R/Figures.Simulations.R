#Figures and tables from results of simulations
#Figure3: jaccard of all data sets for non inflated models N=3-500
#Supp figure3: jaccard of all data sets for non inflated models N=3,5,10 and 30
#Table 2: Summary of performance measures for simulations N=100
#Supplementary tables 3-7:  Summary of performance measures for simulations N=3,5,10,30,500

#header amb tot
source("D:/Doctorat/Simplex/R/Figures.HeaderScript.R")


########################################################################################
######### Figure 4: Resultats de les simulacions: TOTES ###########
########################################################################################

#Fn per plotar les corbes
#la fn sense les inflated està abaix!!!
#per a fer nomes de 3 a 30 descomentar línia scale_x.... dins la fn  fn
#fn sense les inflated (i comptant que als evals hi ha gamlss per simplex), es diu igual!!!
plot.eval.curves <- function(evals,measure=jaccard,ylim_0=0,ylim_1=1){
  measure <- enquo(measure)
  
  eval.all <- cbind (do.call("rbind",evals),N=as.numeric(rep(names(res.list),each=9)))
  eval.all<- data.frame(eval.all,model.p=rownames(eval.all),stringsAsFactors = F)
  
  #aquí elimino les inflated 
  eval.all <- eval.all[!(eval.all$model.p %in% c("p.binf","p.sinf","p.s.gamlss")),]
  
  eval.all$model.p <- factor(recode(eval.all$model.p,  p.b="beta", p.s="simplex",
                                    p.n="normal", p.l="logistic", 
                                    p.q="quantile",p.limma="limma"),
                             levels=c("beta", "simplex","normal","logistic","quantile","limma"))
  
  p <- ggplot(data=eval.all,aes(x=N, y=!!measure, colour=model.p)) +
    geom_point(size = 1)+ geom_line(size=0.6) +ylim(ylim_0,ylim_1)
  p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
             axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
             axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
   scale_x_continuous(breaks=as.numeric(names(res.list)))
 # scale_x_continuous(breaks=c(3,5,10,30),limits=c(0,30)) #la poso i la trec segons necessitat
}    
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

load(file=file.path(GSE50660_data,"simulated.models.list.RData"))
evals.s <- NULL
evals.b <- NULL
evals.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
}

gse50.s <- plot.eval.curves(evals.s)
gse50.b <- plot.eval.curves(evals.b)
gse50.n <- plot.eval.curves(evals.n)

#
load(file=file.path(GSE116339_data,"simulated.models.list.RData"))
evals.s <- NULL
evals.b <- NULL
evals.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
}
gse11.s <- plot.eval.curves(evals.s)
gse11.b <- plot.eval.curves(evals.b)
gse11.n <- plot.eval.curves(evals.n)

#
load(file=file.path(RRBS216_data,"simulated.models.list.RData"))
evals.s <- NULL
evals.b <- NULL
evals.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
}
rrbs.s <- plot.eval.curves(evals.s)
rrbs.b <- plot.eval.curves(evals.b)
rrbs.n <- plot.eval.curves(evals.n)

#
load(file=file.path(WGBS81_data,"simulated.models.list.RData"))
evals.s <- NULL
evals.b <- NULL
evals.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
}
wgbs.s <- plot.eval.curves(evals.s)
wgbs.b <- plot.eval.curves(evals.b)
wgbs.n <- plot.eval.curves(evals.n)

mylegend<-g_legend(gse50.s) 

all <- grid.arrange(arrangeGrob(gse50.s+ theme(legend.position="none"),gse11.s+ theme(legend.position="none"),rrbs.s+ theme(legend.position="none"), wgbs.s+ theme(legend.position="none"),    
                                gse50.b+ theme(legend.position="none"),gse11.b+ theme(legend.position="none"), rrbs.b+ theme(legend.position="none"), wgbs.b+ theme(legend.position="none"),
                                gse50.n+ theme(legend.position="none"),gse11.n+ theme(legend.position="none"), rrbs.n+ theme(legend.position="none"), wgbs.n+ theme(legend.position="none")+ theme(legend.position="none"),
                                nrow=3,ncol=4),
                    legend=mylegend,nrow=2,heights=c(10, 1))
plot(all)

ggsave(file=file.path(resultsDir,"Fig3.jaccard.simulations.datasets.0.500.pdf"), all, width = 20, height = 16, units = "cm")
#ggsave(file=file.path(resultsDir,"SuppFig3.jaccard.simulations.datasets.0.30.pdf"), all, width = 20, height = 16, units = "cm")

# ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.allwithsimplexgamlss.datasets.0.500.png"), all, width = 20, height = 16, units = "cm")
# ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.noinfl.datasets.0.30.png"), all, width = 20, height = 16, units = "cm")
#ggsave(file=file.path(resultsDir,"SuppFig2.jaccard.simulations.datasets.0.30.png"), all, width = 20, height = 16, units = "cm")
#ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.noinfl.datasets.0.500.png"), all, width = 20, height = 16, units = "cm")
#ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.datasets.0.500.png"), all, width = 20, height = 16, units = "cm")

#################################################
##################### TRUE POSITIVES
# load(file=file.path(GSE50660_data,"simulated.models.list.RData"))
# evals.s <- NULL
# evals.b <- NULL
# evals.n <- NULL
# for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
#   evals.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
#   evals.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
#   evals.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
# }
# 
# gse50.s <- plot.eval.curves(evals.s,measure=TP,ylim_1=100)
# gse50.b <- plot.eval.curves(evals.b,measure=TP,ylim_1=100)
# gse50.n <- plot.eval.curves(evals.n,measure=TP,ylim_1=100)
# 
# #
# load(file=file.path(GSE116339_data,"simulated.models.list.RData"))
# evals.s <- NULL
# evals.b <- NULL
# evals.n <- NULL
# for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
#   evals.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
#   evals.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
#   evals.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
# }
# gse11.s <- plot.eval.curves(evals.s,measure=TP,ylim_1=100)
# gse11.b <- plot.eval.curves(evals.b,measure=TP,ylim_1=100)
# gse11.n <- plot.eval.curves(evals.n,measure=TP,ylim_1=100)
# 
# #
# load(file=file.path(RRBS216_data,"simulated.models.list.RData"))
# evals.s <- NULL
# evals.b <- NULL
# evals.n <- NULL
# for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
#   evals.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
#   evals.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
#   evals.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
# }
# rrbs.s <- plot.eval.curves(evals.s,measure=TP,ylim_1=100)
# rrbs.b <- plot.eval.curves(evals.b,measure=TP,ylim_1=100)
# rrbs.n <- plot.eval.curves(evals.n,measure=TP,ylim_1=100)
# 
# #
# load(file=file.path(WGBS81_data,"simulated.models.list.RData"))
# evals.s <- NULL
# evals.b <- NULL
# evals.n <- NULL
# for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
#   evals.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
#   evals.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
#   evals.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
# }
# wgbs.s <- plot.eval.curves(evals.s,measure=TP,ylim_1=100)
# wgbs.b <- plot.eval.curves(evals.b,measure=TP,ylim_1=100)
# wgbs.n <- plot.eval.curves(evals.n,measure=TP,ylim_1=100)
# 
# mylegend<-g_legend(gse50.s) 
# 
# all <- grid.arrange(arrangeGrob(gse50.s+ theme(legend.position="none"),gse11.s+ theme(legend.position="none"),rrbs.s+ theme(legend.position="none"), wgbs.s+ theme(legend.position="none"),    
#                                 gse50.b+ theme(legend.position="none"),gse11.b+ theme(legend.position="none"), rrbs.b+ theme(legend.position="none"), wgbs.b+ theme(legend.position="none"),
#                                 gse50.n+ theme(legend.position="none"),gse11.n+ theme(legend.position="none"), rrbs.n+ theme(legend.position="none"), wgbs.n+ theme(legend.position="none")+ theme(legend.position="none"),
#                                 nrow=3,ncol=4),
#                     legend=mylegend,nrow=2,heights=c(10, 1))
# plot(all)
# #ggsave(file=file.path(resultsDir,"SuppFig2.jaccard.simulations.datasets.0.30.png"), all, width = 20, height = 16, units = "cm")
# ggsave(file=file.path(resultsDir,"Fig4.TP.simulations.datasets.0.30.png"), all, width = 20, height = 16, units = "cm")

#####################################################################################
########################### TAble 2 ###############################################
#####################################################################################
#primer extrec tota la info 

load(file=file.path(GSE50660_data,"simulated.models.list.RData"))
evals.GSE50.s <- NULL
evals.GSE50.b <- NULL
evals.GSE50.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.GSE50.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.GSE50.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.GSE50.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
}

load(file=file.path(GSE116339_data,"simulated.models.list.RData"))
evals.GSE11.s <- NULL
evals.GSE11.b <- NULL
evals.GSE11.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.GSE11.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.GSE11.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.GSE11.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
}

load(file=file.path(RRBS216_data,"simulated.models.list.RData"))
evals.RRBS216.s <- NULL
evals.RRBS216.b <- NULL
evals.RRBS216.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.RRBS216.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.RRBS216.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.RRBS216.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
}

load(file=file.path(WGBS81_data,"simulated.models.list.RData"))
evals.WGBS81.s <- NULL
evals.WGBS81.b <- NULL
evals.WGBS81.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.WGBS81.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.WGBS81.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.WGBS81.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
}

create.table <- function (n){
  tab.s <-cbind(evals.GSE50.s[[n]],evals.GSE11.s[[n]],
          evals.RRBS216.s[[n]],evals.WGBS81.s[[n]])
  tab.b <-cbind(evals.GSE50.b[[n]],evals.GSE11.b[[n]],
                evals.RRBS216.b[[n]],evals.WGBS81.b[[n]])
  tab.n <-cbind(evals.GSE50.n[[n]],evals.GSE11.n[[n]],
                evals.RRBS216.n[[n]],evals.WGBS81.n[[n]])
  tab <- rbind(tab.s,tab.b,tab.n)
  rownames(tab) <- recode(rownames(tab),  p.b="beta", p.binf="beta.inf",p.s="simplex",
                                              p.sinf="simplex.inf", p.n="normal", p.l="logistic",
                                              p.q="quantile",p.limma="limma",p.s.gamlss="sim.gamlss")
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


