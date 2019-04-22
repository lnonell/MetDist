#Figure 4

workingDir<-"D:/Doctorat/Simplex/Data"
setwd(workingDir)
resultsDir <- file.path("./Summary")

#data sets paths
GSE50660_data<-"D:/Doctorat/Simplex/Data/GSE50660_Smoking"
GSE116339_data<-"D:/Doctorat/Simplex/Data/GSE116339_PBB"
RRBS216_data<-"D:/Doctorat/Simplex/Data/RRBS_216_Tissue_CL"
WGBS81_data<-"D:/Doctorat/Simplex/Data/WGBS_81_Bueprint" #a veure si al final el puc utilizar!

#llibreries i sources pels plots
library(ggplot2)
library(gridExtra)
library(dplyr)
source(file=file.path("D:/Doctorat/Simplex/R","SimulationFunctions.R")) #carrego directament la fn est.betabin.params

#fn per a obtenir la llegenda d'un grafic
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


########################################################################################
######### Figure 4: Resultats de les simulacions: TOTES ###########
########################################################################################

#Fn per plotar les corbes
#la fn sense les inflated està abaix!!!
plot.eval.curves <- function(evals,measure=jaccard,ylim_0=0,ylim_1=1){
  measure <- enquo(measure)
  
  eval.all <- cbind (do.call("rbind",evals),N=as.numeric(rep(names(res.list),each=8)))
  eval.all<- data.frame(eval.all,model.p=rownames(eval.all))
  eval.all$model.p <- factor(recode(eval.all$model.p,  p.b="beta", p.binf="beta.inf",p.s="simplex",
                                    p.sinf="simplex.inf", p.n="normal", p.l="logistic", 
                                    p.q="quantile",p.limma="limma"),
                             levels=c("beta", "beta.inf","simplex","simplex.inf","normal","logistic","quantile","limma"))
  
  p <- ggplot(data=eval.all,aes(x=N, y=!!measure, colour=model.p)) +
    geom_point(size = 1)+ geom_line(size=0.6) +ylim(ylim_0,ylim_1)
  p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
             axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
             axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
  scale_x_continuous(breaks=as.numeric(names(res.list)))
  #scale_x_continuous(breaks=c(3,5,10,30),limits=c(0,30)) #la poso i la trec segons necessitat
}    

#per si fem "zoom"
# p +  theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
# axis.text.x = element_text(size=8,angle=45), axis.text.y = element_text(size=8),
# axis.title = element_text(size=7), legend.key.size = unit(0.8,"line")) +
#   scale_x_continuous(breaks=c(3,5,10,30),limits=c(0,30)) 

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

ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.noinfl.datasets.0.30.png"), all, width = 20, height = 16, units = "cm")
#ggsave(file=file.path(resultsDir,"SuppFig2.jaccard.simulations.datasets.0.30.png"), all, width = 20, height = 16, units = "cm")
ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.noinfl.datasets.0.500.png"), all, width = 20, height = 16, units = "cm")
#ggsave(file=file.path(resultsDir,"Fig4.jaccard.simulations.datasets.0.500.png"), all, width = 20, height = 16, units = "cm")

#################################################
##################### TRUE POSITIVES
load(file=file.path(GSE50660_data,"simulated.models.list.RData"))
evals.s <- NULL
evals.b <- NULL
evals.n <- NULL
for (i in 1:6) { #transposo per a poder-ho posar com interessa, doncs voldré plotar index
  evals.s[[i]] <- t(models.eval(models=res.list[[i]][[1]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.b[[i]] <- t(models.eval(models=res.list[[i]][[2]], dml.r=100,alpha=0.05, adjust=TRUE))
  evals.n[[i]] <- t(models.eval(models=res.list[[i]][[3]], dml.r=100,alpha=0.05, adjust=TRUE))
}

gse50.s <- plot.eval.curves(evals.s,measure=TP,ylim_1=100)
gse50.b <- plot.eval.curves(evals.b,measure=TP,ylim_1=100)
gse50.n <- plot.eval.curves(evals.n,measure=TP,ylim_1=100)

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
gse11.s <- plot.eval.curves(evals.s,measure=TP,ylim_1=100)
gse11.b <- plot.eval.curves(evals.b,measure=TP,ylim_1=100)
gse11.n <- plot.eval.curves(evals.n,measure=TP,ylim_1=100)

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
rrbs.s <- plot.eval.curves(evals.s,measure=TP,ylim_1=100)
rrbs.b <- plot.eval.curves(evals.b,measure=TP,ylim_1=100)
rrbs.n <- plot.eval.curves(evals.n,measure=TP,ylim_1=100)

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
wgbs.s <- plot.eval.curves(evals.s,measure=TP,ylim_1=100)
wgbs.b <- plot.eval.curves(evals.b,measure=TP,ylim_1=100)
wgbs.n <- plot.eval.curves(evals.n,measure=TP,ylim_1=100)

mylegend<-g_legend(gse50.s) 

all <- grid.arrange(arrangeGrob(gse50.s+ theme(legend.position="none"),gse11.s+ theme(legend.position="none"),rrbs.s+ theme(legend.position="none"), wgbs.s+ theme(legend.position="none"),    
                                gse50.b+ theme(legend.position="none"),gse11.b+ theme(legend.position="none"), rrbs.b+ theme(legend.position="none"), wgbs.b+ theme(legend.position="none"),
                                gse50.n+ theme(legend.position="none"),gse11.n+ theme(legend.position="none"), rrbs.n+ theme(legend.position="none"), wgbs.n+ theme(legend.position="none")+ theme(legend.position="none"),
                                nrow=3,ncol=4),
                    legend=mylegend,nrow=2,heights=c(10, 1))
plot(all)
#ggsave(file=file.path(resultsDir,"SuppFig2.jaccard.simulations.datasets.0.30.png"), all, width = 20, height = 16, units = "cm")
ggsave(file=file.path(resultsDir,"Fig4.TP.simulations.datasets.0.30.png"), all, width = 20, height = 16, units = "cm")


#fn sense les inflated
plot.eval.curves <- function(evals,measure=jaccard,ylim_0=0,ylim_1=1){
  measure <- enquo(measure)
  
  eval.all <- cbind (do.call("rbind",evals),N=as.numeric(rep(names(res.list),each=8)))
  eval.all<- data.frame(eval.all,model.p=rownames(eval.all),stringsAsFactors = F)
  
  #aquí elimino les inflated 
  eval.all <- eval.all[!(eval.all$model.p %in% c("p.binf","p.sinf")),]
  
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
  #scale_x_continuous(breaks=c(3,5,10,30),limits=c(0,30)) #la poso i la trec segons necessitat
}    
