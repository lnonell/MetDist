#Figure 2: AIC for all data sets comparing simplex, beta and normal
#I will create the same but adding beta binomial for NGS
#11/05/19 ho actualitzo tot

#header amb tot
source("D:/Doctorat/Simplex/MetDist/R/Figures.HeaderScript.R")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


######################################################################
########### Figure 2: Com ajusten les dades per a cada data set ######
######################################################################
colors <- c("#E69F00", "#56B4E9", "#999999")
colors2 <- c("#E69F00", "#56B4E9","springgreen", "#999999")

#fer fn i per a cada data set grar un plot, que després ajuntarem
aic.plot <- function(best.dist.all,tit="A"){
  ll <- nrow(best.dist.all)
  aic.plot <- data.frame(x=1:ll,aic=best.dist.all[,1], model="simplex")
  aic.plot <- rbind(aic.plot,
                    data.frame(x=1:ll,aic=best.dist.all[,2], model="beta"),
                    data.frame(x=1:ll,aic=best.dist.all[,3], model="normal"))
  
  p <-ggplot(aic.plot) + geom_smooth(aes(x=x, y=aic, color=model)) + xlab("CpG")+ ylab("AIC")+
    theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=7,angle=45), axis.text.y = element_text(size=8),
          axis.title = element_text(size=8), legend.key.size = unit(0.8,"line"))+
    scale_color_manual(values=colors) + ggtitle(tit)
  p
}

aic.bb.plot <- function(best.dist.all,best.dist.all.betabin,tit="A"){
  ll <- nrow(best.dist.all)
  aic.plot <- data.frame(x=1:ll,aic=best.dist.all[,1], model="simplex")
  aic.plot <- rbind(aic.plot,
                    data.frame(x=1:ll,aic=best.dist.all[,2], model="beta"),
                    data.frame(x=1:ll,aic=best.dist.all.betabin[,1], model="beta-binomial"),
                    data.frame(x=1:ll,aic=best.dist.all[,3], model="normal"))
  
  p <-ggplot(aic.plot) +  geom_smooth(aes(x=x, y=aic, color=model)) + xlab("CpG")+ylab("AIC")+
    theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=7,angle=45), axis.text.y = element_text(size=8),
          axis.title = element_text(size=8), legend.key.size = unit(0.8,"line"))+
    scale_color_manual(values=colors2) + ggtitle(tit)
  p
}


load(file=file.path(GSE50660_data,"best.dist.all.RData")) #per si de cas
a <- aic.plot(best.dist.all,tit="A")
a

load(file=file.path(GSE116339_data,"best.dist.all.RData"))
b <- aic.plot(best.dist.all,tit="B")
b

load(file=file.path(RRBS188_data,"best.dist.all.RData"))
c <- aic.plot(best.dist.all,tit="C")
c

load(file=file.path(RRBS188_data,"best.dist.betabin.RData"))
c.bb <- aic.bb.plot(best.dist.all,best.dist.all.betabin,tit="C")
c.bb

load(file=file.path(WGBS81_data,"best.dist.all.RData"))
d <- aic.plot(best.dist.all,tit="D")
d

load(file=file.path(WGBS81_data,"best.dist.betabin.RData"))
d.bb <- aic.bb.plot(best.dist.all,best.dist.all.betabin,tit="D")
d.bb

#per a tenir una única llegenda
mylegend<-g_legend(a)

#png(file=file.path(resultsDir,"Fig2.modelest.datasets.png"), width=480, height=480,res=100) #param res és l'important!!
ga <- grid.arrange(arrangeGrob(a + theme(legend.position="none"), #potser treure això de la legend??
                               b + theme(legend.position="none"),
                               c + theme(legend.position="none"),
                               d + theme(legend.position="none"), 
                               nrow=2),
                   mylegend, nrow=2,heights=c(10, 1))
#dev.off()
ggsave(file=file.path(resultsDir,"Fig2.modelest.datasets.110519.pdf"), ga, width = 14, height = 14, units = "cm")


#i ara la mateixa però amb la betabinomial per NGS
mylegend<-g_legend(d)
ga <- grid.arrange(arrangeGrob(a + theme(legend.position="none"), #potser treure això de la legend??
                               b + theme(legend.position="none"),
                               c.bb + theme(legend.position="none"),
                               d.bb + theme(legend.position="none"), 
                               nrow=2),
                   mylegend, nrow=2,heights=c(10, 1))
#dev.off()
ggsave(file=file.path(resultsDir,"Fig2.modelest.datasets.betabin.110519.pdf"), ga, width = 14, height = 14, units = "cm")
