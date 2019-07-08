#Figure 2: AIC for all data sets comparing simplex, beta and normal
#I will create the same but adding beta binomial for NGS
#11/05/19 ho actualitzo tot
#7//19: supp table 7: resultats AIC i KS pels 4 datasets

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
head(best.dist.all)
a <- aic.plot(best.dist.all,tit="A")
a

load(file=file.path(GSE116339_data,"best.dist.all.RData"))
head(best.dist.all)
b <- aic.plot(best.dist.all,tit="B")
b

load(file=file.path(RRBS188_data,"best.dist.all.RData"))
head(best.dist.all)
c <- aic.plot(best.dist.all,tit="C")
c

load(file=file.path(RRBS188_data,"best.dist.betabin.RData"))
c.bb <- aic.bb.plot(best.dist.all,best.dist.all.betabin,tit="C")
c.bb

load(file=file.path(WGBS81_data,"best.dist.all.RData"))
head(best.dist.all)
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
#els pics del final són els cromosomes sexuals.

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


######### KS results ##############
ks.plot <- function(best.dist.all,tit="A"){
  ll <- nrow(best.dist.all)
  ks.plot <- data.frame(x=1:ll,logp=-log10(best.dist.all[,4]), model="simplex")
  ks.plot <- rbind(ks.plot,
                    data.frame(x=1:ll,logp=-log10(best.dist.all[,5]), model="beta"),
                    data.frame(x=1:ll,logp=-log10(best.dist.all[,6]), model="normal"))
  
  p <-ggplot(ks.plot) + geom_smooth(aes(x=x, y=logp, color=model)) + xlab("CpG")+ ylab("-log(p-value)")+
    theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=7,angle=45), axis.text.y = element_text(size=8),
          axis.title = element_text(size=8), legend.key.size = unit(0.8,"line"))+
    scale_color_manual(values=colors) + ggtitle(tit)
  p
}


load(file=file.path(GSE50660_data,"best.dist.all.RData")) #per si de cas
head(best.dist.all)
gse50 <- t(apply(best.dist.all,1,best.aic.ks))
gse50.best <-rbind(table(gse50[,1]),table(gse50[,2])) #1 aic 2 ks best
a <- ks.plot(best.dist.all,tit="A")
a

load(file=file.path(GSE116339_data,"best.dist.all.RData"))
head(best.dist.all)
gse11 <- t(apply(best.dist.all,1,best.aic.ks))
gse11.best <-rbind(table(gse11[,1]),table(gse11[,2])) #1 aic 2 ks best
b <- ks.plot(best.dist.all,tit="B")
b

load(file=file.path(RRBS188_data,"best.dist.all.RData"))
head(best.dist.all)
RRBS <- t(apply(best.dist.all,1,best.aic.ks))
RRBS.best <-rbind(table(RRBS[,1]),table(RRBS[,2])) #1 aic 2 ks best
c <- ks.plot(best.dist.all,tit="C")
c

load(file=file.path(WGBS81_data,"best.dist.all.RData"))
head(best.dist.all)
WGBS <- t(apply(best.dist.all,1,best.aic.ks))
WGBS.best <-rbind(table(WGBS[,1]),table(WGBS[,2])) #1 aic 2 ks best
d <- ks.plot(best.dist.all,tit="D")
d

mylegend<-g_legend(a)
ga <- grid.arrange(arrangeGrob(a + theme(legend.position="none"), #potser treure això de la legend??
                               b + theme(legend.position="none"),
                               c + theme(legend.position="none"),
                               d + theme(legend.position="none"), 
                               nrow=2),
                   mylegend, nrow=2,heights=c(10, 1))

ggsave(file=file.path(resultsDir,"KS.bestdist.datasets.110519.pdf"), ga, width = 14, height = 14, units = "cm")

#table best
table.best <- as.data.frame(rbind(gse50.best,gse11.best,RRBS.best,WGBS.best))

#arrangements
table.best$Method <-rep(c("AIC","K-S"),4)
table.best$"Data set" <- rep(c("a450k-smoking","aEPIC-PBB","RRBS-ES","WGBS-BLUE"),each=2)

table.best <- table.best[,c(5,4,3,1,2)]
colnames(table.best)[3:5] <-c("simplex","beta","normal")
table.best

write.csv2(table.best, file=file.path(resultsDir,"SuppTable7.table.bestdist.AIC.KS.070719.csv"),row.names = F)
