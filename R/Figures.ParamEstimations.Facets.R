#figures and tables from parameter estimates results
#Suppfig 1: Correlations MOM MLE
#Suggfig 2: Densities for the different parameters
#Suppfig XX: betabi correls
#PROVO Facets en les correls pero no va bé pq les escales son diferents per fila i per columna 
#23/4/19 canvio RRBS216 a RRBS188
#11/05/19 ho actualitzo tot

#header amb tot
source("D:/Doctorat/Simplex/Metdist/R/Figures.HeaderScript.R")

######################################################################
########### Param estimations from all datasets  ###########
######################################################################
# carregar totes les dades  i fer correlacions entre mom i mle i guardar-ho com a supplementary
# a més fer un plot per a totes les mu's i sigmes (simplex) de cada data base
# un altre per a s1 i s2 (betes) de cada data base
# i el mateix per a les normals

#############################################
#Supp Figure 1: Correlacions entre els params
#############################################
load(file=file.path(GSE50660_data,"params.est.all.RData"))
est.params.gse50 <- est.params
rm(est.params)
load(file=file.path(GSE116339_data,"params.est.all.RData"))
est.params.gse11 <- est.params
rm(est.params)
#load(file=file.path(RRBS216_data,"params.est.all.090319.RData")) queda substituida per RRBS188
load(file=file.path(RRBS188_data,"params.est.all.RData"))
est.params.rrbs <- rnb.est.params
rm(rnb.est.params)
colnames(est.params.rrbs)
load(file=file.path(WGBS81_data,"params.est.all.RData"))
est.params.wgbs <- rnb.est.params
rm(rnb.est.params)

# apply(est.params.gse50, 2, range, na.rm=T)
# apply(est.params.gse11, 2, range, na.rm=T)
# apply(est.params.rrbs, 2, range, na.rm=T)
# apply(est.params.wgbs, 2, range, na.rm=T)

est.all <- rbind(data.frame(est.params.gse50,dataset="A450k-smoking"),
                  data.frame(est.params.gse11,dataset="EPIC-PBB"),
                  data.frame(est.params.rrbs,dataset="RRBS-ES"),
                  data.frame(est.params.wgbs,dataset="WGBS-BLUE"))

#hem de passar les columnes a files per blocks
# library(tidyr)
# 
# #primers params (x a les grafiques)
# est.all1 <- est.all %>%  gather(model.p1,val1,c(s.mle.mu,s.mle.sig,b.mle.s1,b.mle.s2,n.mle.m,n.mle.sd))
# #segons params (y a les grafiques)
# est.all2 <- est.all %>%  gather(model.p2,val2,c(s.zoip.mu,s.zoip.sig,b.mom.s1,b.mom.s2,n.mom.m,n.mom.sd))
# 
# est.all4facets <- data.frame(est.all1[,c("dataset","model.p1","val1")],
#                              est.all2[,c("model.p2","val2")])
# est.all4facets$model <- ifelse(est.all4facets$model.p1=="s.mle.mu","simplex mu",
#                                ifelse(est.all4facets$model.p1=="s.mle.sig","simplex sigma",
#                                       ifelse(est.all4facets$model.p1=="b.mle.s1","beta shape1",
#                                              ifelse(est.all4facets$model.p1=="b.mle.s2","beta shape2",
#                                                     ifelse(est.all4facets$model.p1=="n.mle.m","normal mean",
#                                                            ifelse(est.all4facets$model.p1=="n.mle.sd","normal sd",
#                                       NA))))))
# 
# #comprovo
# table(est.all4facets$model.p1,est.all4facets$model.p2) #ara sí
# table(est.all4facets$model.p1,est.all4facets$model) # sí
# 
# est.all4facets$model <- factor(est.all4facets$model,
#                          levels=c("simplex mu", "simplex sigma",
#                                   "beta shape1", "beta shape2",
#                                   "normal mean", "normal sd"))

#sembla que estan a escales molt diferents i per tant no es poden ajuntar!!!
#tampoc veic com posar els colors 
colors <- c("green4", "blue", "red2", "orange")

# #potser cal
# g <- ggplot(data = est.all4facets, aes(x = val1, y = val2)) + 
#   geom_point() +  geom_smooth(method = "lm",color="grey43" )
# h <- g + xlab("")+ylab("")  + facet_grid (model ~ dataset, scales="free")
#     
# ggsave(file=file.path(resultsDir,"Supp1.facets.test.kk.png"), h, width = 16, height = 20, units = "cm")

############## correlacions
plot.cor <- function(df, x, y, color,tit){
  x <- enquo(x)
  y <- enquo(y)
  g <- ggplot(data = df, aes(x = !!x, y = !!y)) + 
    geom_point(color=color) +  geom_smooth(method = "lm",color="grey43" )
  g + labs(title=tit)+xlab("")+ylab("")  +
    theme(plot.title = element_text(size=8,hjust = 0.5),
          axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))
}

# g <- ggplot(data = df, aes(x = s.mle.mu, y =s.zoip.mu)) + 
#   geom_point(color=colors[1]) +
#   geom_smooth(method = "lm",color="grey43" )
# g + labs(title="simplex.mu")+xlab("")+ylab("")  +
#   theme(plot.title = element_text(size=8,hjust = 0.5),
#         axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))

#els simplex correlacionem simplex i zoip, la resta mle i mom
s.c.m.gse50 <- plot.cor(df=est.params.gse50, x=s.mle.mu, y=s.zoip.mu, color=colors[1],tit="simplex mu")
s.c.s.gse50 <- plot.cor(df=est.params.gse50, x=s.mle.sig, y=s.zoip.sig, color=colors[1],tit="simplex sigma")
b.c.m.gse50 <- plot.cor(df=est.params.gse50, x=b.mle.s1, y=b.mom.s1, color=colors[1],tit="beta shape1")
b.c.s.gse50 <- plot.cor(df=est.params.gse50, x=b.mle.s2, y=b.mom.s2, color=colors[1],tit="beta shape2")
n.c.m.gse50 <- plot.cor(df=est.params.gse50, x=n.mle.m, y=n.mom.m, color=colors[1],tit="normal mean")
n.c.s.gse50 <- plot.cor(df=est.params.gse50, x=n.mle.sd, y=n.mom.sd, color=colors[1],tit="normal sd")

s.c.m.gse11 <- plot.cor(df=est.params.gse11, x=s.mle.mu, y=s.zoip.mu, color=colors[2],tit="simplex mu")
s.c.s.gse11 <- plot.cor(df=est.params.gse11, x=s.mle.sig, y=s.zoip.sig, color=colors[2],tit="simplex sigma")
b.c.m.gse11 <- plot.cor(df=est.params.gse11, x=b.mle.s1, y=b.mom.s1, color=colors[2],tit="beta shape1")
b.c.s.gse11 <- plot.cor(df=est.params.gse11, x=b.mle.s2, y=b.mom.s2, color=colors[2],tit="beta shape2")
n.c.m.gse11 <- plot.cor(df=est.params.gse11, x=n.mle.m, y=n.mom.m, color=colors[2],tit="normal mean")
n.c.s.gse11 <- plot.cor(df=est.params.gse11, x=n.mle.sd, y=n.mom.sd, color=colors[2],tit="normal sd")

s.c.m.rrbs <- plot.cor(df=est.params.rrbs, x=s.mle.mu, y=s.zoip.mu, color=colors[3],tit="simplex mu")
s.c.s.rrbs <- plot.cor(df=est.params.rrbs, x=s.mle.sig, y=s.zoip.sig, color=colors[3],tit="simplex sigma")
b.c.m.rrbs <- plot.cor(df=est.params.rrbs, x=b.mle.s1, y=b.mom.s1, color=colors[3],tit="beta shape1")
b.c.s.rrbs <- plot.cor(df=est.params.rrbs, x=b.mle.s2, y=b.mom.s2, color=colors[3],tit="beta shape2")
n.c.m.rrbs <- plot.cor(df=est.params.rrbs, x=n.mle.m, y=n.mom.m, color=colors[3],tit="normal mean")
n.c.s.rrbs <- plot.cor(df=est.params.rrbs, x=n.mle.sd, y=n.mom.sd, color=colors[3],tit="normal sd")

s.c.m.wgbs <- plot.cor(df=est.params.wgbs, x=s.mle.mu, y=s.zoip.mu, color=colors[4],tit="simplex mu")
s.c.s.wgbs <- plot.cor(df=est.params.wgbs, x=s.mle.sig, y=s.zoip.sig, color=colors[4],tit="simplex sigma")
b.c.m.wgbs <- plot.cor(df=est.params.wgbs, x=b.mle.s1, y=b.mom.s1, color=colors[4],tit="beta shape1")
b.c.s.wgbs <- plot.cor(df=est.params.wgbs, x=b.mle.s2, y=b.mom.s2, color=colors[4],tit="beta shape2")
n.c.m.wgbs <- plot.cor(df=est.params.wgbs, x=n.mle.m, y=n.mom.m, color=colors[4],tit="normal mean")
n.c.s.wgbs <- plot.cor(df=est.params.wgbs, x=n.mle.sd, y=n.mom.sd, color=colors[4],tit="normal sd")

#ajuntem: ho guardare en png pq pdf pot ser una bogeria
all <- grid.arrange(arrangeGrob(s.c.m.gse50,s.c.m.gse11, s.c.m.rrbs, s.c.m.wgbs,    
                                s.c.s.gse50,s.c.s.gse11, s.c.s.rrbs, s.c.s.wgbs,  
                                b.c.m.gse50,b.c.m.gse11, b.c.m.rrbs, b.c.m.wgbs,    
                                b.c.s.gse50,b.c.s.gse11, b.c.s.rrbs, b.c.s.wgbs,  
                                n.c.m.gse50,n.c.m.gse11, n.c.m.rrbs, n.c.m.wgbs,    
                                n.c.s.gse50,n.c.s.gse11, n.c.s.rrbs, n.c.s.wgbs,  
                                nrow=6,ncol=4), nrow=1)
#all
ggsave(file=file.path(resultsDir,"SuppFig1.correls.est.params.mle.datasets.110519.png"), all, width = 16, height = 20, units = "cm")
ggsave(file=file.path(resultsDir,"SuppFig1.correls.est.params.mle.datasets.110519.pdf"), all, width = 16, height = 20, units = "cm")

#############################################
#Supp Figure XXX: el mateix per betabinomial
#############################################

load(file=file.path(RRBS188_data,"betabin.params.est.RData"))
est.params.rrbs <- as.data.frame(rnb.betabin.est.params)
rm(rnb.betabin.est.params)
colnames(est.params.rrbs)

load(file=file.path(WGBS81_data,"betabin.params.est.RData"))
est.params.wgbs <- as.data.frame(rnb.betabin.est.params)
rm(rnb.betabin.est.params)

bb.c.m.rrbs <- plot.cor(df=est.params.rrbs, x=bb.mle.s1, y=bb.mom.s1, color=colors[3],tit="beta binomial s1")
bb.c.s.rrbs <- plot.cor(df=est.params.rrbs, x=bb.mle.s2, y=bb.mom.s2, color=colors[3],tit="beta binomial s2")

bb.c.m.wgbs <- plot.cor(df=est.params.wgbs, x=bb.mle.s1, y=bb.mom.s1, color=colors[4],tit="beta binomial s1")
bb.c.s.wgbs <- plot.cor(df=est.params.wgbs, x=bb.mle.s2, y=bb.mom.s2, color=colors[4],tit="beta binomial s2")

#ajuntem: ho guardare en png pq pdf pot ser una bogeria
all.bb <- grid.arrange(arrangeGrob(bb.c.m.rrbs, bb.c.s.rrbs,    
                                   bb.c.m.wgbs, bb.c.s.wgbs,  
                                   nrow=2,ncol=2), nrow=1)
#all
ggsave(file=file.path(resultsDir,"SuppFig4.correls.est.params.betbin.mle.datasets.120519.png"), all.bb, width = 16, height = 20, units = "cm")



#############################################
#Supp Figure 2: Density de les estimacions dels params
#############################################

plot.dens <- function(df, param, ds, log2=F){
  param <- enquo(param)
  ds <- enquo(ds)
  s.s <- ggplot(df, aes(x=!!param, fill=!!ds)) + 
    geom_density(alpha=.4) + theme(legend.position="bottom") + theme(legend.title=element_blank())
  if (log2) s.s + scale_color_manual(values=colors) + scale_fill_manual(values=colors) + scale_x_continuous(trans='log2')
  else s.s + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
}

#test no fill, només linies
# plot.dens <- function(df, param, ds, log2=F){
#   param <- enquo(param)
#   ds <- enquo(ds)
#   s.s <- ggplot(df, aes(x=!!param, color=!!ds)) + 
#     geom_density() 
#   if (log2) s.s + scale_color_manual(values=colors) +  scale_x_continuous(trans='log2')
#   else s.s + scale_color_manual(values=colors) 
# }

#ara ja generem els grafics

#simplex
df1 <- data.frame(ds="a450k-smoking",mu=est.params.gse50$s.mle.mu,sig=est.params.gse50$s.mle.sig)
df2 <- data.frame(ds="aEPIC-PBB",mu=est.params.gse11$s.mle.mu,sig=est.params.gse11$s.mle.sig)
df3 <- data.frame(ds="RRBS-ES",mu=est.params.rrbs$s.mle.mu,sig=est.params.rrbs$s.mle.sig) #tornar a carregar params si hem fet betabin
df4 <- data.frame(ds="WGBS-BLUE",mu=est.params.wgbs$s.mle.mu,sig=est.params.wgbs$s.mle.sig) #tornar a carregar params si hem fet betabin
df <- rbind(df1,df2,df3,df4)

s.m <- plot.dens(df, param=mu,ds=ds, log2 =F)
s.s <- plot.dens(df, param=sig,ds=ds, log2 =T)

#beta
df1 <- data.frame(ds="a450k-smoking",s1=est.params.gse50$b.mle.s1,s2=est.params.gse50$b.mle.s2)
df2 <- data.frame(ds="aEPIC-PBB",s1=est.params.gse11$b.mle.s1,s2=est.params.gse11$b.mle.s2)
df3 <- data.frame(ds="RRBS-ES",s1=est.params.rrbs$b.mle.s1,s2=est.params.rrbs$b.mle.s2)
df4 <- data.frame(ds="WGBS-BLUE",s1=est.params.wgbs$b.mle.s1,s2=est.params.wgbs$b.mle.s2)
df <- rbind(df1,df2,df3,df4)

b.m <- plot.dens(df, param=s1,ds=ds, log2 =T)
b.s <- plot.dens(df, param=s2,ds=ds, log2 =T)

#normal
df1 <- data.frame(ds="a450k-smoking",mu=est.params.gse50$n.mle.m,sd=est.params.gse50$n.mle.sd)
df2 <- data.frame(ds="aEPIC-PBB",mu=est.params.gse11$n.mle.m,sd=est.params.gse11$n.mle.sd)
df3 <- data.frame(ds="RRBS-ES",mu=est.params.rrbs$n.mle.m,sd=est.params.rrbs$n.mle.sd)
df4 <- data.frame(ds="WGBS-BLUE",mu=est.params.wgbs$n.mle.m,sd=est.params.wgbs$n.mle.sd)
df <- rbind(df1,df2,df3,df4)

n.m <- plot.dens(df, param=mu,ds=ds, log2 =F)
n.s <- plot.dens(df, param=sd,ds=ds, log2 =F)

#la llegenda
mylegend<-g_legend(n.s)

#monto la figura
all <- grid.arrange(arrangeGrob(s.m + theme(legend.position="none"), 
                                s.s + theme(legend.position="none"),
                                b.m + theme(legend.position="none"), 
                                b.s + theme(legend.position="none"),
                                n.m + theme(legend.position="none"), 
                                n.s + theme(legend.position="none"),
                                nrow=3,ncol=2),
                    mylegend, nrow=2,heights=c(10, 1))
ggsave(file=file.path(resultsDir,"SuppFig2.est.params.mle.datasets.120519.pdf"), all, width = 14, height = 20, units = "cm")
