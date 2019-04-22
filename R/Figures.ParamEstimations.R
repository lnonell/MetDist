#figures and tables from parameter estimates results
#Suppfig 1: Correlations MOM MLE
#Suppfig 2: Densities for the different parameters
#Suppfig XX: betabi correls

#header amb tot
source("D:/Doctorat/Simplex/R/Figures.HeaderScript.R")

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
load(file=file.path(GSE50660_data,"params.est.all.120119.RData"))
est.params.gse50 <- est.params
rm(est.params)
load(file=file.path(GSE116339_data,"params.est.all.210219.RData"))
est.params.gse11 <- est.params
rm(est.params)
#load(file=file.path(RRBS216_data,"params.est.all.120119.RData"))
#load(file=file.path(RRBS216_data,"params.est.all.020319.RData"))
load(file=file.path(RRBS216_data,"params.est.all.090319.RData"))
est.params.rrbs <- rnb.est.params
rm(rnb.est.params)
colnames(est.params.rrbs)
#load(file=file.path(WGBS81_data,"params.est.all.RData"))
#load(file=file.path(WGBS81_data,"params.est.all.020319.RData"))
load(file=file.path(WGBS81_data,"params.est.all.090319.RData"))

est.params.wgbs <- rnb.est.params
rm(rnb.est.params)

apply(est.params.gse50, 2, range, na.rm=T)
#        s.mle.mu s.mle.sig  s.zoip.mu s.zoip.sig    b.mom.s1   b.mom.s2    b.mle.s1    b.mle.s2 binf.mle.s1 binf.mle.s2
# [1,] 0.03675831 0.2413557 0.03675583  0.2413354   0.2108814   0.188386   0.5140319   0.5199418    1.133821   0.8942635
# [2,] 0.96840259 6.0728849 0.96839279  6.0717148 163.6719375 163.913653 188.3485588 189.6722137   57.708526  33.7265609
#       n.mom.m   n.mom.sd    n.mle.m   n.mle.sd
# [1,] 0.03117263 0.03000025 0.03117263 0.02996791
# [2,] 0.97070543 0.32087425 0.97070543 0.32052829
apply(est.params.gse11, 2, range, na.rm=T)
#         s.mle.mu  s.mle.sig  s.zoip.mu s.zoip.sig    b.mom.s1    b.mom.s2    b.mle.s1    b.mle.s2 binf.mle.s1
# [1,] 0.03923843  0.2415039 0.03923187  0.2414835   0.2117655   0.1443674   0.5029186   0.3629927    0.413979
# [2,] 0.96207305 11.1077817 0.96207895 11.1122240 163.8193951 163.5935903 167.3989935 165.7837980   59.933499
#     binf.mle.s2    n.mom.m   n.mom.sd    n.mle.m   n.mle.sd
# [1,]   0.4155871 0.03858037 0.03000018 0.03858037 0.02997808
# [2,]  67.3928258 0.96306584 0.42159840 0.96306584 0.42128783
apply(est.params.rrbs, 2, range, na.rm=T)
#         s.mle.mu s.mle.sig s.zoip.mu s.zoip.sig  b.mom.s1  b.mom.s2    b.mle.s1     b.mle.s2  binf.mle.s1  binf.mle.s2
# [1,] 0.002024828  0.970169     0.001    0.00000 -0.248997 -0.248997   0.1288789    0.1288737 2.524122e-01 4.339461e-01
# [2,] 0.996060266 63.633742     0.999   63.15065       Inf       Inf 131.9927638 2242.5682692 2.050425e+08 4.198953e+07
#       n.mom.m  n.mom.sd      n.mle.m    n.mle.sd
# [1,]       0 0.0000000 0.0002666667 0.002293953
# [2,]       1 0.7071068 0.9972222222 0.500000000

apply(est.params.wgbs, 2, range, na.rm=T)
#       s.mle.mu s.mle.sig  s.zoip.mu s.zoip.sig  b.mom.s1  b.mom.s2  b.mle.s1  b.mle.s2  binf.mle.s1  binf.mle.s2
# [1,] 0.04960212  1.218596 0.04956221   1.218596 -0.248997 -0.248997 0.1289177 0.1289263 2.845620e-01 3.114245e-01
# [2,] 0.94967629 63.273728 0.95043659  63.150645  3.052946  3.054790 6.8885565 7.2263283 8.239343e+06 4.122168e+06
#         n.mom.m  n.mom.sd    n.mle.m  n.mle.sd
# [1,] 0.04166667 0.2000001 0.04166667 0.1416667
# [2,] 0.95833333 0.7071068 0.95833333 0.5000000

#sembla que estan a escales molt diferents
colors <- c("green4", "blue", "red2", "orange")

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
ggsave(file=file.path(resultsDir,"SuppFig1.correls.est.params.mle.datasets.png"), all, width = 16, height = 20, units = "cm")
ggsave(file=file.path(resultsDir,"SuppFig1.correls.est.params.mle.datasets.pdf"), all, width = 16, height = 20, units = "cm")


#############################################
#Supp Figure XXX: el mateix per betabinomial
#############################################

load(file=file.path(RRBS216_data,"betabin.params.est.RData"))
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
ggsave(file=file.path(resultsDir,"SuppFig1.correls.est.params.betbin.mle.datasets.png"), all.bb, width = 16, height = 20, units = "cm")




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
df3 <- data.frame(ds="RRBS-ENC",mu=est.params.rrbs$s.mle.mu,sig=est.params.rrbs$s.mle.sig)
df4 <- data.frame(ds="WGBS-BLUE",mu=est.params.wgbs$s.mle.mu,sig=est.params.wgbs$s.mle.sig)
df <- rbind(df1,df2,df3,df4)

s.m <- plot.dens(df, param=mu,ds=ds, log2 =F)
s.s <- plot.dens(df, param=sig,ds=ds, log2 =T)

#beta
df1 <- data.frame(ds="a450k-smoking",s1=est.params.gse50$b.mle.s1,s2=est.params.gse50$b.mle.s2)
df2 <- data.frame(ds="aEPIC-PBB",s1=est.params.gse11$b.mle.s1,s2=est.params.gse11$b.mle.s2)
df3 <- data.frame(ds="RRBS-ENC",s1=est.params.rrbs$b.mle.s1,s2=est.params.rrbs$b.mle.s2)
df4 <- data.frame(ds="WGBS-BLUE",s1=est.params.wgbs$b.mle.s1,s2=est.params.wgbs$b.mle.s2)
df <- rbind(df1,df2,df3,df4)

b.m <- plot.dens(df, param=s1,ds=ds, log2 =T)
b.s <- plot.dens(df, param=s2,ds=ds, log2 =T)

#normal
df1 <- data.frame(ds="a450k-smoking",mu=est.params.gse50$n.mle.m,sd=est.params.gse50$n.mle.sd)
df2 <- data.frame(ds="aEPIC-PBB",mu=est.params.gse11$n.mle.m,sd=est.params.gse11$n.mle.sd)
df3 <- data.frame(ds="RRBS-ENC",mu=est.params.rrbs$n.mle.m,sd=est.params.rrbs$n.mle.sd)
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
ggsave(file=file.path(resultsDir,"SuppFig2.est.params.mle.datasets.pdf"), all, width = 14, height = 20, units = "cm")
