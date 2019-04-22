#13/1/19 Resum de tot el trobat
#figures per l'article
#per a cada conjunt de dades he fet les mateixes anàlisis a part de l'anàlisi específica del data set

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


######################################################################
########################## Figure X: Cpg's de cada data set ##########
######################################################################
#com que surt un xurro, potser millor fer la dist global: FER pq es penja, objectes massa grans
#FER
######### GSE50660
load(file=file.path(GSE50660_data,"beta.filtered.RData"))
a <- qplot(as.vector(beta), geom="histogram") 
rm(beta)
######### GSE116339
load(file=file.path(GSE116339_data,"PM.filtered.RData"))
b <- qplot(as.vector(PM.f), geom="histogram") 
rm(beta)
######### RRBS216
load(file=file.path(RRBS216_data,"meth.betas.RData"))
c <- qplot(as.vector(rnb.meth), geom="histogram") #aquí hi ha missings
rm(rnb.meth)
######### WGBS81? hauria de ser el d

png(file=file.path(resultsDir,"Fig1.png"), width=480, height=480,res=100) #param res és l'important!!
grid.arrange(arrangeGrob(a + theme(legend.position="none"), #potser treure això de la legend??
                         b + theme(legend.position="none"),
                         c + theme(legend.position="none"),
                         d + theme(legend.position="none"),
                         nrow=2),
             mylegend, nrow=2,heights=c(10, 1))
dev.off()




########################################################################################
######### Figure 4: Resultats de les simulacions: Mesures d'avaluacio 100 100###########
########################################################################################
##plot del jaccard index per a cada conjunt de simulacions de cada data set.
#Tenim 100 i 10 de moment, serien 2 gràfics

#FUNCIó PER A FER-HO PER A CADA DATA SET
eval.plot <- function(sim.simplex=cpgs.simplex.models,
                      sim.beta=cpgs.beta.models,sim.normal=cpgs.normal.models,
                      tit="A"){
  #carregar source i paquet
  #pel recode, millor que car que de vegades costa d'instal·lar
  #ho he de trasposar tot per a poder fer el gràfic
  
  s.eval <- models.eval(models=sim.simplex, dml.r=100,alpha=0.05, adjust=TRUE)
  b.eval <- models.eval(models=sim.beta, dml.r=100,alpha=0.05, adjust=TRUE)
  n.eval <- models.eval(models=sim.normal, dml.r=100,alpha=0.05, adjust=TRUE)
  
  eval.all <- data.frame(t(s.eval), model.p=rownames(t(s.eval)),sim.model="simplex")
  eval.all <- rbind(eval.all,
                data.frame(t(b.eval), model.p=rownames(t(b.eval)),sim.model="beta"),
                data.frame(t(n.eval), model.p=rownames(t(n.eval)),sim.model="normal"))
  
  eval.all$model.p <- factor(recode(eval.all$model.p,  p.b="beta", p.binf="beta.inf",p.s="simplex",
                    p.sinf="simplex.inf", p.n="normal", p.l="logistic", p.q="quantile"),
                    levels=c("beta", "beta.inf","simplex","simplex.inf","normal","logistic","quantile"))
  #pdf? 
  p <- ggplot(data=eval.all, aes(x=model.p, y=jaccard)) +
    geom_bar(aes(fill = sim.model),stat = "identity",position = "dodge")+
    ylim(0,1) +
    theme(legend.position="bottom",legend.title = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=7,angle=45), axis.text.y = element_text(size=8),
          axis.title = element_text(size=7), legend.key.size = unit(0.8,"line"))+ggtitle(tit)
  p
  #potser caldrà orientar l'eix de les x
}

#carrego resultats de les simulacions
load(file=file.path(GSE50660_data,"simulated.models.p.5percdml.2000.100.100.RData"))
a <- eval.plot(tit="A")
a

load(file=file.path(GSE116339_data,"simulated.models.p.5percdml.2000.100.100.RData")) 
b <- eval.plot(tit="B")
b

load(file=file.path(RRBS216_data,"simulated.models.p.5percdml.2000.100.100.RData"))
c <- eval.plot(tit="C")
c

#dades WGBS81?
load(file=file.path(WGBS81_data,"simulated.models.p.5percdml.2000.100.100.RData"))
d <- eval.plot(tit="D")
d

mylegend<-g_legend(a) #fn g_lengend mes a munt

png(file=file.path(resultsDir,"Fig3.simul.100.100.jaccard.png"), width=480, height=480,res=100) #param res és l'important!!
ga <- grid.arrange(arrangeGrob(a + theme(legend.position="none"), #potser treure això de la legend??
                         b + theme(legend.position="none"),
                         c + theme(legend.position="none"),
                         d + theme(legend.position="none"), #si finalment afegeixo WGBS81 
                         nrow=2),
             mylegend, nrow=2,heights=c(10, 1))
dev.off()
ggsave(file=file.path(resultsDir,"Fig4.simul.100.100.jaccard.pdf"), ga, width = 14, height = 14, units = "cm")



#i ara què faig amb aquestes tables 
########################################################################################
######### Supp Figure 3: Resultats de les simulacions: Mesures d'avaluacio 10 10########
########################################################################################
#la mateixa fn per plotar: eval plot, però li he de passar uns altres objectes

load(file=file.path(GSE50660_data,"simulated.models.p.5percdml.2000.10.10.RData"))
a <- eval.plot(sim.simplex=cpgs.simplex.10.models,
               sim.beta=cpgs.beta.10.models,sim.normal=cpgs.normal.10.models,
               tit="A")
a

load(file=file.path(GSE116339_data,"simulated.models.p.5percdml.2000.10.10.RData")) 
b <- eval.plot(sim.simplex=cpgs.simplex.10.models,
               sim.beta=cpgs.beta.10.models,sim.normal=cpgs.normal.10.models,
               tit="B")
b

load(file=file.path(RRBS216_data,"simulated.models.p.5percdml.2000.10.10.RData"))
c <- eval.plot(sim.simplex=cpgs.simplex.10.models,
               sim.beta=cpgs.beta.10.models,sim.normal=cpgs.normal.10.models,
               tit="C")
c

#dades WGBS81?
load(file=file.path(WGBS81_data,"simulated.models.p.5percdml.2000.10.10.RData"))
d <- eval.plot(tit="D")
d

mylegend<-g_legend(a) #fn g_lengend mes a munt

png(file=file.path(resultsDir,"FigX.simul.10.10.jaccard.png"), width=480, height=480,res=100) #param res és l'important!!
ga <- grid.arrange(arrangeGrob(a + theme(legend.position="none"), #potser treure això de la legend??
                         b + theme(legend.position="none"),
                         c + theme(legend.position="none"),
                         d + theme(legend.position="none"), #si finalment afegeixo WGBS81 
                         nrow=2),
             mylegend, nrow=2,heights=c(10, 1))
dev.off()
ggsave(file=file.path(resultsDir,"SuppFig3.simul.10.10.jaccard.pdf"), ga, width = 14, height = 14, units = "cm")


########################################################################################
######### Figure 5?: Comparacions de cada dataset: Venn diagrams dels significatius #####
########################################################################################
#no tinc clar si potser ho hauria de fer sobre els resultats...potser millor, em sembla que 
# library(devtools)
# install_github("machalen/VennPlots")
library(VennPlots)

#no hi ha les dades guardades
# load(file=file.path(GSE50660_data,"simulated.models.p.5percdml.2000.100.100.RData"))
# 
# cpgs.simplex.models.adj.p <- apply(cpgs.simplex.models,2,p.adjust)
# head(cpgs.simplex.models.adj.p)
# 
# #faig venn de beta.inf, simplex.inf, normal, logistic and quantile
# l1 <- rownames(cpgs.simplex.models.adj.p[!is.na(cpgs.simplex.models.adj.p[,3]<0.05) & 
#                                            cpgs.simplex.models.adj.p[,1]<0.05,])
# l2 <- rownames(cpgs.simplex.models.adj.p[!is.na(cpgs.simplex.models.adj.p[,4]<0.05) & 
#                                            cpgs.simplex.models.adj.p[,1]<0.05,])
# l3 <- rownames(cpgs.simplex.models.adj.p[!is.na(cpgs.simplex.models.adj.p[,5]<0.05) & 
#                                            cpgs.simplex.models.adj.p[,1]<0.05,])
# l4 <- rownames(cpgs.simplex.models.adj.p[!is.na(cpgs.simplex.models.adj.p[,6]<0.05) & 
#                                            cpgs.simplex.models.adj.p[,1]<0.05,])
# l5 <- rownames(cpgs.simplex.models.adj.p[!is.na(cpgs.simplex.models.adj.p[,7]<0.05) & 
#                                            cpgs.simplex.models.adj.p[,1]<0.05,])
# l.names=c("simplex infl","beta infl","normal","logistic","quantile")
# 
# Venn5D(list.1=l1,list.2=l2,list.3=l3,list.4=l4,list.5=l5,listNames=l.names,
#          filename=file.path(resultsDir,"Venn.GSE50.100.100.pdf"), CatCex=0.8, CatDist=rep(0.1, 5))
  
Venn5D <- function(list.1=l1,list.2=l2,list.3=l3,list.4=l4,list.5=l5,listNames=l.names,
                   filename=file.path(resultsDir,"Venn.GSE50.100.100"),CatCex=0.8, CatDist=rep(0.1, 5)){
  
 # require(Vennerable)
  require(VennDiagram)
  require(colorfulVennPlot) #Per generar plots amb colors diferents
  require(RColorBrewer)
  
  cols <- c(brewer.pal(8,"Pastel1"), brewer.pal(8,"Pastel2"))  #Fixar els colors per als venn diagrams amb 4 condicions
  cols <- cols[c(8,2,3,15,5,6,7,1,9,10,11,12,13,14,4,16)]
   #remove na's just in case
   list.1 <- list.1[!is.na(list.1)]
   list.2 <- list.2[!is.na(list.2)]
   list.3 <- list.3[!is.na(list.3)]
   list.4 <- list.4[!is.na(list.4)]
   list.5 <- list.5[!is.na(list.5)]
  
   list.venn<-list(list.1,list.2,list.3,list.4,list.5)
   names(list.venn)<-listNames
   
   pdf(filename)
   venn.plot <- draw.quintuple.venn(
     area1 = length(list.1),
     area2 = length(list.2),
     area3 = length(list.3),
     area4 = length(list.4),
     area5 = length(list.5),
     n12 = length(intersect(list.1,list.2)),
     n13 = length(intersect(list.1,list.3)),
     n14 = length(intersect(list.1,list.4)),
     n15 = length(intersect(list.1,list.5)),
     n23 = length(intersect(list.2,list.3)),
     n24 = length(intersect(list.2,list.4)),
     n25 = length(intersect(list.2,list.5)),
     n34 = length(intersect(list.3,list.4)),
     n35 = length(intersect(list.3,list.5)),
     n45 = length(intersect(list.4,list.5)),
     n123 = length(intersect(list.1,intersect(list.2,list.3))),
     n124 = length(intersect(list.1,intersect(list.2,list.4))),
     n125 = length(intersect(list.1,intersect(list.2,list.5))),
     n134 = length(intersect(list.1,intersect(list.3,list.4))),
     n135 = length(intersect(list.1,intersect(list.3,list.5))),
     n145 = length(intersect(list.1,intersect(list.4,list.5))),
     n234 = length(intersect(list.2,intersect(list.3,list.4))),
     n235 = length(intersect(list.2,intersect(list.3,list.5))),
     n245 = length(intersect(list.2,intersect(list.4,list.5))),
     n345 = length(intersect(list.3,intersect(list.4,list.5))),
     n1234 = length(intersect(list.1,intersect(list.2,intersect(list.3,list.4)))),
     n1235 = length(intersect(list.1,intersect(list.2,intersect(list.3,list.5)))),
     n1245 = length(intersect(list.1,intersect(list.2,intersect(list.4,list.5)))),
     n1345 = length(intersect(list.1,intersect(list.3,intersect(list.4,list.5)))),
     n2345 = length(intersect(list.2,intersect(list.3,intersect(list.4,list.5)))),
     n12345 = length(intersect(list.1,intersect(list.2,intersect(list.3,intersect(list.4,list.5))))),
     category = listNames,
     fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
     cat.cex = CatCex,
     cat.dist = CatDist,
     margin = 0.05,
     cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
             1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
     ind = TRUE)
 
   dev.off()          
               
}
  

venn5D_betamodels <-function(beta.models.adj.p,file=file.path(resultsDir,"Venn.GSE50.smoking.pdf")){
  l1 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p$p.sinf<0.05) & beta.models.adj.p$p.sinf<0.05,])
  l2 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p$p.binf<0.05) & beta.models.adj.p$p.binf<0.05,])
  l3 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p$p.n<0.05) & beta.models.adj.p$p.n<0.05,])
  l4 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p$p.l<0.05) & beta.models.adj.p$p.l<0.05,])
  l5 <- rownames(beta.models.adj.p[!is.na(beta.models.adj.p$p.q<0.05) & beta.models.adj.p$p.q<0.05,])
  l.names=c("simplex infl","beta infl","normal","logistic","quantile")
  #fn definida més a munt
  Venn5D(list.1=l1,list.2=l2,list.3=l3,list.4=l4,list.5=l5,listNames=l.names,
         filename=file, CatCex=0.8, CatDist=rep(0.1, 5))
}

#GSE50660
load(file=file.path(GSE50660_data,file="betadata.models.RData"))
beta.models.adj.p <- as.data.frame(apply(beta.models,2,p.adjust))
venn5D_betamodels(beta.models.adj.p,file=file.path(resultsDir,"Venn.GSE50.smoking.pdf"))

#GSE116339
load(file=file.path(GSE116339_data,file="betadata.models.RData"))
beta.models.adj.p <- as.data.frame(apply(beta.models,2,p.adjust))
venn5D_betamodels(beta.models.adj.p,file=file.path(resultsDir,"Venn.GSE116339.sex.pdf"))

#RRBS216
load(file=file.path(RRBS216_data,file="rnb.meth2comp.lineage.models.RData")) 
#ccarrega objecte rnb.meth2comp.lineage.models
#li poso el mateix nom i així ho aprofito tot!
beta.models.adj.p <- as.data.frame(apply(rnb.meth2comp.lineage.models,2,p.adjust))
venn5D_betamodels(beta.models.adj.p,file=file.path(resultsDir,"Venn.RRBS216.lineage.pdf"))
#upps, 0 coincidències


