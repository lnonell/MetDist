workingDir<-"D:/Doctorat/Simplex/Data"
setwd(workingDir)
resultsDir <- file.path("./Summary")

#data sets paths
GSE50660_data<-"D:/Doctorat/Simplex/Data/GSE50660_Smoking"
GSE116339_data<-"D:/Doctorat/Simplex/Data/GSE116339_PBB"
RRBS216_data<-"D:/Doctorat/Simplex/Data/RRBS_216_Tissue_CL"
RRBS188_data<-"D:/Doctorat/Simplex/Data/RRBS_188"
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