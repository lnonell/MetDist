#Simulations generations
#30/11/18

#To simulate data for testing purposes
#1. Generate simulated distribution for N CpGs under a beta, simplex and normal in a realistic scenario in two groups of size S1 and S2 (S1 and S2 >30)
#2. Generate simulated distribution for N CpGs under a beta, simplex and normal in a realistic scenario in two groups of size S1 and S2 (S1 and S2 >10 and <30)
#3. Generate simulated distribution for N CpGs under a beta, simplex and normal in a realistic scenario in two groups of size S1 and S2 (S1 and S2 <10 )

library(simplexreg)

#he de revisar rang a i b, això és el que vaig fer servir per a fer les simulacions inicials
a=seq(0.1,0.9,by=0.1)
b=seq(1,25)

#3. Ns grans
#3.1 per a cada N he d'escollir aleat una a i una b i generar una mostra de 100 individuus
#    assumeixo que hi haurà diferències quan provinguin de dues dist diferents
#    N1=100 CpGs
#3.2 agafo 200 individuus i genero una dist per a N2=100 CpGs, aquí no hi hauria d'haver diferències

#calculation of a & b having mean and var #https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance

#3.2
#Mostres
sam.n=200
#CpGs
cpg.n=100

set.seed(123)
d.sim <- matrix(NA, ncol=sam.n, nrow=cpg.n)
colnames(d.sim) <- paste0("s",1:sam.n)
rownames(d.sim) <- paste0("c",1:cpg.n)

a <- runif(sam.n,0,1) #aquest param és mu
b <- runif(sam.n,1,25) #sigma, dp fa sigma^2 
range(a) #ok
range(b) #ok

#mirar si hi ha alguna relació amb a i b??
#NO: hem de simular per cpg, per tant per fila i no per columna: aquesta matriu no té diferències
for (i in 1:cpg.n) d.sim[i,] <- rsimplex(200,mu=a[i],sig=b[i])
 
#ara hem de fer el mateix però per als 100 primers pacients una dist i pels segons una altra
d.sim.diff <- matrix(NA, ncol=sam.n, nrow=cpg.n)
colnames(d.sim.diff) <- paste0("s",1:sam.n)
rownames(d.sim.diff) <- paste0("c",1:cpg.n)

for (i in 1:cpg.n) {
  d.sim.diff[i,1:100] <- rsimplex(100,mu=a[i],sig=b[i]) #les primeres 100 mostres els mateixos params
  d.sim.diff[i,101:200] <- rsimplex(100,mu=a[100+i],sig=b[100+i]) #CHECK
}


#ja tinc els simplex, falten les betes i les normals

hist(rsimplex(10000,2,10),breaks=seq(0,1,by=0.001),freq=FALSE,main=paste("a=",aa," b=",bb,sep=""))

######### VELL : SIMULACIONS pel treball
hist(rsim(10000,0.5,16),breaks=seq(0,1,by=0.001),freq=FALSE)
pdf("Simplex simulations.020416.pdf")
for (aa in a){
  for (bb in b){
    hist(rsim(10000,aa,bb),breaks=seq(0,1,by=0.001),freq=FALSE,main=paste("a=",aa," b=",bb,sep=""))
  }
}
dev.off()

aa=0.4
bb=5
hist(rsim(10000,aa,bb),freq=FALSE,main=paste("a=",aa," b=",bb,sep=""))

vals=c(0.1,0.5,1,2,3,4,5,8,10)

pdf("beta simulations.pdf")
for (a in vals){
  for (b in vals){
    hist(rbeta(10000,a,b),breaks=seq(0,1,by=0.001),freq=FALSE,main=paste("a=",a," b=",b,sep=""))
  }
}
dev.off()
