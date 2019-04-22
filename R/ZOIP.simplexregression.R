library(ZOIP)


datas <- data.frame(yi = runif(1000,min=0,max=1))
range(datas$yi)
 0.0007627939 0.9985508111
mod <- RM.ZOIP(formula.mu = yi ~ 1, data = datas, family = "Original")
mod

yi <- data.frame(yi = rZOIP(n = 1000, mu = 0.6, sigma = 0.2,
                            p0 = 0.03, p1 = 0.05, family = "R-S")) #beta Rigby-Stasinopoulos
hist(yi$yi)
yi2 <- data.frame(yi = rZOIP(n = 1000, mu = 0.6, sigma = 0.2,
                            p0 = 0.03, p1 = 0.05, family = "R-S")) #beta Ferrari Cribari-Neto
hist(yi2$yi)
yi3 <- data.frame(yi = rZOIP(n = 1000, mu = 0.6, sigma = 0.2,
                             p0 = 0.03, p1 = 0.05, family = "Simplex")) #Simplex dist
hist(yi3$yi)


# La funcion RM.ZOIP estima los par´ametros de una distribucion ZOIP vıa maxima
# verosimilitud utilizando el optimizador deseado (nlminb,optim).
mod <- RM.ZOIP(formula.mu = yi ~ 1, data = yi, family = "Original")
mod #OK

mod <- RM.ZOIP(formula.mu = yi ~ 1, data = yi3, family = "Simplex")
mod #OK uses MLE 


n<-1000
x1<-stats::runif(n)
x2<-stats::runif(n)

datas <- data.frame(yi = runif(1000,min=0,max=1))
range(datas$yi)
mod <- RM.ZOIP(formula.mu = yi ~ 1, data = datas, family = "Original")

datas <- data.frame(yi = runif(1000,min=0,max=1),x1)

mod <- RM.ZOIP(formula.mu = yi ~ x1, data = datas, family = "Original", link=c("log","identity","identity","identity")) #no xuta, sembla que la nova versió conté més params
mod

#afegeixo una var condició a veure si puc fer una regressió
yi <- data.frame(yi = runif(1000,min=0,max=1),xi = c(rep(0,500),rep(1,500)))
mod <- RM.ZOIP(formula.mu = yi ~ xi, data = yi, family = "Original", link=c("log","logit"))
mod



datas<-data.frame(x1,x2)

mod <- RM.ZOIP(formula.mu = x1 ~ x2, data = datas, family = "Simplex",  link="logit")
mod

# 
# formula.mu=yi~x1
# formula.sigma=~x1+x2
# formula.p0=~1
# formula.p1=~x2

mod <- RM.ZOIP(formula.mu=yi~x1,formula.sigma=~x1+x2, formula.p0=~1,formula.p1=~x2, data = yi, 
               link=c('log','log','identity','logit'),
               family='Original')

mod

#beta
yi <- data.frame(yi = as.vector(RSet.betas[1:10,])) #no xuta
yi <- data.frame(yi = runif(1000,min=0,max=1))

mod.b <- RM.ZOIP(formula.mu = yi ~ 1, data= yi, family = "R-S",optimizer="optim")
mod.b


mod.b <- RM.ZOIP(formula.mu = yi ~ 1, data= yi, family = "R-S")
mod.b


mod.s <- RM.ZOIP(formula.mu = yi ~ 1, data= yi, family = "Simplex")
mod.s

y2i <- data.frame(yi = rZOIP(n = 1000, mu = 0.6, sigma = 0.2,
                            p0 = 0.03, p1 = 0.05, family = "Simplex"))
mod2 <- RM.ZOIP(formula.mu = yi ~ 1, data= y2i, family = "Simplex",optimizer="optim")
mod2

###############DADES REALS DE BETES #################
#Provo amb les dades de la presentació (D:\Doctorat\1s year\Projecte defence\Presentation\Presentation_LN_070616_rev150818.Rmd)
#Agafo les dades de minfiData, ho passo a data.frame

dades <- as.data.frame(dades)

#sembla que el nom de les vars no li agraden

colnames(dades) <- paste0("var",1:6)
head(dades)

mod2 <- RM.ZOIP(formula.mu = var1 ~ 1, data=dades, family = "Simplex",optimizer="optim")
#ara sí: 28/10/18
#error de sempre: Error in dZOIP(x = y, mu = mu, sigma = sigma, p0 = p0, p1 = p1, family = family,  : 
#p0+p1 must be between 0 and 1

summary(mod2)


#provo el default del param optimizer
mod2 <- RM.ZOIP(formula.mu = var1 ~ 1, data=dades, family = "Simplex",optimizer="nlminb")
#igual Error in dZOIP(x = y, mu = mu, sigma = sigma, p0 = p0, p1 = p1, family = family,  : 
#p0+p1 must be between 0 and 1

#hi ha Nas, faig subsetting
dades1 <- dades[1:250000,]
range(dades1$var1) # 0 1
mod2 <- RM.ZOIP(formula.mu = var1 ~ 1, data=dades1, family = "Simplex",optimizer="nlminb")
#Error in dZOIP(x = y, mu = mu, sigma = sigma, p0 = p0, p1 = p1, family = family,  : 
#p0+p1 must be between 0 and 1 

#ara provo d'eliminar els 0's i 1s
which(dades1$var1==0)
#15623 118110
which(dades1$var1==1)
#54697 113083 126388

dades2<-dades1[1:15000,]
mod2 <- RM.ZOIP(formula.mu = var1 ~ 1, data=dades2, family = "Simplex",optimizer="nlminb")
#doncs segueix igual...
#Error in dZOIP(x = y, mu = mu, sigma = sigma, p0 = p0, p1 = p1, family = family,  : 
#                 p0+p1 must be between 0 and 1 
range(dades2$var1) #0.002012441 0.995980875

#provo amb mixed effects
mod2 <- RMM.ZOIP(formula.mu = var1 ~ 1, data=dades2, family = "Simplex",optimizer="nlminb")

