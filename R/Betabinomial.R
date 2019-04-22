#Beta binomial, a veure si m'aclaro

#explicat aquí:https://math.stackexchange.com/questions/185352/finding-alpha-and-beta-of-beta-binomial-model-via-method-of-moments


#ara l'exemple de la wiki 
#de la wiki https://en.wikipedia.org/wiki/Beta-binomial_distribution
M <- 0:12
F <-c(3,24,104,286,670,1033,1343,1112,829,478,181,45,7)
N=max(M)


#fitdistr(M, densfun=dbetabinom.ab, start=list(shape1=alpha,shape2=beta), size=12) #del paquet MASS
fitdist(M, distr=dbetabinom.ab, start=list(size=12,shape1=alpha,shape2=beta),lower=c(0,0)) #del paquet fitdistrplus
fitdist(p, distr=dbeta, start=list(shape1=-0.06228748,shape2=-0.3169064),lower=c(0,0))


m1=sum(M*F)/sum(F)
#6.23
m2=sum((M^2*F))/sum(F)
#42.3094
n=12

alpha <-  (n*m1-m2)/(n*(m2/m1-m1-1)+m1)
beta <-  ((n-m1)*(n-(m2/m1)))/(n*((m2/m1)-m1-1)+m1)
alpha
#34.13502
beta
#31.60849

#exemple de la wiki!!!

# s1=34.1350
# s2=31.6085
# 
# mu=s1/(s1+s2)
# disp= 1/(s1+s2)
# 
# mitjana=mu*n #6.23 ara si
# variancia= n*mu*(1-mu)
# 


n=100
# prob = 0.5
# 
# a <-rbetabinom(1000, n, prob = prob)
# hist(a)
# 
# #intento estimar params
# fitdistr(a, densfun = dbetabinom,start=list(prob=0.5))
# 
# fitdist(a, distr=dbetabinom,start=list(size=n,prob=0.5))
# fitdist(a, distr=dbetabinom.ab)
# fitdist(a, distr=dbetabinom.ab, start=list(size=n, shape1=1,shape2=2))
# fitdist(a, distr=dsimplex, start=list(mu=0.5,sig=2),optim.method="Nelder-Mead")
# fitdist(vector, densfun = "beta",  start = list(shape1 = 1, shape2 = 10))
a <- rnorm(100)
fitdist(a, distr=dnorm)

N <- 9; xx <- 0:N; s1 <- 2; s2 <- 3
a <- rbetabinom.ab(100, size=N, shape1=s1, shape2=s2)
a <-c(9,7,5,3,2,6,2,4,2,3,3,4,2,4,3,5,1,1,1,3,2,5,7,5,1,3,6,1,4,1,1,7,2,3,6,3,5,6,1,6,5,6,1,2,4,7,
2,1,1,7,5,2,2,2,5,4,8,1,0,1,1,2,7,5,1,3,7,2,4,7,2,6,9,1,5,5,5,9,0,1,5,3,2,4,5,2,7,4,3,5,5,2,
1,1,1,2,0,5,3,4)
hist(a)
fitdist(a, distr=dbetabinom.ab, start=list(size=N, shape1=1,shape2=1),lower=c(0,0))
#dona igual si poso shapes 2 i 3

#a partir d'aqui intento estimar, de la wiki

n=N
#en aquest cas no sabem de quin total, podem assumir que és 1 o el max, que es 9, dona el mateix!!
F=rep(1,length(a))
F=rep(N,length(a))
M=a
#llavors si que 
m1=sum(M*F)/sum(F)
#3.57
m2=sum((M^2*F))/sum(F)
#17.85
n<-N
alpha <-  (n*m1-m2)/(n*(m2/m1-m1-1)+m1)
beta <-  ((n-m1)*(n-(m2/m1)))/(n*((m2/m1)-m1-1)+m1)

alpha
#[1] 1.919355
beta
#[1] 2.919355
#si!!! corresponen a s1 i s2

#miro si puc reproduir grafics wiki
N <- 10; xx <- 0:N; s1 <- 0.2; s2 <- 0.25
a <- rbetabinom.ab(100, size=N, shape1=s1, shape2=s2)
hist(a)

fitdist(a, distr=dbetabinom.ab, start=list(size=N, shape1=1,shape2=1),lower=c(0,0))
#           estimate Std. Error
# size   10.0000000         NA
# shape1  0.3031255         NA
# shape2  0.2028078         NA

#si!!!