#funció de la variança de la simplex

vu <-function (u) {u^3*(1-u)^3}

x <-runif(1000,0,1)
y <- vu(x)

plot(x,y) #sembla una normal!!

#prova
vu2 <-function (u) {u^2*(1-u)^2}
z <-vu2(x)

plot(x,z)
#semblant, la clau està en les cues...
