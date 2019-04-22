#6/4/19
#Betabin regression
#several options
# aod 
# HRQoL
# bbmle

library(aod) #analysis of overdispersed data

data(orob2)
fm1 <- betabin(cbind(y, n - y) ~ seed, ~ 1, data = orob2)
fm2 <- betabin(cbind(y, n - y) ~ seed + root, ~ 1, data = orob2)
fm3 <- betabin(cbind(y, n - y) ~ seed * root, ~ 1, data = orob2)
# show the model
fm1; fm2; fm3

AIC(fm1, fm2, fm3)
summary(AIC(fm1, fm2, fm3), which = "AICc")
# Wald test for root effect
wald.test(b = coef(fm3), Sigma = vcov(fm3), Terms = 3:4)
# likelihood ratio test for root effect
anova(fm1, fm3)
# model predictions
New <- expand.grid(seed = levels(orob2$seed),
                   root = levels(orob2$root))
data.frame(New, predict(fm3, New, se = TRUE, type = "response"))
# Djallonke sheep data
data(dja)
betabin(cbind(y, n - y) ~ group, ~ 1, dja)
# heterogeneous phi
betabin(cbind(y, n - y) ~ group, ~ group, dja,
        control = list(maxit = 1000))
# phi fixed to zero in group TREAT
betabin(cbind(y, n - y) ~ group, ~ group, dja,
        fixpar = list(4, 0))
# glim without overdispersion
summary(glm(cbind(y, n - y) ~ group,
            family = binomial, data = dja))
# phi fixed to zero in both groups
betabin(cbind(y, n - y) ~ group, ~ group, dja,
        fixpar = list(c(3, 4), c(0, 0))) 


data(dja)
betabin(cbind(y, n - y) ~ group, ~ 1, dja)
# heterogeneous phi
betabin(cbind(y, n - y) ~ group, ~ group, dja,
        control = list(maxit = 1000))
# phi fixed to zero in group TREAT
betabin(cbind(y, n - y) ~ group, ~ group, dja,
        fixpar = list(4, 0))
# glim without overdispersion
summary(glm(cbind(y, n - y) ~ group,
            family = binomial, data = dja))
# phi fixed to zero in both groups
betabin(cbind(y, n - y) ~ group, ~ group, dja,
        fixpar = list(c(3, 4), c(0, 0))) 


#overdispersion param is 1/1+alpha+beta: 
#https://stats.stackexchange.com/questions/24795/types-of-dispersion-parameter-for-binomial-data

#generate data
set.seed(1)
alpha = 1
beta = 1
n = 40
x <- rbeta(1000, alpha, beta)
y <- qbinom(runif(1000), n, x)

#modeling
mb <- betabin(cbind(y,n-y)~1, random=~1, data=as.data.frame(list(y=y)))

#which shows that the used dispersion is the inverse of the sum 1+α+β with α and β the coefficients of the beta-distribution.
#dispersion param
mb@param[2]   
# phi.(Intercept) 
# 0.340159 
1/(alpha+beta+1)
# [1] 0.3333333

s.mb <- summary(mb)
s.mb@Phi["phi.(Intercept)","Pr(> z)"] #0


# https://rpubs.com/cakapourani/beta-binomial: test for overdispersion Tarone's z-statistic

library(VGAM)
library(dplyr)
library(data.table)
library(ggplot2)

set.seed(5)
M <- 1000
alt_hyp = null_hyp <- vector("numeric", length = M)
for (i in 1:M) {
  # Total number of cells
  C <- rbinom(1, 70, runif(1, min = 0.4, max = 0.7))
  # Total number of CpGs for each cell
  n <- rbinom(C, 80, runif(1, min = 0.3, max = 0.8))
  
  ##---
  # Compute Tarone's Z statistic for overdispersion Beta Binomial
  ##---
  # Generate synthetic data with \mu = 0.8 and \rho = 0.05
  m_alt <- VGAM::rbetabinom(length(n), n, prob = 0.8, rho = 0.05)
  p_hat = sum(m_alt) / sum(n)
  S = sum( (m_alt - n * p_hat)^2 / (p_hat * (1 - p_hat)) )
  Z_score = (S - sum(n)) / sqrt(2 * sum(n * (n - 1)))
  alt_hyp[i] <- Z_score 
  
  ##---
  # Compute Tarone's Z statistic for Binomial
  ##---
  # Generate synthetic data with \mu = 0.8 and \rho = 0
  m_null <- rbetabinom(length(n), n, prob = 0.8, rho = 0)
  # Compute Tarone's Z statistic
  p_hat = sum(m_null) / sum(n)
  S = sum( (m_null - n * p_hat)^2 / (p_hat * (1 - p_hat)) )
  Z_score = (S - sum(n)) / sqrt(2 * sum(n * (n - 1)))
  null_hyp[i] <- Z_score
}
# Melt object for ggplot2
dt <- list(Null = null_hyp, Alternative = alt_hyp) %>% 
  melt %>% as.data.table %>% 
  setnames(c("value", "L1"), c("value", "Hypothesis"))
# Create histogram plot
p <- ggplot(dt, aes(x = value, color = Hypothesis, fill = Hypothesis)) +
  geom_histogram(aes(y = ..density..), alpha = 0.3, 
                 position = "dodge", bins = 100) +
  geom_vline(xintercept = 0.2, linetype = "dashed", size = 0.7) + 
  stat_function(fun = dnorm, colour = "cornflowerblue", size = 1, 
                arg = list(mean = 0, sd = 1)) + 
  scale_x_continuous(limits = c(-3,11.45)) +
  labs(title = "Tarone's Z score distribution",x = "Z score",y = "Density") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") 
print(p)
