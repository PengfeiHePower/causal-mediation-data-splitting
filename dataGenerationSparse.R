# generate data with different n and p.

#setwd('Documents/stat research/simulation')
### simulation of DS
#### the sample size n=500, the number of mediators p=500, active mediators p1=50, number of pre-treatment m=20
#### m=gamma_0+gamma_1 t_i+gamma_2 x+ epsilon_m
#### y=beta_0 + beta_1 t_i + beta_2 m + beta_3 x+epsilon

##### Sample t, x
# treat: sample from bernoulli with probability p=0.5, n by 1
# pretreat: sample from N(0, I), n by m

##### Model mediator
# M=gamma.0+gamma.1*treat+gamma.2*pretreat+epsilon.m, n by p
# gamma.0: sample from N(0, theta.gamma0), theta_1 fixed, p by 1
# gamma.1: sample from N(0, theta.gamma1), p by 1
# gamma.2: each column sample from N(0, Sigma.g2), p by m

##### Model outcome
# y=beta.0 + beta.1*treat + beta.2*M + beta.3*pretreat+epsilon.y
# first sample p1 out of p as S1
# beta.0: sample from N(0, theta.beta2), 1
# beta.1: sample from N(0, theta.beta3), 1
# beta.2: for i in S_1, sample from N(0, \delta\sqrt{p/n}), other 0, p by 1
# beta.3: sample from N(0, Sigma.beta3), m by 1

# Most important variable:
# observation: treat, pretreat
# coefficient: gamma.0, gamma.1, gamma.2; beta.0, beta.1, beta.2, beta.3
# error: epsilon.y



######### data generating
#set.seed(1)
library(MASS)
library(stats)

require("getopt", quietly=TRUE)

spec = matrix(c(
    "Sample", "n", 1, "integer",
    "Mediator", "p", 1, "integer"
), byrow=TRUE, ncol=4)

opt = getopt(spec);

p=opt$Mediator
n=opt$Sample
p1=50
m = 20
cat('Sample size:', n, '\n')
cat('Mediator size:', p, '\n')
cat('Active mediator size:', p1, '\n')
cat('Pretreat size:', m, '\n')

Sigma.pretreat = 2 * diag(m)
mu.pretreat = rep(0,m)
# design matrix
treat = rbinom(n, 1, 0.5)
treat = matrix(treat, ncol = 1) # n by 1
pretreat = mvrnorm(n, mu.pretreat, Sigma.pretreat) #X[1,] select a column, n by m

######### Model mediator
theta.gamma0 = 1
theta.gamma1 = 1
mu.gamma2 = rep(0, m)
Sigma.gamma2 = 3 * diag(m)

intercept.M = matrix(rep(1, n), ncol = 1, nrow = n)
gamma.0 = matrix(nrow = p, ncol = 1) # p by 1
gamma.1 = matrix(nrow = p, ncol = 1) # p by 1
gamma.0[,1] = rnorm(p,0,theta.gamma0)
gamma.1[,1] = rnorm(p,0,theta.gamma1)
gamma.2 = mvrnorm(p, mu.gamma2, Sigma.gamma2) # p by m

mu.M = rep(0, p)
Sigma.M = diag(p)
epsilon.M = mvrnorm(n, mu.M, Sigma.M)

M = intercept.M %*% t(gamma.0) + treat %*% t(gamma.1) + pretreat %*% t(gamma.2) + epsilon.M



######### Model outcome
### select active mediator
S1 = sample(1:p, p1)
S0 = setdiff(1:p, S1)
delta = 5 #signal level
beta.2 = matrix(rep(0, p), ncol = 1, nrow = p) # p by 1
for(i in S1){
  beta.2[i] = rnorm(1, 0, delta*sqrt(log(p)/n))
}
theta.beta0 = 1
theta.beta1 = 1
mu.beta3 = rep(0, m)
Sigma.beta3 = diag(m)
intercept.Y = matrix(rep(1, n), ncol = 1) # n by 1
beta.0 = rnorm(1, 0, theta.beta0) # 1 by 1
beta.1 = rnorm(1, 0, theta.beta1) # 1 by 1
beta.3 = mvrnorm(1, mu.beta3, Sigma.beta3) # m by 1
epsilon.y = rnorm(n, 0, 1) # n by 1
epsilon.Y = matrix(epsilon.y, ncol = 1)

Y = intercept.Y * beta.0 + treat * beta.1 + M %*% beta.2 + pretreat %*% beta.3 + epsilon.y


## save data
save.image(file = paste('data/dataSparse_n', as.character(n), '_p', as.character(p), '.RData', sep=''))