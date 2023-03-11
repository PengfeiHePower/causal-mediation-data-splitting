setwd('Documents/stat research/simulation')
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
library(MASS)
library(stats)
n = 500
p = 500
p1 = 50
m = 20

Sigma.pretreat = diag(m)
mu.pretreat = rep(0,m)
# design matrix
treat = rbinom(n, 1, 0.5)
treat = matrix(treat, ncol = 1) # n by 1
pretreat = mvrnorm(n, mu.pretreat, Sigma.pretreat) #X[1,] select a column, n by m

######### Model mediator
theta.gamma0 = 1
theta.gamma1 = 1
mu.gamma2 = rep(0, m)
Sigma.gamma2 = diag(m)

intercept.M = matrix(rep(1, n), ncol = 1)
gamma.0 = matrix(nrow = p, ncol = 1) # p by 1
gamma.1 = matrix(nrow = p, ncol = 1) # p by 1
gamma.0[,1] = rnorm(p,0,theta.gamma0)
gamma.1[,1] = rnorm(p,0,theta.gamma1)
gamma.2 = mvrnorm(p, mu.gamma2, Sigma.gamma2) # p by m


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
epsilon.Y = rnorm(n, 0, 1) # n by 1
epsilon.Y = matrix(epsilon.Y, ncol = 1)


## save data
save.image(file = 'data/meanModel.RData')

