# Generate sythetic data with Identity matrix
## Must return: treat, pretreat, M, Y, S1, S0

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

load(file = 'data/meanModel.RData')
mu.M = rep(0, p)
Sigma.M = diag(p)
epsilon.M = mvrnorm(n, mu.M, Sigma.M)

M = intercept.M %*% t(gamma.0) + treat %*% t(gamma.1) + pretreat %*% t(gamma.2) + epsilon.M
Y = intercept.Y * beta.0 + treat * beta.1 + M %*% beta.2 + pretreat %*% beta.3 + epsilon.Y
# save data
save.image(file = 'data/dataID.RData')