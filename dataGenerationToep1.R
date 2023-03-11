setwd('Documents/stat research/simulation') # if needed
### simulation of DS
#### the sample size n=500, the number of mediators p=500, active mediators p1=50, number of pre-treatment m=20
#### m=gamma_0+gamma_1 t_i+gamma_2 x+ epsilon_m
#### y=beta_0 + beta_1 t_i + beta_2 m + beta_3 x+epsilon

##### Sample t, x
# t_i: sample from bernoulli with probability p=0.5
# x_i: sample from N(0, s), s is a Toeplitz matrix, let's fix it.

##### Model mediator
# m=gamma_0+gamma_1 t_i+gamma_2 x+ epsilon_m, n by p
# gamma_0: sample from N(0, \theta_0), theta_1 fixed, p by 1
# gamma_1: sample from N(0, \theta_1), p by 1
# gamma_2: each column sample from N(0, Sigma.g2), p by m
# epsilon_m: sample from (0, Sigma.em), Sigma.em have different choice starting from I, n by p

##### Model outcome
# y=beta_0 + beta_1 t_i + beta_2 m + beta_3 x+epsilon_y
# first sample 50 out of 500 as S_1
# beta_0: sample from N(0, \theta_2), 1
# beta_1: sample from N(0, \theta_3), 1
# beta_2: for i in S_1, sample from N(0, \delta\sqrt{p/n}), other 0, p by 1
# beta_3: sample from N(0, Sigma.b3), m by 1
# epsilon: sample from N(0, I), n by 1

######### data generating
#set.seed(1)
library(MASS)
library(stats)
n = 500
p = 500
p1 = 50
m = 20
#Sigma.x = toeplitz(c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,0,0,0,0,0,0,0,0,0))
Sigma.x = diag(m)
mu.x = rep(0,20)
# design matrix
t = rbinom(n, 1, 0.5)
t = matrix(t, ncol = 1, nrow = n) # n by 1
X = mvrnorm(n, mu.x, Sigma.x) #X[1,] select a column, n by m

######### Model mediator
theta0 = 5
theta1 = 8
mu.g2 = rep(0, m)
#Sigma.g2 = toeplitz(c(1,0.5,0.4,0.3,0.2,0.1,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0))
Sigma.g2 = diag(m)
mu.em = rep(0,p)
##################################################
bandit = c(1, rep(0.1, 10), rep(0, 500-11))
#Sigma.em = toeplitz(bandit) #identity matrix for mediator
Sigma.em = diag(p)
##################################################
intercept.M = matrix(rep(1, n), ncol = 1, nrow = n)
gamma.0 = matrix(nrow = p, ncol = 1) # p by 1
gamma.1 = matrix(nrow = p, ncol = 1) # p by 1
gamma.0[,1] = rnorm(p,0,theta0)
gamma.1[,1] = rnorm(p,0,theta1)
gamma.2 = mvrnorm(p, mu.g2, Sigma.g2) # p by m
epsilon.m = mvrnorm(n, mu.em, Sigma.em) # n by p

M = intercept.M %*% t(gamma.0) + t %*% t(gamma.1) + X %*% t(gamma.2) + epsilon.m



######### Model outcome
### select active mediator
S1 = sample(1:p, p1)
delta = 5
beta.2 = matrix(rep(0, p), ncol = 1, nrow = p) # p by 1
for(i in S1){
  beta.2[i] = rnorm(1, 0, delta*sqrt(log(p)/n))
}
theta2 = 5
theta3 = 10
mu.b3 = rep(0, m)
#Sigma.b3 = toeplitz(c(1,0.5,0.4,0.3,0.2,0.1,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0))
Sigma.b3 = diag(m)
intercept.Y = matrix(rep(1, n), ncol = 1) # n by 1
beta.0 = rnorm(1, 0, theta2) # 1 by 1
beta.1 = rnorm(1, 0, theta3) # 1 by 1
beta.3 = mvrnorm(1, mu.b3, Sigma.b3) # m by 1
epsilon.y = rnorm(n, 0, 1) # n by 1
epsilon.y = matrix(epsilon.y, ncol = 1)

Y = intercept.Y * beta.0 + t * beta.1 + M %*% beta.2 + X %*% beta.3 + epsilon.y

## save data
save.image(file = 'data/SimToeplitz1.RData')

