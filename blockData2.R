######### data generating
#set.seed(1)
library(MASS)
library(stats)



  
  delta = 5
  p=500
  n=500
  p1 = 50
  m = 20
  print(delta)
  
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
  
  intercept.M = matrix(rep(1, n), ncol = 1, nrow = n)
  gamma.0 = matrix(nrow = p, ncol = 1) # p by 1
  gamma.1 = matrix(nrow = p, ncol = 1) # p by 1
  gamma.0[,1] = rnorm(p,0,theta.gamma0)
  gamma.1[,1] = rnorm(p,0,theta.gamma1)
  gamma.2 = mvrnorm(p, mu.gamma2, Sigma.gamma2) # p by m
  
  ## mediator block covariance matrix

  m2.1 = toeplitz(c(1,0.1,0.1))
  m2.2 = toeplitz(c(1,0.1,0.1,0.1))
  Sigma.M = diag(500)
  ind2.1 = seq(1,200,3)
  ind2.2 = seq(201,400,4)
  for(i in ind2.1){
    Sigma.M[i:(i+2), i:(i+2)] = m2.1
  }
  for(i in ind2.2){
    Sigma.M[i:(i+3), i:(i+3)] = m2.2
  }
  mu.M = rep(0, p)
  epsilon.M = mvrnorm(n, mu.M, Sigma.M)

  
  
  
  
  
  ######### Model outcome
  ### select active mediator
  S1 = sample(1:p, p1)
  S0 = setdiff(1:p, S1) #signal level
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
  epsilon.y = matrix(epsilon.y, ncol = 1)
  
  M = intercept.M %*% t(gamma.0) + treat %*% t(gamma.1) + pretreat %*% t(gamma.2) + epsilon.M
  Y = intercept.Y * beta.0 + treat * beta.1 + M %*% beta.2 + pretreat %*% beta.3 + epsilon.y
  
  
  ## save data
  save.image(file = 'data/dataBlock2.RData')


