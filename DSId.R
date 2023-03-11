##### data structure
##### Model mediator
# M=gamma.0+gamma.1*treat+gamma.2*pretreat+epsilon.m, n by p

##### Model outcome
# y=beta.0 + beta.1*treat + beta.2*M + beta.3*pretreat+epsilon.y

# Most important variable:
# observation: treat, pretreat
# coefficient: gamma.0, gamma.1, gamma.2; beta.0, beta.1, beta.2, beta.3
# error: epsilon.y

setwd('Documents/stat research/simulation')
load(file = 'data/meanModel.RData')
library(MASS)
library(stats)
library(glmnet)

## generate different covariance matrix for Mediator error
mu.M = rep(0, p)
Sigma.M.Id = diag(p)
epsilon.M.Id = mvrnorm(n, mu.M, Sigma.M.Id)

response = function(epsilon.M){
  M = intercept.M %*% t(gamma.0) + treat %*% t(gamma.1) + pretreat %*% t(gamma.2) + epsilon.M
  Y = intercept.Y * beta.0 + treat * beta.1 + M %*% beta.2 + pretreat %*% beta.3 + epsilon.Y
  out = list(M=M, Y=Y)
  return(out)
}

f = function(u,v){
  return(u + v)
}

fdr = function(t, M){
  up = sum(M<(-1*t))
  down = max(1, sum(M>t))
  return(up/down)
}

response.Id = response(epsilon.M.Id)

fdr.DS.Id = c()
power.DS.Id = c()
######## DS ############################
## split data into two equal parts
D = c(1:n)

for(reps in 1:50){
D1 = sample(D, size = n/2)
D2 = setdiff(D, D1)

t1 = treat[D1]
t2 = treat[D2]
t1.m = matrix(t1, ncol=1)
t2.m = matrix(t2, ncol=1)
inter.1 = matrix(rep(1,n/2), ncol = 1) #intercept
inter.full = matrix(rep(1,n), ncol = 1)
X1 = pretreat[D1,]
X2 = pretreat[D2,]
M1 = response.Id$M[D1,]
M2 = response.Id$M[D2,]
Y1 = response.Id$Y[D1,]
Y2 = response.Id$Y[D2,]

# Phase 1: Lasso on outcome model Y
X.full1 = cbind(inter.1, t1.m, X1, M1) #full design matrix
p.fac1 = rep(1, dim(X.full1)[2])
p.fac1[1:22] = 0

cv.lasso.1 = cv.glmnet(X.full1, Y1, penalty.factor = p.fac1)
coef.lasso1.min = coef(cv.lasso.1$glmnet.fit, s = cv.lasso.1$lambda.min, exact = F)

ind.lasso.raw = which(coef.lasso1.min!=0)
# estimation of beta1 we need
beta1.e1 = coef.lasso1.min[ind.lasso.raw[ind.lasso.raw>23]]

ind.lasso = ind.lasso.raw[ind.lasso.raw>23]
ind.lasso = ind.lasso - 23#index of active mediators, p0 mediators are selected
p1.hat = length(ind.lasso)

# Phase2: Do regression on D2 and ind.lasso for gamma and beta
# Mediator model
M.S1 = response.Id$M[, ind.lasso]
X.M = cbind(inter.full, pretreat, treat)
model.M2 = lm(M.S1~X.M)
gamma.e = solve(t(X.M) %*% X.M) %*% t(X.M) %*% M.S1

# estimation of gamma1 we need
gamma.e1 = gamma.e[dim(gamma.e)[1],]

# Outcome model
M2.S1 = M2[, ind.lasso]
X.Y2 = cbind(inter.1, t2.m, X2, M2.S1)
model.Y2 = lm(Y2~X.Y2)
beta2.e = solve(t(X.Y2) %*% X.Y2) %*% t(X.Y2) %*% Y2
# estimation of beta1 we need
beta2.e1 = beta2.e[(23:length(beta2.e))]

# generate mirror statistics
min.1 = c()
min.2 = c()
for(i in 1:length(gamma.e1)){
  min.1 = c(min.1, sign(gamma.e1[i] * beta1.e1[i]) * min(abs(gamma.e1[i]), abs(beta1.e1[i])))
  min.2 = c(min.2, sign(gamma.e1[i] * beta2.e1[i]) * min(abs(gamma.e1[i]), abs(beta2.e1[i])))
}


Mirror = c()
for(i in 1:length(min.1)){
  Mirror = c(Mirror, sign(min.1[i]*min.2[i])*f(abs(min.1[i]), abs(min.2[i])))
}


Mirror.abs = abs(Mirror)
fdr.M = c()
for(i in Mirror.abs){
  fdr.M = c(fdr.M, fdr(i, Mirror))
}
cutoffs = Mirror.abs[fdr.M<=0.05]
cutoff = min(cutoffs)

S1.final = ind.lasso[which(Mirror > cutoff)]#final set of selected values.
S0.final = setdiff(1:p, S1.final)

# evaluate
fdr.eval = length(intersect(S0,S1.final))/length(S1.final) # 0.0421
power.eval = length(intersect(S1, S1.final))/length(S1.final) # 0.85

fdr.DS.Id = c(fdr.DS.Id, fdr.eval)
power.DS.Id = c(power.DS.Id, power.eval)
}

save(fdr.DS.Id, power.DS.Id, file = 'results/DS_Id.RData')

