##### data structure
# t: treatment, n by 1
# X: pre-treatment covariates, n by m

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

# Only change Sigma.em starting from identity, other parameters remain same.
library(MASS)
library(stats)
library(glmnet)

#set arguments
require("getopt", quietly=TRUE)

spec = matrix(c(
    "Filename", "f", 1, "character",
    "Fdr", "q", 1, "numeric",
    "Savename", "s", 1, "character",
    "MdsTime", "m", 1, "integer"
), byrow=TRUE, ncol=4)

opt = getopt(spec);
fdr_level = opt$Fdr
savedfile = paste('results/', opt$Savename, sep='')
mds.times = opt$MdsTime
cat("File loaded:", paste('data/', opt$Filename, sep=''), "\n")
cat("Fdr level:", fdr_level, "\n")
cat("Mds Time:", mds.times, "\n")

# load data
load(file = paste('data/', opt$Filename, sep=''))

# important functions
fdr = function(t, M){
  up = sum(M<(-1*t))
  down = max(1, sum(M>t))
  return(up/down)
}

f = function(u,v){
  return(u * v)
}

find.cutoff = function(q, score){
  cutoff = 0
  sum.score = 0
  for(i in 1:length(score)){
    sum.score = sum.score+score[i]
    if(sum.score<=q){
      cutoff = score[i]
    }
  }
  return(cutoff)
}



fdr.MDS.Id = c()
power.MDS.Id = c()

D = c(1:n)

for(reps in 1:50){
med.score = rep(0, p)

for(i in 1:mds.times){
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
  M1 = M[D1,]
  M2 = M[D2,]
  Y1 = Y[D1,]
  Y2 = Y[D2,]
  
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
  M.S1 = M[, ind.lasso]
  X.M = cbind(inter.full, pretreat, treat)
  gamma.e = solve(t(X.M) %*% X.M) %*% t(X.M) %*% M.S1
  
  # estimation of gamma1 we need
  gamma.e1 = gamma.e[dim(gamma.e)[1],]
  
  # Outcome model
  M2.S1 = M2[, ind.lasso]
  X.Y2 = cbind(inter.1, t2.m, X2, M2.S1)
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
  cutoffs = Mirror.abs[fdr.M<=fdr_level]
  cutoff = min(cutoffs)
  
  S1.ds = ind.lasso[which(Mirror > cutoff)]
  for(i in 1:p){
    if(i %in% S1.ds){
      med.score[i] = med.score[i] + 1/max(1, length(S1.ds))
    }
  }
}

med.score = med.score/mds.times

sort.med.score = sort(med.score)
cutoff = find.cutoff(fdr_level, sort.med.score)
S1.mds = which(med.score>cutoff)
S0.mds = setdiff(1:p, S1.mds)

fdr.mds = length(intersect(S0,S1.mds))/length(S1.mds) #0.033
power.mds = length(intersect(S0, S0.mds))/length(S0.mds) #0.955

fdr.MDS.Id = c(fdr.MDS.Id, fdr.mds)
power.MDS.Id = c(power.MDS.Id, power.mds)
}

cat('FDR:', mean(fdr.MDS.Id))
cat('Power:', mean(power.MDS.Id))

# save(fdr.MDS.Id, power.MDS.Id, file = savedfile)
# cat("Save to:", savedfile, "\n")