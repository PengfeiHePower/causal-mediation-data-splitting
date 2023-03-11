toepW.v = c(2,5,10,20,30,50,100)
toepR.v = c(0.6,0.5,0.4,0.3,0.2,0.1)
load(file = 'data/meanModel.RData')

cond.toep = matrix(rep(0,length(toepW.v)*length(toepR.v)), nrow=length(toepW.v))
for(i in 1:length(toepW.v)){
  toepW = toepW.v[i]
  for(j in 1:length(toepR.v)){
    toepR = toepR.v[j]
    row1 = rep(0,500)
    row1[1] = 1
    exps = seq(1,toepW)
    row1[2:(toepW+1)] = toepR^exps
    Sigma.M = toeplitz(row1)
    cond.toep[i,j] = kappa(Sigma.M[S1,S1])
  }
}


fdr.toep.100 = matrix(rep(0,length(toepW.v)*length(toepR.v)), nrow=length(toepW.v))
power.toep.100 = matrix(rep(0,length(toepW.v)*length(toepR.v)), nrow=length(toepW.v))

for(i in 1:length(toepW.v)){
  for(j in 1:length(toepR.v)){
    fdr = read.csv(file=paste('results/fdr/covariance/toep/toep_100_',toepW.v[i],toepR.v[j]*10,'.csv',sep=''))
    power = read.csv(file=paste('results/power/covariance/toep/toep_100_',toepW.v[i],toepR.v[j]*10,'.csv',sep=''))
    fdr.toep.100[i,j] = mean(fdr$x[fdr$x!=0])
    power.toep.100[i,j] = mean(power$x)
  }
}

fdr.ds.toep = matrix(rep(0,length(toepW.v)*length(toepR.v)), nrow=length(toepW.v))
power.ds.toep = matrix(rep(0,length(toepW.v)*length(toepR.v)), nrow=length(toepW.v))

for(i in 1:length(toepW.v)){
  for(j in 1:length(toepR.v)){
    fdr = read.csv(file=paste('results/fdr/covariance/toep/DStoep',toepW.v[i],toepR.v[j]*10,'.csv',sep=''))
    power = read.csv(file=paste('results/power/covariance/toep/DStoep',toepW.v[i],toepR.v[j]*10,'.csv',sep=''))
    fdr.ds.toep[i,j] = mean(fdr$x[fdr$x!=0])
    power.ds.toep[i,j] = mean(power$x)
  }
}

cond.toep = matrix(rep(0,length(toepW.v)*length(toepR.v)), nrow=length(toepW.v))
for(i in 1:length(toepW.v)){
  for(j in 1:length(toepR.v)){
    row1 = rep(0,500)
    row1[1] = 1
    exps = seq(1,toepW.v[i])
    row1[2:(toepW.v[i]+1)] = toepR.v[j]^exps
    Sigma.M = toeplitz(row1)
    cond.toep[i,j] = kappa(Sigma.M)
  }
}

fdr.toep.1 = matrix(rep(0,length(toepW.v)*length(toepR.v)), nrow=length(toepW.v))
power.toep.1 = matrix(rep(0,length(toepW.v)*length(toepR.v)), nrow=length(toepW.v))

for(i in 1:length(toepW.v)){
  for(j in 1:length(toepR.v)){
    fdr = read.csv(file=paste('results/fdr/covariance/toep/toep_1_',toepW.v[i],toepR.v[j]*10,'.csv',sep=''))
    power = read.csv(file=paste('results/power/covariance/toep/toep_1_',toepW.v[i],toepR.v[j]*10,'.csv',sep=''))
    fdr.toep.1[i,j] = mean(fdr$x[fdr$x!=0])
    power.toep.1[i,j] = mean(power$x)
  }
}


