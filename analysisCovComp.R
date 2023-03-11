compR.v = c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)


cond.comp=c()
for(i in 1:length(compR.v)){
  p=500
  row1 = rep(0,p)
  row1[1] = 1
  row1[2:p] = compR.v[i]
  Sigma.M = toeplitz(row1)
  cond.comp = c(cond.comp, kappa(Sigma.M))
}



fdr.ds.comp = c()
power.ds.comp = c()
for(i in 1:length(compR.v)){
  fdr = read.csv(paste('results/fdr/covariance/compound/DScomp', compR.v[i]*10, '.csv', sep=''))
  power = read.csv(paste('results/power/covariance/compound/DScomp', compR.v[i]*10, '.csv', sep=''))
  fdr.ds.comp = c(fdr.ds.comp, mean(fdr$x))
  power.ds.comp = c(power.ds.comp, mean(power$x))
}

plot(compR.v, fdr.ds.comp, pch=19,ylim=c(0,0.4), main='DS fdr for compound stmmetric covariance', xlab='correlation', ylab='fdr')
abline(0.05,0,col='red')
plot(compR.v, power.ds.comp, pch=19, main='DS power for compound stmmetric covariance', xlab='correlation', ylab='power')

plot(cond.comp, fdr.ds.comp, pch=19, ylim=c(0,0.4), main='DS fdr for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'fdr')
abline(0.05,0,col='red')
plot(cond.comp, power.ds.comp, pch=19, main='DS power for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'power')

fdr.mds.comp.100 = c()
power.mds.comp.100 = c()
for(i in 1:length(compR.v)){
  fdr = read.csv(paste('results/fdr/covariance/compound/comp_100_', compR.v[i]*10, '.csv', sep=''))
  power = read.csv(paste('results/power/covariance/compound/comp_100_', compR.v[i]*10, '.csv', sep=''))
  fdr.mds.comp.100 = c(fdr.mds.comp.100, mean(fdr$x))
  power.mds.comp.100 = c(power.mds.comp.100, mean(power$x))
}

plot(compR.v, fdr.mds.comp.100, pch=19,ylim=c(0,0.4), main='MDS(100) fdr for compound stmmetric covariance', xlab='correlation', ylab='fdr')
abline(0.05,0,col='red')
plot(compR.v, power.mds.comp.100, pch=19, main='MDS(100) power for compound stmmetric covariance', xlab='correlation', ylab='power')

plot(cond.comp, fdr.mds.comp.100, pch=19, ylim=c(0,0.4), main='MDS(100) fdr for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'fdr')
abline(0.05,0,col='red')
plot(cond.comp, power.mds.comp.100, pch=19, main='MDS(100) power for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'power')

fdr.mds.comp.200 = c()
power.mds.comp.200 = c()
for(i in 1:length(compR.v)){
  fdr = read.csv(paste('results/fdr/covariance/compound/comp_200_', compR.v[i]*10, '.csv', sep=''))
  power = read.csv(paste('results/power/covariance/compound/comp_200_', compR.v[i]*10, '.csv', sep=''))
  fdr.mds.comp.200 = c(fdr.mds.comp.200, mean(fdr$x))
  power.mds.comp.200 = c(power.mds.comp.200, mean(power$x))
}

plot(compR.v, fdr.mds.comp.200, pch=19,ylim=c(0,0.4), main='MDS(200) fdr for compound stmmetric covariance', xlab='correlation', ylab='fdr')
abline(0.05,0,col='red')
plot(compR.v, power.mds.comp.200, pch=19, main='MDS(200) power for compound stmmetric covariance', xlab='correlation', ylab='power')

plot(cond.comp, fdr.mds.comp.200, pch=19, ylim=c(0,0.4), main='MDS(200) fdr for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'fdr')
abline(0.05,0,col='red')
plot(cond.comp, power.mds.comp.200, pch=19, main='MDS(200) power for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'power')

fdr.mds.comp.1 = c()
power.mds.comp.1 = c()
for(i in 1:length(compR.v)){
  fdr = read.csv(paste('results/fdr/covariance/compound/comp_1_', compR.v[i]*10, '.csv', sep=''))
  power = read.csv(paste('results/power/covariance/compound/comp_1_', compR.v[i]*10, '.csv', sep=''))
  fdr.mds.comp.1 = c(fdr.mds.comp.1, mean(fdr$x))
  power.mds.comp.1 = c(power.mds.comp.1, mean(power$x))
}

plot(compR.v, fdr.mds.comp.1, pch=19,ylim=c(0,0.4), main='MDS(1) fdr for compound stmmetric covariance', xlab='correlation', ylab='fdr')
abline(0.05,0,col='red')
plot(compR.v, power.mds.comp.1, pch=19, main='MDS(1) power for compound stmmetric covariance', xlab='correlation', ylab='power')

plot(cond.comp, fdr.mds.comp.1, pch=19, ylim=c(0,0.4), main='MDS(1) fdr for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'fdr')
abline(0.05,0,col='red')
plot(cond.comp, power.mds.comp.1, pch=19, main='MDS(1) power for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'power')

fdr.ds.comp.1 = c()
power.ds.comp.1 = c()
for(i in 1:length(compR.v)){
  fdr = read.csv(paste('results/fdr/covariance/compound/comp_ds_', compR.v[i]*10, '.csv', sep=''))
  power = read.csv(paste('results/power/covariance/compound/comp_ds_', compR.v[i]*10, '.csv', sep=''))
  fdr.ds.comp.1 = c(fdr.ds.comp.1, mean(fdr$x))
  power.ds.comp.1 = c(power.ds.comp.1, mean(power$x))
}

plot(compR.v, fdr.ds.comp.1, pch=19,ylim=c(0,0.4), main='DS(1) fdr for compound stmmetric covariance', xlab='correlation', ylab='fdr')
abline(0.05,0,col='red')
plot(compR.v, power.ds.comp.1, pch=19, main='DS(1) power for compound stmmetric covariance', xlab='correlation', ylab='power')

plot(cond.comp, fdr.ds.comp.1, pch=19, ylim=c(0,0.4), main='DS(1) fdr for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'fdr')
abline(0.05,0,col='red')
plot(cond.comp, power.ds.comp.1, pch=19, main='DS(1) power for compound stmmetric covariance',
     xlab = 'condition number', ylab = 'power')
