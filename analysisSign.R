sign.v = c(1,3,5,7,9,13)
fdr.sign = c()
power.sign = c()
for(i in 1:length(sign.v)){
  fdr = read.csv(paste('results/fdr/signal/sign', sign.v[i], '.csv', sep=''))
  power = read.csv(paste('results/power/signal/sign', sign.v[i], '.csv', sep=''))
  fdr.sign = c(fdr.sign, mean(fdr$x))
  power.sign = c(power.sign, mean(power$x))
}
plot(sign.v, fdr.sign, pch=19, ylim=c(0,0.15), 
     main='MDS fdr vs signal strength', xlab='signal strength', ylab='fdr')
abline(0.05,0,col='red')

plot(sign.v, power.sign, pch=19, 
     main='MDS power vs signal strength', xlab='signal strength', ylab='power')

fdr.sign.100 = c()
power.sign.100 = c()
for(i in 1:length(sign.v)){
  fdr = read.csv(paste('results/fdr/signal/sign100_', sign.v[i], '.csv', sep=''))
  power = read.csv(paste('results/power/signal/sign100_', sign.v[i], '.csv', sep=''))
  fdr.sign.100 = c(fdr.sign.100, mean(fdr$x))
  power.sign.100 = c(power.sign.100, mean(power$x))
}
plot(sign.v, fdr.sign.100, pch=19, ylim=c(0,0.15), 
     main='MDS(100) fdr vs signal strength', xlab='signal strength', ylab='fdr')
abline(0.05,0,col='red')

plot(sign.v, power.sign.100, pch=19, 
     main='MDS(100) power vs signal strength', xlab='signal strength', ylab='power')

fdr.ds.sign = c()
power.ds.sign = c()
for(i in 1:length(sign.v)){
  fdr = read.csv(paste('results/fdr/signal/DSsign', sign.v[i], '.csv', sep=''))
  power = read.csv(paste('results/power/signal/DSsign', sign.v[i], '.csv', sep=''))
  fdr.ds.sign = c(fdr.ds.sign, mean(fdr$x))
  power.ds.sign = c(power.ds.sign, mean(power$x))
}

plot(sign.v, fdr.ds.sign, pch=19, ylim=c(0,0.15), 
     main='DS fdr vs signal strength', xlab='signal strength', ylab='fdr')
abline(0.05,0,col='red')

plot(sign.v, power.ds.sign, pch=19, 
     main='DS power vs signal strength', xlab='signal strength', ylab='power')
