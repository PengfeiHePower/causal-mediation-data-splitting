ind = 1:4
fdr.mds.block.50 = c()
power.mds.block.50 = c()
for(i in ind){
  fdr = read.csv(paste('results/fdr/covariance/block/block', i, '.csv', sep=''))
  power = read.csv(paste('results/power/covariance/block/block', i, '.csv', sep=''))
  fdr.mds.block.50 = c(fdr.mds.block.50, mean(fdr$x[fdr$x!=0]))
  power.mds.block.50 = c(power.mds.block.50, mean(power$x))
}

fdr.ds.block = c()
power.ds.block = c()
for(i in ind){
  fdr = read.csv(paste('results/fdr/covariance/block/DSblock', i, '.csv', sep=''))
  power = read.csv(paste('results/power/covariance/block/DSblock', i, '.csv', sep=''))
  fdr.ds.block = c(fdr.ds.block, mean(fdr$x[fdr$x!=0]))
  power.ds.block = c(power.ds.block, mean(power$x))
}

fdr.mds.block.100 = c()
power.mds.block.100 = c()
for(i in ind){
  fdr = read.csv(paste('results/fdr/covariance/block/block100_', i, '.csv', sep=''))
  power = read.csv(paste('results/power/covariance/block/block100_', i, '.csv', sep=''))
  fdr.mds.block.100 = c(fdr.mds.block.100, mean(fdr$x[fdr$x!=0]))
  power.mds.block.100 = c(power.mds.block.100, mean(power$x))
}
