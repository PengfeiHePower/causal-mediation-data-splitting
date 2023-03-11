for(filename in list.files('results')){
  if(!(filename %in% c('fdr', 'power'))){
    load(paste('results',filename, sep = '/'))
  }
}

avg.fdr.DS.ID = mean(fdr.DS.Id)
avg.fdr.MDS.ID.1 = mean(fdr.MDS.Id.1)
avg.fdr.MDS.ID.2 = mean(fdr.MDS.Id.2)
avg.fdr.MDS.ID.4 = mean(fdr.MDS.Id.4)
avg.fdr.MDS.ID.6 = mean(fdr.MDS.Id.6)
avg.fdr.MDS.ID.8 = mean(fdr.MDS.Id.8)
avg.fdr.MDS.ID.10 = mean(fdr.MDS.Id.10)
avg.fdr.MDS.ID.20= mean(fdr.MDS.Id.20)
avg.fdr.MDS.ID.50 = mean(fdr.MDS.Id.50)
avg.fdr.MDS.ID.100 = mean(fdr.MDS.Id.100)
avg.fdr.MDS.ID.150 = mean(fdr.MDS.Id.150)
avg.fdr.MDS.ID.200 = mean(fdr.MDS.Id.200)
avg.fdr.MDS.ID.250 = mean(fdr.MDS.Id.250)

avg.power.DS.ID = mean(power.DS.Id)
avg.power.MDS.ID.1 = mean(power.MDS.Id.1)
avg.power.MDS.ID.2 = mean(power.MDS.Id.2)
avg.power.MDS.ID.4 = mean(power.MDS.Id.4)
avg.power.MDS.ID.6 = mean(power.MDS.Id.6)
avg.power.MDS.ID.8 = mean(power.MDS.Id.8)
avg.power.MDS.ID.10 = mean(power.MDS.Id.10)
avg.power.MDS.ID.20 = mean(power.MDS.Id.20)
avg.power.MDS.ID.50 = mean(power.MDS.Id.50)
avg.power.MDS.ID.100 = mean(power.MDS.Id.100)
avg.power.MDS.ID.150 = mean(power.MDS.Id.150)
avg.power.MDS.ID.200 = mean(power.MDS.Id.200)
avg.power.MDS.ID.250 = mean(power.MDS.Id.250)


avg.fdr.MDS = c(avg.fdr.MDS.ID.1, avg.fdr.MDS.ID.2, avg.fdr.MDS.ID.4, avg.fdr.MDS.ID.6, avg.fdr.MDS.ID.8, 
                avg.fdr.MDS.ID.10, avg.fdr.MDS.ID.20, avg.fdr.MDS.ID.50, avg.fdr.MDS.ID.100,
                avg.fdr.MDS.ID.150, avg.fdr.MDS.ID.200, avg.fdr.MDS.ID.250)

avg.power.MDS = c(avg.power.MDS.ID.1, avg.power.MDS.ID.2, avg.power.MDS.ID.4, avg.power.MDS.ID.6, avg.power.MDS.ID.8, 
                  avg.power.MDS.ID.10, avg.power.MDS.ID.20, avg.power.MDS.ID.50, avg.power.MDS.ID.100,
                  avg.power.MDS.ID.150, avg.power.MDS.ID.200, avg.power.MDS.ID.250)
iters = c(1,2,4,6,8,10,20,50,100,150,200,250)
plot(iters, avg.fdr.MDS, pch=19, ylim=c(0,0.1), main = 'fdr from different iters')
abline(avg.fdr.DS.ID, b=0, col='red')
plot(iters, avg.power.MDS, pch=19, type='b', main = 'power for different iters', ylim=c(0.92,1))
abline(avg.power.DS.ID, b=0, col='red')

