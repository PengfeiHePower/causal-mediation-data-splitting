library(MASS)
library(stats)
library(glmnet)

# dataBlock1
m1 = toeplitz(c(1,0.1,0.1,0))
M1 = diag(500)
ind1 = seq(1,200,4)
for(i in ind1){
  M1[i:(i+3), i:(i+3)] = m1
}
kappa(M1)

# dataBlock2
m2.1 = toeplitz(c(1,0.1,0.1))
m2.2 = toeplitz(c(1,0.1,0.1,0.1))
M2 = diag(500)
ind2.1 = seq(1,200,3)
ind2.2 = seq(201,400,4)
for(i in ind2.1){
  M2[i:(i+2), i:(i+2)] = m2.1
}
for(i in ind2.2){
  M2[i:(i+3), i:(i+3)] = m2.2
}
kappa(M2)

# dataBlock3
m3.1 = toeplitz(c(1,0.1,0.1))
m3.5 = toeplitz(c(1,0.1,0.1,0.1,0.1,0.1,0.1))
M3 = diag(500)
ind3.1 = seq(1,300,3)
ind3.5 = seq(301,400,7)
for(i in ind3.1){
  M3[i:(i+2),i:(i+2)] = m3.1
}
for(i in ind3.5){
  M3[i:(i+6),i:(i+6)] = m3.5
}
kappa(M3)


# dataBlock3
m4.1 = toeplitz(c(1,0.1,0.1))
M4 = diag(500)
ind4.1 = sample(1:450, 100)
for(i in ind4.1){
  M4[i:(i+2),i:(i+2)] = m4.1
}
kappa(M4)
