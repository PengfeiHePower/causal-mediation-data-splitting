setwd("/mnt/ufs18/home-016/hepengf1/Documents/stat-research/causal-mediation-data-splitting")

library(MASS)
library(stats)
library(glmnet)
library(qvalue)
library(qqman)
library(locfdr)
# library(DACT)

nonnullPropEst <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation

  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)

  epsest=NULL

  for (j in 1:length(tt)) {

    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi

    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(max(epsest))
}

EfronCorrect = function(pval){
  z = stats::qnorm(1-pval)
  res <- locfdr(z,nulltype = 1)
  mean.emp = res$fp0["mlest","delta"]
  sd.emp = res$fp0["mlest","sigma"]
  pval.emp = stats::pnorm(z,mean = mean.emp,sd = sd.emp,lower.tail = F)
  return(pval.emp)
}

JCCorrect = function(pval){
  z = stats::qnorm(pval,lower.tail = F)
  res= nullParaEst(z)
  pval.JC = stats::pnorm(z,mean = res$mu,sd = res$s,lower.tail = F)
  return(pval.JC)
}

nullParaEst<-function (x,gamma=0.1)
{
 # x is a vector of z-values
 # gamma is a parameter, default is 0.1
 # output the estimated mean and standard deviation

 n = length(x)
 t = c(1:1000)/200

 gan    = n^(-gamma)
 that   = 0
 shat   = 0
 uhat   = 0
 epshat = 0

 phiplus   = rep(1,1000)
 phiminus  = rep(1,1000)
 dphiplus  = rep(1,1000)
 dphiminus = rep(1,1000)
 phi       = rep(1,1000)
 dphi      = rep(1,1000)

 for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
 }

 ind = min(c(1:1000)[(phi - gan) <= 0])
 tt = t[ind]
 a  = phiplus[ind]
 b  = phiminus[ind]
 da = dphiplus[ind]
 db = dphiminus[ind]
 c  = phi[ind]

 that   = tt
 shat   = -(a*da + b*db)/(tt*c*c)
 shat   = sqrt(shat)
 uhat   = -(da*b - db*a)/(c*c)
 epshat = 1 - c*exp((tt*shat)^2/2)

 return(musigma=list(mu=uhat,s=shat))
}

DACT = function(Z_a,Z_b, p_a, p_b,correction=NULL){
#   Z_a = stats::qnorm(p_a,lower.tail = F)
#   Z_b = stats::qnorm(p_b,lower.tail = F)
  pi0a = 1 - nonnullPropEst(Z_a,0,1)
  pi0b = 1 - nonnullPropEst(Z_b,0,1)
  #pi0a = locfdr::locfdr(Z_a,nulltype = 0)$fp0[5,3]
  #pi0b = locfdr::locfdr(Z_b,nulltype = 0)$fp0[5,3]
  if(pi0a > 1){
    pi0a = 1
  }
  if(pi0b >1){
    pi0b = 1
  }
  p.mat = cbind(p_a,p_b)
  p3 = (apply(p.mat,1,max))^2
  wg1 = pi0a*(1-pi0b)
  wg2 = (1-pi0a)*pi0b
  wg3 = pi0a*pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1,wg2,wg3)/wg.sum
  p_dact = wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3
  if(correction == "Efron"){
  p_dact = EfronCorrect(p_dact)
  }
  if(correction == "JC"){
    p_dact = JCCorrect(p_dact)
  }
  return(p_dact)
}


require("getopt", quietly=TRUE)

spec = matrix(c(
    "Filename", "f", 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec);

filename = opt$Filename
cat('Data location:', filename, '\n')

load(file=filename)

fdr_level = 0.05

Zbeta.ob = c() #for m model
Zgamma.ob = c() #for y model
## first obtain Z-score for beta_1 and gamma_1
for(i in 1:p){
    mediator = M[,i]
    modelm = lm(mediator~treat+pretreat)
    modely = lm(Y~treat+mediator+pretreat)
    summary.m = summary(modelm)
    summary.y = summary(modely)
    Zbeta.ob = c(Zbeta.ob, summary.m$coefficients[2, "t value"])
    Zgamma.ob = c(Zgamma.ob, summary.y$coefficients[3, "t value"])
}

# pbeta.ob <- pnorm(abs(Zbeta.ob), lower.tail = F) * 2
# pgamma.ob <- pnorm(abs(Zgamma.ob), lower.tail = F) * 2
pbeta.ob <- pnorm(Zbeta.ob, lower.tail = F)
pgamma.ob <- pnorm(Zgamma.ob, lower.tail = F)

# Z_a = stats::qnorm(pbeta.ob,lower.tail = F)
# Z_b = stats::qnorm(pgamma.ob,lower.tail = F)
# pi0a = 1 - nonnullPropEst(Z_a ,0,1)
# pi0b = 1 - nonnullPropEst(Z_b,0,1)

res = DACT(Z_a=Zbeta.ob,Z_b=Zgamma.ob,p_a =pbeta.ob, p_b =pgamma.ob,  correction="JC")

S1.dact = which(res<=fdr_level)
S0.dact = setdiff(1:p, S1.dact)

fdr.dact = length(intersect(S0,S1.dact))/max(1,length(S1.dact))
power.dact = length(intersect(S0, S0.dact))/length(S0.dact)
cat('Fdr:', fdr.dact, "\n")
cat('Power:', power.dact, "\n")
