##### data structure
##### Model mediator
# M=gamma.0+gamma.1*treat+gamma.2*pretreat+epsilon.m, n by p

##### Model outcome
# y=beta.0 + beta.1*treat + beta.2*M + beta.3*pretreat+epsilon.y

# Most important variable:
# observation: treat, pretreat
# coefficient: gamma.0, gamma.1, gamma.2; beta.0, beta.1, beta.2, beta.3
# error: epsilon.y
setwd("/mnt/ufs18/home-016/hepengf1/Documents/stat-research/causal-mediation-data-splitting")

library(MASS)
library(stats)
library(glmnet)
library(HIMA)
library(stringr)

require("getopt", quietly=TRUE)

spec = matrix(c(
    "Filename", "f", 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec);

filename = opt$Filename
cat('Data location:', filename, '\n')

load(file=filename)

pheno.df = data.frame(treat = treat)
colnames(pretreat) = c('pre1','pre2','pre3','pre4','pre5','pre6','pre7','pre8',
'pre9','pre10','pre11','pre12','pre13','pre14','pre15','pre16','pre17','pre18','pre19','pre20')
pretreat.df = data.frame(pretreat)
pheno.df = cbind(pheno.df, pretreat.df)
pheno.df$Y = Y

hima.fit = hima2(Y~treat+pre1+pre2+pre3+pre4+pre5+pre6+pre7+pre8+pre9+pre10+pre11+pre12+
pre13+pre14+pre15+pre16+pre17+pre18+pre19+pre20,
data.pheno=pheno.df,
data.M=M,
outcome.family = "gaussian",
mediator.family = "gaussian",
penalty = "SCAD",
scale = FALSE)

#hima.fit
selected = rownames(hima.fit)
S1.hima = c()
for(name in selected){
  S1.hima = c(S1.hima, as.numeric(str_extract(name, "\\d+")))
}
S0.hima = setdiff(1:p, S1.hima)
fdr.hima = length(intersect(S0,S1.hima))/max(1,length(S1.hima))
power.hima = length(intersect(S0, S0.hima))/length(S0.hima)
cat('Fdr:', fdr.hima, "\n")
cat('Power:', power.hima, "\n")