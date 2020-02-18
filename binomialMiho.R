# binomialMiho.R

# dfData = read.csv('dataExternal/NTD methylation_2probes and covariates for Umar.csv', header=T)
# 
library(lme4)

fit.2 = glm(as.matrix(dfData[,lCn[[10]]]) ~ 1 + dfData.cov$Diets , family='binomial')
# 
# # load some test data
# data("cbpp")
# dfData = cbpp
# str(dfData)
# ## test in stan
# library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# stanDso = rstan::stan_model(file='binomialRegressionRandomEffects.stan')
# 
# lStanData = list(Ntotal=nrow(dfData),
#                  Nlevels1 = nlevels(dfData$herd),
#                  NfactorMap1 = as.numeric(dfData$herd),
#                  y=dfData$incidence,
#                  Ntrials=dfData$size)
# 
# fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('intercept', 'nCoefFactor1', 'mu',
#                                                                          'sigmaFactor1'),
#                     cores=4)
# print(fit.stan, c('intercept', 'nCoefFactor1', 'sigmaFactor1'), digits=3)
# 
# fit.1 = glmer(cbind(incidence, size - incidence) ~ 1 + (1 | herd), data = dfData, family = binomial)
# summary(fit.1)
# 
# # compare the random effects
# r = ranef(fit.1)
# s = extract(fit.stan)$nCoefFactor1
# dim(s)
# plot(r$herd[,1], colMeans(s))
# cor(r$herd[,1], colMeans(s))
# 
# ## 3 random effects
# dfData = read.csv('dataExternal/NTD methylation_2probes and covariates for Umar.csv', header=T)
# str(dfData)
# library(lme4)
# 
# fit.1 = glmer(cbind(probe2_Gtf3a_methylated, probe2_Gtf3a_unmethylated) ~ 1 + (1|diet) +
#                 (1 | sex) + (1 | genotype), data=dfData, family='binomial')
# 
# summary(fit.1)
# 
# #stanDso = rstan::stan_model(file='binomialRegression3RandomEffects.stan')
# 
# lStanData = list(Ntotal=nrow(dfData),
#                  Nlevels1 = nlevels(dfData$diet),
#                  NfactorMap1 = as.numeric(dfData$diet),
#                  Nlevels2 = nlevels(dfData$sex),
#                  NfactorMap2 = as.numeric(dfData$sex),
#                  Nlevels3 = nlevels(dfData$genotype),
#                  NfactorMap3 = as.numeric(dfData$genotype),
#                  y=dfData$probe2_Gtf3a_methylated,
#                  Ntrials=dfData$probe2_Gtf3a_methylated + dfData$probe2_Gtf3a_unmethylated)
# 
# fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('intercept', 'nCoefFactor1',
#                                                                          'nCoefFactor2', 
#                                                                          'nCoefFactor3',
#                                                                          'sigmaFactor1',
#                                                                          'sigmaFactor2',
#                                                                          'sigmaFactor3'),
#                     cores=4)
# print(fit.stan, c('intercept', 'nCoefFactor1',  
#                   'nCoefFactor2', 'nCoefFactor3',
#                   'sigmaFactor1', 'sigmaFactor2', 'sigmaFactor3'), digits=3)
# 
# traceplot(fit.stan, c('intercept', 'nCoefFactor1',  
#                   'nCoefFactor2', 'nCoefFactor3',
#                   'sigmaFactor1', 'sigmaFactor2', 'sigmaFactor3'))
# 
# 
# summary(fit.1)
# 
# 
# 
# # compare the random effects
# r = ranef(fit.1)
# s = extract(fit.stan)$nCoefFactor1
# dim(s)
# plot(r$diet[,1], colMeans(s))
# cor(r$diet[,1], colMeans(s))


###################################################### start here
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegressionRandomEffects.stan')
dfData = read.csv('dataExternal/S34_D1spMF_D3spMF_prom_1x_meth_unmeth_FEB2020.csv', header=T)

## format the data 
dfData.cov = dfData[,c(1:5)]
dfData.orig = dfData
cn = colnames(dfData.orig[,-c(1:5)])
i = gsub('methylated|unmethylated', cn, replacement = '')
i = factor(i)
lCn = split(cn, f = i)
###################################### function to fit model
modelFunction = function(lData){
  fs = tryCatch(sampling(stanDso, data=lData, iter=1000, chains=4, pars=c('intercept', 'nCoefFactor1'),
                                                                           #'sigmaFactor1'),
                      cores=4), error=function(e) NULL)
  return(fs)
}

## format the lStan data input before calling function
lStanData.1 = list(Ntotal=nrow(dfData.cov),
                 Nlevels1 = nlevels(dfData.cov$Diets),
                 NfactorMap1 = as.numeric(dfData.cov$Diets),
                 y=dfData[,lCn[['prom_68215_']][1]],
                 Ntrials=dfData[,lCn[['prom_68215_']][1]]+dfData[,lCn[['prom_68215_']][2]])

fit.stan.1 = modelFunction(lStanData.1)


##################################### extract results
## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

levels(dfData.cov$Diets)

extractResult = function(fit.stan, base='HFA', deflection='LFA', lev=levels(dfData.cov$Diets)){
  if (is.null(fit.stan)) return(NULL)
  mCoef = extract(fit.stan)$nCoefFactor1
  colnames(mCoef) = lev
  dif = getDifference(mCoef[,deflection], mCoef[,base])
  r = data.frame(coef.base=mean(mCoef[,base]), 
                 coef.deflection=mean(mCoef[,deflection]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
}

lResults = vector(mode = 'list', length=length(lCn))
# extractResult(fit.stan.1)
for (i in 1188:length(lCn)) {
  x = lCn[[i]]
  ## format the lStan data input before calling function
  lStanData.x = list(Ntotal=nrow(dfData.cov),
                     Nlevels1 = nlevels(dfData.cov$Diets),
                     NfactorMap1 = as.numeric(dfData.cov$Diets),
                     y=dfData[,x[1]],
                     Ntrials=dfData[,x[1]]+dfData[,x[2]])
  
  fit.stan.x = modelFunction(lStanData.x)
  lResults[[i]] = extractResult(fit.stan.x)
  rm(fit.stan.x)
  gc(reset = T)
  save(lResults, file='results/miho_2.rds')
  cat(paste(i, ' ......'))
}


# lResults = lapply(lCn, function(x){
#   ## format the lStan data input before calling function
#   lStanData.x = list(Ntotal=nrow(dfData.cov),
#                      Nlevels1 = nlevels(dfData.cov$Diets),
#                      NfactorMap1 = as.numeric(dfData.cov$Diets),
#                      y=dfData[,x[1]],
#                      Ntrials=dfData[,x[1]]+dfData[,x[2]])
#   
#   fit.stan.x = modelFunction(lStanData.x)
#   ret = extractResult(fit.stan.x)
#   rm(fit.stan.x)
#   gc(reset = T)
#   return(ret)
# })
# 
# save(lResults, file='results/miho.rds')
