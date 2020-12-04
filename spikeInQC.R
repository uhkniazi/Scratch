# Name: spikeInQC.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 4/12/2020
# Desc: testing functions for QC of RNA-Seq spikein


library(zebrafishRNASeq)
data("zfGenes")

mData = as.matrix(zfGenes)
head(mData)
fGroup = gl(2, k = 3, labels = c('Ctl', 'Trt'))

f1 = rowMeans(mData) < 3
f2 = apply(mData, 1, function(x) length(x[x>5]) > 2)
table(f1); table(f2)
table(names(f1[f1]) %in% names(f2[!f2]))

## choose filter 2
mData = mData[f2,]
dim(mData)

## filter the spike ins
i = grep('ercc', x = rownames(mData), ignore.case = T)
mErcc = mData[i,]
mData = mData[-i,]
dim(mData); dim(mErcc)

####### proportion of reads mapping to erccs
iD = colSums(mData)
iE = colSums(mErcc)
iE = iE/(iD+iE)

barplot(iE)

plot.ercc.proportions = function(ivTotal, ivErcc, fGroupings, ...){
  pErcc = ivErcc/ivTotal * 100
  c = rainbow(nlevels(fGroupings))
  barplot(pErcc[order(fGroupings)], col=c[as.numeric(fGroupings)[order(fGroupings)]],
          ...)
  abline(h = mean(pErcc), lty=2)
}

plot.ercc.proportions(colSums(rbind(mData, mErcc)), colSums(mErcc),
                      fGroupings = fGroup,
                      ylim=c(0, 10), 
                      ylab='Percentage of total reads mapping to ERCC')


library(MASS)
library(car)
iE = logit(iE)
t.test(iE ~ fGroup)

############ model and ma plot of data

## log posterior function
library(LearnBayes)
library(numDeriv)

mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=10000), method='Nelder-Mead', data=data)
  # calculate hessian
  fit$hessian = (hessian(logpost, fit$par, data=data))
  colnames(fit$hessian) = names(mode)
  rownames(fit$hessian) = names(mode)
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  stuff = list(mode = mode, var = h, converge = fit$convergence == 
                 0)
  return(stuff)
}

mylogpost = function(theta, data){
  ## parameters to track/estimate
  betas = theta[grep('betas', names(theta))] # vector of betas i.e. regression coefficients for population
  sigmaPop = exp(theta['sigmaPop']) # population level sd
  ## data
  resp = data$resp # response variable
  mModMatrix = data$mModMatrix # design matrix 
  
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # write the priors and likelihood 
  lp = dnorm(betas[1], 0, 2, log=T) + sum(dnorm(betas[-1], 0, 2, log=T)) + 
    dexp(sigmaPop, 1/28, log=T)
  lik = sum(dnorm(resp,mean = iFitted, sd=sigmaPop, log=T))
  val = lik + lp
  return(val)
}

## create data for modelling
m = mData[sample(1:nrow(mData), 100, F), ]
m = log(rbind(m, mErcc)+1)
#m = m[, 1:2]
head(m)
f = gl(2, 100, labels = c('C1', 'C3'))
f=fGroup

lData = list(resp=c(temp1, temp2), mModMatrix=model.matrix(c(temp1, temp2) ~  f - 1))
start = c('sigmaPop'=log(1), 'betas'=rep(0, times=ncol(lData$mModMatrix)))

mylogpost(start, lData)
mylaplace(mylogpost, start, lData)
