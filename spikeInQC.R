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
# iD = colSums(mData)
# iE = colSums(mErcc)
# iE = iE/(iD+iE)
# 
# barplot(iE)

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


# library(MASS)
# library(car)
# iE = logit(iE)
# t.test(iE ~ fGroup)

plot.ercc.MD = function(mLogCounts, fGroupings, ...){
  if (nlevels(fGroupings) != 2) stop('grouping factor should be 2 levels only')
  mResults = apply(mLogCounts, 1, function(x){
    tapply(x, fGroupings, mean)
  })
  mResults = t(mResults)
  ## difference
  d = mResults[,2] - mResults[,1]
  ## choose colour
  c = rep('grey', length.out=length(d))
  c[grep('ercc', names(d), ignore.case = T)] = 'blue'
  plot(rowMeans(mLogCounts), d, pch=20, col=c, xlab='mean log(count+1)', 
       ylab='log difference', ...)
  abline(h = 0, lty=2)
  i = grep('ercc', rownames(mLogCounts), ignore.case = T)
  lines(lowess(rowMeans(mLogCounts)[-i], d[-i]), col='black')
  lines(lowess(rowMeans(mLogCounts)[i], d[i]), col='blue')
}

## prepare input data and select groupings
m = log(rbind(mData, mErcc)+1)
head(m)
m = m[,c(1,3)]
head(m)
f = gl(2, k = 1, labels = c('C1', 'C5'))

plot.ercc.MD(m, f, ylim=c(-5, 5), main='MD plot')

############ model and ma plot of data

# ## log posterior function
# library(LearnBayes)
# library(numDeriv)
# 
# mylaplace = function (logpost, mode, data) 
# {
#   options(warn = -1)
#   fit = optim(mode, logpost, gr = NULL,  
#               control = list(fnscale = -1, maxit=10000), method='Nelder-Mead', data=data)
#   # calculate hessian
#   fit$hessian = (hessian(logpost, fit$par, data=data))
#   colnames(fit$hessian) = names(mode)
#   rownames(fit$hessian) = names(mode)
#   options(warn = 0)
#   mode = fit$par
#   h = -solve(fit$hessian)
#   stuff = list(mode = mode, var = h, converge = fit$convergence == 
#                  0)
#   return(stuff)
# }
# 
# mylogpost = function(theta, data){
#   ## parameters to track/estimate
#   betas = theta[grep('betas', names(theta))] # vector of betas i.e. regression coefficients for population
#   sigmaPop = exp(theta['sigmaPop']) # population level sd
#   ## data
#   resp = data$resp # response variable
#   mModMatrix = data$mModMatrix # design matrix 
#   
#   # calculate fitted value
#   iFitted = mModMatrix %*% betas
#   # write the priors and likelihood 
#   lp = dcauchy(betas[1], 0, 5, log=T) + sum(dcauchy(betas[-1], 0, 5, log=T)) + 
#     dexp(sigmaPop, 1/28, log=T)
#   lik = sum(dnorm(resp,mean = iFitted, sd=sigmaPop, log=T))
#   val = lik + lp
#   return(val)
# }

## create data for modelling
# m = mData[sample(1:nrow(mData), 100, F), ]
# m = log(rbind(m, mErcc)+1)
# m = m[, 1:2]
# head(m)
# f = gl(2, 1, labels = c('C1', 'C3'))
# f = fGroup
# lData = list(resp=m[1,], mModMatrix=model.matrix(m[1,] ~  f - 1))
# start = c('sigmaPop'=log(1), 'betas'=rep(0, times=ncol(lData$mModMatrix)))

# mylogpost(start, lData)
# fit.1 = mylaplace(mylogpost, start, lData)
# 
# getDifference = function(ivData, ivBaseline){
#   stopifnot(length(ivData) == length(ivBaseline))
#   # get the difference vector
#   d = ivData - ivBaseline
#   # get the z value
#   z = mean(d)/sd(d)
#   # get 2 sided p-value
#   p = pnorm(-abs(mean(d)/sd(d)))*2
#   return(list(d=mean(d), z=z, p=p))
# }
# 
# tpar = list(m=fit.1$mode, var=fit.1$var*2, df=4)
# ## get a sample directly and using sir (sampling importance resampling with a t proposal density)
# s = sir(mylogpost, tpar, 5000, lData)
# colnames(s)[-1] = levels(f)
# 
# getDifference(s[,3], s[,2])

############# put all this into a single function
modelFunction = function(start, lData, cBase, cDeflection){
  
  ### internal functions
########################  
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
    lp = dcauchy(betas[1], 0, 5, log=T) + sum(dcauchy(betas[-1], 0, 5, log=T)) + 
      dexp(sigmaPop, 1/28, log=T)
    lik = sum(dnorm(resp,mean = iFitted, sd=sigmaPop, log=T))
    val = lik + lp
    return(val)
  }
  
  getDifference = function(ivData, ivBaseline){
    stopifnot(length(ivData) == length(ivBaseline))
    # get the difference vector
    d = ivData - ivBaseline
    # get the z value
    z = mean(d)/sd(d)
    # get 2 sided p-value
    p = pnorm(-abs(mean(d)/sd(d)))*2
    return(list(d=mean(d), z=z, p=p))
  }
  
###############
  fit.1 = mylaplace(mylogpost, start, lData)
  tpar = list(m=fit.1$mode, var=fit.1$var*2, df=4)
  ## get a sample directly and using sir (sampling importance resampling with a t proposal density)
  s = sir(mylogpost, tpar, 5000, lData)
  colnames(s)[-1] = levels(f)
  return(getDifference(s[,cDeflection], s[,cBase]))
}
dim(mData)
m = mData#[sample(1:nrow(mData), 1000, F), ]
m = log(rbind(m, mErcc)+1)
head(m)
m = m[, c(1,3)]
head(m)
dim(m)
f = gl(2, 1, labels = c('C1', 'C5'))
levels(f)
#f = fGroup
#f = factor(c('C1', 'C2', 'C2'))
# lData = list(resp=m[2,], mModMatrix=model.matrix(m[2,] ~  f - 1))
# start = c('sigmaPop'=log(1), 'betas'=rep(0, times=ncol(lData$mModMatrix)))
# levels(f)
# summary(lm(lData$resp ~ f))
# modelFunction(start, lData, 'C1', 'C2')
mResults = apply(m, 1, function(x){
  tapply(x, f, mean)
})
mResults = t(mResults)
dim(mResults)
head(mResults)
d = mResults[,2] - mResults[,1]
c = rep('grey', length.out=length(d))
c[grep('ercc', names(d), ignore.case = T)] = 'blue'
identical(rownames(m), names(d))
plot(rowMeans(m), d, pch=20, col=c, xlab='mean log(count+1)', 
     ylab='log difference', ylim=c(-6, 6))
abline(h = 0, lty=2)
i = grep('ercc', rownames(m), ignore.case = T)
lines(lowess(rowMeans(m)[-i], d[-i]), col='black')
lines(lowess(rowMeans(m)[i], d[i]), col='blue')

##########
lResult = apply(m, 1, function(x){
  lData = list(resp=x, mModMatrix=model.matrix(x ~ f - 1))
  start = c('sigmaPop'=log(1), 'betas'=rep(0, times=ncol(lData$mModMatrix)))
  unlist(tryCatch(modelFunction(start, lData, 'C1', 'C2'), error=function(e) NULL))
})

f = sapply(lResult, is.null)
table(f)
lResult = lResult[!f]

dfResult = data.frame(do.call(rbind, lResult))

hist(dfResult$p)

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p < p.cut] = 'red'
  col[grep('ercc', rownames(dat), ignore.case = T)] = 'blue'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$d, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC', ylim=c(-2, 2), cex=0.6)
  #abline(h = 0, lty=1)
  i = grep('ercc', rownames(dat), ignore.case = T)
  lines(lowess(m[i], dat$d[i]), col='blue')
  lines(lowess(m[-i], dat$d[-i]), col='black')
}

plotMeanFC(rowMeans(m[rownames(dfResult),]), dfResult, 0.05, 'test')
