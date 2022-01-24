# File: preterm.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: preliminary analysis of preterm data
# Date: 21/01/2022

load('dataExternal/pretermdata/PretermData.rda')

tapply(CGAdelivery, factor(Ccohort), summary)
tapply(CGAsampling, factor(Ccohort), summary)
tapply(CGSamDel, factor(Ccohort), summary)
xtabs(~ Ccohort + Csite)
## prepare data to use in the analysis
# col names are very long, replace with shorter names
# dfKey = data.frame(key = paste('Met', 1:ncol(CMetabolomics), sep=':'), name=colnames(CMetabolomics))
# colnames(CMetabolomics) = dfKey$key

m = log((CTranscriptomics)+1)
table(is.finite(m))
table(is.na(m))
table(is.na(rowSums(m)))
table(is.na(colSums(m)))
# i = which(is.na(colSums(m)))
rm(m)
dfData = data.frame((CTranscriptomics))
table(complete.cases(dfData))
dim(dfData)
## remove variables with 0 sd i.e. not changing 
s = apply(dfData, 2, sd)
summary(s)
s = which(s == 0)
length(s)
dfData = dfData[,-s]
dim(dfData)
lData.train = list(data=dfData, covariates=factor(Ccohort))

p.vals = lapply(1:ncol(lData.train$data), function(x){
  df = data.frame(y=lData.train$data[,x], d=lData.train$covariates, c=factor(Csite))
  f = lm(y ~ d + c, data=df)
  s = summary(f)$coefficients
  return(s['d1', 4])
})

names(p.vals) = colnames(lData.train$data)
dfPvals = do.call(rbind, p.vals)
dfPvals = cbind(dfPvals, p.adjust(dfPvals[,1], method = 'BH'))
colnames(dfPvals) = c('pvalue', 'p.adj')
dfPvals = data.frame(dfPvals)
f = which(dfPvals$pvalue < 0.01)
length(f)
cvTopVariables.lm = rownames(dfPvals)[f]


if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

########################## perform a random forest step
dfData = lData.train$data[,cvTopVariables.lm]
fGroups = lData.train$covariates

oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)

plot.var.selection(oVar.r)

######################## Stan section for binomial regression approach
dfData = lData.train$data[,cvTopVariables.lm]
dim(dfData)
dfData$fGroups = lData.train$covariates
table(dfData$fGroups)
rm(fGroups)
levels(dfData$fGroups)
lData = list(resp=ifelse(dfData$fGroups == '1', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rethinking)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
# initf = function(chain_id = 1) {
#   list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
# }


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4,# init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

# save(fit.stan, file='temp/fit.stan.binom_preterm_met.rds')

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest
mCoef = extract(fit.stan)$betas2
dim(mCoef)
## get the intercept
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]

## coeftab object 
ct.1 = coeftab(fit.stan)
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
rownames(ct.1@se)[i[-1]] = colnames(mCoef)
plot(ct.1, pars=colnames(mCoef))

m = abs(colMeans(mCoef))
m = sort(m, decreasing = T)

l2 = barplot(m, 
             las=2, xaxt='n', col='grey', main='Top Variables')
axis(1, at = l2, labels = names(m), tick = F, las=2, cex.axis=0.7 )

### predictive performance of the model
## binomial prediction
mypred = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  #iFitted = plogis(iFitted)
  return(iFitted)
}


mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
library(lattice)
## get the predicted values
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', 
       ylab='Predicted Probability of Pre-term',
       data=dfData)


# ## find correlated variables
dim(dfData)
mData = as.matrix(dfData[,-87])
length(as.vector(mData))
mCor = cor(mData, use="na.or.complete")
library(caret)
image(mCor)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.7, names=T)
data.frame(n)
# sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
# 
n = sapply(n, function(x) {
  rownames(mCor)[(abs(mCor[,x]) >= 0.7)]
})

cvTopVariables.rf = rownames(CVariableSelection.RandomForest.getVariables(oVar.r))[1:20]
cvTopVariables.bin = names(m)[1:20]
table(cvTopVariables.bin %in% cvTopVariables.rf)
cvTopVariables = unique(c(cvTopVariables.rf, cvTopVariables.bin))

## subset selection
dfData = lData.train$data[,cvTopVariables]
fGroups = lData.train$covariates
oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)
plot.var.selection(oVar.sub)
table(fGroups)
log(40)
# select 4 variables
cvVar = CVariableSelection.ReduceModel.getMinModel(oVar.sub, size = 3)

table(cvVar %in% cvTopVariables.bin)
table(cvVar %in% cvTopVariables.rf)

# cross validation
oCV.lda = CCrossValidation.LDA(dfData[,cvVar], dfData[,cvVar], fGroups, fGroups, level.predict = '1',
                               boot.num = 100, k.fold = 10) 

plot.cv.performance(oCV.lda)

