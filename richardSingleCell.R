# richardSingleCell.R
# some DE analysis of a count matrix

## load count matrix
mCounts = read.csv(file.choose(), header=T)
mCounts = mCounts[!duplicated(mCounts$X), ]
rownames(mCounts) = mCounts$X
mCounts = mCounts[,-1]
mCounts = as.matrix(mCounts)
colnames(mCounts) = gsub('X', '', colnames(mCounts))
# load the covariates
dfCovariates = read.csv(file.choose(), header=T)
identical(colnames(mCounts), as.character(dfCovariates$SampleID))
dfCovariates$patientId = factor(dfCovariates$patientId)

mData = mCounts
dim(mData)

dim(na.omit(mData))

summary(mData)
# drop the samples where average across rows is less than 1
i = rowMeans(mData)
table( i < 1)
mData = mData[!(i< 1),]
dim(mData)

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData) | inData < 1 ] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

table(ivProb < 0.4)
mData = mData[!(ivProb < 0.4), ]
dim(mData)
## load CDiagnostics and test
## compare the normalised and raw data
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(log(mData+0.5), 'Normalised')
# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
str(dfCovariates)
fBatch = dfCovariates$patientId
fBatch2 = dfCovariates$treatment

#pdf('results/matrixClustering.pdf')

## check using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.1, fBatch2, legend.pos = 'topright', axis.label.cex = 0.7)

par(mfrow=c(1,2))
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7, legend.pos = 'topright')
plot.mean.summary(oDiag.1, fBatch2, axis.label.cex = 0.7, legend.pos = 'topright')

par(mfrow=c(1,2))
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.1, fBatch2, axis.label.cex = 0.7)

par(mfrow=c(1,2))
plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)

par(mfrow=c(1,1))
plot.PCA(oDiag.1, fBatch)
plot.PCA(oDiag.1, fBatch, csLabels = dfCovariates$treatment, legend.pos = 'topleft')

plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)
plot.dendogram(oDiag.1, fBatch2, labels_cex = 0.7)

###############################################################################################
################ clusters and covariates explaining variance
plot(oDiag.1@lData$PCA$sdev)
mPC = oDiag.1@lData$PCA$x[,1:4]
i = which(mPC[,1] > 150)
fNewBatch = rep(1, times=length(dfCovariates$SampleID))
fNewBatch[i] = 2
fNewBatch = factor(fNewBatch)

## check if batch assigned correctly
plot.PCA(oDiag.1, fNewBatch)
plot.dendogram(oDiag.1, fNewBatch, labels_cex=2)
## try a linear mixed effect model to account for varince
library(lme4)
dfData = data.frame(mPC)
dfData = stack(dfData)
dfData$fBatch = dfCovariates$treatment
dfData$fAdjust1 = dfCovariates$patientId
dfData$fAdjust2 = fNewBatch

dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj1 = factor(dfData$fAdjust1:dfData$ind)
dfData$Coef.adj2 = factor(dfData$fAdjust2:dfData$ind)

str(dfData)
plot(density(dfData$values))

fit.lme1 = lmer(values ~ 1  + (1 | Coef) + (1 | Coef.adj1) + (1 | Coef.adj2), data=dfData)
summary(fit.lme1)

fit.lme2 = lmer(values ~ 1  + (1 | Coef) + (1 | Coef.adj1), data=dfData)
summary(fit.lme2)

fit.lme3 = lmer(values ~ 1  + (1 | Coef), data=dfData)
summary(fit.lme3)

anova(fit.lme1, fit.lme2, fit.lme3)

plot((fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme1)), resid(fit.lme1)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme1)))

plot((fitted(fit.lme2)), resid(fit.lme2), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme2)), resid(fit.lme2)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme2)))
lines(density(fitted(fit.lme1)), col=2)

plot((fitted(fit.lme3)), resid(fit.lme3), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme3)), resid(fit.lme3)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme2)))
lines(density(fitted(fit.lme1)), col=2)
lines(density(fitted(fit.lme3)), col=3)

## fit model with stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='tResponse3RandomEffectsNoFixed.stan')

### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj1),
                 Nclusters3=nlevels(dfData$Coef.adj2),
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj1),
                 NgroupMap3=as.numeric(dfData$Coef.adj2),
                 y=dfData$values, 
                 gammaShape=0.5, gammaRate=1e-4)

fit.stan = sampling(stanDso, data=lStanData, iter=10000, chains=4,
                    pars=c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaRan3',
                           'nu', 'sigmaPop', 'mu',
                           'rGroupsJitter1', 'rGroupsJitter2', 'rGroupsJitter3'),
                    cores=4)#, init=initf, control=list(adapt_delta=0.99, max_treedepth = 12))

print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaRan3', 'sigmaPop', 'nu'), digits=3)
traceplot(fit.stan, 'betas')
traceplot(fit.stan, 'sigmaRan1')
traceplot(fit.stan, 'sigmaRan2')
traceplot(fit.stan, 'sigmaRan3')
m = cbind(extract(fit.stan)$sigmaRan1, extract(fit.stan)$sigmaRan2, extract(fit.stan)$sigmaRan3) 
dim(m)
m = log(m)
colnames(m) = c('Treatment', 'Patient', 'Adjustment')
pairs(m, pch=20, cex=0.5, col='grey')
library(lattice)
df = stack(data.frame(m))
histogram(~ values | ind, data=df, xlab='Log SD')


hist(dfData$values, prob=T)
plot(density(dfData$values))
mFitted = extract(fit.stan)$mu

apply(mFitted[sample(1:nrow(mFitted), size = 100), ], 1, function(x) lines(density(x)))

m = colMeans(mFitted)
r = dfData$values - m
plot(m, r, pch=20)
lines(lowess(m, r), col=2)

dfCovariates$fAdjustment = fNewBatch
write.csv(dfCovariates, file='dataExternal/richardSingleCell/covariates_2.csv')

###################################################################################### 
### second part of analysis fitting model 
mData = round(mData, 0)
mData.bk = mData
#mData = mData.bk[rownames(mData.bk) %in% cvGenes, ]
dim(mData)
dfData = data.frame(t(mData))
dfData = stack(dfData)
dim(dfData)
str(dfCovariates)
dfData$fBatch = factor(dfCovariates$treatment)
dfData$fAdjust1 = factor(dfCovariates$patientId)
dfData$fAdjust2 = factor(dfCovariates$fAdjustment)
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj1 = factor(dfData$fAdjust1:dfData$ind)
dfData$Coef.adj2 = factor(dfData$fAdjust2:dfData$ind)
dim(dfData)
dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef, dfData$Coef.adj1, dfData$Coef.adj2), ]
str(dfData)

# # ## setup the model
# library(lme4)
# fit.lme1 = glmer.nb(values ~ 1  + (1 | Coef) + (1 | Coef.adj1) + (1 | Coef.adj2), data=dfData)
# summary(fit.lme1)
# fit.lme2 = glmer.nb(values ~ 1 + (1 | Coef) + (1 | Coef.adj1), data=dfData)
# summary(fit.lme2)
# fit.lme3 = glmer.nb(values ~ 1 + (1 | Coef) + (1 | Coef.adj2), data=dfData)
# summary(fit.lme3)
# fit.lme4 = glmer.nb(values ~ 1 + (1 | Coef), data=dfData)
# summary(fit.lme4)
# 
# anova(fit.lme1, fit.lme2, fit.lme3, fit.lme4)
# 
# ran = ranef(fit.lme1, condVar=F)
# 
# plot((fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
# lines(lowess((fitted(fit.lme1)), resid(fit.lme1)), col=2)
# plot(density(dfData$values))
# lines(density(fitted(fit.lme1)))
# lines(density(fitted(fit.lme2)))
# lines(density(fitted(fit.lme3)))
# lines(density(fitted(fit.lme4)))
# 
# plot(resid(fit.lme1), dfData$fBatch)
# plot(resid(fit.lme1), dfData$fAdjust1)
# plot(resid(fit.lme1), dfData$fAdjust2)
# plot(resid(fit.lme1), dfData$ind)

# utility function to calculate gamma distribution used as a prior for scale
gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(log(dfData$values+0.5)), 2*sd(log(dfData$values+0.5)))

stanDso = rstan::stan_model(file='nbinomResp3RandomEffectsMultipleScales.stan')
## subset the data to get the second level of nested parameters
## this is done to avoid loops in the stan script to map the scale parameters
## of each ind/gene to the respective set of coefficients for jitters
d = dfData[!duplicated(dfData$Coef), ]

lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj1),
                 Nclusters3=nlevels(dfData$Coef.adj2),
                 NScaleBatches1 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj1),
                 NgroupMap3=as.numeric(dfData$Coef.adj2),
                 NBatchMap1=as.numeric(d$ind), # this is where we use the second level mapping
                 y=dfData$values, 
                 gammaShape=l$shape, gammaRate=l$rate)


# 
# fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=6,
#                     pars=c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaRan3',
#                            'iSize', #'mu',
#                            'rGroupsJitter1'),# 'rGroupsJitter2', 'rGroupsJitter3'),
#                     cores=6)#, control=list(adapt_delta=0.99, max_treedepth = 12))#, init=initf, 
# save(fit.stan, file='Temp/fit.stan.10Jan.rds')

print(fit.stan, c('betas', 'sigmaRan1[1]', 'sigmaRan2', 'sigmaRan3', 'iSize'), digits=3)
traceplot(fit.stan, 'betas')
traceplot(fit.stan, c('sigmaRan2', 'sigmaRan3'))
traceplot(fit.stan, c('sigmaRan1[1]'))
## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# ## get the intercept at population level
#iIntercept = as.numeric(extract(fit.stan)$intercept)
## add the intercept to each random effect variable, to get the full coefficient
#mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifferenceSim = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = min(c(sum(d <= 0)/length(d), sum(d >= 0)/length(d))) * 2
  #p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='MCSF', deflection='MUC1-ST') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifferenceSim(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
# ## some p-values are essentially zero, add a random jitter to it for plotting
# p = -1*log10(dfResults$pvalue+1e-200)
# p = runif(nrow(dfResults), p, p+50)
# p = 10^(-1*p)
# dfResults$P.Value = p
head(rownames(dfResults))
dfResults$SYMBOL = rownames(dfResults)
# ## produce the plots
# range(dfResults$logFC)
# f_plotVolcano(dfResults, 'MUC1-ST vs MCSF', 0.01, fc.lim=range(dfResults$logFC))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(log(m), dfResults, 0.01, '')
table(dfResults$pvalue < 0.05)
table(dfResults$adj.P.Val < 0.01)

# make sure names match as - is replaced by . in factor names
table(dfResults$SYMBOL %in% rownames(mData))
s = gsub('\\.', '-', dfResults$SYMBOL)
table(s %in% rownames(mData))
si = which(!(s %in% rownames(mData)))
s[si] = gsub('-', '.', s[si])
table(s %in% rownames(mData))
dfResults$SYMBOL = s

## add variance term
mVar = extract(fit.stan)$sigmaRan1
dim(mVar)
colnames(mVar) = levels(dfData$ind)

table(colnames(mVar) %in% rownames(mData))
s = gsub('\\.', '-', colnames(mVar))
table(s %in% rownames(mData))
si = which(!(s %in% rownames(mData)))
s[si] = gsub('-', '.', s[si])
table(s %in% rownames(mData))
colnames(mVar) = s
table(colnames(mVar) %in% dfResults$SYMBOL)
i = match(dfResults$SYMBOL, colnames(mVar))
mVar = mVar[,i]
# sanity check
identical(colnames(mVar), dfResults$SYMBOL)
dfResults$Variance = colMeans(mVar)

## save the results 
write.csv(dfResults, file='Temp/DEAnalysisMUC1-STvsMCSF.xls')

# load old results
dfResults.old = read.csv(file.choose(), header=T, stringsAsFactors = F)
#rownames(dfResults.old) = dfResults.old$Gene

table(dfResults$SYMBOL %in% dfResults.old$Gene)
dfResults.sub = dfResults[(dfResults$SYMBOL %in% dfResults.old$Gene), ]
# put in same order
i = match(dfResults.sub$SYMBOL, dfResults.old$Gene)
dfResults.old = dfResults.old[i,]
identical(dfResults.sub$SYMBOL, dfResults.old$Gene)

## plot fc
plot(dfResults.sub$logFC, log(dfResults.old$Ratio), pch=20, xlab='NB GLM FC', ylab='Older FC', sub='Missing 429 Genes in Older', main='Common 12744 Genes')
plot(dfResults.sub$pvalue, dfResults.old$P.value, pch=20)
plot(dfResults.sub$adj.P.Val, dfResults.old$FDR, pch=20, xlab='NB GLM FDR', ylab='Older FDR', sub='Missing 429 Genes in Older', main='Common 12744 Genes')
# coplot(dfResults.sub$adj.P.Val ~ dfResults.old$FDR | abs(dfResults.sub$logFC), pch=20,
#        ylab='NB GLM', xlab='Older')

coplot(dfResults.sub$logFC ~ log(dfResults.old$Ratio) | dfResults.sub$Variance, pch=20,
       ylab='NB GLM', xlab='Older')

# which genes are missing
table(dfResults$SYMBOL %in% dfResults.old$Gene)
dfResults.missing = dfResults[!(dfResults$SYMBOL %in% dfResults.old$Gene), ]

par(mfrow=c(2,2))
hist(dfResults.sub$logFC)
hist(dfResults.sub$Variance)

hist(dfResults.missing$logFC)
hist(dfResults.missing$Variance)

hist(dfResults.sub$pvalue)
hist(dfResults.missing$pvalue)

table(dfResults.sub$adj.P.Val < 0.01) # 17%
table(dfResults.missing$adj.P.Val < 0.01) # 25%

# compare the 2 datasets
oDiag.sub = CDiagnosticPlots(log(mData[rownames(mData) %in% dfResults.sub$SYMBOL, ]+0.5), 'Common')
oDiag.missing = CDiagnosticPlots(log(mData[rownames(mData) %in% dfResults.missing$SYMBOL, ]+0.5), 'Missing')
# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
str(dfCovariates)
fBatch = dfCovariates$patientId
fBatch2 = dfCovariates$treatment

#pdf('results/matrixClustering.pdf')

## check using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.sub, fBatch2, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.missing, fBatch2, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.sub, fBatch2, axis.label.cex = 0.7, legend.pos = 'topright')
plot.mean.summary(oDiag.missing, fBatch2, axis.label.cex = 0.7, legend.pos = 'topright')

plot.sigma.summary(oDiag.sub, fBatch2, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.missing, fBatch2, axis.label.cex = 0.7)

plot.missing.summary(oDiag.sub, fBatch2, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.missing, fBatch2, axis.label.cex = 0.7, cex.main=1)
## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.sub)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.sub = CDiagnosticPlotsSetParameters(oDiag.sub, l)
oDiag.missing = CDiagnosticPlotsSetParameters(oDiag.missing, l)

plot.PCA(oDiag.sub, fBatch2, legend.pos = 'right')
plot.PCA(oDiag.missing, fBatch2, legend.pos = 'right')

plot.dendogram(oDiag.sub, fBatch2, labels_cex = 1)
plot.dendogram(oDiag.missing, fBatch2, labels_cex = 1)

## order by fold changes
dfResults.sub.order = dfResults.sub[order(abs(dfResults.sub$logFC), decreasing = T),]
dfResults.order = dfResults[order(abs(dfResults$logFC), decreasing = T),]
dfResults.old.order = dfResults.old[order(abs(log(dfResults.old$Ratio)), decreasing = T),]
table(dfResults.sub.order$SYMBOL[1:100] %in% dfResults.old.order$Gene[1:100])
table(dfResults.order$SYMBOL[1:100] %in% dfResults.old.order$Gene[1:100])

## save the results 
write.csv(dfResults.order, file='Temp/DEAnalysisMUC1-STvsMCSF_ordered.xls')
write.csv(dfResults.missing, file='Temp/DEAnalysisMUC1-STvsMCSF_missing.xls')
write.csv(dfResults.old.order, file='Temp/old_ordered.xls')

## looking at top left region of plot
dfFDR = data.frame(sym=dfResults.sub$SYMBOL, mine=dfResults.sub$pvalue, their=dfResults.old$P.value, my.fc=dfResults.sub$logFC, their.fc=log(dfResults.old$Ratio))
str(dfFDR)
head(dfFDR)
dfFDR.topleft = dfFDR[dfFDR$mine < 0.01 & dfFDR$their > 0.8, ]
dfFDR.bottomright = dfFDR[dfFDR$mine > 0.8 & dfFDR$their < 0.1, ]

dfFC.topright = dfFDR[dfFDR$my.fc > 5, ]


