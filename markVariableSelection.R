# some testing of variable selection on data from mark

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

dfData = read.csv(file.choose(), header = T, na.strings = c('na', 'NA', 'NaN'))
dfData = t(dfData)
dfData = data.frame(dfData)
cn = t(dfData[1,])
colnames(dfData) = cn[,1]
dfData = dfData[-1,]

i = which(dfData$`Class ID` %in% c('ACLF', 'AD', 'SC'))
length(i)
dfData = dfData[i,]
dim(dfData)
dfData = droplevels.data.frame(dfData)

write.csv(dfData, '~/Downloads/FILES_FOR_BRC_ECR/LNEG LCMS WITH METADATA Transposed mass spec.csv', row.names = F)

dfData = read.csv('~/Downloads/FILES_FOR_BRC_ECR/LNEG LCMS WITH METADATA Transposed mass spec.csv', header = T, na.strings = c('na', 'NA', 'NaN'))
fGroups = factor(dfData$X90.day)

dfData = dfData[,-c(1:94)]
gc(reset = T)

dfData = na.omit(dfData)
mDat = log(dfData+1)
mDat = t(mDat)

url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(mDat, 'lneg')

# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
# e.g. in this case it is different lanes/machines
fBatch = fGroups

## compare the 2 methods using various plots
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)

plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)
plot.PCA(oDiag.1, fBatch)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)

# calculate the coef of var for each gene
cv = apply(mDat, 1, function(x) sd(x)/abs(mean(x)))
# check cv
summary(cv)

# cut data into groups based on quantiles of cv
cut.pts = quantile(cv, probs = 0:10/10)
groups = cut(cv, breaks = cut.pts, include.lowest = T)
groups = cut(cv, breaks = cut.pts, include.lowest = T, labels = 0:9)
iMean = apply(mDat, 1, mean)
iVar = apply(mDat, 1, var)
coplot((cv) ~ iMean | groups)
coplot(iVar ~ iMean | groups)
# choose genes with small cv
# f = cv < 0.2
# choosing groups from quantile 0 to 40
mDat = mDat[groups %in% c(0, 1, 2, 3),]

# select a subset of genes that show differential expression
p.t = apply(mDat, 1, function(x) t.test(x ~ fGroups)$p.value)
p.t.adj = p.adjust(p.t, 'BH')
p.t.adj = sort(p.t.adj, decreasing = F)
t = names(p.t.adj[1:200])
mDat.sub = mDat[rownames(mDat) %in% t,]

# set variables
dfData.bk = dfData
dfData = dfData[,rownames(mDat.sub)]
dim(dfData)
dfData = na.omit(dfData)
dim(dfData)

## perform nested random forest
## adjust boot.num as desired
oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)
# plot the top 20 genes based on importance scort with 95% confidence interval for standard error
plot.var.selection(oVar.r)
# get the variables
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:10]

# use the top 30 genes to find top combinations of genes
dfData = dfData[,colnames(dfData) %in% cvTopGenes]

## look at colinear variables
## find correlated variables
m = NULL;

for (i in 1:ncol(dfData)){
  m = cbind(m, dfData[,i])
}
colnames(m) = colnames(dfData)
mCor = cor(m, use="na.or.complete")
library(caret)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.7, names=T)
data.frame(n)
sapply(n, function(x) {
  (abs(mCor[,x]) >= 0.7)
})

cvKeep = c('X1091.89_623.234')
n = n[!(n%in% cvKeep)]
i = which(colnames(dfData) %in% n)
cn = colnames(dfData)[-i]

dfData.bk2 = dfData
dfData = dfData[,cn]

oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 50)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# print variable combinations
for (i in 1:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

## 10 fold nested cross validation with various variable combinations
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 1:4){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = data.frame(dfData[,cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  
  dfData.test = dfData.train
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = fGroups,
                             train.groups = fGroups, level.predict = '0', boot.num = 500)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}