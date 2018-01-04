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

dfData = dfData[,1:93]
gc(reset = T)

str(dfData)
write.csv(dfData, '~/Downloads/FILES_FOR_BRC_ECR/LNEG LCMS WITH METADATA Transposed.csv', row.names = F)

dfData = read.csv('~/Downloads/FILES_FOR_BRC_ECR/LNEG LCMS WITH METADATA Transposed.csv', header = T, na.strings = c('na', 'NA', 'NaN'))

s = sapply(1:ncol(dfData), function(x) sum(is.na(dfData[,x])))
s
names(s) = colnames(dfData)
s


# set variables
dfData.bk = dfData
dfData = dfData[,-c(1:3)]
dfData = na.omit(dfData)
dim(dfData)
fGroups = factor(dfData$X90.day)
dfData = dfData[,-(which(colnames(dfData) %in% c('X90.day', 'X30.day', 'X1.Year')))]

## perform nested random forest
## adjust boot.num as desired
oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 50)
# plot the top 20 genes based on importance scort with 95% confidence interval for standard error
plot.var.selection(oVar.r)
# get the variables
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:30]

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

cvKeep = c('MELD', 'Neutrophil', 'MAP')
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