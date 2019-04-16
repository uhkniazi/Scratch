# downloading some geo data for mansoor

library(GEOquery)
library(Biobase)
library(lumi)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(annotate)
library(limma)

# load the data, clean and create factors
dir.create('dataExternal/sepsis', showWarnings = F)
gse =  getGEO('GSE54514', GSEMatrix = T, destdir = 'dataExternal/sepsis/', AnnotGPL = T, getGPL = T)
oExp = gse$GSE54514_series_matrix.txt.gz

# add lumi nuIDs 
oExp.lumi = addNuID2lumi(oExp, lib.mapping = 'lumiHumanIDMapping' )
# remove any NA data
exprs(oExp.lumi) = na.omit(exprs(oExp.lumi))

# get the grouping factor
fSamples = as.character(pData(oExp.lumi)$source_name_ch1)
# comparing sepsis_nonsurvivor over days
# remove control and survivors
i = grep('Control|sepsis_survivor', fSamples)
oExp.lumi = oExp.lumi[,-i]
dfSamples = pData(oExp.lumi)

fSamples = as.character(dfSamples$source_name_ch1)
# get the covariates
fDays = gsub('[\\w+ ,]+Day (\\d+).+', '\\1', x = fSamples, perl = T)
fGender = as.character(dfSamples$characteristics_ch1.4)
fGender = gsub('\\gender: (\\w)', '\\1', fGender)
fAge = as.character(dfSamples$characteristics_ch1.5)
iAge = as.numeric(gsub('age \\(years\\): (\\d+)', '\\1', fAge))
fAge = cut(iAge, quantile(iAge, 0:5/5), include.lowest = T, labels = paste0('Age', 1:5))

# normalize the data
lumi.n = lumiN(oExp.lumi, method = 'rsn')

library(downloader)

## check data matrix
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(exprs(lumi.n), 'lumi norm')

fBatch = fAge
fBatch = factor(fGender)
fBatch = factor(fDays)
levels(fBatch)
## compare the 2 methods using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

## if all the vector 0 for PCA
oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)

plot.PCA(oDiag.1.2, fBatch, legend.pos = 'topright')

plot.dendogram(oDiag.1.2, fBatch, labels_cex = 0.8, cex.main=0.7)

fDays = factor(fDays)
fGender = factor(fGender)

cvSym = getSYMBOL(rownames(lumi.n), 'lumiHumanAll.db')
table(is.na(cvSym))
dim(lumi.n)
lumi.n = lumi.n[!is.na(cvSym),]
dim(lumi.n)

mDat = exprs(lumi.n)
### lmer library and function
library(lme4)

f_get.lme.sd = function(x){
  f = lmer(x ~ 1 + iAge + (1 | fDays) + (1 | fGender))
  f2 = summary(f)
  return(as.data.frame(f2$varcor)[,5])
}

mVar = apply(mDat, 1, f_get.lme.sd)
rownames(mVar) = c('fDays', 'fGender', 'residual')
mVar = t(mVar)
identical(rownames(mVar), rownames(mDat))

dfResults = data.frame(round(mDat,3), round(mVar,3), symbol=getSYMBOL(rownames(mDat), 'lumiHumanAll.db'))
write.csv(dfResults, file='dataExternal/sepsis/dfResults.csv')

## get top 2k
dfResults = dfResults[order(dfResults$fDays, decreasing = T), ]

mDat = mDat[rownames(dfResults)[1:2000],]
dim(mDat)
rownames(mDat) = as.character(dfResults$symbol[1:2000])
lData.sepsis = list(data=mDat, grouping=fDays, adjust.1=iAge, adjust.2=fGender)
save(lData.sepsis, file='dataExternal/sepsis/lData.sepsis.rds')
