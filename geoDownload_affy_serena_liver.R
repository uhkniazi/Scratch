# Name: geoDownload_affy_serena_liver.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 29/4/2022
# Desc: Download the relevant datasets, normalise and perform EDA

source('header.R')

library(GEOquery)
library(downloader)

## open the soft format and raw data
gse =  getGEO(filename = 'dataExternal/serena/GSE84954_series_matrix.txt.gz')

# ## read raw data CEL files and normalize
# library(affy)
# setwd('dataExternal/serena/')
# untar('GSE84954_RAW.tar')
# oData = ReadAffy()
# setwd(gcsWD)
# # normalize the data
# x.affy = rma(oData)

# get the samples from the expression set object
dfSamples = pData(gse)
dim(exprs(gse))
identical(rownames(dfSamples), colnames(exprs(gse)))

# # col names of the affy expression matrix
# cn = colnames(exprs(x.affy))
# # remove the last .CEL.gz from the names
# cn = gsub('^(GSM\\d+).CEL\\.gz', replacement = '\\1', cn, perl = T)

# order the sample names according to cel files
# table(cn %in% rownames(dfSamples))
# i = match(cn, rownames(dfSamples))
# dfSamples = dfSamples[i,]
# ## sanity check
# identical(rownames(dfSamples), cn)

# ## complete the expression set object by adding this sample information
# pData(x.affy) = dfSamples
# colnames(exprs(x.affy)) = cn
# # sanity check
# identical(colnames(x.affy), rownames(pData(x.affy)))

# ## annotation data
# fData(x.affy) = fData(gse)
# ## store the grouping factors in affy object
# x.affy$fCondition = dfSamples$group1
# x.affy$fGender = dfSamples$group2
# x.affy$age = as.numeric(as.character(dfSamples$group3))

## check normalisation 
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

mCounts = exprs(gse)
range(mCounts)

oDiag.1 = CDiagnosticPlots(log(mCounts), 'Soft Log format')
oDiag.2 = CDiagnosticPlots(mCounts, 'Soft format')
colnames(dfSamples)
fDisease = factor(gse$`disease:ch1`)
fSubject = factor(gse$`subjectid:ch1`)
fTissue = factor(gse$`tissue:ch1`)
levels(fDisease); levels(fSubject); levels(fTissue)

xtabs( ~ fDisease + fSubject + fTissue)
fBatch = fTissue
levels(fBatch)
## compare the 2 methods using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.PCA(oDiag.2, fBatch, cex.main=1)

plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

## this should give an error as scaling can't be done
## if all the vector 0 for PCA
oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)
oDiag.2.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)

plot.PCA(oDiag.1.2, fBatch, legend.pos = 'topright')
plot.PCA(oDiag.2.2, fBatch, legend.pos = 'topright')

plot.dendogram(oDiag.1.2, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2.2, fBatch, labels_cex = 0.8, cex.main=0.7)

par(mfrow=c(1,1))
plot.PCA(oDiag.1.2, fBatch, legend.pos = 'topleft')
plot.PCA(oDiag.2.2, fBatch, legend.pos = 'bottomright')

### use the data from each tissue separately and in log format
iLiver = which(fTissue == 'liver')
gse.Liver = gse[, iLiver]

### repeat the analysis on only one tissue
mCounts = log(exprs(gse.Liver))
range(mCounts)

oDiag.L = CDiagnosticPlots(mCounts, 'Liver: log Soft Format')
fSubject = factor(gse.Liver$`subjectid:ch1`)
fDisease = factor(gse.Liver$`disease:ch1`)
levels(fDisease); levels(fSubject)

xtabs( ~ fDisease + fSubject)
fBatch = fDisease
levels(fBatch)

l = CDiagnosticPlotsGetParameters(oDiag.L)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.L = CDiagnosticPlotsSetParameters(oDiag.L, l)
plot.PCA(oDiag.L, fBatch, csLabels = as.character(fSubject), cex=1.5)

### muscle tissue
### use the data from each tissue separately and in log format
iMuscle = grep('muscle', gse$`tissue:ch1`)
gse.Muscle = gse[, iMuscle]

### repeat the analysis on only one tissue
mCounts = log(exprs(gse.Muscle))
range(mCounts)

oDiag.M = CDiagnosticPlots(mCounts, 'Muscle: log Soft Format')
fSubject = factor(gse.Muscle$`subjectid:ch1`)
fDisease = factor(gse.Muscle$`disease:ch1`)
levels(fDisease); levels(fSubject)

xtabs( ~ fDisease + fSubject)
fBatch = fDisease
levels(fBatch)

l = CDiagnosticPlotsGetParameters(oDiag.M)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.M = CDiagnosticPlotsSetParameters(oDiag.M, l)
plot.PCA(oDiag.M, fBatch, csLabels = as.character(fSubject), cex=1.5)

### fat tissue
### use the data from each tissue separately and in log format
iFat = grep('fat', gse$`tissue:ch1`)
gse.Fat = gse[, iFat]

### repeat the analysis on only one tissue
mCounts = log(exprs(gse.Fat))
range(mCounts)

oDiag.F = CDiagnosticPlots(mCounts, 'Fat: log Soft Format')
fSubject = factor(gse.Fat$`subjectid:ch1`)
fDisease = factor(gse.Fat$`disease:ch1`)
levels(fDisease); levels(fSubject)

xtabs( ~ fDisease + fSubject)
fBatch = fDisease
levels(fBatch)

l = CDiagnosticPlotsGetParameters(oDiag.F)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.F = CDiagnosticPlotsSetParameters(oDiag.F, l)
plot.PCA(oDiag.F, fBatch, csLabels = as.character(fSubject), cex=1.5)


