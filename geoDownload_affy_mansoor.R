# downloading some geo data for mansoor

library(GEOquery)
library(downloader)

## open the soft format and raw data
gse =  getGEO(filename = 'dataExternal/mansoor/GSE47908_series_matrix.txt.gz')

## read raw data CEL files and normalize
library(affy)
setwd('dataExternal/mansoor/')
untar('GSE47908_RAW.tar')
oData = ReadAffy()
setwd('../../')
# normalize the data
x.affy = rma(oData)

# get the samples from the expression set object
dfSamples = pData(gse)

# col names of the affy expression matrix
cn = colnames(exprs(x.affy))
# remove the last .CEL.gz from the names
cn = gsub('^(GSM\\d+)_.+.CEL\\.gz', replacement = '\\1', cn, perl = T)

# order the sample names according to cel files
table(cn %in% rownames(dfSamples))
i = match(cn, rownames(dfSamples))
dfSamples = dfSamples[i,]
## sanity check
identical(rownames(dfSamples), cn)

## complete the expression set object by adding this sample information
pData(x.affy) = dfSamples
colnames(exprs(x.affy)) = cn
# sanity check
identical(colnames(x.affy), rownames(pData(x.affy)))

## annotation data
fData(x.affy) = fData(gse)
str(dfSamples)

## check normalisation 
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(log(exprs(gse)+1), 'Soft format')
oDiag.2 = CDiagnosticPlots(exprs(x.affy), 'Raw normalised')

fBatch = as.character(dfSamples$characteristics_ch1.1)
fBatch = gsub('disease state: ', '', fBatch)
fBatch[fBatch == 'ulcerative colitis-associated dysplasia'] = 'UC'
fBatch = factor(fBatch)
levels(fBatch)
## compare the 2 methods using various plots
par(mfrow=c(1,1))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

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

### save the x.affy object and create metafile entry
getwd()
n = make.names(paste('GSE47908.rds'))
n2 = paste0('dataExternal/mansoor/', n)
save(x.affy, file=n2)
