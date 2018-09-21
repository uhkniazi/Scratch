# geoDownload.R
# downloading some geo data sets for testing

library(GEOquery)
library(downloader)

## open the soft format and raw data
gse = getGEO('GSE114651')

lFiles = list.files('~/Downloads/temp/GSE114651_RAW/')

df = lapply(lFiles, function(x) read.csv(x, header=T,
                                         sep=' ', row.names=1))
df = do.call(cbind, df)

df = as.matrix(df)
class(gse$GSE114651_series_matrix.txt.gz)
exprs(gse$GSE114651_series_matrix.txt.gz) = df
p = pData(gse$GSE114651_series_matrix.txt.gz)

oExp = ExpressionSet(df)
pData(oExp) = p

cn = colnames(exprs(oExp))
dfSamples = pData(oExp)
# order the sample names according to pheno data table
table(cn %in% as.character(dfSamples$title))
i = match(cn, as.character(dfSamples$title))
dfSamples = dfSamples[i,]
## sanity check
identical(as.character(dfSamples$title), cn)
rownames(dfSamples) = dfSamples$title
## complete the expression set object by adding this sample information
pData(oExp) = dfSamples
# sanity check
identical(colnames(oExp), rownames(pData(oExp)))
save(oExp, file='oExp_raw_GSE114651.rds')
###########################################################

# create database entries
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and metafile table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe MetaFile;'))$Field[-1]

## create the entry for samples
cSampleCol
g_pid = 9
g_did = 16
title = paste(as.character(dfSamples$title), rownames(dfSamples))
description = paste(as.character(dfSamples$source_name_ch1), as.character(dfSamples$characteristics_ch1))
group1 = gsub('fibroblasts_(\\w+)_patient\\d+', '\\1', as.character(dfSamples$title))
group2 = gsub('(female|male).+', '\\1', as.character(dfSamples$characteristics_ch1))
group3 = gsub('\\D+', '', as.character(dfSamples$characteristics_ch1))

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=title, description=description, group1=group1, 
                       group2=group2, group3=group3)
# write this data to the database
rownames(dfSamples) = NULL

# # write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)

## annotation data
fData(x.affy) = fData(gse)
## store the grouping factors in affy object
x.affy$fCondition = dfSamples$group1
x.affy$fGender = dfSamples$group2
x.affy$age = as.numeric(as.character(dfSamples$group3))

## check normalisation 
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(log(exprs(gse)+1), 'Soft format')
oDiag.2 = CDiagnosticPlots(exprs(x.affy), 'Raw normalised')

fBatch = x.affy$fCondition

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

### save the x.affy object and create metafile entry
getwd()
n = make.names(paste('GSE7669 Synovial fibroblasts, RA versus OA.rds'))
n2 = paste0('~/Data/MetaData/', n)
save(x.affy, file=n2)

## access the mysql database to create entry for this data object
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='GSE7669 Synovial fibroblasts, RA versus OA, normalised object')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)
