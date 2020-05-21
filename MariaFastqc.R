# Name: MariaFastqc.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 19/5/2020
# Desc: fastqc for some files

library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CFastqQuality/experimental/CFastqQuality.R'
download(url, 'CFastqQuality.R')

# load the required packages
source('CFastqQuality.R')
# delete the file after source
unlink('CFastqQuality.R')

dfMetaData = read.csv('dataExternal/AO-3370-S331-S369_sample_attributes.txt', 
                      header=T, sep='\t')
## format the metadata
s = strsplit(as.character(dfMetaData$fastq), ',')
class(s)
s[[1]]
df1 = dfMetaData
df1$fastq = sapply(s, function(x) return(x[1]))
head(df1)
df1$read_direction = '1'
df2 = dfMetaData
df2$fastq = sapply(s, function(x) return(x[2]))
df2$read_direction = '2'
dfMetaData = rbind(df1, df2)

dfMetaData$read_direction = factor(dfMetaData$read_direction)
str(dfMetaData)
dfMetaData$fastq = gsub(' ', '', dfMetaData$fastq)
dfMetaData$fastq = gsub('1\\.', '1_subsampled.', dfMetaData$fastq)
## load fastq file names
paths = dir('dataExternal/Partek_AO-3370-S331-S369_Subsample_reads_Subsampled_reads/',
            pattern='*.fastq.gz', full=F)

## check if file names are identical to the names from meta data sheet
table(paths %in% as.character(dfMetaData$fastq))

## set working directory to directory with fastq files
cCurDir = getwd()
setwd('dataExternal/Partek_AO-3370-S331-S369_Subsample_reads_Subsampled_reads/')

paths = as.character(dfMetaData$fastq)

# set factors required for analysis
fReadDirection = factor(dfMetaData$read_direction)
cNames = paste0(as.character(dfMetaData$Sample.name), '_', as.character(dfMetaData$read_direction))
# any additional factors that need including in the metadata list
lMetaData = list(meta=dfMetaData)

## create object, go have a cup of coffee while this runs
ob = CFastqQualityBatch(paths, cNames, fReadDirection, lMetaData)
setwd(cCurDir)
save(ob, file='maria.rds')

## get read counts in millions
iGetReadCount(ob)
barplot.readcount(ob)
plot.alphabetcycle(ob)
plot.qualitycycle(ob)

######### some additional diagnostic plots on the data matrix
### some diagnostic plots
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

## extract the base quality matrix 
mBatch = mGetReadQualityByCycle(ob)
dim(mBatch)
mBatch[1:10, 1:4]

## creat an object of diagnostics class to make plots
oDiag = CDiagnosticPlots(mBatch, 'Base Quality')

## turning off automatic jitters
## we set jitter to FALSE for PCA, otherwise, in sparse matrix a random jitter is added to avoid divisions by zeros
l = CDiagnosticPlotsGetParameters(oDiag)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag = CDiagnosticPlotsSetParameters(oDiag, l)

fBatch = ob@fReadDirection
str(ob@lMeta$meta)
fBatch = factor(ob@lMeta$meta$Lane)
## try some various factors to make the plots of low dimensional summaries
plot.mean.summary(oDiag, fBatch)
plot.sigma.summary(oDiag, fBatch)
boxplot.median.summary(oDiag, fBatch)
plot.PCA(oDiag, fBatch, csLabels = '')
plot.dendogram(oDiag, fBatch, labels_cex = 0.8)

## looking at alphabets 
## change direction and alphabet i.e. base as required
i = grep('1', ob@fReadDirection)

lAlphabets = lapply(i, function(x){
  m = t(mGetAlphabetByCycle(ob@lData[[x]]))
  m = m[,c('A', 'T', 'G', 'C')]
  r = rowSums(m)
  m = sweep(m, 1, r, '/')
  return(m)
})

mAlphabet = do.call(cbind, lapply(lAlphabets, function(x) return(x[,'A'])))
dim(mAlphabet)
colnames(mAlphabet) = (as.character(ob@lMeta$meta$Sample.name[ob@lMeta$meta$read_direction == '1']))
oDiag.2 = CDiagnosticPlots(mAlphabet, 'forward base A')

## turning off automatic jitters
## we set jitter to FALSE for PCA, otherwise, in sparse matrix a random jitter is added to avoid divisions by zeros
l = CDiagnosticPlotsGetParameters(oDiag.2)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)

str(ob@lMeta$meta)
fBatch = factor(ob@lMeta$meta$Lane[ob@lMeta$meta$read_direction == '1'])
# f2 = factor(ob@lMeta$meta$Gender)
# fBatch = fBatch:f2
## try some various factors to make the plots of low dimensional summaries
plot.PCA(oDiag.2, fBatch)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8)
