# quick rocr and summary data for rob

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

dfData = read.csv(file.choose(), header=T)
dfData = dfData[dfData$Condition != 'TB',]
dfData = droplevels.data.frame(dfData)
str(dfData)

fGroups = dfData$Condition
dfData = dfData[,-5]

dim(dfData)
table(fGroups)

# create the cross validation object
oCV = CCrossValidation.LDA(dfData, dfData, fGroups, fGroups, level.predict = 'SA',
                           boot.num = 100, k.fold = 10) 

plot.cv.performance(oCV)

## get a table of the cutoffs for the cross validation
dfCVCutoffs = getCutoffTprFprCrossValidation(oCV)[,-4]
dfCVCutoffs = dfCVCutoffs[,-4]


## get without cross validation
p = oCV@oPerf.val

dfCutoffWithoutCV = data.frame(cutoff=p@alpha.values[[1]], tpr=p@y.values[[1]], fpr=p@x.values[[1]])

pdf('Temp/rob_plots.pdf')
par(family='Helvetica')
plot.cv.performance(oCV)
dev.off(dev.cur())

write.csv(dfCVCutoffs, file='Temp/rob1.csv')
write.csv(dfCutoffWithoutCV, file='Temp/rob2.csv')