# leoniGeneComparisons.R

table(grepl(pattern = genesIn[1], week.3$gene, ignore.case = T))

s = sapply(seq_along(genesIn), function(x) grep(pattern = paste0("^", genesIn[x], "$"), week.3$gene, ignore.case = T))
week.3.sub = week.3[unlist(s),]

s = sapply(seq_along(genesIn), function(x) grep(pattern = paste0("^", genesIn[x], "$"), week.8$gene, ignore.case = T))
week.8.sub = week.8[unlist(s),]

s = sapply(seq_along(genesIn), function(x) grep(pattern = paste0("^", genesIn[x], "$"), week.11$gene, ignore.case = T))
week.11.sub = week.11[unlist(s),]

identical(week.3.sub$gene, week.8.sub$gene)
identical(week.3.sub$gene, week.11.sub$gene)

w3 = data.frame(wt=week.3.sub$value_1, tg=week.3.sub$value_2)
w3 = stack(w3)
w3$gene = week.3.sub$gene
w3$time = 3

w8 = data.frame(wt=week.8.sub$value_1, tg=week.8.sub$value_2)
w8 = stack(w8)
w8$gene = week.8.sub$gene
w8$time = 8

w11 = data.frame(wt=week.11.sub$value_1, tg=week.11.sub$value_2)
w11 = stack(w11)
w11$gene = week.11.sub$gene
w11$time = 11

dfData = rbind(w3, w8, w11)
dfData$values = log(dfData$values+0.5)
library(lattice)
xyplot(values ~ time | gene, data=dfData, groups=ind, pch=20, col=1:2,
       key=list(columns=2, text=list(levels(dfData$ind)), points=list(pch=20, col=1:2)),
       type=c('g', 'r', 'p'), ylab='Log Average', xlab='Time Weeks')

dfExport = rbind(week.3.sub, week.8.sub, week.11.sub)
write.csv(dfExport, file='dataExternal/leoni/selectedGenes.csv')
