dfData = read.csv('Temp/DeProteinsJens.csv', header=F, stringsAsFactors = F)
str(dfData)
library(GO.db)
columns(GO.db)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
dfEntrezID = AnnotationDbi::select(org.Hs.eg.db, keys = dfData$V1, columns = c('SYMBOL', 'ENTREZID', 'GO'), keytype = 'SYMBOL')
columns(org.Hs.eg.db)
columns(GO.db)
dfGO = AnnotationDbi::select(GO.db, dfEntrezID$GO, columns=c('DEFINITION', 'TERM'), keytype = 'GOID')
identical(dfEntrezID$GO, dfGO$GOID)
dfEntrezID$TERM = dfGO$TERM
dfEntrezID$DEFINITION = dfGO$DEFINITION
write.csv(dfEntrezID, file='Temp/GoAnnodations.csv')
library(GOstats)
# get the universe of genes with go terms
univ = keys(org.Hs.eg.db, 'ENTREZID')
dfUniv = AnnotationDbi::select(org.Hs.eg.db, keys = univ, columns = c('GO'), keytype = 'ENTREZID')
dim(dfUniv)
dfUniv = na.omit(dfUniv)
dim(dfUniv)
univ = unique(dfUniv$ENTREZID)
length(univ)
## make hypergeometric test object for each type, CC, BP and MF
params = new('GOHyperGParams', geneIds=unique(dfEntrezID$ENTREZID),
             annotation='org.Hs.eg.db',
             universeGeneIds=univ,
             ontology='BP',
             pvalueCutoff= 0.01,
             conditional=FALSE,
             testDirection='over')

oGOStat = hyperGTest(params) 
# get pvalues
ivPGO = pvalues(oGOStat)
# fdr
ivPGO.adj = p.adjust(ivPGO, 'BH')

cvSigGO.ID = names(ivPGO.adj[ivPGO.adj < 0.01])
dfSigGO = summary(oGOStat)

dfSigGO.1FDR = dfSigGO[dfSigGO$GOBPID %in% cvSigGO.ID,]
write.csv(dfSigGO, 'Temp/goreport_BiologicalProcess.csv')
htmlReport(oGOStat, 'Temp/goreport_BiologicalProcess.html')

## molecular function
params = new('GOHyperGParams', geneIds=unique(dfEntrezID$ENTREZID),
             annotation='org.Hs.eg.db',
             universeGeneIds=univ,
             ontology='MF',
             pvalueCutoff= 0.01,
             conditional=FALSE,
             testDirection='over')

oGOStat = hyperGTest(params) 
# get pvalues
ivPGO = pvalues(oGOStat)
# fdr
ivPGO.adj = p.adjust(ivPGO, 'BH')

cvSigGO.ID = names(ivPGO.adj[ivPGO.adj < 0.01])
dfSigGO = summary(oGOStat)

dfSigGO.1FDR = dfSigGO[dfSigGO$GOBPID %in% cvSigGO.ID,]
write.csv(dfSigGO, 'Temp/goreport_MolecularFunction.csv')
htmlReport(oGOStat, 'Temp/goreport_MolecularFunction.html')

## CC
params = new('GOHyperGParams', geneIds=unique(dfEntrezID$ENTREZID),
             annotation='org.Hs.eg.db',
             universeGeneIds=univ,
             ontology='CC',
             pvalueCutoff= 0.01,
             conditional=FALSE,
             testDirection='over')

oGOStat = hyperGTest(params) 
# get pvalues
ivPGO = pvalues(oGOStat)
# fdr
ivPGO.adj = p.adjust(ivPGO, 'BH')

cvSigGO.ID = names(ivPGO.adj[ivPGO.adj < 0.01])
dfSigGO = summary(oGOStat)

dfSigGO.1FDR = dfSigGO[dfSigGO$GOBPID %in% cvSigGO.ID,]
write.csv(dfSigGO, 'Temp/goreport_CellularCompartment.csv')
htmlReport(oGOStat, 'Temp/goreport_CellularCompartment.html')
