# checks on bam file with coverage issues

## load the transcript db objects
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicAlignments)
library(rtracklayer)

setwd('dataExternal/')
csFile = list.files('.', pattern = '*.bam$', recursive = T)

# get the exons into GRangesList object
oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')

# for each of these bam file lists do the counting
# for each of these bam file lists do the counting
lCounts = lapply(csFile, function(bfl){
  ## create a bamfiles list object
  oBamFiles = BamFileList(bfl, index=paste0(bfl, '.bai'))
  return(assays(summarizeOverlaps(oGRLgenes, oBamFiles, ignore.strand = F, singleEnd=F))$counts)
})

save(lCounts, file='../Temp/lCounts.rds')

mCounts = do.call(cbind, lCounts)

library(org.Hs.eg.db)
columns(org.Hs.eg.db)

##### 2 - find position of gene of interest
df = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(mCounts), keytype = 'ENTREZID', columns = 'SYMBOL')

identical(rownames(mCounts), df$ENTREZID)
df = cbind(df, mCounts)
head(df)

write.csv(df, file='../Temp/counts_all.csv')
