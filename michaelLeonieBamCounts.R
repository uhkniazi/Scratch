# checks on bam file with coverage issues

## load the transcript db objects
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicAlignments)
library(rtracklayer)

csFile = file.choose()

# get the exons into GRangesList object
oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')

# for each of these bam file lists do the counting
oBamFiles = BamFileList(csFile, index=paste0(csFile, '.bai'))
mCounts = (assays(summarizeOverlaps(oGRLgenes, oBamFiles, ignore.strand = F, singleEnd=F))$counts)
save(mCounts, file='Temp/mCounts.rds')

library(org.Hs.eg.db)
columns(org.Hs.eg.db)

##### 2 - find position of gene of interest
df = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(mCounts), keytype = 'ENTREZID', columns = 'SYMBOL')

identical(rownames(mCounts), df$ENTREZID)
df = cbind(df, mCounts)
head(df)
