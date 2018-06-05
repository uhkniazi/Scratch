# checks on bam file with coverage issues

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CBamQuality/master/CBamQuality.R'
download(url, 'CBamQuality.R')

# load the required packages
source('CBamQuality.R')
# delete the file after source
unlink('CBamQuality.R')

csFile = file.choose()

#### logic
## 1- load the list of annotated genes
## 2- find the position of gene of interest
## 3- read the bam file over this gene
## 4- show coverage over this gene

###### 1 - load the list of annotated genes
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
oGRgene = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
length(oGRgene)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

##### 2 - find position of gene of interest
df = AnnotationDbi::select(org.Hs.eg.db, keys = 'FCGR3A', columns = 'ENTREZID', keytype = 'SYMBOL')
oGRgene.sel = oGRgene['2214']

##### 3 - read the bam file and load the specific gene position
bf = BamFile(csFile)
seqinfo(bf)

which = oGRgene.sel
flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE, 
                   hasUnmappedMate = FALSE, isMinusStrand = NA, isMateMinusStrand = NA,
                   isFirstMateRead = NA, isSecondMateRead = NA, isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
                   isDuplicate = FALSE)

param = ScanBamParam(flag=flag, what = scanBamWhat(), which=which)
# read the GAlignments object
oGA = readGAlignmentPairs(bf, param=param)

# load the trainscript for this gene
# get the exons into GRangesList object
oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')
oGRLgenes = oGRLgenes['2214']

#### 4- show coverage over this gene
# get the coverage over this transcript in the bam file
cov = coverageByTranscript(oGA, oGRLgenes, ignore.strand=FALSE)

# this will be an rle object
cov = cov$`2214`
cov

# 
# # display the coverage over each exon
# oViews = Views(cov, start = start(reduce(unlist(oGRLgenes))), end = end(reduce(unlist(oGRLgenes))))

## number of reads aligned 
cov2 = assays(summarizeOverlaps(oGRLgenes, bf, ignore.strand = F, singleEnd=F,
                          param=param))$counts


