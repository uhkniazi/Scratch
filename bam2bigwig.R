# bam2bigwig.R
# testing conversion of bam to bigwig files

(!require(Rsamtools) || !require(GenomicAlignments))
library(rtracklayer)

# load a sample bamfile
csFile = system.file('extdata', 'ex1.bam', package='Rsamtools')

# load the bam file as GAlignment object
oGABam = readGAlignments(csFile)

cov1 = coverage(oGABam)

oGRBam = as(oGABam, 'GRanges')

cov2 = coverage(oGRBam)

## check resize
ga = oGABam[1:5]
gr = oGRBam[1:5]

gr.50 = resize(gr, 50)
cov.50 = coverage(gr.50)
cov.gr = coverage(gr)
# resizing can mess things up
as.vector(cov.50$seq1)[1:60]
as.vector(cov.gr$seq1)[1:60]

## as is and normalised
tr = length(oGABam)
cov1
cov1.rpm = lapply(cov1, function(x) signif(10^6*x/tr,3))
cov1.rpm = as(cov1.rpm, 'SimpleRleList')

export.bw(cov1, con='test.bw')
