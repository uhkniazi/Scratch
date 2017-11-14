library(GenomicRanges)
r1 = IRanges(c(1, 10. 20). width=c(5, 12, 10))
r1 = IRanges(c(1, 10. 20). width=c(5, 12, 10)))
r1 = IRanges(start=c(1, 10. 20), width=c(5, 12, 10))
r1 = IRanges(start=c(1, 10, 20), width=c(5, 12, 10))
r1
r2 = IRanges(start=c(5, 15, 25), width=c(5, 12, 10))
r2
union(r1, r2)
r1
r2
intersect(r1, r2)
reduce(r1, r2)
reduce(r1)
rl = IRangesList(r1, r2)
rl
reduce(rl)
union(rl)
IRanges::union(rl, rl)
reduce(r2)
rl
union(r1, r2, r3)
union(r1, r2, r1)
union(r1, r2)
reduce(r1)
reduce(unlist(rl))

