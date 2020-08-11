# Name: anton_simulation.R
# Date: 11/8/2020
# Auth: umar.niazi@kcl.ac.uk
# Desc: simulation to calculate overlaps

### a list of genes/characters
cvSeed = as.character(1:16000)

## random sample without replacement
selGenes = function(prop){
  return(sample(cvSeed, size = prop*length(cvSeed), replace = F))
}

## proportions, change as required
p1 = 20/100
p2 = 4/100
p3 = 7/100

## simulate one experiment and intersection
simulateOne = function(){
  l1 = selGenes(p1)
  l2 = selGenes(p2)
  l3 = selGenes(p3)
  return(length(Reduce(intersect, list(l1, l2, l3))))
}

## repeat experiment X number of times
ivSim = replicate(1000, simulateOne())

## get a P value for observed proportion in your data w.r.t. simulated 
## observations
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}

## change as required
iObs = 15

## figure and 2 sided p-value
hist(ivSim, main='simulated intersections', 
     xlab='')
points(iObs, 0, pch=20, col=2, cex=2)
getPValue(ivSim, iObs)
