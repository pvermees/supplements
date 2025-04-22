# to be merged with helper.R and/or main.R

rm(list=ls())
graphics.off()

source("helper.R")
source("polygon_overlap.R")

json <- IsoplotR:::fromJSON(file='../json/results.json')
roiss <- IsoplotR:::fromJSON(file='../json/roiss.json')
Pieter <- IsoplotR:::fromJSON(file='../json/Pieter.json')
Andy <- IsoplotR:::fromJSON(file='../json/Andy.json')

PieterResults <- parseJSON(Pieter,defaultrois=roiss)
AndyResults <- parseJSON(Andy)
grains <- intersect(names(PieterResults),names(AndyResults))

count <- function(xyP,xyA,i=0,cutoff=20){
    nP <- nrow(xyP)
    nA <- nrow(xyA)
    xy <- rbind(xyP,xyA)
    d <- dist(xy)
    dPA <- as.matrix(d)[1:nP,(nP+1):(nP+nA)]
    j <- which.min(dPA)
    if (dPA[j]<cutoff & nP>1 & nA>1){
        rP <- row(dPA)[j]
        rA <- col(dPA)[j]
        i <- count(xyP=xyP[-rP,,drop=FALSE],
                   xyA=xyA[-rA,,drop=FALSE],
                   i=i+1,cutoff=cutoff)
    }
    i
}

N00 <- c()
for (grain in grains){
    xyP <- PieterResults[[grain]]$counts
    xyA <- AndyResults[[grain]]$counts
    N00[grain] <- count(xyP,xyA)
}
