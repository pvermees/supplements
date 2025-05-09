rm(list=ls())
graphics.off()

library(IsoplotR)
library(provenance)
library(ggplot2)
library(jpeg)
library(grid)
library(gridExtra)

source("helper.R")
source("polygon_overlap.R")

json <- IsoplotR:::fromJSON(file='../json/results.json')
roiss <- IsoplotR:::fromJSON(file='../json/roiss.json')
Pieter <- IsoplotR:::fromJSON(file='../json/Pieter.json')
Andy <- IsoplotR:::fromJSON(file='../json/Andy.json')

PieterResults <- parseJSON(Pieter,defaultrois=roiss)
AndyResults <- parseJSON(Andy)
grains <- intersect(names(PieterResults),names(AndyResults))

#### 1. calculate the density ratios for PV and AC ####
PAsPA <- NULL
total_counts_in_A0 <- c(N00=0,nP=0,nA=0)
total_counts <- c(N1=0,N2=0)
total_area <- c(A0=0,A1=0,A2=0)
for (i in seq_along(grains)){
    grain <- grains[i]
    ROIP <- PieterResults[[grain]]$ROI
    ROIA <- AndyResults[[grain]]$ROI
    ROI0 <- polygon_intersection(ROIP,ROIA)
    A1 <- polygon_area(ROIP)
    A2 <- polygon_area(ROIA)
    A0 <- polygon_area(ROI0)
    total_area <- total_area + c(A0,A1,A2)
    xyP <- PieterResults[[grain]]$counts
    xyA <- AndyResults[[grain]]$counts
    xyPinA0 <- xyP[apply(xyP,1,point_in_polygon,ROI0),]
    xyAinA0 <- xyA[apply(xyA,1,point_in_polygon,ROI0),]
    N1 <- nrow(xyP)
    N2 <- nrow(xyA)
    total_counts <- total_counts + c(N1,N2)
    counts_in_A0 <- getN00(xyPinA0,xyAinA0)
    N00 <- counts_in_A0['N00']
    total_counts_in_A0 <- total_counts_in_A0 + counts_in_A0
    PA <- (N1*A2)/(N2*A1)
    sPA <- PA * sqrt(1/N1+1/N2-2*N00/(N1*N2))
    PAsPA <- rbind(PAsPA,c(PA,sPA))
}
colnames(PAsPA) <- c('PA','sPA')

#### 2. Plot all the superimposed ROIs and counts for PV and AC ####
pdf(file='../output/AvProis.pdf',width=6,height=6,onefile=TRUE)
op <- par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
for (grain in grains){
    plotROIs(PieterResults[[grain]],AndyResults[[grain]])
    mtext(text=paste("grain",grain),line=-1)
}
par(op)
dev.off()

#### 3. Compare PV and AC's results for grain 4648 ####
pdf(file='../output/AvP.pdf',width=9,height=5)
op <- par(mar=c(4,4,3,1),mfrow=c(1,2),mgp=c(2.5,1,0))
plotROIs(PieterResults[['4648']],AndyResults[['4648']])
legend('topleft','a)',bty='n')
PAradial(PAsPA,cex=1.0,spacing=1.2)
legend('topleft','b)',bty='n')
par(op)
dev.off()

#### 4. Plot all crowdsourcing results on 2D histograms ####
allgrains <- list()
for (grain in json){
    allgrains[[as.character(grain$grain)]] <- grain
}
pdf(file="../output/crowdsourcing.pdf",onefile=TRUE,width=10,height=5)
for (grain in grains){
    grid.arrange(
        plotimage(idir='../screenshots',grain=grain), 
        plotcounts(allgrains[[grain]],roiss) + labs(title=paste("grain",grain)),
        nrow = 1,
        widths = c(1, 1)
    )
}
dev.off()

#### 5. Show crowdsourcing results for two selected grains ####
grain1 <- '4649'
grain2 <- '4673'
pdf(file="../output/4649vs4673.pdf",onefile=FALSE)
op <- par(mar=c(4,4,0,1))
image1 <- plotimage(idir='../screenshots',
                    grain=allgrains[[grain1]]$grain) +
    labs(tag = "a)") +
    theme(legend.position="right",
          legend.direction="vertical",
          legend.key.height=unit(1.0, "cm"))
counts1 <- plotcounts(allgrains[[grain1]],roiss) +
    labs(tag = "c)") +
    theme(legend.position="bottom",
          legend.direction="horizontal",
          legend.key.width=unit(1.0, "cm"))
image2 <- plotimage(idir='../screenshots',
                    grain=allgrains[[grain2]]$grain) +
    labs(tag = "b)") +
    theme(legend.position="right",
          legend.direction="vertical",
          legend.key.height=unit(1.0, "cm"))
counts2 <- plotcounts(allgrains[[grain2]],roiss) +
    labs(tag = "d)") +
    theme(legend.position="right",
          legend.direction="vertical",
          legend.key.height=unit(1.0, "cm"))
grid.arrange(
    image1,counts1,
    image2,counts2,
    nrow = 2,ncol = 2
)
par(op)
dev.off()

#### 6. Create a boxplot of all the crowdsourcing results ####
####    and a radial plot of the two grains in part 5.    ####
pdf(file='../output/radialcrowd.pdf',width=10,height=5)
op <- par(mfrow=c(1,2))
tab <- list2table(allgrains[grains])
# remove Andy Carter, who used a different count area:
dat <- tab[tab[,'worker'] != 292,]
ip <- par(mar=c(4,4,1,1))
boxplot(Ns ~ index,data=dat,col=NA)
legend('topleft',legend="a)",bty='n')
par(ip)
hasgrain1 <- which(dat[,'grain'] == grain1)
hasgrain2 <- which(dat[,'grain'] == grain2)
userswithgrains12 <- intersect(dat[hasgrain1,'worker'],dat[hasgrain2,'worker'])
inum <- which(dat[hasgrain1,'worker'] %in%  userswithgrains12)
iden <- which(dat[hasgrain2,'worker'] %in%  userswithgrains12)
num <- dat[hasgrain1[inum],'Ns']
den <- dat[hasgrain2[iden],'Ns']
numden <- cbind(num,den)
index1 <- which(unique(dat[,'grain']) %in% grain1)
index2 <- which(unique(dat[,'grain']) %in% grain2)
colnames(numden) <- c(paste0('Ns(',index1,')'),
                      paste0('Ns(',index2,')'))
counts <- provenance:::as.counts(numden)
bg <- rep('white',length(userswithgrains12))
bg[userswithgrains12 == 2] <- 'blue'
provenance:::radialplot.counts(counts,bg=bg)
legend('topleft',legend="b)",bty='n')
par(op)
dev.off()
