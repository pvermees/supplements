rm(list=ls())
graphics.off()

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
PAsPA1 <- PAsPA2 <- PAsPA3 <- NULL
for (i in seq_along(grains)){
    grain <- grains[i]
    PieterNsA <- grain2NsA(PieterResults[[grain]])
    AndyNsA <- grain2NsA(AndyResults[[grain]])
    A0 <- polygon_overlap_area(PieterResults[[grain]]$ROI,
                               AndyResults[[grain]]$ROI)
    # optimistic scenario: Rex Galbraith formula
    PA <- (PieterNsA[1]*AndyNsA[2])/(PieterNsA[2]*AndyNsA[1])
    # pessimistic scenario: assuming independence
    N1 <- PieterNsA[1]
    N2 <- AndyNsA[1]
    sPA1 <- PA * sqrt(1/N1 + 1/N2)
    PAsPA1 <- rbind(PAsPA1,c(PA,sPA1))
    rho <- (N1+N2)/(PieterNsA[2]+AndyNsA[2])
    sPA2 <- PA * sqrt(1/(rho*PieterNsA[2]) +
                      1/(rho*AndyNsA[2]) -
                      2*A0/(rho*PieterNsA[2]*AndyNsA[2]))
    PAsPA2 <- rbind(PAsPA2,c(PA,sPA2))
    # realistic scenario: Rex Galbraith formula 2
    xyP <- PieterResults[[grain]]$counts
    xyA <- AndyResults[[grain]]$counts
    N00 <- count(xyP,xyA)
    sPA3 <- PA * sqrt(1/N1+1/N2-2*N00/(N1*N2))
    PAsPA3 <- rbind(PAsPA3,c(PA,sPA3))
}
colnames(PAsPA1) <- colnames(PAsPA2) <- colnames(PAsPA3) <- c('PA','sPA')

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
pdf(file='../output/AvP.pdf',width=9,height=3)
op <- par(mar=c(4,4,3,1),mfrow=c(2,2),mgp=c(2.5,1,0))
plotROIs(PieterResults[['4648']],AndyResults[['4648']])
legend('topleft','a)',bty='n')
PAradial(PAsPA1,cex=0.7,spacing=1.2)
legend('topleft','b)',bty='n')
PAradial(PAsPA2,cex=0.7,spacing=1.2)
legend('topleft','c)',bty='n')
PAradial(PAsPA3,cex=0.7,spacing=1.2)
legend('topleft','d)',bty='n')
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
