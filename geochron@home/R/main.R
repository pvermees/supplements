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

results <- read.csv('../data/results.csv')
json <- IsoplotR:::fromJSON(file='../data/results.json')
roiss <- IsoplotR:::fromJSON(file='../data/roiss.json')
Pieter <- IsoplotR:::fromJSON(file='../data/Pieter.json')
Andy <- IsoplotR:::fromJSON(file='../data/Andy.json')

PieterResults <- parseJSON(Pieter,defaultrois=roiss)
AndyResults <- parseJSON(Andy)
grains <- intersect(names(PieterResults),names(AndyResults))

#### 1. calculate the density ratios for PV and AC ####
total_counts_in_A0 <- c(N00=0,nP=0,nA=0)
total_counts <- c(N1=0,N2=0)
total_area <- c(A0=0,A1=0,A2=0)
px2um2 <- mean(results[,'area_pixels']/results[,'area_mm2'])/1e6
ng <- length(grains)
PAsPA <- matrix(NA,nrow=ng,ncol=2,
                dimnames=list(grains,c('PA','sPA')))
NA120 <- matrix(NA,nrow=ng,ncol=7,
                dimnames=list(1:ng,c('grain','N1','A1','N2','A2','N0','A0')))
rhoP <- rhoA <- w <- r <- rep(NA,ng)
for (i in seq_along(grains)){
    grain <- grains[i]
    ROIP <- PieterResults[[grain]]$ROI
    ROIA <- AndyResults[[grain]]$ROI
    ROI0 <- polygon_intersection(ROIP,ROIA)
    A1 <- polygon_area(ROIP,px2um2=px2um2)
    A2 <- polygon_area(ROIA,px2um2=px2um2)
    A0 <- polygon_area(ROI0,px2um2=px2um2)
    total_area <- total_area + c(A0,A1,A2)
    xyP <- PieterResults[[grain]]$counts
    xyA <- AndyResults[[grain]]$counts
    xyPinA0 <- xyP[apply(xyP,1,point_in_polygon,ROI0),]
    xyAinA0 <- xyA[apply(xyA,1,point_in_polygon,ROI0),]
    w[i] <- nrow(xyPinA0)
    r[i] <- nrow(xyAinA0)/nrow(xyPinA0)
    N1 <- nrow(xyP)
    N2 <- nrow(xyA)
    total_counts <- total_counts + c(N1,N2)
    counts_in_A0 <- getN00(xyPinA0,xyAinA0)
    N00 <- counts_in_A0['N00']
    total_counts_in_A0 <- total_counts_in_A0 + counts_in_A0
    rhoP[i] <- 100*N1/A1 # x 1e-6 cm-2
    rhoA[i] <- 100*N2/A2 # x 1e-6 cm-2
    PAsPA[i,'PA'] <- rhoP[i]/rhoA[i]
    PAsPA[i,'sPA'] <- PAsPA[i,'PA'] * sqrt(1/N1+1/N2-2*N00/(N1*N2))
    NA120[i,] <- c(grain,N1,A1,N2,A2,N00,A0)
}
rp <- sum(w*r)/sum(w)
sp <- sqrt(sum(w*(r-rp)^2/(ng-1)))

write.table(NA120,file='../output/PA.csv',sep=',')

#### 2. Plot all the superimposed ROIs and counts for PV and AC ####
pdf(file='../output/AvProis.pdf',width=6,height=6,onefile=TRUE)
op <- par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
for (grain in grains){
    plotROIs(PieterResults[[grain]],AndyResults[[grain]])
    mtext(text=paste("grain",grain),line=-1)
}
par(op)
dev.off()

#### 3. Compare PV and AC's results for grain 3 (a.k.a. 4648) ####
pdf(file='../output/AvP.pdf',width=15,height=5)
op <- par(mar=c(4,4,.5,.5),mfrow=c(1,3),mgp=c(2.5,1,0),cex=1.0)
plotROIs(PieterResults[['4648']],AndyResults[['4648']])
legend('topleft','a)',bty='n')
plot(rhoA,rhoP,type='n',
     xlab=expression(hat(rho)[AC]*"(cm"^-2%*%"10"^-6*")"),
     ylab=expression(hat(rho)[PV]*"(cm"^-2%*%"10"^-6*")"),
     bty='n',log='xy',xlim=c(0.1,2),ylim=c(0.1,2))
abline(a=0,b=1)
points(rhoA,rhoP,type='p')
legend('topleft','b)',bty='n')
mtext(text='log ratios PV/AC',line=-1)
mtext(text=paste0('mean = ',signif(mean(log(PAsPA[,'PA'])),3)),line=-2)
mtext(text=paste0('s.d = ',signif(sd(log(PAsPA[,'PA'])),3)),line=-3)
mtext(text=paste0('n = ',nrow(PAsPA)),line=-4)
par(mar=c(4,4,3,1))
PAradial(PAsPA,cex=1.0,spacing=1.2)
legend('topleft','c)',bty='n')
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
grain1 <- list(number='4685',index=23) # list(number='4682',index=39)
grain2 <- list(number='4674',index=25) # list(number='4680',index=43)
pdf(file="../output/23vs25.pdf",onefile=FALSE)
op <- par(mar=c(4,4,0,1))
image1 <- plotimage(idir='../screenshots',
                    grain=allgrains[[grain1$number]]$grain) +
    labs(tag = grain1$index) +
    theme(legend.position="right",
          legend.direction="vertical",
          legend.key.height=unit(1.0, "cm"))
counts1 <- plotcounts(allgrains[[grain1$number]],roiss) +
    theme(legend.position="bottom",
          legend.direction="horizontal",
          legend.key.width=unit(1.0, "cm"))
image2 <- plotimage(idir='../screenshots',
                    grain=allgrains[[grain2$number]]$grain) +
    labs(tag = grain2$index) +
    theme(legend.position="right",
          legend.direction="vertical",
          legend.key.height=unit(1.0, "cm"))
counts2 <- plotcounts(allgrains[[grain2$number]],roiss) +
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

trustworthy_results <- clean_results(results)

#### 6. Scatter plot of crowdsourcing results ####
pdf(file='../output/crowd-correlation.pdf',width=5,height=5)
op <- par(mar=c(3,3,1,1),mgp=c(2,1,0),xpd=NA,bty='n')
lim <- c(0,70)
compare_grains(trustworthy_results,grain_x=grain2$index,grain_y=grain1$index,xlim=lim,ylim=lim)
lines(x=lim,y=lim,col='grey50',lwd=1.5)
par(op)
dev.off()

# generalised linear fit to crowd-sourced data
fit <- glm(count ~ user_id + index, data=trustworthy_results, family="poisson")
mswd <- summary(fit)$deviance/summary(fit)$df.residual

# summary tables
lst <- crowdtable(trustworthy_results)
counts2latex(lst,destination='../output/crowdtable.txt')
counts2latex(lst,destination='../output/shortcrowdtable.txt',short=TRUE)

# stripchart
pdf(file='../output/stripchart.pdf',width=4,height=6,pointsize=7.5)
op <- par(mar=c(4,3,0,9))
grouped_list <- split(trustworthy_results$count, trustworthy_results$index)
ordered_list <- grouped_list[colnames(lst$tab)]
stripchart(ordered_list,method='stack',pch=16,cex=0.65,
           offset=1/5,frame.plot=FALSE,axes=FALSE)
ns <- length(ordered_list)
axis(side=1)
mtext(expression('Number of tracks counted'~'(N'[s]*')'),side=1,line=2.5)
axis(side=2,lwd=0,at=1:ns,labels=names(ordered_list),line=-1,las=1)
mtext("Grain ID in order of PV's counts",side=2,line=1.5)
lines(x=c(0,0),y=c(1,ns))
axis(side=4,lwd=0,at=1:ns,labels=lst$n,line=1,las=1)
axis(side=4,lwd=0,at=1:ns,labels=format(round(lst$ybar,1),nsmall=1),line=3,las=1)
axis(side=4,lwd=0,at=1:ns,labels=format(round(lst$sy,1),nsmall=1),line=5.5,las=1)
axis(side=4,lwd=0,at=ns+1,labels="n",line=1,las=1)
axis(side=4,lwd=0,at=ns+1,labels="mean",line=3,las=1)
axis(side=4,lwd=0,at=ns+1,labels="s.d.",line=5.5,las=1)
medians <- apply(lst$tab,2,median,na.rm=TRUE)
points(x=medians,y=(1:ns)-0.2,pch=24,bg='yellow',cex=0.8)
points(x=lst$tab[1,],y=(1:ns)-0.2,pch=24,bg='blue',cex=0.8)
par(op)
dev.off()
