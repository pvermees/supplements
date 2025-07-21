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
PAsPA <- NULL
total_counts_in_A0 <- c(N00=0,nP=0,nA=0)
total_counts <- c(N1=0,N2=0)
total_area <- c(A0=0,A1=0,A2=0)
px2um2 <- mean(results[,'area_pixels']/results[,'area_mm2'])/1e6
PAsPA <- matrix(NA,nrow=length(grains),ncol=2)
colnames(PAsPA) <- c('PA','sPA')
rhoP <- rhoA <- rep(NA,length(grains))
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
}

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
pdf(file='../output/AvP.pdf',width=15,height=5)
op <- par(mar=c(4,4,.5,.5),mfrow=c(1,3),mgp=c(2.5,1,0),cex=1.0)
plotROIs(PieterResults[['4648']],AndyResults[['4648']])
legend('topleft','a)',bty='n')
plot(rhoA,rhoP,type='n',
     xlab=expression(hat(rho)[AC]*"(cm"^2%*%"10"^-6*")"),
     ylab=expression(hat(rho)[PV]*"(cm"^2%*%"10"^-6*")"),
     bty='n',log='xy',xlim=c(0.1,2),ylim=c(0.1,2))
abline(a=0,b=1)
points(rhoA,rhoP,type='p')
legend('topleft','b)',bty='n')
mtext(text='log ratios PV/AC',line=-1)
mtext(text=paste0('mean = ',signif(mean(log(PAsPA[,'PA'])),3)),line=-2)
mtext(text=paste0('st.dev = ',signif(sd(log(PAsPA[,'PA'])),3)),line=-3)
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

#### 6. Boxplot, scatter plot and radial plot of crowdsourcing results ####
grain1 <- 23
grain2 <- 25
trustworthy_results <- clean_results(results)
pdf(file='../output/radialcrowd.pdf',width=7,height=5)
layout(rbind(c(1,1,1,2,2),
             c(1,1,1,3,3)))
p1 <- par(mar=c(3,3,1,8),mgp=c(2,1,0))
# box plot
grouped_list <- split(trustworthy_results$count, trustworthy_results$index)
boxplot(grouped_list,horizontal=TRUE,
        xlab='Ns',ylab='grain number',
        xaxt='n',las=2,col=NA)
add_admin_count_to_boxplot(trustworthy_results,grouped_list)
axis(1)
legend('topleft',legend='a)',bty='n',cex=1.5,adj=c(1,0))
# table
add_table(grouped_list)
# scatter plot
p2 <- par(mar=c(3,2.5,1,1),xpd=NA,bty='n')
compare_grains(trustworthy_results,grain1,grain2)
legend('topleft',legend='b)',bty='n',cex=1.0,adj=c(1,0))
# radial plot
radialcrowd(trustworthy_results,grain1,grain2,from=0.4,to=3.0,t0=1)
legend('topleft',legend='c)',bty='n',cex=1.0,adj=c(1,0))
par(p2)
par(p1)
dev.off()

# generalised linear fit to crowd-sourced data
fit <- glm(count ~ index + user_id, data=trustworthy_results, family="poisson")
mswd <- summary(fit)$deviance/summary(fit)$df.residual
