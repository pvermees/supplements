rm(list=ls())
library(ggplot2)
library(jpeg)
library(grid)
library(gridExtra)
source("polygon_overlap.R")
setwd('/home/pvermees/Dropbox/fissiontracks/geochron@home/paper/code')
graphics.off()

parseJSON <- function(json,defaultrois){
    out <- list()
    for (analysis in json){
        grain <- as.character(analysis$grain)
        out[[grain]] <- list()
        regions <- analysis$regions[[1]]
        if (is.null(regions)){
            for (item in defaultrois){
                if (item$grain_id == analysis$grain){
                    roi <- unlist(item$regions[[1]]$vertices)
                }
            }
        } else {
            roi <- unlist(regions)
        }
        out[[grain]]$ROI <- matrix(roi,ncol=2,byrow=TRUE)
        colnames(out[[grain]]$ROI) <- c('x','y')
        cnts <- unlist(analysis$grainpoints)
        out[[grain]]$counts <- cbind(x=as.numeric(cnts[names(cnts) == 'x_pixels']),
                                     y=as.numeric(cnts[names(cnts) == 'y_pixels']))
        rownames(out[[grain]]$counts) <- NULL
        out[[grain]]$id <- analysis$id
    }
    out
}

# Plots ROIs for PV and AC
plotROIs <- function(Pgrain,Agrain){
    plot(rbind(Pgrain$ROI,Agrain$ROI),type='n',bty='n',xlab='x',ylab='y',asp=1)
    polygon(Pgrain$ROI,border='blue')
    polygon(Agrain$ROI,border='red')
    points(Pgrain$counts,pch=22,bg='blue')
    points(Agrain$counts,pch=21,bg='red')
    legend("top",legend=c("PV","AC"),
           pch=c(22,21),pt.bg=c('blue','red'),
           xpd=NA,bty='n',horiz=TRUE,inset=-0.05)
}

# estimate track count and counting area from list of track coordinates and vertices
grain2NsA <- function(lst,pix2mm2=1){
    Ns <- nrow(lst$counts)
    x <- lst$ROI[,'x']
    y <- lst$ROI[,'y']
    pix2 <- 0.5 * abs(sum(x * c(y[-1], y[1])) - sum(y * c(x[-1], x[1])))
    A <- pix2mm2 * pix2
    c(Ns,A)
}

PAradial <- function(PAsPA,cex=1.0,spacing=1.0){
    fit <- IsoplotR::radialplot(PAsPA,title=FALSE,bg='white',
                                xlab='precision',z0=1)
    tst <- IsoplotR:::roundit(fit$age)
    maintit <- substitute(
        paste("central ",rho[PV],"/",rho[AC],"-ratio" == a%+-%b, " (", n == c, ")"),
        list(a = tst[1], b = tst[2], c = nrow(PAsPA))
    )
    mswdtit <- substitute(
        paste(MSWD == a, ", ", p == b),
        list(
            a = signif(fit$mswd, 2),
            b = signif(fit$p.value, 3)
        )
    )
    dsd <- IsoplotR:::roundit(100*fit$disp)
    disptit <- substitute(
        paste("dispersion" == a%+-%b, "%"),
        list(a = dsd[1], b = dsd[2])
    )
    mtext(text=maintit,line=spacing,cex=cex)
    mtext(text=mswdtit,line=0,cex=cex)
    if (fit$p.value<0.05){
        mtext(text=disptit,line=-spacing,cex=cex)
    }
}

# gets ROI from rois.json
getroi <- function(grain_id,rois){
    for (roi in rois){
        if (roi$grain_id == grain_id){
            if (length(roi$regions)>1){
                warning('This script cannot handle multi-part ROIs.')
            }
            region <- roi$regions[[1]]
            shift <- region$shift
            vertices <- region$vertices
            nvertices <- length(vertices)
            coordinates <- unlist(vertices)
            xs <- coordinates[seq(from=1,to=2*nvertices-1,by=2)]
            ys <- coordinates[seq(from=2,to=2*nvertices,by=2)]
            lat_scaled <- (roi$image_height-ys)/roi$image_width
            lon_scaled <- xs/roi$image_width
            return(cbind(c(lon_scaled,lon_scaled[1]),
                         c(lat_scaled,lat_scaled[1])))
        }
    }
    return(lon=c(0,1,1,0,0),lat=c(0,0,1,1,0))
}

plotimage <- function(grain){
    img <- readJPEG(paste0(grain,".jpg"))
    out <- ggplot() +
        annotation_custom(
            rasterGrob(img, interpolate = TRUE),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
        ) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5))
    out
}

plotcounts <- function(grain,rois){
    ofname <- paste0(grain,'.pdf')
    roi <- getroi(grain$grain,rois)
    colnames(roi) <- c('lon','lat')
    title <- paste0('grain=',grain$grain)
    xy = data.frame(lon=c(),lat=c())
    for (result in grain$results){
        ncounts <- length(result$latlngs)
        latlngs <- unlist(result$latlngs)
        if (ncounts>0){
            lat <- latlngs[seq(from=1,to=2*ncounts-1,by=2)]
            lon <- latlngs[seq(from=2,to=2*ncounts,by=2)]
            xy <- rbind(xy,cbind(lon,lat))
        }
    }
    
    dx <- diff(range(roi[,'lon']))
    out <- ggplot(xy, aes(x = lon, y = lat)) +
        geom_hex(binwidth = rep(dx,2)/30) +
        coord_fixed() + 
        scale_fill_gradient(low = "white", high = "black") + 
        geom_polygon(data = as.data.frame(roi),
                     mapping = aes(x = lon, y = lat),
                     fill = NA,
                     linewidth = 1,
                     color = "red") +
        theme(
            panel.background = element_rect(fill='transparent'),
            axis.ticks = element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            legend.key = element_rect(colour = 'black',
                                      fill = 'pink',
                                      linewidth = 0.5,
                                      linetype='dashed')
        )
    out
}

# create a results.csv-like table from a list of raw counts and vertices
list2table <- function(lst){
    index <- grains <- workers <- Ns <- A <- NULL
    for (i in seq_along(lst)){
        grain <- lst[[i]]
        for (result in grain$results){
            index <- append(index,i)
            grains <- append(grains,grain$grain)
            workers <- append(workers,result$worker$id)
            Ns <- append(Ns,result$result)
            A <- append(A,grain$area_mm2)
        }
    }
    data.frame(index=index,grain=grains,worker=workers,Ns=Ns,A=A)
}

json <- IsoplotR:::fromJSON(file='results.json')
roiss <- IsoplotR:::fromJSON(file='roiss.json')
Pieter <- IsoplotR:::fromJSON(file='Pieter.json')
Andy <- IsoplotR:::fromJSON(file='Andy.json')

PieterResults <- parseJSON(Pieter,defaultrois=roiss)
AndyResults <- parseJSON(Andy)
grains <- intersect(names(PieterResults),names(AndyResults))

#### 1. generate Markdown table of links to results for PV and AC ####
md <- paste0("| i | grain | PV (N<sub>s</sub>,μm<sup>2</sup>) | AC (N<sub>s</sub>,μm<sup>2</sup>) |\n",
             "|---|-------|-----------------------------------|-----------------------------------|\n")
url <- "https://isoplotr.es.ucl.ac.uk/geochron@home/ftc/result"
PAsPA1 <- PAsPA2 <- NULL
for (i in seq_along(grains)){
    grain <- grains[i]
    PieterNsA <- grain2NsA(PieterResults[[grain]])
    AndyNsA <- grain2NsA(AndyResults[[grain]])
    A0 <- polygon_overlap_area(PieterResults[[grain]]$ROI,
                               AndyResults[[grain]]$ROI)
    md <- paste0(md,
                 "| ", i, " | ", grain , 
                 "| ( [",PieterNsA[1],",",round(PieterNsA[2]),"](",
                 file.path(url,PieterResults[[grain]]$id),
                 ") ) ",
                 "| ( [",AndyNsA[1],",",round(AndyNsA[2]),"](",
                 file.path(url,AndyResults[[grain]]$id),
                 ") )|\n")
    # optimistic scenario: Rex Galbraith formula
    PA <- (PieterNsA[1]*AndyNsA[2])/(PieterNsA[2]*AndyNsA[1])
    # pessimistic scenario: assuming independence
    sPA1 <- PA * sqrt(1/PieterNsA[1] + 1/AndyNsA[1])
    PAsPA1 <- rbind(PAsPA1,c(PA,sPA1))
    rho <- (PieterNsA[1]+AndyNsA[1])/(PieterNsA[2]+AndyNsA[2])
    sPA2 <- PA * sqrt(1/(rho*PieterNsA[2]) +
                      1/(rho*AndyNsA[2]) -
                      2*A0/(rho*PieterNsA[2]*AndyNsA[2]))
    PAsPA2 <- rbind(PAsPA2,c(PA,sPA2))
}
colnames(PAsPA1) <- colnames(PAsPA2) <- c('PA','sPA')
cat(md,file="PVvAC.md")

#### 2. Plot all the superimposed ROIs and counts for PV and AC ####
pdf(file='../figures/AvProis.pdf',width=6,height=6,onefile=TRUE)
op <- par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
for (grain in grains){
    plotROIs(PieterResults[[grain]],AndyResults[[grain]])
    mtext(text=paste("grain",grain),line=-1)
}
par(op)
dev.off()

#### 3. Compare PV and AC's results for grain 4648 ####
pdf(file='../figures/AvP.pdf',width=9,height=3)
op <- par(mar=c(4,4,3,1),mfrow=c(1,3),mgp=c(2.5,1,0))
plotROIs(PieterResults[['4648']],AndyResults[['4648']])
legend('topleft','a)',bty='n')
PAradial(PAsPA1,cex=0.7,spacing=1.2)
legend('topleft','b)',bty='n')
PAradial(PAsPA2,cex=0.7,spacing=1.2)
legend('topleft','c)',bty='n')
par(op)
dev.off()

#### 4. Plot all crowdsourcing results on 2D histograms ####
allgrains <- list()
for (grain in json){
    allgrains[[as.character(grain$grain)]] <- grain
}
pdf(file="../figures/crowdsourcing.pdf",onefile=TRUE,width=10,height=5)
for (grain in grains){
    grid.arrange(
        plotimage(grain), 
        plotcounts(allgrains[[grain]],roiss) + labs(title=paste("grain",grain)),
        nrow = 1,
        widths = c(1, 1)
    )
}
dev.off()

#### 5. Show crowdsourcing results for two selected grains ####
grain1 <- '4649'
grain2 <- '4673'
pdf(file="../figures/4649vs4673raw.pdf",onefile=FALSE)
op <- par(mar=c(4,4,0,1))
image1 <- plotimage(allgrains[[grain1]]$grain) +
    labs(tag = "a)") +
    theme(legend.position="right",
          legend.direction="vertical",
          legend.key.height=unit(1.0, "cm"))
counts1 <- plotcounts(allgrains[[grain1]],roiss) +
    labs(tag = "c)") +
    theme(legend.position="bottom",
          legend.direction="horizontal",
          legend.key.width=unit(1.0, "cm"))
image2 <- plotimage(allgrains[[grain2]]$grain) +
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
pdf(file='../figures/radialcrowd.pdf',width=10,height=5)
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
