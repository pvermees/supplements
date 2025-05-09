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
        x <- as.numeric(cnts[names(cnts) == 'x_pixels'])
        y <- as.numeric(cnts[names(cnts) == 'y_pixels'])
        out[[grain]]$counts <- cbind(x=x,y=y)
        rownames(out[[grain]]$counts) <- NULL
        out[[grain]]$id <- analysis$id
    }
    out
}

# Plots ROIs for PV and AC
plotROIs <- function(Pgrain,Agrain){
    plot(rbind(Pgrain$ROI,Agrain$ROI),
         type='n',bty='n',xlab='x',ylab='y',asp=1)
    polygon(Pgrain$ROI,border='blue')
    polygon(Agrain$ROI,border='red')
    points(Pgrain$counts,pch=22,bg='blue')
    points(Agrain$counts,pch=21,bg='red')
    legend("top",legend=c("PV","AC"),
           pch=c(22,21),pt.bg=c('blue','red'),
           xpd=NA,bty='n',horiz=TRUE,inset=-0.05)
}

# estimate track count and counting area from a
# list of track coordinates and vertices
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
        paste("central ",rho[PV],"/",rho[AC],
              "-ratio" == a%+-%b,
              " (", n == c, ")"),
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

plotimage <- function(idir='',grain){
    img <- readJPEG(file.path(idir,paste0(grain,".jpg")))
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

# counts tracks identified by analysts P and A within a cutoff distance
getN00 <- function(xyP,xyA,i=0,cutoff=20){
    nP <- nrow(xyP)
    nA <- nrow(xyA)
    xy <- rbind(xyP,xyA)
    d <- dist(xy)
    dPA <- as.matrix(d)[1:nP,(nP+1):(nP+nA)]
    j <- which.min(dPA)
    if (dPA[j]<cutoff & nP>1 & nA>1){
        rP <- row(dPA)[j]
        rA <- col(dPA)[j]
        iPA <- getN00(xyP=xyP[-rP,,drop=FALSE],
                      xyA=xyA[-rA,,drop=FALSE],
                      i=i+1,cutoff=cutoff)
    } else {
        iPA <- c(N00=i,nP=nP,nA=nA)
    }
    return(iPA)
}
