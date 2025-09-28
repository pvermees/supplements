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
    invisible(fit)
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
    } else if (dPA[j]<cutoff){
        iPA <- c(N00=i+1,nP=nP,nA=nA)
    } else {
        iPA <- c(N00=i,nP=nP,nA=nA)
    }
    return(iPA)
}

clean_results <- function(results,ACgrains){
    time_obj <- sapply(results$create_date,
                       FUN=as.POSIXct,
                       format = "%Y-%m-%d %H:%M:%OS",
                       tz = "UTC")
    timely_submission <- time_obj < as.POSIXct("2025-02-11 00:00:00.00+00:00",
                                               format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
    is_trustworthy <- !(results$user_id %in% c(237,238)) # non UCL students
    uses_default_roi <- !(results$user_id %in% c(292,2)) # AC and PV
    is_admin <- results$user_id == 1
    AC_ids <- unique(results$index[results$user_id==292])
    good_grains <- results$index %in% AC_ids
    good_counts <- (timely_submission & is_trustworthy & uses_default_roi) | is_admin
    keep <- good_grains & good_counts
    out <- results[keep,]
    out$index <- as.factor(out$index)
    out$user_id <- as.factor(out$user_id)
    out
}

add_admin_count_to_boxplot <- function(trustworthy_results,grouped_list){
    admin_results <- trustworthy_results[trustworthy_results$user_id == 1,]
    ACids <- as.numeric(names(grouped_list))
    admin_counts <- admin_results$count[match(ACids,admin_results$index)]
    points(admin_counts,1:length(ACids),pch=21,bg='blue')
}

compare_grains <- function(results,grain1,grain2,plot=TRUE,...){
    results1 <- results[results$index == grain1,]
    results2 <- results[results$index == grain2,]
    users <- intersect(results1$user_id,
                       results2$user_id)
    i1 <- match(users,results1$user_id)
    i2 <- match(users,results2$user_id)
    out <- list(name=paste0(grain1,' vs ',grain2))
    class(out) <- 'counts'
    out$x <- cbind(results1$count[i1],
                   results2$count[i2])
    colnames(out$x) <- c(paste0('Ns(',grain1,')'),
                         paste0('Ns(',grain2,')'))
    ns <- nrow(out$x)
    if (plot){
        bg <- rep(NA,ns)
        iadmin <- which(users==1)
        bg[iadmin] <- 'blue'
        duplo <- duplicated(out$x)
        X <- out$x
        X[duplo,] <- jitter(X[duplo,])
        plot(X,
             xlab=paste0('grain ',grain1),
             ylab=paste0('grain ',grain2),
             pch=21,bg=bg,...)
        rho <- format(round(cor(out$x)[1,2],2),nsmall=2)
        mtext(line=0,text=paste0('correlation=',rho),cex=1.0)
        PVratio <- out$x[iadmin,1]/out$x[iadmin,2]
        mtext(line=-2,
              text=paste0("PV's ratio=",
                          format(round(PVratio,2),nsmall=2)),
              cex=1.0)
        pooled_ratio <- sum(results1$count[i1])/sum(results2$count[i2])
        mtext(line=-3,
              text=paste0("pooled ratio=",
                          format(round(pooled_ratio,2),nsmall=2)),
              cex=1.0)
    }
    out
}

radialcrowd <- function(results,grain1,grain2,from=NA,to=NA,t0=NA,...){
    radialdat <- compare_grains(results,grain1,grain2,plot=FALSE)
    duplicates <- duplicated(radialdat)
    radialdat$x[duplicates,] <- jitter(radialdat$x[duplicates,])
    provenance:::radialplot.counts(radialdat,title=FALSE,from=from,to=to,t0=t0,...)
    pooledratio <- sum(radialdat$x[,1])/sum(radialdat$x[,2])
    admincount1 <- results$count[results$user_id==1 & results$index==grain1]
    admincount2 <- results$count[results$user_id==1 & results$index==grain2]
    adminratio <- admincount1/admincount2
    graphics::mtext(paste0('pooled ratio = ',signif(pooledratio,3)),line=-1,cex=0.7)
    graphics::mtext(paste0('PV ratio = ',signif(adminratio,3)),line=-2,cex=0.7)
    zs_admin <- provenance:::x2zs(tail(radialdat$x,1),from=from,to=to,t0=t0)
    IsoplotR:::plot_radial_points(zs_admin,pch=21,bg='blue')
    invisible(radialdat)
}

add_table <- function(grouped_list,sigdig=1,nsmall=1,x=c(70,75,85)){
    y <- length(grouped_list)+1
    text(x=x[1],y=y,labels='n',xpd=NA,pos=4)
    text(x=x[2],y=y,labels='mean',xpd=NA,pos=4)
    text(x=x[3],y=y,labels='s.d',xpd=NA,pos=4)
    for (i in seq_along(grouped_list)){
        text(x=x[1],y=i,labels=length(grouped_list[[i]]),xpd=NA,pos=4)
        text(x=x[2],y=i,labels=format(round(mean(grouped_list[[i]]),sigdig),nsmall=nsmall),xpd=NA,pos=4)
        text(x=x[3],y=i,labels=format(round(sd(grouped_list[[i]]),sigdig),nsmall=nsmall),xpd=NA,pos=4)
    }
}

crowdtable <- function(trustworthy_results){
    count_table <- reshape(
        trustworthy_results[,c('user_id','index','count')],
        direction = "wide",
        idvar = "user_id",
        timevar = "index",
        sep = ""
    )
    area_table <- reshape(
        trustworthy_results[,c('user_id','index','area_pixels')],
        direction = "wide",
        idvar = "user_id",
        timevar = "index",
        sep = ""
    )
    order_table <- function(unordered_table,first){
        res <- data.matrix(subset(unordered_table,select=-user_id))
        rownames(res) <- unordered_table$user_id
        colnames(res) <- substring(colnames(res),first=first)
        res[order(as.integer(rownames(res))),
            order(as.integer(colnames(res)))]
    }
    ordered_counts <- order_table(count_table,first=6) # 6 characters in "count."
    ng <- rp <- sp <- rep(NA,nrow(ordered_counts))
    for (i in 1:nrow(ordered_counts)){
        good <- !is.na(ordered_counts[i,])
        ng[i] <- sum(good)
        w <- ordered_counts[1,good]
        r <- ordered_counts[i,good]/ordered_counts[1,good]
        rp[i] <- sum(w*r)/sum(w)
        sp[i] <- sqrt(sum(w*(r-rp[i])^2)/(ng[i]-1))
    }
    n <- colSums(!is.na(ordered_counts))
    ybar <- colMeans(ordered_counts,na.rm=TRUE)
    sy <- apply(ordered_counts,2,sd,na.rm=TRUE)
    roworder <- c(1,1+order(rp[-1],decreasing=TRUE))
    colorder <- order(ordered_counts[1,])
    list(tab=ordered_counts[roworder,colorder],
         ng=ng[roworder],rp=rp[roworder],sp=sp[roworder],
         n=n[colorder],ybar=ybar[colorder],sy=sy[colorder])
}

counts2latex <- function(lst,destination,short=FALSE){
    ng <- lst$ng
    rp <- sprintf('%.2f',lst$rp)
    sp <- sprintf('%.2f',lst$sp)
    short_rp <- gsub('0\\.', '.', rp)
    short_sp <- gsub('0\\.', '.', sp)
    n <- lst$n
    ybar <- sprintf('%.1f',lst$ybar)
    sy <- sprintf('%.1f',lst$sy)
    short_ybar <- gsub('0\\.', '.', ybar)
    short_sy <- gsub('0\\.', '.', sy)
    medians <- apply(lst$tab[-1,],2,median,na.rm=TRUE)
    r_m <- sprintf('%.2f',medians/lst$tab[1,])
    short_rm <- gsub('0\\.', '.',r_m)
    out <- rbind(c('user','~',colnames(lst$tab),'~','$n_g$','$r_p$','$s_p$'),
                 c('~','~',rep(NA,ncol(lst$tab)),'~','~','~','~'),
                 c('PV',NA,lst$tab[1,],NA,ng[1],short_rp[1],short_sp[1]),
                 NA,
                 cbind(rownames(lst$tab)[-1],NA,lst$tab[-1,],NA,
                       ng[-1],short_rp[-1],short_sp[-1]),
                 c('$r_m$','~',short_rm,'~','~','~','~'))
    if (short){
        out <- out[-(15:(nrow(out)-11)),]
        out[15,] <- '~'
        out[15,1] <- "$\\vdots$"
        out <- out[-c(2,4),]
        out <- out[,-c(2,ncol(out)-3)]
    }
    write2tabular(out,medians,destination,short)
}

write2tabular <- function(out,medians,destination,short=FALSE){
    nr <- nrow(out)
    nc <- ncol(out)
    file_conn <- file(destination, open = "w")
    cat('\\begin{tabular}{',file=file_conn)
    if (short){
        fmt <- c('c|',rep('c@{~}',nc-5),'c|',rep('c@{~}',3))
    } else {
        fmt <- c('c@{~}','c@{~~~}',rep('c@{~}',nc-6),'c@{~~~}',rep('c@{~}',3))
    }
    cat(paste0(fmt,collapse=''),file=file_conn)
    cat('}\n',file=file_conn)
    cat(paste0('\\multicolumn{',nc,'}{c}{grain ID}\\cr\n'),file=file_conn)
    for (i in 1:nr){
        values <- out[i,]
        values[is.na(values)] <- ''
        cat(paste0(values[1:2],collapse=' & '),file = file_conn)
        cat(' & ',file = file_conn)
        for (j in 3:(nc-2)){
            needle <- round(as.numeric(values[j]),digits=1)
            haystack <- round(medians[j-2]+c(-0.5,0,0.5),digits=1)
            if (any(haystack == needle,na.rm=TRUE) && i>2 && i<(nr-2) && !short){
                cat('\\bf ',file = file_conn)
            }
            cat(values[j],file = file_conn)
            cat(' & ',file = file_conn)
        }
        cat(paste0(tail(values,2),collapse=' & '),file = file_conn)
        if (i%in%c(1,2,nr-1)){
            cat(' \\cr\\hline\n',file = file_conn)
        } else {
            cat(' \\cr\n',file = file_conn)
        }
    }
    cat('\\end{tabular}\n',file=file_conn)
    close(file_conn)    
}
