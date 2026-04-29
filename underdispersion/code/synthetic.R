source("helper.R")
library(IsoplotR)

intro <- function(){

    pdf(file='../figures/synthetic-ecdf.pdf',width=4,height=4)
    
    set.seed(1)
    
    op <- par(mar=c(3.5,3,0.2,0.2),mgp=c(2,1,0),cex=0.8)
    plot(x=c(0,1),y=c(0,1),type='l',xlab='p',ylab='Fn(p)',bty='n',asp=1)
    text(x=0.13,y=0.1,labels='best case scenario (1:1)',pos=4,srt=45,col='black')
    arrows(x0=0.65,x1=0.35,y0=0.75,y1=0.75,length=0.05)
    text(x=0.55,y=0.85,labels='genuine')
    text(x=0.55,y=0.8,labels='overdispersion')
    text(x=0.45,y=0.7,labels='or Type I')
    text(x=0.45,y=0.65,labels='cherry picking')
    arrows(x0=0.825,x1=0.9,y0=0.775,y1=0.7,length=0.05)
    text(x=0.8,y=0.65,labels="'too good",xpd=NA,pos=4)
    text(x=0.8,y=0.6,labels="to be true'",xpd=NA,pos=4,offset=2/3)
    arrows(x0=0.65,x1=0.95,y0=0.4,y1=0.4,length=0.05)
    text(x=0.8,y=0.45,labels='outlier removal')
    text(x=0.8,y=0.35,labels='Type II cherry picking')
    nr <- 100
    x <- rchisq(nr,df=10,ncp=7)
    p1 <- 1-pchisq(x,df=10)
    p2 <- 1-c(pchisq(x[1:round(nr/2)],df=10),
              pchisq(x[round(nr/2+1):nr],df=26))
    xp <- rep(knots(ecdf(p2)),each=2)
    yp <- ecdf(p2)(xp)
    forbidden <- (yp < xp)
    nf <- sum(forbidden)
    xf <- xp[forbidden][-1]
    yf <- yp[forbidden][-nf]
    yf0 <- yp[which(forbidden)[1]-1]
    polygon(x=c(yf0,xf[1],xf,1),
            y=c(yf0,yf0,yf,1),col=rgb(1,0,0,0.5),border=NA)
    lines(ecdf(p1),verticals=TRUE,pch=NA,col='blue',col.01line = NA)
    lines(ecdf(p2),verticals=TRUE,pch=NA,col='red',col.01line = NA)
    legend('bottomright',legend='(a)',bty='n',cex=1.2,inset=1/20)
    par(op)

    dev.off()
}

fissiontrack_cherries <- function(seed=1,ng=20){

    set.seed(seed)

    pick_cherries_FT <- function(samp,alpha=0.05){
        ng <- length(samp)
        cherries <- samp
        cherries$x <- samp$x[1,,drop=FALSE]
        nchanges <- r <- rj <- 0
        avgs <- rep(samp$x[1,1]/samp$x[1,2],ng)
        for (j in 2:ng){
            change <- FALSE
            cherries$x <- rbind(cherries$x,samp$x[j,])
            avg <- central(cherries)
            while (avg$p.value < alpha){
                change <- TRUE
                r <- sum(cherries$x[1:(j-1),1])/sum(cherries$x[1:(j-1),2])
                rj <- cherries$x[j,1]/cherries$x[j,2]
                if (rj > r){
                    cherries$x[j,1] <- cherries$x[j,1] - 1
                } else {
                    cherries$x[j,1] <- cherries$x[j,1] + 1
                }
                avg <- central(cherries)
            }
            if (change) nchanges <- nchanges + 1
            avgs[j] <- sum(cherries$x[1:j,1])/sum(cherries$x[1:j,2])
        }
        list(cherries=cherries,nchanges=nchanges,avgs=avgs)
    }

    N <- 100
    p <- pc <- rep(NA,N)
    samps <- cherries <- avgs <- list()
    nchanges <- rep(0,N)
    for (i in 1:N){
        samps[[i]] <- FTsamp(pop=2,ng=ng,tt=100,sigma=0.2)
        p[i] <- central(samps[[i]])$p.value
        ch <- pick_cherries_FT(samps[[i]])
        cherries[[i]] <- ch$cherries
        avgs[[i]] <- ch$avgs
        nchanges[i] <- ch$nchanges
        pc[i] <- central(cherries[[i]])$p.value
    }

    pdf(file='../figures/FT.pdf',width=4.5,height=4.5)
    
    op <- par(mgp=c(2,1,0),mar=c(3.5,3.5,.5,.5))
    i <- which.max(nchanges)
    r <- samps[[i]]$x[,1]/samps[[i]]$x[,2]
    rc <- cherries[[i]]$x[,1]/cherries[[i]]$x[,2]
    plot(1:ng,r,pch=21,bg='black',xlab='j',
         ylab=expression('N'[s]*'/N'[i]))
    points(1:ng,rc,pch=22,bg='black')
    arrows(x0=1:ng,x1=1:ng,y0=r,y1=rc,length=0.15)
    lines(1:ng,avgs[[i]])

    par(op)

    dev.off()

    #plot(ecdf(pc))

    return(cherries)
}

pick_cherries_wtdmean <- function(dat,alpha=NA,MSWD=NA){
    original <- dat
    ng <- nrow(dat)
    if (is.na(alpha) && is.na(MSWD)){
        alpha <- 0.05
    } else if (is.na(alpha)){
        df <- ng - 1
        X2 <- MSWD * df
        alpha <- 1-pchisq(X2,df=df)
    }
    for (i in 1:ng){
        wtdm <- weightedmean(dat,detect.outliers=FALSE,plot=FALSE)
        if (wtdm$p.value < alpha){
            worst <- which.max(abs(dat[,1]-wtdm$mean[1])/dat[,1])
            dat <- dat[-worst,]
        } else {
            break
        }
    }
    outliers <- which(!original[,1] %in% dat[,1])
    return(list(outliers=outliers,p.value=wtdm$p.value))
}

wtdmean_cherries <- function(seed=1,ng=20){

    set.seed(seed)
    pi1 <- 0.3
    pi2 <- 0.5
    pi3 <- 0.2
    mu <- 100
    s1 <- 20
    s2 <- 5
    s3 <- 10
    err <- 1 + exp(rnorm(n=ng))/50
    n1 <- rbinom(n=1,size=ng,prob=pi1)
    n2 <- ifelse(n1<ng,rbinom(n=1,size=ng-n1,prob=pi2/(pi2+pi3)),0)
    n3 <- ng-n1-n2
    t1 <- mu + abs(rnorm(n=n1,mean=0,sd=s1)) + rnorm(n=n1,mean=0,sd=err[1:n1])
    t2 <- mu + rnorm(n=n2,mean=0,sd=s2) + rnorm(n=n2,mean=0,sd=err[(n1+1):(n1+n2)])
    t3 <- mu - abs(rnorm(n=n3,mean=0,sd=s3)) + rnorm(n=n3,mean=0,sd=err[(n1+n2+1):ng])
    wtdmean_data <- cbind(x=c(t1,t2,t3),sx=err)
    pdf(file='../figures/wtdmean.pdf',width=9,height=4.5)
    op <- par(mfrow=c(1,2),mgp=c(2,1,0),mar=c(3.5,3.5,3.5,0.5))
    weightedmean(wtdmean_data,ranked=TRUE,detect.outliers=FALSE)
    ch <- pick_cherries_wtdmean(wtdmean_data,alpha=0.05)
    weightedmean(wtdmean_data,ranked=TRUE,omit=ch$outliers)
    par(op)
    dev.off()

}

isochron_cherries <- function(seed=1,ng=20){
    set.seed(seed)

    ng <- 20
    Ar4036 <- 298.5
    Ar36 <- runif(ng,min=200,max=1000)
    sAr36 <- sqrt(Ar36)
    cV36 <- sAr36/Ar36
    Ar39 <- exp(rnorm(ng,mean=log(10000),sd=1))
    sAr39 <- sqrt(Ar39)
    cV39 <- sAr39/Ar39
    Ar40 <- Ar36 * Ar4036 + Ar39*10
    sAr40 <- sqrt(Ar40)
    cV40 <- sAr40/Ar40
    x <- Ar39/Ar40
    sx <- x * sqrt(cV39^2 + cV40^2)
    y <- Ar36/Ar40
    sy <- y * sqrt(cV36^2 + cV40^2)
    rho <- (cV40^2)/(sqrt(cV40^2 + cV39^2) + sqrt(cV40^2 + cV36^2))
    sxy <- rho*sx*sy
    XY <- cbind(rep(NA,ng),rep(NA,ng))
    for (i in 1:ng){
        E <- rbind(c(sx[i]^2 + (x[i]/20)^2,sxy[i]),
                   c(sxy[i],sy[i]^2 + (y[i]/10)^2))
        XY[i,] <- MASS::mvrnorm(1,mu=c(x[i],y[i]),Sigma=E)
    }
    yd <- cbind(X=XY[,1],sX=sx,Y=XY[,2],sY=sy,rho=rho)

    pick_cherries_isochron <- function(yd,alpha=0.05){
        omit <- c()
        ng <- nrow(yd)
        good <- 1:ng
        omit <- c()
        fit <- york(yd)
        pval <- fit$p.value
        while (TRUE){
            if (pval>alpha){
                break
            } else {
                pvalvec <- rep(NA,length(good))
                for (i in seq_along(good)){
                    fit <- york(yd[-c(omit,good[i]),])
                    pvalvec[i] <- fit$p.value
                }
                pval <- max(pvalvec)
                bad <- good[which.max(pvalvec)]
                omit <- c(omit,bad)
                good <- good[-which.max(pvalvec)]
            }
        }
        omit
    }

    pdf(file='../figures/isochron.pdf',width=9,height=4.5)

    op <- par(mfrow=c(1,2),mgp=c(2,1,0),mar=c(3.5,3.5,2.5,0.5))
    fit1 <- york(yd)
    xlim <- c(0,-fit1$a[1]/fit1$b[1])
    ylim <- c(0,fit1$a[1]*1.15)
    scatterplot(yd,fit=fit1,xlim=xlim,ylim=ylim,
                xlab=expression(""^36*"Ar/"^40*"Ar"),
                ylab=expression(""^39*"Ar/"^40*"Ar"))
    mtext(line=1,paste0("MSWD = ",signif(fit1$mswd,2),
                        ", p-value = ",signif(fit1$p.value,2),
                        " (",nrow(yd),"/",nrow(yd),")"))
    omit <- pick_cherries_isochron(yd,alpha=0.05)
    fit2 <- york(yd[-omit,])
    scatterplot(yd,fit=fit2,xlim=xlim,ylim=ylim,
                omit=omit,omit.stroke="black",
                xlab=expression(""^36*"Ar/"^40*"Ar"),
                ylab=expression(""^39*"Ar/"^40*"Ar"))
    mtext(line=1,paste0("MSWD = ",signif(fit2$mswd,2),
                        ", p-value = ",signif(fit2$p.value,2),
                        " (",nrow(yd)-length(omit),"/",nrow(yd),")"))
    par(op)

    dev.off()

}

concordia_cherries <- function(seed=6){

    set.seed(seed)

    selection2concordia <- function(X,Y){
        mX <- mean(X)
        mY <- mean(Y)
        EXY <- cov(cbind(X,Y))/length(X)
        selection <- matrix(c(mean(X),sqrt(EXY[1,1]),
                              mean(Y),sqrt(EXY[2,2]),
                              cov2cor(EXY)[1,2]),
                            nrow=1)
        dat <- as.UPb(selection)
        tt <- tryCatch({
            IsoplotR:::concordia_age(dat)
        }, error = function(e){
            NA
        })
        tt
    }

    pick_cherries_concordia <- function(X,Y,alpha=0.05){
        ns <- length(X)
        best_window <- c(1,2)
        best_p_value <- 0
        for (i in 1:(ns-2)){
            for (j in (i+1):ns){
                tt <- selection2concordia(X[i:j],Y[i:j])
                if (all(is.na(tt))){
                # skip
                } else {
                    pval <- tt$p.value['concordance']
                    if (pval>alpha && (j-i)>diff(best_window)){
                        best_window <- c(i,j)
                        best_p_value <- pval
                    }
                }
            }
        }
        list(best_window=best_window,best_p_value=best_p_value)
    }

    tc <- 1050
    stc <- 2
    wc <- age2ratio(tc,stc,ratio="Wetherill")
    U85 <- settings('iratio','U238U235')[1]
    nblk <- 20
    nsig <- 80
    dwell38 <- 1
    dwell06 <- 4
    dwell07 <- 20
    down <- exp(-(1:nsig)/(2*nsig))
    true38sig <- 10000*down
    true06sig <- wc$x[2]*true38sig*dwell06
    true07sig <- wc$x[1]*true38sig*dwell07/U85
    brown38<- cumsum(rnorm(nsig)*true38sig/50)
    brown06 <- cumsum(rnorm(nsig)*true06sig/50)
    brown07 <- cumsum(rnorm(nsig)*true07sig/50)
    # background
    blk38 <- 10
    blk06 <- 5
    blk07 <- 2
    sig38 <- rpois(nblk+nsig,lambda=blk38)
    sig06 <- rpois(nblk+nsig,lambda=blk06)
    sig07 <- rpois(nblk+nsig,lambda=blk07)
    # signal
    isig <- (nblk+1):(nblk+nsig)
    sig38[isig] <- rpois(nsig,lambda=true38sig) + brown38
    sig06[isig] <- rpois(nsig,lambda=true06sig) + brown06
    sig07[isig] <- rpois(nsig,lambda=true07sig) + brown07

    # get average ratios
    buffer <- 5
    meas38 <- (sig38-blk38)[isig[buffer:(nsig-buffer)]]
    meas35 <- meas38/U85
    meas06 <- ((sig06-blk06)/dwell06)[isig[5:(nsig-5)]]
    meas07 <- ((sig07-blk07)/dwell07)[isig[5:(nsig-5)]]
    X <- meas07/meas35
    Y <- meas06/meas38
    E <- cov(cbind(X,Y))/length(X)
    ell <- ellipse(x=mean(X),y=mean(Y),covmat=E)
    fitall <- selection2concordia(X,Y)

    # cherry picking
    cherries <- pick_cherries_concordia(X,Y)
    selection <- cherries$best_window[1]:cherries$best_window[2]
    E_cherries <- cov(cbind(X,Y)[selection,])/diff(cherries$best_window)
    ell_cherries <- ellipse(x=mean(X[selection]),y=mean(Y[selection]),
                            covmat=E_cherries)
    fitcherries <- selection2concordia(X[selection],Y[selection])
    pall <- fitcherries$p.value['concordance']
    
    # plot
    pdf(file='../figures/concordia.pdf',width=9,height=4.5)
    op <- par(mfrow=c(1,2),mgp=c(2,1,0),mar=c(3.5,3.5,3,1))
    tt <- 1:(nblk+nsig)
    matplot(x=cbind(tt,tt,tt),y=cbind(sig06,sig07,sig38),col='black',
            type='b',pch=c(0,2,5),lty=1,xlab="t (s)",ylab="signal (counts)")
    abline(v=c(nblk+buffer,nblk+nsig-buffer),lty=1)
    abline(v=c(nblk+buffer+min(selection),nblk+buffer+max(selection)),lty=2)
    usr <- par('usr')
    arrows(x0=nblk+buffer,
           x1=nblk+nsig-buffer,
           y0=usr[4]*1.1,
           y1=usr[4]*1.1,code=3,length=0.1,xpd=NA,col='blue')
    arrows(x0=nblk+buffer+min(selection),
           x1=nblk+buffer+max(selection),
           y0=usr[4]*1.05,
           y1=usr[4]*1.05,code=3,length=0.1,xpd=NA,col='red')
    text(x=nblk+buffer,
         y=usr[4]*1.1,
         labels='Window 1',pos=2,xpd=NA,col='blue')
    text(x=nblk+buffer+min(selection),
         y=usr[4]*1.05,
         labels='Window 2',pos=2,xpd=NA,col='red')
    concordia(tlim=c(700,1200))
    points(x=X,y=Y,pch=21)
    polygon(ell,lwd=2,col=rgb(0,0,1,0.5))
    polygon(ell_cherries,lwd=2,col=rgb(1,0,0,0.5))
    mtext(bquote("Window 1: MSWD =" ~ .(signif(fitall$mswd["concordance"], 2)) *
                     ", p(" * chi^2 * ") =" ~ .(signif(fitall$p.value["concordance"], 2))), 
          side = 3, line = 1)
    mtext(bquote("Window 2: MSWD =" ~ .(signif(fitcherries$mswd["concordance"], 2)) *
                     ", p(" * chi^2 * ") =" ~ .(signif(fitcherries$p.value["concordance"], 2))), 
          side = 3, line = 0)
    par(op)
    dev.off()

}

synthetic_ecdfs <- function(ns=100,seed=1){

    get_meta_data <- function(ns=100,maxdisp=0){
        out <- list()
        ttrue <- 100
        sigma <- runif(ns)*maxdisp
        err <- 2
        ng <- runif(ns,min=5,max=50)
        for (i in 1:ns){
            epsilon <- rnorm(ng[i],mean=0,sd=err)
            tmeas <- exp(rnorm(ng[i],mean=log(ttrue),sd=sigma[i])) + epsilon
            out[[i]] <- cbind(x=tmeas,sx=err)
        }
        out
    }

    get_pval_ecdf <- function(datlist,alpha=NA,MSWD=NA){
        ns <- length(datlist)
        pval <- rep(NA,ns)
        for (i in 1:ns){
            dat <- datlist[[i]]
            cdat <- pick_cherries_wtdmean(dat,alpha=alpha,MSWD=MSWD)
            if (length(cdat$outliers)>0){
                fit <- weightedmean(dat[-cdat$outliers,],
                                    detect.outliers=FALSE,plot=FALSE)
            } else {
                fit <- weightedmean(dat,detect.outliers=FALSE,plot=FALSE)
            }
            pval[i] <- fit$p.value
        }
        pval
    }
    
    set.seed(seed)

    pdf(file='../figures/synthetic_ecdfs.pdf',width=9,height=4.5)
    
    op <- par(mfrow=c(1,2),mgp=c(2,1,0),mar=c(3.5,3.5,0.5,0.5))

    non_dispersed_data <- get_meta_data(ns=ns,maxdisp=0)
    pval_unpicked <- get_pval_ecdf(non_dispersed_data,alpha=0)
    pval_alpha_picked <- get_pval_ecdf(non_dispersed_data,alpha=0.05)
    pval_MSWD_picked <- get_pval_ecdf(non_dispersed_data,MSWD=1)

    plot(ecdf(pval_unpicked),verticals=TRUE,
         pch=NA,main='',xlab='p',ylab='Fn(p)')
    lines(ecdf(pval_alpha_picked),verticals=TRUE,pch=NA,col='blue')
    lines(ecdf(pval_MSWD_picked),verticals=TRUE,pch=NA,col='red')
    lines(x=c(0,1),y=c(0,1))

    legend('topleft','a)',bty='n',adj=c(1,2))

    dispersed_data <- get_meta_data(ns=ns,maxdisp=0.03)
    pval_dispersed_unpicked <- get_pval_ecdf(dispersed_data,alpha=0)
    pval_dispersed_alpha_picked <- get_pval_ecdf(dispersed_data,alpha=0.05)
    pval_dispersed_MSWD_picked <- get_pval_ecdf(dispersed_data,MSWD=1)

    plot(ecdf(pval_dispersed_unpicked),verticals=TRUE,
         pch=NA,main='',xlab='p',ylab='Fn(p)')
    lines(ecdf(pval_dispersed_alpha_picked),verticals=TRUE,pch=NA,col='blue')
    lines(ecdf(pval_dispersed_MSWD_picked),verticals=TRUE,pch=NA,col='red')
    pval_mixed <- c(sample(x=pval_dispersed_unpicked,size=ns/2),
                    #sample(x=pval_dispersed_alpha_picked,size=ns/3),
                    sample(x=pval_dispersed_MSWD_picked,size=ns/2))
    lines(ecdf(pval_mixed),verticals=TRUE,pch=NA,col='orange')
    lines(x=c(0,1),y=c(0,1))

    legend('topleft','b)',bty='n',adj=c(1,2))
    
    par(op)

    dev.off()
}

detrital_cherries <- function(){
    dat <- read.csv("../data/DLNR.csv")
    UPb <- as.UPb(as.matrix(dat[,-1]),format=2)
    tt <- age(UPb,discordance=list(option='p',before=TRUE))
    dat <- read.csv("../data/DLNRcherries.csv")
    cherries <- as.UPb(as.matrix(dat[,-1]),format=2)
    ttc <- age(cherries,discordance=list(option='p',before=TRUE))
    pdf(file='../figures/detritals.pdf',width=9,height=3.5)
    op <- par(mfrow=c(1,3),mgp=c(2,1,0),mar=c(3.5,3.5,1,0.5))
    concordia(UPb,tlim=c(0,3000),ellipse.fill="#0000FF80")
    axis(side=3)
    concordia(cherries,tlim=c(0,3000),ellipse.fill="#FF000080")
    axis(side=4)
    plot(ecdf(tt[,'p[conc]']),verticals=TRUE,pch=NA,
         col="#0000FF80",main=NA,xlab='p',ylab='Fn(p)')
    lines(ecdf(ttc[,'p[conc]']),verticals=TRUE,pch=NA,col="#FF000080")
    legend('topleft',legend=c('default selection windows',
                              'cherry-picked'),
           col=c("#0000FF80","#FF000080"),bty='n',lty=rep(1,2))
    par(op)
    dev.off()
}

type1cherries <- function(){
    nr <- 3
    nc <- 6
    pdf(file='../figures/type1.pdf',width=10,height=10*nr/nc)
    m <- rbind(
        c(1,rep(2,nc),3),
        c(1,4+(1:nc),3),
        c(1,4+((nc+1):(2*nc)),3),
        c(1,4+((2*nc+1):(3*nc)),3),
        c(1,rep(4,nc),3)
    )
    op <- par(mar=rep(0,4),mgp=c(1.4,0.5,0))
    layout(m,widths=c(0.03,rep(0.97/nc,nc),0.001),
           heights=c(0.02,rep(0.9/nr,nr),0.08))#, layout.show(n=nr*nc+3)
    set.seed(5)
    n <- 10
    plot.new()
    plot.new()
    plot.new()
    plot.new()
    for (i in 1:nr){
        for (j in 1:nc){
            X <- runif(10)
            Y <- runif(10)
            plot(X,Y,xlim=c(0,1.0),ylim=c(0,1.3),asp=1,xaxt='n',yaxt='n',type='n')
            fit <- lm(Y ~ X)
            smry <- summary(fit)
            r2 <- signif(smry$r.squared,2)
            pval <- signif(smry$coefficients[2,4],2)
            mtext(bquote('r'^2*'='*.(r2)*', p='*.(pval)),line=-1.6,cex=0.8)
            if (pval < 0.05){
                lines(c(0,1),predict(fit,newdata=data.frame(X=c(0,1))))
            }
            points(X,Y,pch=21,bg='white')
            if (i==nr){
                axis(side=1,at=c(0,0.5,1))
                mtext('x',side=1,line=1.5)
            }
            if (j==1){
                axis(side=2,at=c(0,0.5,1))
                mtext('y',side=2,line=1.5)
            }
        }
    }
    par(op)
    dev.off()
}

intro()

FTcherries <- fissiontrack_cherries()

UPbcherries <- wtdmean_cherries()

ArArCherries <- isochron_cherries(seed=26)

concordiaCherries <- concordia_cherries()

synthetic_ecdfs(n=200,seed=3)

detrital_cherries()

type1cherries()
