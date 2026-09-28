read_data <- function(fn='Desch2.csv'){
    out <- dat <- read.csv(fn,check.names=FALSE)
    out[,3:4] <- dat[,3:4] * 1e-7
    out[,5:6] <- dat[,5:6] * 1e-6
    out[,7:8] <- dat[,7:8] * 1e-5
    out
}

get_pars <- function(nuclide){
    if (nuclide=='Al26'){
        n <- 26
        t12 <- 0.717
        ylab <- expression('ln('^26*'Al/'^27*'Mg)')
        ycol <- '(Al26/Al27)_0'
        yerrcol <- '2se (Al26/Al27)_0'
        ySSlab <- expression(omega * '('^26 * 'Al/'^27 * 'Mg)'[ss])
    } else if (nuclide=='Mn53'){
        n <- 53
        t12 <- 3.7
        ylab <- expression('ln('^53*'Mn/'^55*'Mn)')
        ycol='(Mn53/Mn55)_0'
        yerrcol='2se (Mn53/Mn55)_0'
        ySSlab <- expression(omega * '('^53 * 'Mn/'^55 * 'Mn)'[ss])
    } else if (nuclide=='Hf182'){
        n <- 182
        t12 <- 8.9
        ylab <- expression('ln('^182*'Hf/'^180*'Hf)')
        ycol='(Hf182/Hf180)_0'
        yerrcol='2se (Hf182/Hf180)_0'
        ySSlab <- expression(omega * '('^182 * 'Hf/'^180 * 'Hf)'[ss])
    }
    list(n=n,
         lambda=log(2)/t12,
         ylab=ylab,
         ycol=ycol,
         yerrcol=yerrcol,
         ySSlab=ySSlab)
}

LL <- function(lw,X,vX,Y,vY){
    w2 <- exp(lw)^2
    x <- (Y*w2+X*vY+Y*vX)/(w2+vY+vX)
    LL <- ((Y-x)^2/vY+(X-x)^2/(w2+vX))/2+log(abs(vY)*abs(w2+vX))/2
    sum(LL)
}

get_w <- function(X,sX,Y,sY){
    fit <- optimise(f=LL,
                    lower=min(log(sX)),
                    upper=log(diff(range(X))),
                    X=X,vX=sX^2,Y=Y,vY=sY^2)
    lw <- fit$minimum
    H <- optimHess(par=lw,fn=LL,X=X,vX=sX^2,Y=Y,vY=sY^2)
    c(w=exp(lw),relerr=sqrt(solve(H)))
}

data2york <- function(dat,Pbcol,Pberrcol,ycol,yerrcol,lambda){
    y_vs_Pb <- dat[,c(Pbcol,Pberrcol,ycol,yerrcol)]
    complete <- complete.cases(y_vs_Pb)
    y_vs_Pb_complete <- data.matrix(y_vs_Pb[complete,])
    groups <- dat[complete,'Type']
    ty <- log(y_vs_Pb_complete[,3])/lambda
    sty <- y_vs_Pb_complete[,4]/(2*y_vs_Pb_complete[,3]*lambda)
    cbind(X=y_vs_Pb_complete[,1],
          sX=y_vs_Pb_complete[,2]/2,
          Y=ty,
          sY=sty,
          rXY=0)
}

statchron <- function(dat,
                      Pbcol='tPbPb (original)',
                      Pberrcol='2se (tPbPb original)',
                      nuclide='Al26'){

    pars <- get_pars(nuclide)

    raw <- data2york(dat=dat,Pbcol=Pbcol,Pberrcol=Pberrcol,
                     ycol=pars$ycol,yerrcol=pars$yerrcol,
                     lambda=pars$lambda)

    fit <- IsoplotR:::MLyork(raw,anchor=c(2,1),model=1)

    scaled <- data.frame(tPbPb=raw[,'X'],stPbPb=raw[,'sX'],
                         ty=raw[,'Y']-fit$a[1],sty=raw[,'sY'])
    
    wy <- get_w(X=scaled$ty,sX=scaled$sty,Y=scaled$tPbPb,sY=scaled$stPbPb)

    plot(x=scaled$tPbPb,y=scaled$ty,pch=16,
         xlim=rev(range(c(scaled$tPbPb+2*scaled$stPbPb,
                          scaled$tPbPb-2*scaled$stPbPb))),
         ylim=rev(range(c(scaled$ty+2*scaled$sty,scaled$ty-2*scaled$sty))),
         xlab=expression('t('^207*'Pb/'^206*'Pb)'),
         ylab=bquote(paste(.(signif(-fit$a[1],5)) + 
                           .(pars$ylab[[1]]),""[0]*"/"*lambda[.(pars$n)])),
         asp=1,bty='n')
    mtext(text=bquote('MSWD' == .(signif(fit$mswd,3)) *
                      ", p("*chi^2 *")="* .(signif(fit$p.value,2))),
          line=1.0,cex=0.8)
    mtext(text=bquote(.(pars$ySSlab[[1]]) == .(signif(wy[1]*pars$lambda,2))
                      %+-% .(signif(100*wy[2],2)) * "%"),
          line=-0.5,cex=0.8)
    abline(a=0,b=1)
    arrows(x0=scaled$tPbPb,x1=scaled$tPbPb,
           y0=scaled$ty-2*scaled$sty,
           y1=scaled$ty+2*scaled$sty,code=3,angle=90,len=0.05)
    arrows(x0=scaled$tPbPb-2*scaled$stPbPb,
           x1=scaled$tPbPb+2*scaled$stPbPb,
           y0=scaled$ty,y1=scaled$ty,code=3,angle=90,len=0.05)

    invisible(scaled)
}

Xcalibration <- function(prefix="Desch2",
                         Pbcol='tPbPb (original)',
                         Pberrcol='2se (tPbPb, original)',
                         ofile=paste0(prefix,'.pdf')){
    ifile <- paste0(prefix,'.csv')
    dat <- read_data(fn=ifile)
    pdf(file=ofile,width=9,height=3.5)
    op <- par(mfrow=c(1,3),mgp=c(2.5,1,0),mar=c(4,4,3,1),oma=c(0,0,2,0))
    yd_Al <- statchron(dat,Pbcol=Pbcol,Pberrcol=Pberrcol)
    legend('topleft',legend='a)',bty='n',xpd=NA,cex=1.2,adj=c(2,-3))
    yd_Mn <- statchron(dat,nuclide='Mn53',Pbcol=Pbcol,Pberrcol=Pberrcol)
    legend('topleft',legend='b)',bty='n',xpd=NA,cex=1.2,adj=c(2,-3))
    yd_Hf <- statchron(dat,nuclide='Hf182',Pbcol=Pbcol,Pberrcol=Pberrcol)
    legend('topleft',legend='c)',bty='n',xpd=NA,cex=1.2,adj=c(2,-3))
    yd <- rbind(yd_Al,yd_Mn,yd_Hf)
    w_x <- get_w(X=yd$tPbPb,sX=yd$stPbPb,Y=yd$ty,sY=yd$sty)
    mtext(bquote(omega[t(Pb)] ==
                 .(signif(w_x[1], 2)) ~ 'Myr'
                 %+-% .(signif(100*w_x[2],2)) * "%"), 
          side=3, line=0, outer=TRUE, cex=0.8)
    par(op)
    dev.off()
}

Xcalibration("Desch2",Pbcol='tPbPb (original)',Pberrcol='2se (tPbPb original)',ofile='Desch2_original.pdf')
Xcalibration("Desch2",Pbcol='tPbPb (adjusted)',Pberrcol='2se (tPbPb adjusted)',ofile='Desch2_adjusted.pdf')
