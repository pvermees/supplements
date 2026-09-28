library(IsoplotR)
setwd('/home/pvermees/Dropbox/Amelin')

read_data <- function(fn='ESSchronometry.csv'){
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

# get horizontal overdispersion
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

Xcalibration <- function(prefix="ESSchronometry",
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

DeschCalculator <- function(){
    dat <- read.csv('Desch1_table2.csv',check.names=TRUE)[1:7,]
    x <- dat[,2]
    sx <- dat[,3]/2
    y <- dat[,4]
    sy <- dat[,5]/2
    m <- (x/sx^2 + y/sy^2)/(1/sx^2 + 1/sy^2)
    X2 <- ((x-m)/sx)^2 + ((y-m)/sy)^2
    MSWD <- sum(X2)/length(x)
    MSWD
}

DeschSimulator <- function(){
    n <- 1000 # number of iterations    
    A <- 10 # number of meteorites
    s26 <- 1 # uncertainty of t26
    s53 <- 0.5 # uncertainty of t53
    s182 <- 1.5 # uncertainty of t182
    sPb <- 2 # uncertainty of tPb
    X2 <- rep(NA,n) # initialise
    for (i in 1:n){ # loop through the iterations
        # random values from standard normal distributions:
        t26 <- rnorm(A,sd=s26)
        t53 <- rnorm(A,sd=s53)
        t182 <- rnorm(A,sd=s182)
        tPb <- rnorm(A,sd=sPb)
        # weighted mean:
        dt <- (t26/s26^2 + t53/s53^2 + t182/s182^2 + tPb/sPb^2)/
            (1/s26^2 + 1/s53^2 + 1/s182^2 + 1/sPb^2)
        # chi-squared statistic:
        X2[i] <- sum(((t26-dt)/s26)^2 + ((t53-dt)/s53)^2 +
                     ((t182-dt)/s182)^2 + ((tPb-dt)/sPb)^2)
    }
    hist(X2,probability=TRUE,main='',xlab=expression(chi^2))
    usr <- par('usr')
    x2 <- seq(from=usr[1],to=usr[2],length.out=50)
    lines(x2,dchisq(x2,df=3*A),xpd=NA)
    dev.copy2pdf(file='X2.pdf')
}

Ycalibration <- function(dat,
                         nuclideX="Al26",
                         nuclideY="Mn53"){
    parsX <- get_pars(nuclideX)
    parsY <- get_pars(nuclideY)
    lX <- parsX$lambda
    lY <- parsY$lambda
    xlab <- bquote(.(parsX$ylab[[1]])[0] * "/" * lambda[.(parsX$n)])
    ylab <- bquote(.(parsY$ylab[[1]])[0] * "/" * lambda[.(parsY$n)])
    yd <- data.frame(X=log(dat[,1])/lX,sX=dat[,2]/(2*dat[,1]*lX),
                     Y=log(dat[,3])/lY,sY=dat[,4]/(2*dat[,3]*lY),
                     rXY=0)
    yd <- yd[complete.cases(yd),]
    yfit <- IsoplotR:::MLyork(yd,anchor=c(2,1),model=1)
    scatterplot(yd,show.ellipses=2,xlab=xlab,ylab=ylab,asp=1,
                xlim=range(c(yd[,1]-2*yd[,2],yd[,1]+2*yd[,2])),
                ylim=range(c(yd[,3]-2*yd[,4],yd[,3]+2*yd[,4])))
    abline(a=yfit$a[1],b=1)
    mtext(bquote('MSWD' == .(signif(yfit$mswd,2)) ~
                     ', p('*chi^2*')' == .(signif(yfit$p.value,2))),
          line=0.5,cex=0.8)
}

extinct <- function(){
    pdf(file='extinct.pdf',width=9,height=3.5)
    dat <- read_data(fn='ESSchronometry.csv')
    op <- par(mfrow=c(1,3),mgp=c(2.2,1,0),mar=c(3.5,3.5,2,0.5))
    Ycalibration(dat[,3:6],nuclideX="Al26",nuclideY="Mn53")
    Ycalibration(dat[,c(3,4,7,8)],nuclideX="Al26",nuclideY="Hf182")
    Ycalibration(dat[,5:8],nuclideX="Mn53",nuclideY="Hf182")
    par(op)
    dev.off()
}

#Xcalibration()
Xcalibration("Desch2",Pbcol='tPbPb (original)',Pberrcol='2se (tPbPb original)',ofile='Desch2_original.pdf')
Xcalibration("Desch2",Pbcol='tPbPb (adjusted)',Pberrcol='2se (tPbPb adjusted)',ofile='Desch2_adjusted.pdf')
#DeschCalculator()
#DeschSimulator()
#extinct()
