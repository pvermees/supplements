plot_prior <- function(a,x=seq(from=a$m,to=a$M,length.out=50)){
    y <- IsoplotR:::prior(x,a,log=FALSE)
    plot(x,y,type='b',
         xlab=expression('[4/8]'[i]*''),
         ylab='prior probability')
}

plot_likelihood <- function(a,
                            x=seq(from=a$x-3*a$sx,
                                  to=a$x+3*a$sx,
                                  length.out=30),
                            ...){
    y <- dnorm(x,mean=a$x,sd=a$sx,log=FALSE)
    plot(x,y,type='l',
         xlab=expression('[4/8]'[m]*''),
         ylab='likelihood',...)
}

plot_dens <- function(x,y,xlab=NA,...){
    nn <- length(x)
    dx <- diff(x)
    dx <- (c(dx[1]-(dx[2]-dx[1]),dx) +
           c(dx,dx[nn-1]+(dx[nn-1]-dx[nn-2])))/2
    good <- (dx!=0)
    DX <- dx[good]
    X <- x[good]
    if (missing(y)) dens <- abs(1/(nn*DX))
    else dens <- (y[good]/DX)/sum(y[good]/DX)
    plot(X,dens,type='b',xlab=xlab,ylim=c(0,1.2*max(dens)),...)
    rug(X)
}

MCplot <- function(UPb,U48m,tt,U48i){
    op <- par(mfrow=c(2,2),mgp=c(1.5,0.5,0),mar=c(3,3,4,0.5))
    good <- (U48i>0)
    col <- rep('black',length(tt))
    col[good] <- 'white'
    isochron(UPb,plot=TRUE)
    legend('topright','a)',bty='n',xpd=NA)
    plot_dens(U48m,xlab=expression('['^234*'U/'^238*'U]'[m]),
              pch=21,bg=col,bty='n',xpd=NA)
    legend('topleft','b)',bty='n',xpd=NA)
    plot_dens(U48i[good],xlab=expression('['^234*'U/'^238*'U]'[i]),
              pch=21,bg='white',bty='n')
    legend('topleft','c)',bty='n',xpd=NA)
    plot_dens(tt[good],xlab='t (Ma)',
              pch=21,bg='white',bty='n')
    legend('topright','d)',bty='n',xpd=NA)
    par(op)
}
