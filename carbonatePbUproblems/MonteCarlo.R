getquantiles <- function(nn,option='norm',mean=0,sd=1,min=0,max=20){
    if (identical(option,'norm')){
        ptile <- seq(from=1/nn,to=(nn-1)/nn,length.out=nn)
        out <- qnorm(p=ptile,mean=mean,sd=sd)
    } else if (identical(option,'uniform')){
        out <- seq(from=min,to=max,length.out=nn)
    } else {
        stop('Invalid getquantiles option')
    }
}

quickfit <- function(UPb){
    misfit <- function(tt,d,yfit){
        McL <- mclean(tt,d)
        x <- 1/McL$Pb206U238
        y <- McL$Pb207Pb206
        log(abs(y-yfit$a[1]-yfit$b[1]*x)/sqrt(1+yfit$b[1]^2))
    }
    yfit <- york(UPb$x)
    fit <- optimize(misfit,lower=0,upper=10,d=UPb$d,yfit=yfit,
                    tol = .Machine$double.eps^0.5)
    fit
}

MC <- function(UPb,nn=50){
    U48 <- UPb$d$U48$x
    sU48 <- UPb$d$U48$sx
    U48m <- getquantiles(nn,option='norm',mean=U48,sd=sU48)
    tt <- U48i <- rep(NA,nn)
    for (i in 1:nn){
        dat <- UPb
        dat$d <- diseq(U48=list(x=U48m[i],option=2),
                       ThU=list(x=0,option=1),
                       RaU=list(x=0,option=1),
                       PaU=list(x=0,option=1))
        tt[i] <- quickfit(dat)$minimum
        McL <- mclean(tt[i],dat$d)
        U48i[i] <- McL$U48i
    }
    MCplot(UPb,U48m,tt,U48i)
}
