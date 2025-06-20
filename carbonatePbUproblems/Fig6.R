op <- par(mfrow=c(1,2),mgp=c(2,1,0),mar=c(3,3.5,0,0.75))

U48merr <- 0.002 # precision of the 234U/238U activity ratio measurements, from Walker (2006)

# panel a)

nt <- 150
tt <- seq(from=0.1,to=3.5,length.out=nt)
maxerr <- 100 # y axis limit
plot(range(tt),c(-25,maxerr),type='n',xlab=expression("t"[c]*"(Ma)"),
     ylab=expression("Magnitude of the correction ([t"[r]*"-t"[c]*"]/t"[r]*", %)"),bty='n',log='')
legend('topright','a)',bty='n')
U48i <- c(0.1,0.5,1,2,4,8,16)
ThU <- list(x=0,sx=0,option=1)
RaU <- list(x=0,sx=0,option=1)
for (i in seq_along(U48i)){
    d <- diseq(U48=list(x=U48i[i],sx=0,option=1),ThU=ThU,RaU=RaU)
    terr <- traw <- U48m <- truncated <- rep(NA,nt)
    for (j in seq_along(tt)){
        McL <- mclean(tt[j],d=d)
        traw[j] <- age(McL$Pb206U238,method="U238-Pb206")[1]
        terr[j] <- 100*(traw[j]-tt[j])/traw[j]
    }
    lines(tt,terr)
    last <- rev(which(is.finite(terr)))[1]
    text(tt[last],terr[last],U48i[i],pos=4,xpd=NA)
}
label <- expression("["^{234}*"U/"^{238}*"U]"[i])
text(tt[last]-0.1,terr[last]+5,label,pos=3,xpd=NA)

# panel b)

maxerr <- 140 # y axis limit
plot(range(tt),c(1,maxerr),type='n',xlab=expression("t"[c]*"(Ma)"),
     ylab=expression("Uncertainty of the correction ([t"[u]*"-t"[l]*"]/t"[c]*", %)"),bty='n',log='y')
legend('topright','b)',bty='n')
U48i <- c(0.1,0.25,0.5,1,2,4,8,16)
for (i in seq_along(U48i)){
    d <- diseq(U48=list(x=U48i[i],sx=0,option=1),ThU=ThU,RaU=RaU)
    terr <- tu <- tl <- U48m <- truncated <- rep(NA,nt)
    for (j in seq_along(tt)){
        McL <- mclean(tt[j],d=d)
        U48m[j] <- ifelse(U48i[i]==1,1.000055,McL$U48)
        Pb6U8 <- McL$Pb206U238
        dl <- diseq(U48=list(x=U48m[j]+U48merr,sx=0,option=2),
                    ThU=ThU,RaU=RaU)
        du <- diseq(U48=list(x=U48m[j]-U48merr,sx=0,option=2),
                    ThU=ThU,RaU=RaU)
        tl[j] <- age(Pb6U8,method="U238-Pb206",d=dl)[1]
        tu[j] <- age(Pb6U8,method="U238-Pb206",d=du)[1]
        McLmin <- mclean(tl[j],d=dl)
        McLmax <- mclean(tu[j],d=du)
        truncated[j] <- (McLmin$truncated | McLmax$truncated)
        if (!truncated[j]){
            terr[j] <- 100*(tu[j]-tl[j])/tt[j]
        }
    }
    lines(tt,terr)
    last <- rev(which(is.finite(terr)))[1]
    if (i==4){
        label <- expression("["^{234}*"U/"^{238}*"U]"[i]*"=1")
    } else {
        label <- U48i[i]
    }
    text(tt[last],terr[last],label,pos=2)

}

# worst case scenario envelope
dm <- diseq(U48=list(x=1,sx=0,option=1),ThU=ThU,RaU=RaU)
dM <- diseq(U48=list(x=12,sx=0,option=1),ThU=ThU,RaU=RaU)
envelope <- NA*tt
for (j in seq_along(tt)){
    tm <- age(age2ratio(tt=tt[j],st=0,ratio="Pb206U238",d=dm)[1],
              method="U238-Pb206")[1]
    tM <- age(age2ratio(tt=tt[j],st=0,ratio="Pb206U238",d=dM)[1],
              method="U238-Pb206")[1]
    envelope[j] <- 100*(tM-tm)/tt[j]
}
lines(tt,envelope,lty=2,xpd=NA)

dev.copy2pdf(file='output/Fig6.pdf',width=9,height=4.5)

par(op)
