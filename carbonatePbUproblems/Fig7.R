tab <- data.matrix(read.csv('data/Vaks.csv',check.names=FALSE,row.names=1))
ns <- nrow(tab)

UPb119 <- tab[,c('U8Pb6','s[U8Pb6]','Pb86','s[Pb86]','rho[U8Pb6,Pb86]')]
UPb1210 <- tab[,c('U5Pb7','s[U5Pb7]','Pb87','s[Pb87]','rho[U5Pb7,Pb87]')]
l4 <- settings('lambda','U234')[1]*1000
l8 <- settings('lambda','U238')[1]
U48 <- tab[,c('U48','s[U48]')]*l4/l8
Pb86c <- tab[,c('Pb86c','s[Pb86c]')]
Pb87c <- tab[,c('Pb87c','s[Pb87c]')]

tU8Pb6 <- tU5Pb7 <- matrix(NA,ns,4)
colnames(tU8Pb6) <- colnames(tU5Pb7) <- c('t','ll','ul','uncorr')
U48i <- matrix(NA,ns,3)
colnames(U48i) <- c('U48i','ll','ul')
rownames(tU8Pb6) <- rownames(tU5Pb7) <- rownames(U48i) <- rownames(tab)

## uncomment to use two permil uncertainty for the U234/U238 activity ratio
#    U48[,'s[U48]'] <- 0.002

for (i in 1:ns){
    print(i)
    settings('iratio','Pb206Pb208',1/Pb86c[i,'Pb86c'])
    settings('iratio','Pb207Pb208',1/Pb87c[i,'Pb87c'])
    try({    
        d <- diseq(ThU=list(x=0,sx=0,option=1))
        dat <- as.UPb(UPb119[i,,drop=FALSE],format=119,d=d)
        fit <- isochron(x=dat,model=1,anchor=1,plot=FALSE)
        tU8Pb6[i,'uncorr'] <- fit$age[1]
        if (is.na(U48[i,'U48'])){
            tU8Pb6[i,'t'] <- fit$age[1]
            tU8Pb6[i,'ll'] <- exp(log(fit$age[1])-2*fit$age[2]/fit$age[1])
            tU8Pb6[i,'ul'] <- exp(log(fit$age[1])+2*fit$age[2]/fit$age[1])
        } else {
            d <- diseq(U48=list(x=U48[i,'U48'],sx=U48[i,'s[U48]'],option=2),
                       ThU=list(x=0,sx=0,option=1))
            dat <- as.UPb(UPb119[i,,drop=FALSE],format=119,d=d)
            fit <- isochron(x=dat,model=1,anchor=1,plot=FALSE)
            tU8Pb6[i,'t'] <- fit$par['t']
            tU8Pb6[i,c('ll','ul')] <- IsoplotR:::bayesci(XL=fit$posterior$t)
            U48i[i,'U48i'] <- fit$par['U48i']
            U48i[i,c('ll','ul')] <- IsoplotR:::bayesci(XL=fit$posterior$U48i)
        }
        d <- diseq(PaU=list(x=0,sx=0,option=1))
        dat <- as.UPb(UPb1210[i,,drop=FALSE],format=1210,d=d)
        fit <- isochron(x=dat,model=1,anchor=1,plot=FALSE)
        tU5Pb7[i,'uncorr'] <- fit$age[1]
        tU5Pb7[i,'ll'] <- exp(log(fit$age[1])-2*fit$age[2]/fit$age[1])
        tU5Pb7[i,'ul'] <- exp(log(fit$age[1])+2*fit$age[2]/fit$age[1])
    })
}

gap <- 0.15
op <- par(mfrow=c(2,1),mar=c(3,3,0.25,0.25),mgp=c(1.75,1,0))
plot(x=c(0,ns),y=c(0,max(tU8Pb6[,'uncorr'],na.rm=TRUE)),
     type='n',xlab='n',ylab='t (Ma)')
points(x=(1:ns)-gap,y=tU8Pb6[,'uncorr'],cex=0.5)
arrows(x0=(1:ns)-gap,y0=tU8Pb6[,'ll'],
       x1=(1:ns)-gap,y1=tU8Pb6[,'ul'],
       code=3,length=0.025,angle=90)
arrows(x0=gap+(1:ns),y0=tU5Pb7[,'ll'],
       x1=gap+(1:ns),y1=tU5Pb7[,'ul'],
       code=3,length=0.025,angle=90,col='blue')
text(x=5,y=1.9,labels='5')
arrows(x0=5,y0=1.5,x1=5,y1=1.75,code=1,length=0.05)
text(x=41,y=1.0,labels='41')
arrows(x0=41,y0=1.35,x1=41,y1=1.1,code=1,length=0.05)
legend('topleft','a)',bty='n',adj=c(1,0.5))
plot(x=c(0,ns),y=c(0,max(U48i[,'ul'],na.rm=TRUE)),
     type='n',xlab='n',ylab=expression(''^'234'*'U/'^'238'*'U'))
points(x=1:ns,y=U48[,'U48'],cex=0.5)
arrows(x0=1:ns,y0=U48i[,'ll'],x1=1:ns,y1=U48i[,'ul'],
       code=3,length=0.025,angle=90)
legend('topleft','b)',bty='n',adj=c(1,0.5))
par(op)

dev.copy2pdf(file='output/Fig7.pdf',width=7,height=6)
