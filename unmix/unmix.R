setwd('/home/pvermees/git/supplements/unmix')
library(IsoplotR)
rm(list=ls())
graphics.off()
set.seed(1)

synthetic_data <- function(m=c(0,100),
                           s=100,
                           p=1/2,
                           n=100){
    AA <- rnorm(n=n,mean=m[1],sd=s)
    BB <- rnorm(n=n,mean=m[2],sd=s)
    nA <- round(n*p)
    nB <- round(n*(1-p))
    if (TRUE){
        CC <- c(sample(AA,size=nA),sample(BB,size=nB))
    } else {
        CC <- c(rnorm(n=nA,mean=m[1],sd=s),
                rnorm(n=nB,mean=m[2],sd=s))
    }
    out <- list(A=AA,B=BB,C=CC)
    class(out) <- append(class(out),'detritals')
    out
}


pdf(file='algorithm.pdf',width=8,height=5)
m <- c(100,200)
s <- 30
p <- 1/3
layout(mat=rbind(c(1,3,5,0),c(2,4,6,0)),
       widths=c(5,5,5,1))
op <- par(mar=c(3,3,0.5,0),mgp=c(2,1,0),cex=1)
mx <- 0
Mx <- 300
x <- mx:Mx
dA <- dnorm(mx:Mx,mean=m[1],sd=s)
dB <- dnorm(mx:Mx,mean=m[2],sd=s)
dC <- p*dA + (1-p)*dB
matplot(x=cbind(x,x),y=cbind(dA,dB),type='l',bty='n',
        xlab='t',ylab='f(t)',lty=c(1,1),
        col=c('red','blue'),lwd=2)
usr <- par('usr')
arrows(x0=m[1],x1=m[2],y0=usr[4],y1=usr[4],
       code=3,length=0.05,angle=45,xpd=NA)
text(x=mean(m),y=usr[4],labels=expression(Delta(t)),
     pos=1,cex=0.8)
arrows(x0=m[1]-s,x1=m[1]+s,
       y0=dnorm(m[1]-s,m[1],s),
       y1=dnorm(m[1]+s,m[1],s),
       code=3,length=0.05,angle=45,xpd=NA)
text(x=m[1],y=dnorm(m[1]+s,m[1],s),
     labels=expression(sigma),cex=0.8,
     pos=3,offset=0.2)
arrows(x0=m[2]-s,x1=m[2]+s,
       y0=dnorm(m[2]-s,m[2],s),
       y1=dnorm(m[2]+s,m[2],s),
       code=3,length=0.05,angle=45,xpd=NA)
text(x=m[2],y=dnorm(m[2]+s,m[2],s),
     labels=expression(sigma),cex=0.8,
     pos=3,offset=0.2)
legend('topleft',legend='a)',bty='n',adj=c(1,0))
plot(x=x,y=dC,type='l',bty='n',
     xlab='t',ylab='f(t)',lwd=2)
legend('topleft',legend='b)',bty='n',adj=c(1,0))
matplot(x=cbind(x,x),y=cbind(cumsum(dA),cumsum(dB)),
        type='l',bty='n',xlab='t',ylab='Fn(t)',
        lty=c(1,1),col=c('red','blue'),lwd=2)
legend('topleft',legend='c)',bty='n',adj=c(1,0))
plot(x=x,y=cumsum(dC),type='l',bty='n',
     xlab='t',ylab='Fn(t)',lwd=2)
lines(x=c(mx,Mx),y=rep(p,2),lty=2)
arrows(x0=Mx,x1=Mx,y0=0,y1=p,code=3,length=0.05,angle=45)
arrows(x0=Mx,x1=Mx,y0=p,y1=1,code=3,length=0.05,angle=45)
text(x=0.95*c(Mx,Mx),y=c(p/2,p+(1-p)/2),
     labels=c(expression(W[1]),
              expression(W[2] == 1-W[1])),
     pos=NULL,srt=90,cex=0.8)
legend('topleft',legend='d)',bty='n',adj=c(1,0))
DZ <- synthetic_data(m=m,s=s,p=p,n=100)
kde_A <- kde(DZ$A,from=mx,to=Mx,plot=FALSE)
kde_B <- kde(DZ$B,from=mx,to=Mx,plot=FALSE)
kde_C <- kde(DZ$C,from=mx,to=Mx,plot=FALSE)
i <- 36
ti <- sort(DZ$C)[i]
plot(x=c(mx,Mx),y=c(0,1),type='n',bty='n',
     xlab='t',ylab='Fn(t)',xlim=c(mx,Mx),lwd=2)
ecdfA <- ecdf(DZ$A)
ecdfB <- ecdf(DZ$B)
ecdfC <- ecdf(DZ$C)
lines(ecdfA,verticals=TRUE,pch=NA,
      col.01line=NA,col='red',lwd=2)
lines(ecdfB,verticals=TRUE,pch=NA,
      col.01line=NA,col='blue',lwd=2)
x1 <- grconvertX(x=ti,from="user",to="ndc")
y1A <- grconvertY(y=ecdfA(ti),from="user",to="ndc")
y1B <- grconvertY(y=ecdfB(ti),from="user",to="ndc")
x2 <- grconvertX(x=Mx*1.1,from="user",to="ndc")
text(x=Mx*1.1,y=ecdfA(ti),
     labels=expression(X["i,1"]),xpd=NA,pos=4)
text(x=Mx*1.1,y=ecdfB(ti),
     labels=expression(X["i,2"]),xpd=NA,pos=4)
legend('topleft',legend='e)',bty='n',adj=c(1,0))
plot(ecdfC,verticals=TRUE,pch=NA,bty='n',
     main='',xlab='t',ylab='Fn(t)',col.01line=NA,
     xlim=c(mx,Mx),lwd=2)
text(x=Mx*1.1,y=ecdfC(ti),
     labels=expression(Y[i]),xpd=NA,pos=4)
text(x=ti,y=0.05,labels=expression(t[i]),
     xpd=NA,pos=4,offset=0.2)
x0 <- grconvertX(x=ti,from="user",to="ndc")
y0 <- grconvertY(y=0,from="user",to="ndc")
y1C <- grconvertY(y=ecdfC(ti),from="user",to="ndc")
rug(DZ$C)
legend('topleft',legend='f)',bty='n',adj=c(1,0))
par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
segments(x0=c(x0,x1,x1,x0),
         x1=c(x1,x2,x2,x2),
         y0=c(y0,y1A,y1B,y1C),
         y1=c(y1A,y1A,y1B,y1C))
par(op)
dev.off()

pdf(file='bootstrap.pdf',width=4,height=4)
op <- par(mar=c(3,3,0.5,0.5),mgp=c(2,1,0))
dt <- seq(from=0,to=3,length.out=10)
ll <- ul <- NULL
p <- 1/3
for (i in seq_along(dt)){
    DZ <- synthetic_data(m=c(0,dt[i]),s=1,p=p)
    fit <- unmix(DZ,source_names=c('A','B'),
                 boot=TRUE,plot=FALSE)
    ll <- c(ll,fit$ll[1])
    ul <- c(ul,fit$ul[1])
}
plot(x=range(dt),y=c(0,1),type='n',
     xlab=expression(Delta(t)),
     ylab=expression(W[1]),bty='n')
arrows(x0=dt,x1=dt,y0=ll,y1=ul,
       code=3,length=0.05,angle=90)
abline(h=p,lty=2)
par(op)
dev.off()

pdf(file='Sickmann.pdf',width=5,height=7)
DZ <- read.data(x='Sickmann.csv',method='detritals')
unmix(DZ,source_names=c('P1','P2','P3'),xlim=c(50,300))
dev.off()
