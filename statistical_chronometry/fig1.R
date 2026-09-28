ta <- 4567
tb <- seq(from=ta-1,to=1,length.out=50)
U <- 1/137.818
l8 <- 0.000155125
sl8 <- 0.000000083
l5 <- 0.00098485
sl5 <- 0.00000067
rho <- 0
vRa <- vRb <- 0
vdt <- vRb/((U*l8*(exp(l5*tb)-1)*exp(l8*tb))/(exp(l8*tb)-1)^2-(U*l5*exp(l5*tb))/(exp(l8*tb)-1))^2+vRa/((U*l8*(exp(l5*ta)-1)*exp(l8*ta))/(exp(l8*ta)-1)^2-(U*l5*exp(l5*ta))/(exp(l8*ta)-1))^2+((U*tb*(exp(l5*tb)-1)*exp(l8*tb))/((exp(l8*tb)-1)^2*((U*l8*(exp(l5*tb)-1)*exp(l8*tb))/(exp(l8*tb)-1)^2-(U*l5*exp(l5*tb))/(exp(l8*tb)-1)))-(U*ta*(exp(l5*ta)-1)*exp(l8*ta))/((exp(l8*ta)-1)^2*((U*l8*(exp(l5*ta)-1)*exp(l8*ta))/(exp(l8*ta)-1)^2-(U*l5*exp(l5*ta))/(exp(l8*ta)-1))))*(sl8^2*((U*tb*(exp(l5*tb)-1)*exp(l8*tb))/((exp(l8*tb)-1)^2*((U*l8*(exp(l5*tb)-1)*exp(l8*tb))/(exp(l8*tb)-1)^2-(U*l5*exp(l5*tb))/(exp(l8*tb)-1)))-(U*ta*(exp(l5*ta)-1)*exp(l8*ta))/((exp(l8*ta)-1)^2*((U*l8*(exp(l5*ta)-1)*exp(l8*ta))/(exp(l8*ta)-1)^2-(U*l5*exp(l5*ta))/(exp(l8*ta)-1))))+rho*sl5*sl8*((U*ta*exp(l5*ta))/((exp(l8*ta)-1)*((U*l8*(exp(l5*ta)-1)*exp(l8*ta))/(exp(l8*ta)-1)^2-(U*l5*exp(l5*ta))/(exp(l8*ta)-1)))-(U*tb*exp(l5*tb))/((exp(l8*tb)-1)*((U*l8*(exp(l5*tb)-1)*exp(l8*tb))/(exp(l8*tb)-1)^2-(U*l5*exp(l5*tb))/(exp(l8*tb)-1)))))+((U*ta*exp(l5*ta))/((exp(l8*ta)-1)*((U*l8*(exp(l5*ta)-1)*exp(l8*ta))/(exp(l8*ta)-1)^2-(U*l5*exp(l5*ta))/(exp(l8*ta)-1)))-(U*tb*exp(l5*tb))/((exp(l8*tb)-1)*((U*l8*(exp(l5*tb)-1)*exp(l8*tb))/(exp(l8*tb)-1)^2-(U*l5*exp(l5*tb))/(exp(l8*tb)-1))))*(rho*sl5*sl8*((U*tb*(exp(l5*tb)-1)*exp(l8*tb))/((exp(l8*tb)-1)^2*((U*l8*(exp(l5*tb)-1)*exp(l8*tb))/(exp(l8*tb)-1)^2-(U*l5*exp(l5*tb))/(exp(l8*tb)-1)))-(U*ta*(exp(l5*ta)-1)*exp(l8*ta))/((exp(l8*ta)-1)^2*((U*l8*(exp(l5*ta)-1)*exp(l8*ta))/(exp(l8*ta)-1)^2-(U*l5*exp(l5*ta))/(exp(l8*ta)-1))))+sl5^2*((U*ta*exp(l5*ta))/((exp(l8*ta)-1)*((U*l8*(exp(l5*ta)-1)*exp(l8*ta))/(exp(l8*ta)-1)^2-(U*l5*exp(l5*ta))/(exp(l8*ta)-1)))-(U*tb*exp(l5*tb))/((exp(l8*tb)-1)*((U*l8*(exp(l5*tb)-1)*exp(l8*tb))/(exp(l8*tb)-1)^2-(U*l5*exp(l5*tb))/(exp(l8*tb)-1)))))

vR <- 0
tt <- tb
vta <- vR/((U*l8*(exp(l5*tt)-1)*exp(l8*tt))/(exp(l8*tt)-1)^2-(U*l5*exp(l5*tt))/(exp(l8*tt)-1))^2-(U*tt*(exp(l5*tt)-1)*exp(l8*tt)*((U*rho*sl5*sl8*tt*exp(l5*tt))/((exp(l8*tt)-1)*((U*l8*(exp(l5*tt)-1)*exp(l8*tt))/(exp(l8*tt)-1)^2-(U*l5*exp(l5*tt))/(exp(l8*tt)-1)))-(U*sl8^2*tt*(exp(l5*tt)-1)*exp(l8*tt))/((exp(l8*tt)-1)^2*((U*l8*(exp(l5*tt)-1)*exp(l8*tt))/(exp(l8*tt)-1)^2-(U*l5*exp(l5*tt))/(exp(l8*tt)-1)))))/((exp(l8*tt)-1)^2*((U*l8*(exp(l5*tt)-1)*exp(l8*tt))/(exp(l8*tt)-1)^2-(U*l5*exp(l5*tt))/(exp(l8*tt)-1)))+(U*tt*exp(l5*tt)*((U*sl5^2*tt*exp(l5*tt))/((exp(l8*tt)-1)*((U*l8*(exp(l5*tt)-1)*exp(l8*tt))/(exp(l8*tt)-1)^2-(U*l5*exp(l5*tt))/(exp(l8*tt)-1)))-(U*rho*sl5*sl8*tt*(exp(l5*tt)-1)*exp(l8*tt))/((exp(l8*tt)-1)^2*((U*l8*(exp(l5*tt)-1)*exp(l8*tt))/(exp(l8*tt)-1)^2-(U*l5*exp(l5*tt))/(exp(l8*tt)-1)))))/((exp(l8*tt)-1)*((U*l8*(exp(l5*tt)-1)*exp(l8*tt))/(exp(l8*tt)-1)^2-(U*l5*exp(l5*tt))/(exp(l8*tt)-1)))

sdt <- sqrt(vdt)
sta <- sqrt(vta)
par(mgp=c(2,0.75,0),mar=c(3.5,3.5,0.5,3.5))
plot(x=tb,y=sdt,
     type='l',bty='n',lwd=2,
     xaxt='n',yaxt='n',
     xlab=NA,ylab=NA,
     xlim=c(0,ta),
     ylim=c(0,4.6),col='red')
par(fg='red')
xticks <- c(0,1000,2000,3000,4000,ta)
labels <- as.character(xticks)
labels[length(labels)] <- expression(t[SS])
axis(side=1,at=xticks,labels=labels)
mtext(side=1,line=2,text=expression(t~"[Ma]"))
axis(side=2,col.lab='red',col.axis='red')
mtext(side=2,line=2,text=expression(sigma(lambda~"|"~t[SS]-t)~"[Myr]"))
lines(x=tt,y=sta,lty=1,lwd=2,col='blue')
par(fg='blue')
axis(side=4,col.lab='blue',col.axis='blue')
mtext(side=4,line=2,text=expression(sigma(lambda~"|"~t)~"[Myr]"))

dev.copy2pdf(file='lambda_err.pdf',width=4.5,height=4.5)
