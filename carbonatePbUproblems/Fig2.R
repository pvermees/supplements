source("AV03.R")

MC(AV03)

op <- par(fig=c(0.0,0.5,0.5,1.0),new=TRUE,mar=rep(0,4),cex=0.6)
plot.window(xlim=c(0,1),ylim=c(0,1))
polygon(x=c(0.13,0.17,0.17,0.13),
        y=c(0.79,0.79,0.82,0.82),xpd=NA)
par(op)
op <- par(fig=c(0.25,0.48,0.70,0.90),new=TRUE,mar=rep(0,4),cex=0.6)
IsoplotR::scatterplot(AV03$x,fit=york(AV03$x))
par(op)

dev.copy2pdf(file='output/Fig2.pdf',width=6,height=7)
