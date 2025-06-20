source ("AV03.R")

op <- par(mfrow=c(2,2),mgp=c(1.75,0.5,0),mar=c(3,3,0.5,0.5))
ludwig(AV03,type=1,plot=TRUE,add=TRUE)
ludwig(AV03eq,type=1,plot=TRUE,add=TRUE)
par(op)

dev.copy2pdf(file='output/Fig5.pdf',width=6,height=6)
