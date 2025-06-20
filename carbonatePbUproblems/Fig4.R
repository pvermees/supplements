source ("Corchia.R")

op <- par(mfrow=c(2,2),mgp=c(1.75,0.5,0),mar=c(3,3,0.5,0.5))
plot_prior(Corchia$d$U48)
plot_likelihood(Corchia$d$U48)
ludwig(Corchia,type=1,plot=TRUE,add=TRUE)
par(op)

dev.copy2pdf(file='output/Fig4.pdf')
