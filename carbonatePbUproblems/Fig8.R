source('ASH15.R')

op <- par(mfcol=c(2,3),mgp=c(1.75,0.5,0),mar=c(3,3,3,0.5))
plot_prior(ASH15diseq$d$U48)
plot_likelihood(ASH15diseq$d$U48)
ludwig(ASH15diseq,type=1,plot=TRUE,add=TRUE)
isochron(ASH15diseq,joint=FALSE,type=1,model=3,taxis=TRUE)
isochron(ASH15eq,joint=FALSE,type=2,taxis=TRUE)

dev.copy2pdf(file='output/Fig8.pdf',width=9,height=6)

ASH15uprior <- ASH15diseq
ASH15uprior$d$U48[c('m','M','sd')] <- c(0,20,100)
fit <- ludwig(ASH15uprior,model=3,type=1,nsteps=100)
tci <- IsoplotR:::bayesci(fit$posterior$t)
t48i <- IsoplotR:::bayesci(fit$posterior$U48i)
message('ASH15 06/38 isochron with uniform prior: ',
        't=',signif(fit$par['t'],4),
        ', ci(t)=[',signif(tci[1],4),',',signif(tci[2],4),']',
        ', [4/8]i=',signif(fit$par['U48i'],4),
        ', ci([4/8]i)=[',signif(t48i[1],4),',',signif(t48i[2],4),']')
