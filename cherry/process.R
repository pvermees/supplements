rm(list=ls())
graphics.off()

## Figure 1 ##

nr <- 3
nc <- 4
m <- rbind(
    c(1,rep(2,nc),3),
    c(1,4+(1:nc),3),
    c(1,4+((nc+1):(2*nc)),3),
    c(1,4+((2*nc+1):(3*nc)),3),
    c(1,rep(4,nc),3)
)
par(mar=rep(0,4),mgp=c(1.4,0.5,0))
marginwidth <- 0.07
marginheight <- 0.09
delta <- 0.001
layout(m,widths=c(marginwidth-delta,rep((1-marginwidth)/nc,nc),delta),
       heights=c(delta,rep((1-marginheight)/nr,nr),marginheight-delta))
set.seed(5)
n <- 10
alpha <- 0.05
plot.new()
legend('topleft',legend='(a)',bty='n',cex=2,xpd=NA,inset=c(-0.85,-0.01))
plot.new()
plot.new()
plot.new()
for (i in 1:nr){
    for (j in 1:nc){
        X <- runif(10)
        Y <- runif(10)
        fit <- lm(Y ~ X)
        smry <- summary(fit)
        r2 <- signif(smry$r.squared,2)
        pval <- signif(smry$coefficients[2,4],2)
        col <- ifelse(pval<alpha,'black','white')
        plot(X,Y,xlim=c(0,1.0),ylim=c(0,1.3),asp=1,
             xaxt='n',yaxt='n',type='n')
        mtext(bquote('r'^2*'='*.(r2)*', p='*.(pval)),line=-1.8,cex=0.8)
        if (pval < alpha){
            lines(c(0,1),predict(fit,newdata=data.frame(X=c(0,1))))
        }
        points(X,Y,pch=21,bg=col)
        if (i==nr){
            axis(side=1,at=c(0,0.5,1))
            mtext('x',side=1,line=1.5)
        }
        if (j==1){
            axis(side=2,at=c(0,0.5,1))
            mtext('y',side=2,line=1.5)
        }
    }
}
dev.copy2pdf(file='output/type1.pdf',width=6,height=4)

## Figure 2 ##

graphics.off()

op <- par(mar=c(3.5,3,0.2,0.2),mgp=c(2,1,0),cex=0.8)
set.seed(20)
mu <- 100
sigma <- 10
ns <- 20
size <- sample(25:100,ns)
err <- 50/sqrt(size)
y <- rnorm(ns,mean=mu,sd=sqrt(sigma^2+err^2))
x <- 1:ns
y0 <- y - err
y1 <- y + err
avg_all <- weighted.mean(y,1/err^2)
mswd_all <- sum(((y-avg_all)/err)^2)/(ns-1)
matplot(c(0,ns),range(c(y0,y1)),type='n',
        bty='o',xlab='x',ylab='y',xpd=NA,cex=0.5)
lines(range(x),rep(avg_all,2))
dy <- (y-avg_all)^2
cherries <- order(dy)[1:10]
nc <- length(cherries)
avg_cherries <- weighted.mean(y[cherries],1/err[cherries]^2)
mswd_cherries <- sum(((y[cherries]-avg_cherries)/err[cherries])^2)/(nc-1)
col <- rep('red',ns)
col[cherries] <- 'blue'
arrows(x0=x,x1=x,y0=y0,y1=y1,len=0.1,angle=90,code=3,col=col)
usr <- par('usr')
mtext(text=paste0('MSWD=',signif(mswd_all,2),
                  ', n=',ns,
                  ', p=',signif(1-pchisq(mswd_all*ns,df=ns-1),2),
                  ' '),
      side=3,line=-2,cex=0.9,
      adj=c(1,1),col='red')
dfcherries <- length(cherries)-1      
mtext(text=paste0('MSWD=',signif(mswd_cherries,2),
                  ', n=',length(cherries),
                  ', p=',signif(1-pchisq(mswd_cherries*dfcherries,
                                         df=dfcherries),2),
                  ' '),
      side=3,line=-3.5,cex=0.9,
      adj=c(1,1),col='blue')
legend('topleft',legend='(b)',bty='n',xpd=NA,cex=1.7,inset=c(-0.14,-0.02))
par(op)
dev.copy2pdf(file='output/type2.pdf',width=6,height=4)

## Figure 3 ##

set.seed(1)
op <- par(mar=c(3.5,3,0.2,0.2),mgp=c(2,1,0),cex=0.8)
plot(x=c(0,1),y=c(0,1),type='l',xlab='p',ylab='Fn(p)',bty='n',asp=1)
text(x=0.13,y=0.1,labels='best case scenario (1:1)',pos=4,srt=45,col='black')
arrows(x0=0.65,x1=0.35,y0=0.75,y1=0.75,length=0.05)
text(x=0.55,y=0.85,labels='genuine')
text(x=0.55,y=0.8,labels='overdispersion')
text(x=0.45,y=0.7,labels='or type-1')
text(x=0.45,y=0.65,labels='cherry picking')
arrows(x0=0.825,x1=0.9,y0=0.775,y1=0.7,length=0.05)
text(x=0.8,y=0.65,labels="'too good",xpd=NA,pos=4)
text(x=0.8,y=0.6,labels="to be true'",xpd=NA,pos=4,offset=2/3)
arrows(x0=0.65,x1=0.95,y0=0.4,y1=0.4,length=0.05)
text(x=0.8,y=0.35,labels='type-2 cherry picking')
nr <- 100
x <- rchisq(nr,df=10,ncp=7)
p1 <- 1-pchisq(x,df=10)
p2 <- 1-c(pchisq(x[1:round(nr/2)],df=10),
          pchisq(x[round(nr/2+1):nr],df=26))
xp <- rep(knots(ecdf(p2)),each=2)
yp <- ecdf(p2)(xp)
forbidden <- (yp < xp)
nf <- sum(forbidden)
xf <- xp[forbidden][-1]
yf <- yp[forbidden][-nf]
yf0 <- yp[which(forbidden)[1]-1]
polygon(x=c(yf0,xf[1],xf,1),
        y=c(yf0,yf0,yf,1),col=rgb(1,0,0,0.5),border=NA)
lines(ecdf(p1),verticals=TRUE,pch=NA,col='blue',col.01line = NA)
lines(ecdf(p2),verticals=TRUE,pch=NA,col='red',col.01line = NA)
legend('bottomright',legend='(a)',bty='n',cex=1.2,inset=1/20)
par(op)

dev.copy2pdf(file='output/synthetic-ecdf.pdf',width=4,height=4)

## Figure 4 ##

graphics.off()

splitter <- function(string){
    split_result <- strsplit(string,"\\),\\(")
    split_result <- lapply(split_result, function(x) {
        gsub("^\\(|\\)$", "", x)
    })
    split_result[[1]]
}

iwano <- read.csv(file='Iwano2018.csv',check.names=FALSE,sep='\t')
pvals_Iwano <- as.numeric(iwano[,'P-value'])

tab <- read.csv(file='p-values_FT.csv',check.names=FALSE)
good <- !is.na(tab[,'p-value'])
pvals_Geochron <- as.numeric(tab[good,'p-value'])/100

tab <- read.delim('Nature.csv',sep='\t',
                  check.names=FALSE,header=TRUE)
ns <- nrow(tab)
titles <- authors <- sname <- MSWD <- method <- nn <- dof <- c()
for (i in 1:ns){
    entry <- tab[i,"(name,MSWD,type,n,df)"]
    title <- tab[i,'Title']
    author <- tab[i,'Authors']
    if (!is.na(entry)){
        split <- splitter(entry)
        for (string in split){
            titles <- append(titles,title)
            authors <- append(authors,author)
            vec <- strsplit(string,",")[[1]]
            sname <- append(sname,vec[1])
            MSWD <- append(MSWD,as.numeric(vec[2]))
            method <- append(method,vec[3])
            nn <- append(nn,as.numeric(vec[4]))
            dof <- append(dof,as.numeric(vec[5]))
        }
    }
}
X2 <- MSWD*dof
pvals_Nature <- pchisq(X2,df=dof,lower.tail=FALSE)

op <- par(mar=c(3.5,3,0.2,0.2),mgp=c(2,1,0),cex=0.8)
plot(x=c(0,1),y=c(0,1),type='l',xlab='p',ylab='Fn(p)',bty='n',asp=1)
lines(ecdf(pvals_Iwano/100),verticals=TRUE,pch=NA,col='blue',col.01line = NA)
lines(ecdf(pvals_Geochron),verticals=TRUE,pch=NA,col='red',col.01line = NA)
lines(ecdf(pvals_Nature),verticals=TRUE,pch=NA,col='black',col.01line = NA)
lines(x=0:1,y=0:1)
text(x=0.65,y=0.25,labels="'forbidden zone'",pos=1,srt=0,col='black')
legend(x=0,y=1,legend=c('Iwano et al. (2018)',
                        'Fission tracks (geochron.org)',
                        'Nature compilation'),
       lty=rep(1,3),pch=NA,col=c('blue','red','black'),
       bty='n',y.intersp=rep(2,3),xpd=NA)
legend('bottomright',legend='(b)',bty='n',cex=1.2,inset=1/20)
par(op)

dev.copy2pdf(file='output/ecdf.pdf',width=4,height=4)
