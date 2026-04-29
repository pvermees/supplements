graphics.off()

fissiontracks <- function(){
    iwano <- read.csv(file='../data/Iwano2018.csv',check.names=FALSE,sep='\t')
    pvals_Iwano <- as.numeric(iwano[,'P-value'])

    tab <- read.csv(file='../data/p-values_FT.csv',check.names=FALSE)
    good <- !is.na(tab[,'p-value'])
    pvals_Geochron <- as.numeric(tab[good,'p-value'])/100
    
    pdf(file='../figures/Iwano.pdf',width=4,height=4)
    op <- par(mar=c(3,3,0.5,0.5),mgp=c(2,1,0))
    colours <- c('#E69F00','#56B4E9')
    plot(x=c(0,1),y=c(0,1),type='l',xlab='p',ylab='Fn(p)',bty='n')
    lines(ecdf(pvals_Iwano/100),verticals=TRUE,main=NA,pch=NA,
          col=colours[1],col.01line = NA)
    lines(ecdf(pvals_Geochron),verticals=TRUE,pch=NA,
          col=colours[2],col.01line = NA)
    legend(x=0,y=1,legend=c('Iwano et al. (2018)',
                            'geochron.org'),
           lty=rep(1,2),pch=NA,col=colours,
           bty='n',y.intersp=rep(1,2),xpd=NA)
    par(op)
    dev.off()

}

Nature <- function(){
    splitter <- function(string){
        split_result <- strsplit(string,"\\),\\(")
        split_result <- lapply(split_result, function(x) {
            gsub("^\\(|\\)$", "", x)
        })
        split_result[[1]]
    }

    tab <- read.delim('../data/Nature.csv',sep='\t',
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
    is.average <- grepl("average",method) | grepl("mean",method)
    is.plateau <- grepl("plateau",method)
    is.isochron <- grepl("isochron",method) |
        grepl("regression",method) | grepl("discordia",method)
    is.ArAr <- grepl("Ar-Ar",method) | grepl("plateau",method)
    is.UPb <- grepl("U-Pb",method) |
        grepl("concordia",method) | grepl("discordia",method)
    X2 <- MSWD*dof
    pvals_Nature <- pchisq(X2,df=dof,lower.tail=FALSE)

    pdf(file='../figures/Nature.pdf',width=8,height=4)
    op <- par(mfrow=c(1,2),mar=c(3,3,0.5,0.5),mgp=c(2,1,0))
    # panel a
    plot(x=c(0,1),y=c(0,1),type='l',xlab='p',ylab='Fn(p)')
    selectors <- c('UPb','ArAr')
    colours <- c('#E69F00','#56B4E9','#000000')
    for (i in seq_along(selectors)){
        selector <- selectors[i]
        p <- pvals_Nature[get(paste0("is.",selector))]
        cdf <- ecdf(p)
        lines(cdf,pch=NA,verticals=TRUE,col=colours[i])
        np <- length(p)
        message(selector,': n=',np,', f=',1-which.min(cdf(sort(p))>sort(p))/np)
    }
    leg <- c(expression(''^206*'Pb/'^238*'U'),
             expression(''^40*'Ar/'^39*'Ar'),
             'all')
    legend(x=0,y=1,legend=leg,lty=rep(1,2),pch=NA,
           col=colours,bty='n',y.intersp=rep(1,2),xpd=NA)
    legend('bottomright','a)',bty='n')
    # panel b
    plot(x=c(0,1),y=c(0,1),type='l',xlab='p',ylab='Fn(p)')
    selectors <- c('average','plateau','isochron')
    colours <- c('#E69F00','#56B4E9','#009E73','#000000')
    for (i in seq_along(selectors)){
        selector <- selectors[i]
        p <- pvals_Nature[get(paste0("is.",selector))]
        cdf <- ecdf(p)
        lines(cdf,pch=NA,verticals=TRUE,col=colours[i])
        np <- length(p)
        message(selector,': n=',np,', f=',1-which.min(cdf(sort(p))>sort(p))/np)
    }
    leg <- c('average','plateau','isochron','all')
    legend(x=0,y=1,legend=leg,lty=rep(1,4),pch=NA,
           col=colours,bty='n',y.intersp=rep(1,4),xpd=NA)
    legend('bottomright','b)',bty='n')
    par(op)
    dev.off()
}

GTS <- function(){
    tab_GTS <- read.csv('../data/GTS.csv')
    X2 <- tab_GTS$MSWD * tab_GTS$df
    good <- !is.na(X2)
    pvals_GTS <- 1-pchisq(X2[good],df=tab_GTS$df[good])

    pdf(file='../figures/GTS.pdf',width=4,height=4)
    op <- par(mar=c(3,3,0.5,0.5),mgp=c(2,1,0))
    plot(x=c(0,1),y=c(0,1),type='l',xlab='p',ylab='Fn(p)',bty='n')
    lines(ecdf(pvals_GTS),verticals=TRUE,pch=NA,col='red')
    is.UPb <- (tab_GTS[,'age_type'] == '206Pb/238U')
    is.ArAr <- (tab_GTS[,'age_type'] == '40Ar/39Ar')
    lines(ecdf(pvals_GTS[is.UPb]),verticals=TRUE,pch=NA,col='#E69F00')
    lines(ecdf(pvals_GTS[is.ArAr]),verticals=TRUE,pch=NA,col='#56B4E9')
    lines(ecdf(pvals_GTS),verticals=TRUE,pch=NA,col='black')
    leg <- c(expression(''^206*'Pb/'^238*'U'),
             expression(''^40*'Ar/'^39*'Ar'),
             'all')
    colours <- c('#E69F00','#56B4E9','#000000')
    legend(x=0,y=1,legend=leg,lty=rep(1,3),pch=NA,
           col=colours,bty='n',y.intersp=rep(1,3),xpd=NA)
    dev.off()
}

fissiontracks()

Nature()

GTS()
