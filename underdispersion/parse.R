rm(list=ls())

library(readxl)
fnames <- list.files('geochron.org',full.names=TRUE)
nf <- length(fnames)
tab <- matrix(NA,nrow=nf,ncol=4)
colnames(tab) <- c('IGSN','dosimeter','irradiation','p-value')
for (i in seq_along(fnames)){
    xlsx_data <- read_excel(fnames[i])
    rnames <- xlsx_data[,1][[1]]
    tab[i,'IGSN'] <- as.character(xlsx_data[which(rnames == 'IGSN')[1],2])
    tab[i,'dosimeter'] <- as.character(xlsx_data[which(rnames == 'Dosimeter Glass')[1],2])
    tab[i,'irradiation'] <- as.character(xlsx_data[which(rnames == 'Irradiation')[1],2])
    tab[i,'p-value'] <- as.character(xlsx_data[which(rnames == 'P (Chi-squred)')[1],2])
    print(tab[i,])
}
write.csv(tab,file='p-values_FT.csv')
