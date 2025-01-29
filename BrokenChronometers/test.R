library(IsoplotR)

pdf(file='test.pdf',onefile=TRUE)

Hoogland <- read.data("Hoogland.csv",method="U-Pb",format=2,ierr=4,
                      d=diseq(U48=list(x=1.00105,sx=0.001,option=2),
                              ThU=list(x=0,sx=0,option=1),
                              RaU=list(x=0,sx=0,option=1),
                              PaU=list(x=0,sx=0,option=1)))
Taung <- read.data("Taung.csv",method="U-Pb",format=5,ierr=3,
                   d=diseq(U48=list(x=1.0117,sx=0.00587,option=2),
                           ThU=list(x=0,sx=0,option=1),
                           RaU=list(x=0,sx=0,option=1),
                           PaU=list(x=0,sx=0,option=1)))
SB1625 <- read.data("SB-1625-22.csv",method="U-Pb",format=8,ierr=2)
SB72 <- read.data("SB-72-8.csv",method="U-Pb",format=8,ierr=2)

isochron(Hoogland,xlim=c(0,2500)) # 2a
isochron(Taung,joint=FALSE,type=1,taxis=TRUE) # 8c
isochron(Taung,joint=FALSE,type=2,taxis=TRUE) # 8d
isochron(SB1625,model=1,joint=FALSE,type=1,taxis=TRUE) # 9a
isochron(SB1625,joint=FALSE,type=2,taxis=TRUE) # 9b
isochron(SB72,joint=FALSE,type=1,taxis=TRUE) # 9c
isochron(SB72,joint=FALSE,type=2,taxis=TRUE) # 9d

ludwig(Hoogland,type=1,plot=TRUE) # Figures 5b and 5c
ludwig(Taung,type=1,plot=TRUE) # Figures 8a and 8b

dev.off()
