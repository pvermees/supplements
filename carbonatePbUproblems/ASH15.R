ThU <- RaU <- PaU <- list(x=0,sx=0,option=1,sd=100)

ASH15diseq <- read.data("data/ASH15K.csv",method="U-Pb",format=5,ierr=4,
                        d=diseq(U48=list(x=0.99925,sx=0.0015/2,option=2,
                                         m=0.5,x0=1.081,M=1.5,sd=0.2),
                                ThU=ThU,RaU=RaU,PaU=PaU))
ASH15eq <- read.data("data/ASH15K.csv",method="U-Pb",format=5,ierr=4,
                     d=diseq(ThU=ThU,PaU=PaU))
