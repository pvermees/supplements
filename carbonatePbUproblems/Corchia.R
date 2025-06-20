ThU <- RaU <- PaU <- list(x=0,sx=0,option=1,sd=100)
U48 <- list(x=0.9512,sx=0.0013/2,option=2,m=0.5,x0=0.75,M=1.0,sd=100)
Corchia <- read.data("data/Corchia.csv",method="U-Pb",format=2,ierr=4,
                     d=diseq(U48=U48,ThU=ThU,RaU=RaU,PaU=PaU))
