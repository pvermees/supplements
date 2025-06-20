ThU <- RaU <- PaU <- list(x=0,sx=0,option=1,sd=100)

U48 <- list(x=1.0046,sx=0.0063/2,option=2,x0=10,sd=100)
AV03 <- read.data("data/AV03.csv",method="U-Pb",format=2,ierr=4,
                  d=diseq(U48=U48,ThU=ThU,RaU=RaU,PaU=PaU))

U48 <- list(x=1,sx=0.0063/2,option=2,x0=10,sd=100)
AV03eq <- read.data("data/AV03.csv",method="U-Pb",format=2,ierr=4,
                    d=diseq(U48=U48,ThU=ThU,RaU=RaU,PaU=PaU))
