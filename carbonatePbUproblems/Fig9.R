source('SB.R')

op <- par(mfrow=c(2,2),mgp=c(1.75,0.5,0),mar=c(3,3,3,0.5))
isochron(SB1625,model=3,joint=FALSE,type=1,taxis=TRUE) # 9a
isochron(SB1625,joint=FALSE,type=2,taxis=TRUE) # 9b
isochron(SB72,joint=FALSE,type=1,taxis=TRUE) # 9c
isochron(SB72,joint=FALSE,type=2,taxis=TRUE) # 9d
par(op)

dev.copy2pdf(file='output/Fig9.pdf')
