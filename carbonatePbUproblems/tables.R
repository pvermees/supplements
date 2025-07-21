d <- diseq(U48=list(x=1.006,option=2,M=100),
           PaU=list(x=0,option=1))
ttrue <- seq(from=0.5,to=3,by=0.5)
t75 <- t68 <- i48 <- NA*ttrue
for (i in seq_along(ttrue)){
    McL <- IsoplotR:::mclean(tt=ttrue[i],d=d)
    i48[i] <- d1$U48$x <- McL$U48i
    t68[i] <- age(x=McL$Pb206U238,method='U238-Pb206')[1]
    t75[i] <- age(x=McL$Pb207U235,method='U235-Pb207')[1]
}

print(
    rbind(
        t = ttrue,
        maxU48i = i48,
        relerr = 100*(t68-ttrue)/ttrue)
)

print(
    rbind(t = ttrue,
          relerr = 100*(ttrue-t75)/ttrue)
)
