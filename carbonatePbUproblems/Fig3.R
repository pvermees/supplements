source("AV03.R")

AV03eq$d$U48$option <- 0
MC(AV03eq)

dev.copy2pdf(file='output/Fig3.pdf',width=6,height=7)
