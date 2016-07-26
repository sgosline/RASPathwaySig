#source("../../bin/TcgaMutationalData.R")

#res<-getMutDataForGene('KRAS',plot=TRUE)
#res<-getMutDataForGene('NRAS',plot=TRUE)
#res<-getMutDataForGene('BRAF',plot=TRUE)

source("../../bin/elasticNetPred.R")
all.a<-sapply(c(0.1,0.01,0.001),doPanDisease)


