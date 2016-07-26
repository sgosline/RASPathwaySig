##goal is to update predictor code to hold out testdata
##what is metric of each individual predictor on test data
source("../../bin/elasticNetPred.R")
library(ggplot2)


plotRes<-function(kvals,pref=''){
    rems=which(apply(kvals,2,function(x){ all(is.na(unlist(x)))||all(is.null(unlist(x)))}))
  nk=kvals[,-rems]
  df<-list()
  for(i in 1:nrow(nk)){
    df$Mutated=c(df$Mutated,rep(rownames(nk)[i],ncol(nk)))
    df$Disease=c(df$Disease,colnames(nk))
    df$AUC=c(df$AUC,unlist(nk[i,]))
  }
  df<-as.data.frame(df)
  png(paste(pref,'AUCValues.png',sep=''))
  g<-ggplot(df)+geom_bar(aes(x=Disease,y=AUC,fill=Mutated),stat='identity',position='dodge')
  g<-g+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  print(g)
  dev.off()
  
}

k5=getDiseaseAUC(k=5)
plotRes(k5,'fiveFoldCV_enet')
#now need to plot...

k10=getDiseaseAUC(k=10)
plotRes(k10,'tenFoldCV_enet')
##now plot these values
##did'nt get to confusion matrix, do this tomorrow
