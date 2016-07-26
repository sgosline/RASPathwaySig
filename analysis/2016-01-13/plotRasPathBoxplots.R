##get predicted value of each set of patients, by mutation

source('../../bin/elasticNetPred.R')
library(pheatmap)
expr.pats<-toPatientId(colnames(alldat))

##TODO: get better list of tumors by disease.
tbd=tumsByDis[which(names(tumsByDis)!='PANCAN')]
dis.inds<-lapply(tbd,function(x) which(expr.pats%in%toPatientId(x)))


##should we expand to include other genes? 

makeDisBoxPlots<-function(mutdata,gene='KRAS',predfrom='COAD'){
  require(ggplot2)
  newdf=subset(mutdata,PredictedFrom==predfrom)
  ntp=names(which(table(subset(newdf,MutationStatus=='MUTANT')$PredictedTo)>5))
  newdf=subset(newdf,PredictedTo%in%ntp)
  png(paste('predicted',gene,'MutantScoresFrom',predfrom,'.png',sep=''))
  p<-ggplot(newdf)+geom_boxplot(aes(PredictedTo,MutationPrediction,color=MutationStatus))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  dev.off()
}

makeGeneBoxPlost<-function(mutdata,predGene,otherdata,secondGene){
  
}

##need to extract single predictions scores for each patient. 
collectPanCancerResponses<-function(mutdata=kras.mut,gene='KRAS',alpha=0.1){
  
  ##get all patients iwth mutations in the gene of interest
  mut.pats=toPatientId(as.character(mutdata$Tumor))
  #values to collect
  predictedFrom=c()#what disease are we predicted from
  predictedTo=c()#what disease are we predicting in
  mutantGene=c()#which gene
  mutationStatus=c() #status of gene in predicted sample
  mutProb=c()
  for(i in names(dis.inds)){
#  pred.mat<-sapply(names(dis.inds),function(i){
    mod.inds=dis.inds[[i]]
    #get expression values
    mod.exp<-cbind(as.character(alldat[,1]),alldat[,mod.inds])
    expr.pats<-toPatientId(colnames(mod.exp)[-1])
    
    ##get mutation values
    mod.mut<-rep('WT',length(mod.inds))
    mod.mut[match(mut.pats,expr.pats)]<-'MUTANT'
    #if(length(which(mod.mut=='MUTANT'))<1){
    
    
    mod.fit=NULL
    try(mod.fit<-model.build(mod.exp,mod.mut,alpha=alpha,doPlot=FALSE))
    if(is.null(mod.fit)){
      res=rep(NA,length(dis.inds))
      names(res)<-paste('TO',names(dis.inds))
      next
    }
#    pred.row<-sapply(names(dis.inds),function(j){
    
    for(j in names(dis.inds)){
      test.inds=dis.inds[[j]]
      test.exp=cbind(as.character(alldat[,1]),alldat[,test.inds])
      test.pats<-toPatientId(colnames(test.exp)[-1])
      
      test.mut=rep('WT',length(test.inds))
      test.mut[match(mut.pats,test.pats)]<-'MUTANT'
      pred=model.pred(mod.fit,test.exp,test.mut,alpha=alpha,doPlot=FALSE)
      print(paste('AUC of',j,'data from',i,'model:',pred$AUC))
      resp=pred$Response
      predictedFrom<-c(predictedFrom,rep(i,nrow(resp)))
      predictedTo=c(predictedTo,rep(j,nrow(resp)))
      mutantGene=c(mutantGene,rep(gene,nrow(resp)))
      mutationStatus=c(mutationStatus,as.character(paste(gene,resp$MutationStatus,sep='')))
      mutProb=c(mutProb,resp$Prediction)
          
    }
  }

  
  #now create epic data frame over everything
  df=data.frame(PredictedFrom=predictedFrom,PredictedTo=predictedTo,MutantGene=mutantGene,MutationStatus=mutationStatus,MutationPrediction=mutProb)
  return(df)
}

#kras_nras.df<-collectPanCancerResponses(rbind(kras.mut,nras.mut),'KRAS_NRAS')
#braf.df<-collectPanCancerResponses(braf.mut,'BRAF')

##nown do boxplos fo each disease
sapply(unique(kras_nras.df$PredictedFrom),function(x){
  makeDisBoxPlots(kras_nras.df,'KRAS_NRAS',x)
  
})
sapply(unique(braf.df$PredictedFrom),function(x){
  makeDisBoxPlots(braf.df,'BRAF',x)
})


  