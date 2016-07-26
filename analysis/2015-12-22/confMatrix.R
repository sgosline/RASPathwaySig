##now do confusion matrix....


source('../../bin/elasticNetPred.R')
library(pheatmap)
expr.pats<-toPatientId(colnames(alldat))

##TODO: get better list of tumors by disease.
tbd=tumsByDis[which(names(tumsByDis)!='PANCAN')]
dis.inds<-lapply(tbd,function(x) which(expr.pats%in%toPatientId(x)))


buildConfMatrix<-function(mutdata=kras.mut,gene='KRAS',alpha=0.1){
    

  ##get all patients iwth mutations in the gene of interest
  mut.pats=toPatientId(as.character(mutdata$Tumor))
  
  #for(i in names(dis.inds)){
  pred.mat<-sapply(names(dis.inds),function(i){
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
      return(res)
    }
    #now with this model, predict across all other disease
    pred.row<-sapply(names(dis.inds),function(j){
    #for(j in names(dis.inds)){
      test.inds=dis.inds[[j]]
      test.exp=cbind(as.character(alldat[,1]),alldat[,test.inds])
      test.pats<-toPatientId(colnames(test.exp)[-1])
      
      test.mut=rep('WT',length(test.inds))
      test.mut[match(mut.pats,test.pats)]<-'MUTANT'
      pred=model.pred(mod.fit,test.exp,test.mut,alpha=alpha,doPlot=FALSE)
      print(paste('AUC of',j,'data from',i,'model:',pred$AUC))
      
      return(pred$AUC)
    })
    names(pred.row)<-paste('TO',names(pred.row))
    return(pred.row)
  })
  colnames(pred.mat)<-paste('FROM',colnames(pred.mat))
  
  na.cols=which(apply(pred.mat,2,function(x) all(is.na(x))))
  if(length(na.cols)>0)
    nmat<-pred.mat[,-na.cols]
  zrows=which(apply(nmat,1,function(x) all(x==0)))
  if(length(zrows)>0)
    nmat<-nmat[-zrows,]
  pheatmap(nmat,cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,main=paste('AUC for',gene,'mutants'),file=paste(gene,'AUCVals.png',sep=''))
  return(nmat)
  }

##now call the function to make the plots
kmat=buildConfMatrix(kras.mut,'KRAS')
bmat=buildConfMatrix(braf.mut,'BRAF')
nmat=buildConfMatrix(nras.mut,'NRAS')
nmat.kmat<-buildConfMatrix(rbind(kras.mut,nras.mut),'KRAS_and_NRAS')

