##begin to do elastic net predictions of kras, braf, nras in TCGA data.

##now submit model to do predictions
library(glmnet)
library(ggplot2)
library(pROC)
library(caret)


#'Build model from data
model.build<-function(exprdata,mut.vec,pref='',alpha=0.1,doPlot=TRUE){

  if(length(which(mut.vec=='MUTANT'))<1){
    print(paste("Less than 1 mutants in patient data for",pref))
    return(NULL)
  }
  exprdata[which(is.na(exprdata),arr.ind=T)]<-0.0
  
  ##create model and identifier predictive genes
  #cvfit<-cv.glmnet(x=t(exprdata[,-1]),y=as.factor(mut.vec),
  cvfit<-cv.glmnet(x=t(exprdata),y=as.factor(mut.vec),
                   family='binomial',type='auc',alpha=alpha)

  if(doPlot){
    png(paste('cvEnetModelFor',pref,'.png',sep=''))
    plot(cvfit)
    dev.off()
  }
  return(cvfit)

}

#'predict model from fit
model.pred<-function(cvfit,exprdata,mut.vec,pref='',alpha=0.1,doPlot=TRUE){

  #extract information and predict
  minlambda=cvfit$lambda.min

  coeffs=coef(cvfit, s = "lambda.min")

  auc.val=0
  pfit<-NULL
  exprdata[which(is.na(exprdata),arr.ind=T)]<-0.0
  #new.exprdata<-apply(exprdata,2,function(x) {x[which(!is.finite(as.numeric(x)))]<-0.0; return(as.numeric(x))})
  #   try(pfit <- predict(cvfit,t(exprdata[,-1]),s=minlambda,type="response"))
  try(pfit <- predict(cvfit,t(exprdata),s=minlambda,type="response"))
  if(is.null(pfit))
	  return(auc.val)

  df<-data.frame(Prediction=pfit[,1],MutationStatus=mut.vec)
  ##plot predicted scores of MT vs WT
  if(doPlot){
    png(paste('mutClassifierFor',pref,'.png',sep=''))
    tstring<-paste('Predictions for\n',gsub('_',' ',pref,fixed=T))
    p<-  ggplot(df)+geom_boxplot(aes(y=Prediction,x=MutationStatus),position='dodge')+ggtitle(tstring)
    print(p)
    dev.off()
  }

  try(auc.val<-auc(roc(MutationStatus~Prediction,df)))

  return(list(Response=df,Coeff=coeffs,AUC=auc.val))
}

#'takes a model and scores each column with the likelihood of being a mutant!
model.score<-function(cvfit,exprdata,pref='',alpha=0.1){
  #extract information and predict
  minlambda=cvfit$lambda.min
  
  coeffs=coef(cvfit, s = "lambda.min")
  
  pred.val=rep(NA,ncol(exprdata))
  names(pred.val)<-colnames(exprdata)
  pfit<-NULL
  exprdata[which(is.na(exprdata),arr.ind=T)]<-0.0
  
  #   try(pfit <- predict(cvfit,t(exprdata[,-1]),s=minlambda,type="response"))
  try(pfit <- predict(cvfit,t(exprdata),s=minlambda,type="response"))
  if(is.null(pfit))
    return(pred.val)
  
  pred.val<-pfit[,1]
  names(pred.val)<-rownames(pfit)
  return(pred.val)
  
#  rerun(pre[,1],MutationStatus=mut.vec)
  
}

buildModelFromData<-function(exprdata,mutdata,pref='',alpha=0.1,doPlot=TRUE,doKfold=10){
  #calculate ground truth in mutation data
  mut.pats=toPatientId(as.character(mutdata$Tumor))
  expr.pats<-toPatientId(colnames(exprdata)[-1])

  mut.vec=rep('WT',length(expr.pats))
  mut.vec[match(mut.pats,expr.pats)]<-'MUTANT'

  #build model
  fit=model.build(exprdata,mut.vec,pref,alpha,doPlot)
  res=model.pred(fit,exprdata,mut.vec,pref,alpha,doPlot)
  df=res$Response
  coeffs=res$Coeff
  #predict
    #also, which genes are nz in model?
  genes<-exprdata[which(coeffs[,1]>0),1]
  gn<-sapply(as.character(genes),function(x) unlist(strsplit(x,split='|',fixed=T))[1])
  write(gn,file=paste('genesInPredictorFor',pref,'.txt',sep=''))
  print(paste('Returning',length(gn),'genes for alpha',alpha))
  return(gn)
}

crossVal<-function(exprdata,mutdata,prefix='',alpha=0.1,k=10){
  #calculate ground truth in mutation data
  mut.pats=toPatientId(as.character(mutdata$Tumor))
  expr.pats<-toPatientId(colnames(exprdata)[-1])

  mut.vec=rep('WT',length(expr.pats))
  mut.vec[match(mut.pats,expr.pats)]<-'MUTANT'
  if(length(which(mut.vec=='MUTANT'))<10){
    print(paste("Less than 10 mutants in patient data for",prefix,'no crossval'))
    return(NA)
  }
  testIdx=createFolds(y=mut.vec,k=k)
  rocs<-sapply(testIdx,function(x){
    y=c(1,1+x)
    test.dat<-exprdata[,y]
    test.mut<-mut.vec[x]

    train.dat<-cbind(exprdata[,1],exprdata[,-y])
    train.mut<-mut.vec[-x]

    mod=model.build(train.dat,train.mut,alpha=alpha,doPlot=FALSE)
    if(is.null(mod))
      return(NA)

    pred=model.pred(mod,test.dat,test.mut,alpha=alpha,doPlot=FALSE)
    auc.val=pred$AUC
    return(auc.val)
  })
  return(mean(rocs,na.rm=T))
}

doPanCancer<-function(){
    ##first get mutational data for each
    source("../../bin/TcgaExpressionData.R")
    source("../../bin/TcgaMutationalData.R")


  ##Now run the test:
  #1 for all three set  s, all mutations
  nras.all.genes<-buildModelFromData(alldat,nras.mut,'nras.all')
  kras.all.genes<-buildModelFromData(alldat,kras.mut,'kras.all')
  braf.all.genes<-buildModelFromData(alldat,braf.mut,'braf.all')

  #2 for all three sets, missense or nonsense only.
  nras.ms.ns<-subset(nras.mut,VariantClassification%in%c('Missense_Mutation','Nonsense_Mutation'))
  nras.genes<-buildModelFromData(alldat,nras.ms.ns,'nras.missense_nonsense')
  kras.ms.ns<-subset(kras.mut,VariantClassification%in%c('Missense_Mutation','Nonsense_Mutation'))
  kras.genes<-buildModelFromData(alldat,kras.ms.ns,'kras.missense_nonsense')
  braf.ms.ns<-subset(braf.mut,VariantClassification%in%c('Missense_Mutation','Nonsense_Mutation'))
  braf.genes<-buildModelFromData(alldat,braf.ms.ns,'braf.missense_nonsense')
}

doPanDisease<-function(alpha=0.01,doPlot=TRUE){
    ##first get mutational data for each
    source("../../bin/TcgaExpressionData.R")
    source("../../bin/TcgaMutationalData.R")


  #3 break down by disease, all mutations
  expr.pats<-toPatientId(colnames(alldat))

  ##TODO: get better list of tumors by disease.
  tbd=tumsByDis[which(names(tumsByDis)!='PANCAN')]
  dis.inds<-lapply(tbd,function(x) which(expr.pats%in%toPatientId(x)))

  res.genes<-sapply(names(dis.inds),function(x){
    inds<-dis.inds[[x]]
    if(length(inds)<50){
      return(list(nras=c(),kras=c(),braf=c()))
    }else{
      red.exp<-cbind(as.character(alldat[,1]),alldat[,inds])
      print(paste(x,'NRAS'))
      nras.all.genes<-buildModelFromData(red.exp,nras.mut,paste(x,'nras.alpha',alpha,sep='.'),alpha,doPlot = doPlot)
      print(paste(x,'KRAS'))
      kras.all.genes<-buildModelFromData(red.exp,kras.mut,paste(x,'kras.alpha',alpha,sep='.'),alpha,doPlot = doPlot)
      print(paste(x,'BRAF'))
      braf.all.genes<-buildModelFromData(red.exp,braf.mut,paste(x,'braf.alpha',alpha,sep='.'),alpha,doPlot = doPlot)
      return(list(nras=nras.all.genes,kras=kras.all.genes,braf=braf.all.genes))
    }
  })
  return(res.genes)
}

getDiseaseAUC<-function(alpha=0.1,k=10){
      source("../../bin/TcgaExpressionData.R")
    source("../../bin/TcgaMutationalData.R")

  #3 break down by disease, all mutations
  expr.pats<-toPatientId(colnames(alldat))

  ##TODO: get better list of tumors by disease.
  tbd=tumsByDis[which(names(tumsByDis)!='PANCAN')]
  dis.inds<-lapply(tbd,function(x) which(expr.pats%in%toPatientId(x)))

  res.auc<-sapply(names(dis.inds),function(x){
    inds<-dis.inds[[x]]
    if(length(inds)<50){
      return(list(nras=c(),kras=c(),braf=c()))
    }else{
      red.exp<-cbind(as.character(alldat[,1]),alldat[,inds])
      print(paste(x,'NRAS'))
      nras<-crossVal(red.exp,nras.mut,paste(x,'nras.alpha',alpha,sep='.'),alpha,k = k)
      print(paste(x,'KRAS'))
      kras<-crossVal(red.exp,kras.mut,paste(x,'kras.alpha',alpha,sep='.'),alpha,k = k)
      print(paste(x,'BRAF'))
      braf<-crossVal(red.exp,braf.mut,paste(x,'braf.alpha',alpha,sep='.'),alpha,k = k)
      return(list(nras=nras,kras=kras,braf=braf))}})

  return(res.auc)
}

##copy confusion matrix code from 2015/12/22
buildConfMatrix<-function(mutdata=kras.mut,gene='KRAS',alpha=0.1){
    source("../../bin/TcgaExpressionData.R")
    source("../../bin/TcgaMutationalData.R")
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

##this is similar to the confusion matrix for each disease, but instead focuses on the
##Average within-disease AUC for various genes
buildCrossGeneMatrix<-function(mutdata,exprdata,genelist=c("KRAS","NRAS","NF1","TP53"),dis.name='GBM',alpha=0.1){
  ##mutdata- list of per-patint

    source("../../bin/TcgaMutationalData.R")
  gene.res<-sapply(genelist,function(gene){
    ##get mutdata
        ##get all patients iwth mutations in the gene of interest
    mut.pats=toPatientId(as.character(mutdata$Tumor))

  mod.inds=dis.inds[[dis.name]]
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
    res=rep(NA,length(genelist))
    names(res)<-paste('TO',genelist)
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
