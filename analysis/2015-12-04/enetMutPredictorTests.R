##begin to do elastic net predictions of kras, braf, nras in TCGA data.

##first get mutational data for each
source("../../bin/TcgaMutationalData.R")

if(!exists('kras.mut'))
  kras.mut<-getMutDataForGene("KRAS",TRUE)
if(!exists('nras.mut'))
  nras.mut<-getMutDataForGene("NRAS",TRUE)
if(!exists("braf.mut"))
  braf.mut<-getMutDataForGene("BRAF",TRUE)

##now submit model to do predictions
library(glmnet)
library(ggplot2)
source("../../bin/TcgaExpressionData.R")

buildModelFromData<-function(exprdata,mutdata,pref=''){
  ##first create vector of mutation data
  ##match patient identifiers and fuse into data frame
  mut.pats=toPatientId(as.character(mutdata$Tumor))
  expr.pats<-toPatientId(colnames(exprdata)[-1])
  
  mut.vec=rep('WT',length(expr.pats))
  mut.vec[match(mut.pats,expr.pats)]<-'MUTANT'
  if(length(which(mut.vec=='MUTANT'))<3)
    return(c())
   ##create model and identifier predictive genes
  cvfit<-cv.glmnet(x=t(exprdata[,-1]),y=as.factor(mut.vec), family='binomial',type='class',alpha=0.1)
  png(paste('cvEnetModelFor',pref,'.png',sep=''))
  plot(cvfit)
  dev.off()
  minlambda=cvfit$lambda.min
  
  coeffs=coef(cvfit, s = "lambda.min")
 
  pfit = predict(cvfit,t(exprdata[,-1]),s=minlambda,type="response")
  
  df<-data.frame(Prediction=pfit[,1],MutationStatus=mut.vec)
  ##plot predicted scores of MT vs WT
  png(paste('mutClassifierFor',pref,'.png',sep=''))
  p<-  ggplot(df)+geom_histogram(aes(x=Prediction,fill=MutationStatus))+scale_y_log10()
  print(p)
  dev.off()
  #also, which genes are nz in model? 
  genes<-exprdata[which(coeffs[,1]>0),1]
  gn<-sapply(as.character(genes),function(x) unlist(strsplit(x,split='|',fixed=T))[1])
  return(gn)
}

doPanCancer<-function(){
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

doPanDisease<-function(){
  #3 break down by disease, all mutations
  expr.pats<-toPatientId(colnames(alldat))

  ##TODO: get better list of tumors by disease.
  tumsByDis=tumsByDis[which(names(tumsByDis)!='pancan12')]
  dis.inds<-lapply(tumsByDis,function(x) which(expr.pats%in%x))

  res.genes<-sapply(names(dis.inds),function(x){
    inds<-dis.inds[[x]]
    if(length(inds)<10){
      return(list(nras=c(),kras=c(),braf=c()))
    }else{
      print(x)
      red.exp<-cbind(as.character(alldat[,1]),alldat[,inds])
      nras.all.genes<-buildModelFromData(red.exp,nras.mut,paste(x,'nras.all',sep='.'))
      kras.all.genes<-buildModelFromData(red.exp,kras.mut,paste(x,'kras.all',sep='.'))
      braf.all.genes<-buildModelFromData(red.exp,braf.mut,paste(x,'braf.all',sep='.'))
      return(list(nras=nras.all.genes,kras=kras.all.genes,braf=braf.all.genes))
    }
  })
  return(res.genes)
}

#4 break down by disease, missense/nonsense only

