
source("../../bin/elasticNetPred.R")
require(parallel)
##patient identifiers from the expression data
expr.pats<-toPatientId(colnames(alldat))

##TODO: get better list of tumors by disease.
dis.inds<-lapply(tumsByDis,function(x) which(expr.pats%in%toPatientId(x)))


##get a list of genes, compare cross-gene predictivity'
crossGenePreds<-function(genelist,cancerType='PANCAN',minPat=10){
  #iterate through the gene list
  df=do.call('rbind',mclapply(genelist,function(g){
    ##get mutation data, including patients with mutation
    mutdata<-getMutDataForGene(g,FALSE,cancerType)
    mut.pats=toPatientId(as.character(mutdata$Tumor))
    genevals=rep(0,length(genelist))
    names(genevals)<-genelist
    if(length(mut.pats)<minPat){
      return(genevals)
    }
          #get expression data
    exprdata<-alldat
    expr.pats<-toPatientId(colnames(exprdata)[-1])
    
    mut.vec=rep('WT',length(expr.pats))
    mut.vec[match(mut.pats,expr.pats)]<-'MUTANT'
    mut.vec=factor(mut.vec,levels=c("WT","MUTANT"))
    #build model
    fit=model.build(exprdata,mut.vec,pref=g,doPlot=FALSE)
    genevals<-mclapply(genelist,function(g2){
      other.muts<-getMutDataForGene(g2,FALSE,cancerType)
      other.mut.pats<-toPatientId(as.character(other.muts$Tumor))
      other.vec<-rep("WT",length(expr.pats))
      other.vec[match(other.mut.pats,expr.pats)]<-'MUTANT'
      other.vec<-factor(other.vec,levels=c("WT","MUTANT"))
      res=model.pred(fit,exprdata,other.vec,pref=paste(g,g2,sep='_to_'),doPlot=F)
      return(res$AUC)
    },mc.cores=10)
    return(genevals)
  },mc.cores=10))
  return(df)
}

genelist=c("RASA1","SPRED1","NF1","TP53","NRAS","KRAS","BRAF","EGFR")
rdf<-crossGenePreds(genelist)
