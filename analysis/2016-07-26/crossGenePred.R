source("../../bin/elasticNetPred.R")
library(pheatmap)

require(parallel)
##patient identifiers from the expression data
expr.pats<-toPatientId(colnames(alldat))

##TODO: get better list of tumors by disease.
dis.inds<-lapply(tumsByDis,function(x) which(expr.pats%in%toPatientId(x)))


#'Get AUC values by predicting from one mutation to another
#' @param genelist list of genes to compare across
#' @param cancerType TCGA cancer abbreviation
#' @param minPat number of patients to require in predictor
crossGenePreds<-function(genelist,cancerType='PANCAN',minPat=10){
  #iterate through the gene list
  df=do.call('rbind',lapply(genelist,function(g){
    ##get mutation data, including patients with mutation
    mutdata<-getMutDataForGene(g,FALSE,cancerType)
    mut.pats=toPatientId(as.character(mutdata$Tumor))
    genevals=rep(0,length(genelist))
    names(genevals)<-genelist
          #get expression data
    exprdata<-alldat
    expr.pats<-toPatientId(colnames(exprdata)[-1])

    mut.vec=rep('WT',length(expr.pats))
    mut.vec[match(mut.pats,expr.pats)]<-'MUTANT'
    mut.vec=factor(mut.vec,levels=c("WT","MUTANT"))

    if(length(which(mut.vec=='MUTANT'))<minPat){
      return(genevals)
    }
    #build model
    fit=model.build(exprdata,mut.vec,pref=g,doPlot=FALSE)
    genevals<-lapply(genelist,function(g2){
      other.muts<-getMutDataForGene(g2,FALSE,cancerType)
      other.mut.pats<-toPatientId(as.character(other.muts$Tumor))


      other.vec<-rep("WT",length(expr.pats))
      other.vec[match(other.mut.pats,expr.pats)]<-'MUTANT'
      other.vec<-factor(other.vec,levels=c("WT","MUTANT"))
      if(length(which(other.vec=='MUTANT'))<2)
          return(0.0)
      res=model.pred(fit,exprdata,other.vec,pref=paste(g,g2,sep='_to_'),doPlot=F)
      return(res$AUC)
    })
    return(genevals)
}))
  rownames(df)<-paste("From",genelist)
  colnames(df)<-paste("To",genelist)

  dmat<-apply(df,2,unlist)
  print(dmat)
  pheatmap(dmat,cellheight=10,cellwidth=10,main=paste(cancerType,'predictor AUC values'),filename=paste(cancerType,'min',minPat,'patientPredictorAUCvals.pdf',sep=''))
  write.table(df,quote=F,file=paste(cancerType,'min',minPat,'patientPredictorAUCvals.txt',sep=''),sep='\t')



  return(df)


}



genelist=c("RASA1","SPRED1","NF1","TP53","NRAS","KRAS","BRAF","EGFR","SHC1","GRB2","MAP2K1","MAP2K","CDK4","RB1","PAK1","SOS1","PTEN","AKT1","PDK1","MTOR")

genelist<-c("KRAS","SPRED1","NF1")

getPredStats<-function(genelist){
    #make the cluster
    cl<-makeCluster(30,outfile='cluster.txt')
    clusterExport(cl,"crossGenePreds")
    clusterEvalQ(cl,source("../../bin/elasticNetPred.R"))
    clusterEvalQ(cl,library(pheatmap))

    res<-do.call('rbind',parLapply(cl,list(names(tumsByDis)),function(ct){
        df<-crossGenePreds(genelist,cancerType=ct)
        print(paste('Finished',ct))
                                        #get offdiagonal predictions
        ndmat<-apply(df,2,unlist)*1-diag(nrow(df))
                                        #now collect mean values
        stats<-c(apply(ndmat,1,function(x) mean(x[x>0])),apply(ndmat,2,function(x) mean(x[x>0])))
        return(stats)
    }))

    rownames(res)<-names(tumsByDis)
    write.table(res,file='pathwayStats.txt',sep='\t')
    stopCluster(cl)
}

res<-getPredStats(genelist)
