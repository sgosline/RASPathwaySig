source("../../bin/elasticNetPred.R")
library(pheatmap)
require(parallel)
source("../../bin/cBioPortalData.R")


#'Get AUC values by predicting from one mutation to another
#' @param genelist list of genes to compare across
#' @param cancerType TCGA cancer abbreviation
#' @param minPat number of patients to require in predictor
crossGenePreds<-function(genelist,mutMatrix,exprMatrix,cancerType='',prefix='',minPat=3){
                                        #iterate through the gene list
    cl<-makeCluster(min(10,length(genelist)),outfile='cluster.txt')
                                        #  clusterExport(cl,"getMutDataForGene")
    #print(dim(exprMatrix))
    #print(dim(mutMatrix))
    clusterExport(cl,"mutMatrix",envir=environment())
    clusterExport(cl,"exprMatrix",envir=environment())
    clusterExport(cl,"minPat",envir=environment())
    clusterExport(cl,"genelist",envir=environment())
    #exporting the whole source is too much, but might work after moving around some source commands
    clusterEvalQ(cl,source("../../bin/elasticNetPred.R"))
    clusterEvalQ(cl,library(pheatmap))

    dlist<-parLapply(cl,as.list(genelist),function(g){
    # dlist<-lapply(genelist,function(g){
        print(paste('Creating predictive model for',g,'across disease to run against',length(genelist),'other genes'))
        ##get mutation data, including patients with mutation
        gr<-which(rownames(mutMatrix)==g)
        genevals=rep(0,length(genelist))
        names(genevals)<-genelist
        if(length(gr)>0){
          gmuts<-which(mutMatrix[gr,])
        }else{
          gmuts<-c()
        }
        print(paste('Found',length(gmuts),'samples with mutated',g))
        mut.vec=rep('WT',ncol(exprMatrix))

                                        #create a factor vector to feed into predictive model
        mut.vec[gmuts]<-'MUTANT'
        mut.vec=factor(mut.vec,levels=c("WT","MUTANT"))

        if(length(which(mut.vec=='MUTANT'))<=minPat){
            print("Not enough mutants here, returning predictions of 0")
            return(genevals)
        }

                                        #build model, currently only elastic net
        fit=model.build(exprMatrix,mut.vec,pref=g,doPlot=FALSE)

        ##now use model to predict
        genevals<-lapply(genelist,function(g2,mutMatrix,fit,exprMatrix){
          or<-which(rownames(mutMatrix)==g2)
          genevals=rep(0,length(genelist))
          names(genevals)<-genelist

          if(length(or)>0){
            othermuts<-which(mutMatrix[or,])
          }else{
            othermuts<-c()
          }
          print(paste('Found',length(othermuts),'samples with mutated',g2))

          other.vec=rep('WT',ncol(exprMatrix))

          #create a factor vector to feed into predictive model
          other.vec[othermuts]<-'MUTANT'
          other.vec=factor(other.vec,levels=c("WT","MUTANT"))

          if(length(which(other.vec=='MUTANT'))<2)
                return(0.0)
          res=model.pred(fit,exprMatrix,other.vec,pref=paste(g,g2,sep='_to_'),doPlot=F)
            return(res$AUC)
        },mutMatrix,fit,exprMatrix)
    names(genevals)<-genelist
    return(genevals)
    })

    df=do.call('rbind',dlist)

   # print(df)
    stopCluster(cl)
    rownames(df)<-paste("From",genelist)
    colnames(df)<-paste("To",genelist)

    dmat<-apply(df,2,unlist)
    print(dmat)
    if(mean(dmat,na.rm=T)!=0){
        pheatmap(dmat,cellheight=10,cellwidth=10,main=paste(ifelse(cancerType=='','All',cancerType),prefix,'predictor AUC values'),filename=paste(cancerType,'min',minPat,'patientPredictorAUCvals.pdf',sep=''))
    }
    write.table(df,quote=F,file=paste(ifelse(cancerType=='','All',cancerType),prefix,'min',minPat,'patientPredictorAUCvals.txt',sep=''),sep='\t')

    return(df)
}

genelist=c("RASA1","SPRED1","NF1","TP53","NRAS","KRAS","BRAF","EGFR","SHC1","GRB2","MAP2K1","MAP2K","CDK4","RB1","PAK1","SOS1","PTEN","AKT1","PDK1","MTOR")

#genelist<-c("KRAS","SPRED1","NF1")


tcga.list<-c("",tcga.cancer.types)
ccle.list<-c("","BREAST","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","LUNG","SKIN","CENTRAL_NERVOUS_SYSTEM","LARGE_INTESTINE","OVARY")

getTcgaPredStats<-function(genelist){
    #TODO still need to get list of cancer types
    dlist<-rev(tcga.list)
    #make the cluster
    res<-do.call('rbind',lapply(dlist,function(ct,genelist){
      mutdata<-getDisMutationData(ct)
      exprdata<-getDisExpressionData(ct)
      compats<-intersect(colnames(mutdata),colnames(exprdata))

      print(paste("Found",length(compats),'patients with expression and mutation data for',ct))
      md<-mutdata[,compats]
      #print(dim(md))
      ed<-exprdata[,compats]
      df<-crossGenePreds(genelist,mutMatrix=md,exprMatrix=ed,cancerType=ct,prefix='tcga')
      print(paste('Finished',ct))
                                        #get offdiagonal predictions
        ndmat<-apply(df,2,unlist)*1-diag(nrow(df))
                                        #now collect mean values
        stats<-c(apply(ndmat,1,function(x) mean(x[x>0])),apply(ndmat,2,function(x) mean(x[x>0])))
        return(stats)
    },genelist))

    rownames(res)<-dlist
    write.table(res,file='TCGApathwayStats.txt',sep='\t')

}


getCclePredStats<-function(genelist){
  #TODO still need to get list of cancer types
  dlist<-ccle.list
                                        #make the cluster
  all.muts<-getCcleMutationData('')
  all.expr<-getCcleExpressionData('')
  res<-do.call('rbind',lapply(dlist,function(ct,genelist,all.muts,all.expr){
#    mutdata<-getCcleMutationData(ct)
                                        #    exprdata<-getCcleExpressionData(ct)
      mutdata<-all.muts
      exprdata<-all.expr

      compats<-intersect(colnames(mutdata),colnames(exprdata))
      if(ct!='')
          compats=compats[grep(ct)]
    df<-crossGenePreds(genelist,mutMatrix=mutdata[,compats],exprMatrix=exprdata[,compats],cancerType=ct,prefix='ccle')
    print(paste('Finished',ct))
    #get offdiagonal predictions
    ndmat<-apply(df,2,unlist)*1-diag(nrow(df))
    #now collect mean values
    stats<-c(apply(ndmat,1,function(x) mean(x[x>0])),apply(ndmat,2,function(x) mean(x[x>0])))
    return(stats)
  },genelist))

  rownames(res)<-dlist
  write.table(res,file='TCGApathwayStats.txt',sep='\t')

}

#res<-crossGenePreds(genelist,'PANCAN')
#res1<-getTcgaPredStats(genelist)
#res2<-getCclePredStats(genelist)
