source("../../bin/elasticNetPred.R")
library(pheatmap)
require(parallel)
source("../../bin/cBioPortalData.R")


tcga.list<<-c("allTcga",tcga.cancer.types)
ccle.list<<-c("allCcle","BREAST","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","LUNG","SKIN","CENTRAL_NERVOUS_SYSTEM","LARGE_INTESTINE","OVARY")



#'Get mutation-like score for each cell line/gene put into function
#' @param genelist list of genes to compare across
#' @param cancerType TCGA cancer abbreviation
#' @param minPat number of patients to require in predictor
crossDataScoresPerGene<-function(gene,datasetList,minPat=3){
                                        #iterate through the gene list
    cl<-makeCluster(min(2,length(datasetList)),outfile='cluster.txt')
  
   # require(gtools)
   # combos<-combinations(datasetList)

    
    #  clusterExport(cl,"getMutDataForGene")
    #print(dim(exprMatrix))
    #print(dim(mutMatrix))
    clusterExport(cl,c("gene","minPat",'datasetList','tcga.list','ccle.list'),envir=environment())
      #exporting the whole source is too much, but might work after moving around some source commands
    clusterEvalQ(cl,source("../../bin/elasticNetPred.R"))
    clusterEvalQ(cl,source("../../bin/cBioPortalData.R"))
    
    clusterEvalQ(cl,library(pheatmap))

    dlist<-parLapply(cl,as.list(datasetList),function(ds){
    # dlist<-lapply(genelist,function(g){
        print(paste('Creating predictive model for',ds,'across for gene',gene,' to run against',length(datasetList),'other datasetse'))
        ##get mutation data, including patients with mutation
        ##get training dataset - expression and mutation
        if(ds%in%tcga.list){
          exprMatrix<-getDisExpressionData(ds,getZscores=T)
          mutMatrix<-getDisMutationData(ds)
        }else{
          mutMatrix<-getCcleMutationData(ds)
          exprMatrix<-getCcleExpressionData(ds,getZscores=T)
        }
          
        gr<-which(rownames(mutMatrix)==g)
        genevals=rep(0,length(genelist))
        names(genevals)<-genelist
        if(length(gr)>0){
          gmuts<-which(mutMatrix[gr,])
        }else{
          gmuts<-c()
        }
        print(paste('Found',length(gmuts),'samples with mutated',gene))
        mut.vec=rep('WT',ncol(exprMatrix))

                                        #create a factor vector to feed into predictive model
        mut.vec[gmuts]<-'MUTANT'
        mut.vec=factor(mut.vec,levels=c("WT","MUTANT"))

        if(length(which(mut.vec=='MUTANT'))<=minPat){
            print("Not enough mutants here, returning predictions of 0")
            return(genevals)
        }

                                        #build model, currently only elastic net

        ##now use model to predict
    genevals<-lapply(datasetList,function(ds2){
      
      
      if(ds2%in%tcga.list){
        testExpr=getDisExpressionData(ds2,getZscores=T)
        testMatrix<-getDisMutationData(ds2)
              }else{
        testMatrix<-getCcleMutationData(ds2)
        testExpr<-getCcleExpressionData(ds2,getZscores=T)
              }
      #make sure we're looking at common patients
      comm.pats<-intersect(colnames(testExpr),colnames(testMatrix))
      testExpr<-testExpr[,comm.pats]
      testMatrix<-testMatrix[,comm.pats]
      
      ##and common genes between train and test
      comm.genes<-intersect(rownames(exprMatrix),rownames(testExpr))
      
      fit=model.build(exprMatrix[comm.genes,],mut.vec,pref=g,doPlot=FALSE)
      
      testExpr<-testExpr[comm.genes,]     
      
      gr2<-which(rownames(testMatrix)==gene)
      genevals=rep(0,length(datasetList))
      names(genevals)<-datasetList
      if(length(gr2)>0){
        gmuts2<-which(testMatrix[gr2,])
      }else{
        gmuts2<-c()
      }
      print(paste('Found',length(gmuts2),'samples with mutated',gene))
      mut.vec2=rep('WT',ncol(testExpr))
      
      #create a factor vector to feed into predictive model
      mut.vec2[gmuts2]<-'MUTANT'
      mut.vec2=factor(mut.vec2,levels=c("WT","MUTANT"))
      if(length(which(mut.vec2=='MUTANT'))<2)
        return(0.0)
      res=model.pred(fit,testExpr,mut.vec2,pref=paste(ds,'to',ds2,'forGene',gene,sep='_'),doPlot=T)
      return(res$AUC)
    })
    
    return(genevals)
    })

    df=do.call('rbind',dlist)
    zrows<-which(apply(df,1,function(x) all(x==0)))
    if(length(zrows)>0)
      df<-df[-zrows,]
    zcols<-which(apply(df,2,function(x) all(x==0)))
    if(length(zcols)>0)
      df<-df[,-zcols]
       # print(df)
    stopCluster(cl)
    colnames(df)<-paste("From",datasetList)
    rownames(df)<-paste("To",datasetList)

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



gene='NF1'
res<-crossDataScoresPerGene(datasetList=ccle.list,gene=gene)    

for(g in genelist){
  ##sample all combinations of datasets - ccle, tcga, to see how each predicts the other.
}
