source("../../bin/elasticNetPred.R")
library(pheatmap)
require(parallel)
source("../../bin/cBioPortalData.R")

#'Cross data preds - evaluate the ability of a predictor to evaluate mutations from one
#'dataset to another
crossDataPreds<-function(genelist,mutMatrix,exprMatrix,testExpr,testMut,p1,p2,minPat=3){
  print(paste('Evaluating gene list preds from',p1,'to',p2))
  cl<-makeCluster(min(2,length(genelist)),outfile='cluster.txt')
  #  clusterExport(cl,"getMutDataForGene")
  #print(dim(exprMatrix))
  #print(dim(mutMatrix))
  
  com.genes<-intersect(rownames(exprMatrix),rownames(testExpr))
  exprMatrix=exprMatrix[com.genes,]
  testExpr=testExpr[com.genes,]
  
  
  clusterExport(cl,c("mutMatrix","exprMatrix","minPat","genelist","testExpr","testMut"),envir=environment())
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
      return(mean(genevals))
    }
    
    #build model, currently only elastic net
    fit=model.build(exprMatrix,mut.vec,pref=g,doPlot=FALSE)
    
    ##now use model to predict in second datasets
    or<-which(rownames(testMut)==g)
    sampvals=rep(0,ncol(testExpr))
    #names(genevals)<-genelist
      
    if(length(or)>0){
      othermuts<-which(testMut[or,])
    }else{
      othermuts<-c()
    }
    print(paste('Found',length(othermuts),'test samples with mutated',g))
      
    other.vec=rep('WT',ncol(testExpr))
      
    #create a factor vector to feed into predictive model
    other.vec[othermuts]<-'MUTANT'
    other.vec=factor(other.vec,levels=c("WT","MUTANT"))
    
    if(length(which(other.vec=='MUTANT'))<2)
      return(0.0)
    
    res=model.pred(fit,as.matrix(testExpr),other.vec,pref=paste(g,'prediction from',p1,'to',p2),doPlot=F)
  #  print(res)
    return(res$AUC)
    })
  names(dlist)<-genelist
  dlist
  
}


#'Get mutation-like score for each cell line/gene put into function
#' @param genelist list of genes to compare across
#' @param cancerType TCGA cancer abbreviation
#' @param minPat number of patients to require in predictor
crossDataScores<-function(genelist,mutMatrix,exprMatrix,testMatrix,cancerType='',prefix='',minPat=3){
                                        #iterate through the gene list
    cl<-makeCluster(min(2,length(genelist)),outfile='cluster.txt')
  
    
    com.genes<-intersect(rownames(exprMatrix),rownames(testMatrix))
    exprMatrix=exprMatrix[com.genes,]
    testMatrix=testMatrix[com.genes,]
    
    #  clusterExport(cl,"getMutDataForGene")
    #print(dim(exprMatrix))
    #print(dim(mutMatrix))
    clusterExport(cl,c("mutMatrix","exprMatrix","minPat","genelist","testMatrix"),envir=environment())
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
    genevals<-model.score(fit,testMatrix)
 #   print(genevals)
#    names(genevals)<-genelist
    return(genevals)
    })

    df=do.call('rbind',dlist)

   # print(df)
    stopCluster(cl)
    colnames(df)<-paste("Sample",colnames(testMatrix))
    rownames(df)<-paste(genelist,'Mutated')

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

if(load){
#test data first
tcga.brca.expr<-getDisExpressionData('brca',getZscores=T)
tcga.brca.mut<-getDisMutationData('brca')
ccle.breast.expr<-getCcleExpressionData("BREAST",getZscores=T)
ccle.breast.mut<-getCcleMutationData("BREAST")

##find the common samples in each.

load<-FALSE


tcga.coad.expr<-getDisExpressionData('coadread',getZscores=T)
tcga.coad.mut<-getDisMutationData("coadread")
ccle.colon.mut<-getCcleMutationData('LARGE_INTESTINE')
ccle.colon.expr<-getCcleExpressionData("LARGE_INTESTINE",getZscores=T)

comm.tcga<-intersect(colnames(tcga.coad.expr),colnames(tcga.coad.mut))
comm.ccle<-intersect(colnames(ccle.colon.expr),colnames(ccle.colon.mut))
}
##first test out the ability to evaluate performance from one disease to another, using AUC
#res<-crossDataPreds(genelist,mutMatrix=tcga.coad.mut[,comm.tcga],exprMatrix=as.matrix(tcga.coad.expr[,comm.tcga]),testExpr=as.matrix(ccle.colon.expr[,comm.ccle]),testMut=ccle.colon.mut[,comm.ccle],'TCGA','CCLE')

##second test out the ability to just score one disease
res<-crossDataScores(genelist,mutMatrix=tcga.coad.mut[,comm.tcga],exprMatrix=as.matrix(tcga.coad.expr[,comm.tcga]),testMatrix=as.matrix(ccle.colon.expr[,comm.ccle]),'TCGA','CCLE')

