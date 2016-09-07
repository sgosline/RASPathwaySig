source("../../bin/elasticNetPred.R")
library(pheatmap)
require(parallel)
source("../../bin/cBioPortalData.R")
library(synapseClient)
synapseLogin()
tcga.list<<-c("allTcga",tcga.cancer.types)
ccle.list<<-c("allCcle","BREAST","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","LUNG","SKIN","CENTRAL_NERVOUS_SYSTEM","LARGE_INTESTINE","OVARY")


#'Get mutation-like score for each cell line/gene put into function
#' @param genelist list of genes to compare across
#' @param cancerType TCGA cancer abbreviation
#' @param minPat number of patients to require in predictor
scoreNFforGene<-function(gene,datasetList,testExpr,mut.vec2,dataset,minPat=3,fullMut=NA,fullExpr=NA){
                                        #iterate through the gene list

    if(is.na(fullMut))
        fullMut<-getDisMutationData('') #getCcleMutationData('')
    if(is.na(fullExpr))
        fullExpr<-getDisExpressionData('',getZscores=TRUE)#getCcleExpressionData('',getZscores=T)

  dlist<-lapply(datasetList,function(ds){
    # dlist<-lapply(genelist,function(g){
        print(paste('Creating predictive model for',ds,'across for gene',gene,' to run against',dataset))
        ##get mutation data, including patients with mutation
        ##get training dataset - expression and mutation
  #  mutMatrix=fullMut
   # exprMatrix=fullExpr
     if(ds=='allTcga'){
       mutMatrix=fullMut
       exprMatrix=fullExpr
   }else{
       samps<-sapply(getSamplesForDisease(ds),function(x) gsub('-','.',x))
       mutMatrix<-fullMut[,intersect(samps,colnames(fullMut))]
       exprMatrix<-fullExpr[,intersect(samps,colnames(fullExpr))]

     }

      gr<-which(rownames(mutMatrix)==gene)
      genevals=rep(0,length(datasetList))
      names(genevals)<-datasetList
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
            return(0)
        }

                                        #build model, currently only elastic net

     # testMatrix<-testMut
      #make sure we're looking at common patients
      comm.pats<-intersect(colnames(testExpr),names(mut.vec2))
      testExpr<-testExpr[,comm.pats]
      mut.vec2<-mut.vec2[comm.pats]

      ##and common genes between train and test
      comm.genes<-intersect(rownames(exprMatrix),rownames(testExpr))

      fit=model.build(exprMatrix[comm.genes,],mut.vec,pref=gene,doPlot=FALSE)

      #testExpr<-testExpr[comm.genes,]


      res=model.pred(fit,testExpr[comm.genes,],mut.vec2,pref=paste(ds,'to',dataset,'forGene',gene,sep='_'),doPlot=T)
      return(res$AUC)
    })
    dlist<-unlist(dlist)
    names(dlist)<-datasetList
    df<-data.frame(Sample=datasetList,AUC=as.numeric(unlist(dlist)))
    print(df)
    require(ggplot2)
    pdf(paste(gene,'mutationsPredictedIn',dataset,'.pdf',sep=''))
    g<-ggplot(df)+geom_point(aes(x=Sample,y=AUC))+ggtitle(paste('Predicting NF1 status in dataset from',gene,'in cell lines'))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(g)
    dev.off()
    #    if(mean(dmat,na.rm=T)!=0){
    #    pheatmap(dmat,cellheight=10,cellwidth=10,main=paste(ifelse(cancerType=='','All',cancerType),prefix,'predictor AUC values'),filename=paste(cancerType,'min',minPat,'patientPredictorAUCvals.pdf',sep=''))
    #}
    #write.table(df,quote=F,file=paste(ifelse(cancerType=='','All',cancerType),prefix,'min',minPat,'patientPredictorAUCvals.txt',sep=''),sep='\t')

    #return(df)
    return(dlist)
}

genelist=c("RASA1","SPRED1","NF1","TP53","NRAS","KRAS","BRAF","EGFR","SHC1","GRB2","MAP2K1","MAP2K","CDK4","RB1","PAK1","SOS1","PTEN","AKT1","PDK1","MTOR")

#genelist=c("RASA1","NF1","TP53","NRAS","KRAS","EGFR","SHC1","CDK4","RB1","PAK1","SOS1","PTEN","AKT1","PDK1","MTOR")

#genelist<-c("KRAS","SPRED1","NF1")



gene='NF1'
#res<-crossDataScoresPerGene(datasetList=ccle.list,gene=gene)
pnfData<-read.table(synGet('syn7124098')@filePath,sep='\t',header=T)
colnames(pnfData)<- gsub('..',' (',gsub('.clone.',' clone)',colnames(pnfData),fixed=T),fixed=T)
phenoData<-read.table(synGet('syn7139168')@filePath,sep='\t',header=T)

mutVec=phenoData$Genotype
names(mutVec)<-rownames(phenoData)

mutVecHetsAsNeg<-rep('MUTANT',length(mutVec))
mutVecHetsAsNeg[which(mutVec=='++')]<-'WT'
mutVecHetsAsNeg<-factor(mutVecHetsAsNeg,levels=c("WT","MUTANT"))
names(mutVecHetsAsNeg)<-rownames(phenoData)

mutVecHetsAsPos<-rep('MUTANT',length(mutVec))
mutVecHetsAsPos[which(mutVec!='--')]<-'WT'
mutVecHetsAsPos<-factor(mutVecHetsAsPos,levels=c("WT","MUTANT"))
names(mutVecHetsAsPos)<-rownames(phenoData)

cl<-makeCluster(min(5,length(genelist)),outfile='pnf_cluster.txt')
##get all data
load('exprData.Rdata')
fullExpr<-exprData
load('mutData.Rdata')
fullMut<-mutData
#  fullMut<-getDisMutationData('') #getCcleMutationData('')
#  fullExpr<-getDisExpressionData('',getZscores=TRUE)#getCcleExpressionData('',getZscores=T)

clusterExport(cl,c("scoreNFforGene","mutVecHetsAsNeg","mutVecHetsAsPos","pnfData",'tcga.list','fullMut','fullExpr'))
clusterEvalQ(cl,source("../../bin/elasticNetPred.R"))
clusterEvalQ(cl,source("../../bin/cBioPortalData.R"))


dlist<-parLapply(cl,as.list(genelist),function(g){
datasetList<-tcga.list
#for(g in genelist){
  ##sample all combinations of datasets - ccle, tcga, to see how each predicts the other.
  res<-scoreNFforGene(g,datasetList,pnfData,mutVecHetsAsNeg,'pnfCellsHetsareMuts',minPat=3,fullMut,fullExpr)
  res2<-scoreNFforGene(g,datasetList,pnfData,mutVecHetsAsPos,'pnfCellsHetsAreWT',minPat=3,fullMut,fullExpr)

  })
stopCluster(cl)
