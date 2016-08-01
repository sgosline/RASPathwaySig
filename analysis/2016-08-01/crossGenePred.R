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
crossGenePreds<-function(genelist,cancerType='PANCAN',minPat=5){
                                        #iterate through the gene list
    cl<-makeCluster(min(30,length(genelist)),outfile='cluster.txt')
                                        #  clusterExport(cl,"getMutDataForGene")
    mutdata<<-getMutStatusByDisease(cancerType)
   # print(dim(mutdata))
    clusterExport(cl,list("toPatientId","mutdata","alldat"))
#    clusterExport(cl,"alldat")
    #exporting the whole source is too much, but might work after moving around some source commands
    clusterEvalQ(cl,source("../../bin/elasticNetPred.R"))
#    clusterExport(cl,"model.build")
#    clusterExport(cl,"model.pred")
    clusterEvalQ(cl,library(pheatmap))



    dlist<-parLapply(cl,as.list(genelist),function(g,mutdata,genelist){
        print(paste('Creating predictive model for',g,'across disease to run against',length(genelist),'other genes'))
        ##get mutation data, including patients with mutation
        gmuts<-subset(mutdata,Gene==g)
        mut.pats=toPatientId(as.character(gmuts$Tumor))
        genevals=rep(0,length(genelist))
        names(genevals)<-genelist
                                        #get expression data
                                        #calculate expression data
        exprdata<-alldat
        expr.pats<-toPatientId(colnames(exprdata)[-1])

                                        #create a factor vector to feed into predictive model
        mut.vec=rep('WT',length(expr.pats))
        mut.vec[match(mut.pats,expr.pats)]<-'MUTANT'
        mut.vec=factor(mut.vec,levels=c("WT","MUTANT"))

        if(length(which(mut.vec=='MUTANT'))>=minPat){
#            print("Not enough mutants here, returning predictions of 0")
#            return(genevals)

                                        #build model, currently only elastic net
        fit=model.build(exprdata,mut.vec,pref=g,doPlot=FALSE)

        ##now use model to predict
        genevals<-lapply(genelist,function(g2,mutdata,fit,exprdata){
            other.muts<-subset(mutdata,Gene==g2)
            print(paste('Found',nrow('other.muts'),'for',g2,'to assess predictions'))
            other.mut.pats<-toPatientId(as.character(other.muts$Tumor))

            other.vec<-rep("WT",length(expr.pats))
            other.vec[match(other.mut.pats,expr.pats)]<-'MUTANT'
            other.vec<-factor(other.vec,levels=c("WT","MUTANT"))
            if(length(which(other.vec=='MUTANT'))<2)
                return(0.0)
            res=model.pred(fit,exprdata,other.vec,pref=paste(g,g2,sep='_to_'),doPlot=F)
            return(res$AUC)
        },mutdata,fit,exprdata)
    }
    return(genevals)
    },mutdata,genelist)

    df=do.call('rbind',dlist)

   # print(df)
    stopCluster(cl)
    rownames(df)<-paste("From",genelist)
    colnames(df)<-paste("To",genelist)

    dmat<-apply(df,2,unlist)
                                        #    print(dmat)
    if(mean(dmat,na.rm=T)!=0){
        pheatmap(dmat,cellheight=10,cellwidth=10,main=paste(cancerType,'predictor AUC values'),filename=paste(cancerType,'min',minPat,'patientPredictorAUCvals.pdf',sep=''))
    }
    write.table(df,quote=F,file=paste(cancerType,'min',minPat,'patientPredictorAUCvals.txt',sep=''),sep='\t')

    return(df)
}

genelist=c("RASA1","SPRED1","NF1","TP53","NRAS","KRAS","BRAF","EGFR","SHC1","GRB2","MAP2K1","MAP2K","CDK4","RB1","PAK1","SOS1","PTEN","AKT1","PDK1","MTOR")

#genelist<-c("KRAS","SPRED1","NF1")

getPredStats<-function(genelist){
    #make the cluster
    dlist<-names(tumsByDis)
 #   dlist<-c("GBM","COAD")
    res<-do.call('rbind',lapply(dlist,function(ct,genelist){
        df<-crossGenePreds(genelist,cancerType=ct)
        print(paste('Finished',ct))
                                        #get offdiagonal predictions
        ndmat<-apply(df,2,unlist)*1-diag(nrow(df))
                                        #now collect mean values
        stats<-c(apply(ndmat,1,function(x) mean(x[x>0])),apply(ndmat,2,function(x) mean(x[x>0])))
        return(stats)
    },genelist))

    rownames(res)<-dlist
    write.table(res,file='pathwayStats.txt',sep='\t')

}

res<-getPredStats(genelist)
