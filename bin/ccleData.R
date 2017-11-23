##script to aquire and process CCLE data
##originally retrieved from https://ctd2.nci.nih.gov/dataPortal/
##ftp://caftpd.nci.nih.gov/pub/dcc_ctd2/TGen/CCLE_RNA-seq_Analysis/

require(synapseClient)
require(data.table)
synapseLogin()


mapEnsToHugo<-function(){
  mfile<-as.data.frame(fread('../../data/gencodeGeneTranscriptMap.csv',sep=',',header=T))
  ensg<-sapply(mfile$target_id,function(x) unlist(strsplit(unlist(strsplit(x,split='|',fixed=T))[2],split='.',fixed=T))[1])
  ret=mfile$gene
  names(ret)<-ensg
  return(ret)
  }

#'get Annotations for cell lines - specifically which disease class they belong to for clustering
getCellLineAnnotations<-function(){
    synid='syn7112975'
}


#'
#'Get CCLE Count data originally from
#'CTDD repository
getCCLEDataCounts<-function(removeDupes=TRUE){
    sid="syn5616077"
    tab<-as.data.frame(fread(synGet(sid)@filePath,sep=',',header=T))
    colnames(tab)[1]<-'EnsGene'

    map=mapEnsToHugo()
    idx=match(tab$EnsGene,names(map))
    tab$GeneSymbol=map[idx]
    if(removeDupes){
      tab=tab[-which(duplicated(tab$GeneSymbol)),]
      rnaMat<-tab[which(!is.na(tab$GeneSymbol)),-c(1,ncol(tab))]
      rownames(rnaMat)<-tab$GeneSymbol[which(!is.na(tab$GeneSymbol))]


    }else{
      rnaMat<-tab[,-1]
      rownames(rnaMat)<-tab[,1]
      rnaMat<-tab
    }

    cell.names=sapply(colnames(rnaMat),function(x) toupper(gsub('.','',x,fixed=T)))
    cell.names[grep("^X",cell.names)]<-sapply(cell.names[grep("^X",cell.names)],function(x) gsub("^X","",x))
    colnames(rnaMat)<-cell.names

    zv=which(apply(rnaMat,1,var)==0)
    if(length(zv)>0)
      rnaMat<-rnaMat[-zv,]

    return(rnaMat)

}
#'
#'Get CCLE TPM data
#'from CTDD repository
getCCLEDataTPM<-function(removeDupes=TRUE){
  sid="syn5616092"
  tab<-as.data.frame(fread(synGet(sid)@filePath,sep=',',header=T))
  colnames(tab)[1]<-'EnsGene'

  map=mapEnsToHugo()
  idx=match(tab$EnsGene,names(map))
  tab$GeneSymbol=map[idx]

  if(removeDupes){
    tab=tab[-which(duplicated(tab$GeneSymbol)),]
    rnaMat<-tab[which(!is.na(tab$GeneSymbol)),-c(1,ncol(tab))]
    rownames(rnaMat)<-tab$GeneSymbol[which(!is.na(tab$GeneSymbol))]

  }else{
    rnaMat<-tab[,-1]
    rownames(rnaMat)<-tab[,1]
    rnaMat<-tab
  }
  cell.names=sapply(colnames(rnaMat),function(x) toupper(gsub('.','',x,fixed=T)))
  cell.names[grep("^X",cell.names)]<-sapply(cell.names[grep("^X",cell.names)],function(x) gsub("^X","",x))
  colnames(rnaMat)<-cell.names

  zv=which(apply(rnaMat,1,var)==0)
  if(length(zv)>0)
    rnaMat<-rnaMat[-zv,]
  return(rnaMat)
  }
