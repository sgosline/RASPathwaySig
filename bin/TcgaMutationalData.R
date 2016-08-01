##mutational data analysis


library(synapseClient)
library(data.table)
library(ggplot2)
library(R.utils)
library(dplyr)
synapseLogin()

mafdir='../../data/somatic_mafs_cleaned'
maffiles=list.files(mafdir)

toPatientId<-function(vec){
  #takes vector of tumor ids and reduces to patient
  sapply(as.character(vec),function(x) paste(unlist(strsplit(x,split='-'))[1:4],collapse='-'))
}


if(!exists('combinedMaf')){
  gpath<-synGet('syn4924181')@filePath
  df1.2.maf<-gsub('.gz','',gpath)

  if(!file.exists(df1.2.maf))
    gunzip(gpath,remove=F)

  combinedMaf<-as.data.frame(fread(df1.2.maf))
  combinedMaf$Patient<-toPatientId(combinedMaf$Tumor_Sample_Barcode)
}

if(!exists('pat.dis')){
  pat.dis<<-synTableQuery('SELECT acronym,samples from syn3281840')
  tumsByDis<<-sapply(unique(pat.dis@values$acronym),function(x) pat.dis@values$samples[which(pat.dis@values$acronym==x)])
  tumsByDis$COADREAD=c(tumsByDis$COAD,tumsByDis$READ)
}


names(maffiles)<-sapply(maffiles,function(x) unlist(strsplit(x,split='_'))[1])

getTumorIdsFromMaf<-function(maffile){
                                        #get a unique set of patients
    disname=unlist(strsplit(basename(maffile),split='_'))[1]
    tab<-as.data.frame(fread(maffile))
    pats<-toPatientId(as.character(unique(tab$Tumor_Sample_Barcode)))
    return(pats)

}


getMutStatusByDisease<-function(disname){
    dpats<-unique(toPatientId(tumsByDis[[toupper(disname)]]))
  #  dpats<-dpats[which(names(dpats)!='PANCAN')]
    print(paste('Found',length(dpats),'samples for',disname))
    maftab<-combinedMaf[which(combinedMaf$Patient%in%dpats),]
    maftab$Disease<-rep(disname,nrow(maftab))
    res<-select(maftab,match(c('Disease',"Hugo_Symbol",'Variant_Classification','PolyPhen','Tumor_Sample_Barcode','HGVSp'),colnames(maftab)))
    colnames(res)<-c('Disease',"Gene",'VariantClassification','Score','Tumor','AAChange')
    return(res)
}

getGeneStatusByDisease<-function(disname,gene='NF1'){
    maftab<-getMutStatusByDisease(disname)
    genetab=maftab[which(maftab$Gene==gene),]
    print(paste('Found',nrow(maftab),'entries for',disname,'and',nrow(genetab),'of those for',gene))
  #  p<-ggplot(genetab)+
  #    geom_bar(aes(Variant_Classification))
    return(genetab)
}

getGeneStatusFromMaf<-function(maffile,gene='NF1'){
  disname=unlist(strsplit(basename(maffile),split='_'))[1]
  tab<-as.data.frame(fread(maffile))
  genetab=tab[which(tab$gene_name==gene),]
  print(paste("Found",nrow(genetab),'entries with',gene,'mutated in',disname))
  p<-ggplot(genetab)+
  geom_bar(aes(Variant_Classification))
  return(genetab)
}

getNF1Data<-function(plot=FALSE){
    getMutDataForGene("NF1",plot)
}

getMutDataForGene<-function(geneName,plot=FALSE,cancerType=NA){
  if(is.na(cancerType)){
    cancerType=setdiff(names(tumsByDis),'PANCAN')
  }else{
    cancerType=c(cancerType)
  }
  df<-do.call('rbind',lapply(cancerType,function(x){
    res<-getGeneStatusByDisease(x,geneName)
   # res$Disease<-rep(x,nrow(res))
    return(res)}))


    #df<-data.frame(Disease=unlist(dis),VariantClassification=unlist(variant),Score=unlist(score),Tumor=unlist(tumor_id),AAChange=unlist(aa_change))
  if(plot){
    png(paste(geneName,'mutationstatus.png',sep=''))
    p<-ggplot(df)+geom_bar(aes(VariantClassification,fill=Disease)) + ggtitle(paste(geneName,"mutations in TCGA")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
    dev.off()
  }
  return(df)
}
