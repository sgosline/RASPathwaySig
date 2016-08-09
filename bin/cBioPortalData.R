##create a script to parse cBIOportal data
#perhaps this will facilitate analysis across the different datasets...
library(cgdsr)
library(data.table)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

all.genes<-unique(fread('../../data/ucsc_kgXref_hg19_2015_10_29.csv')$geneSymbol)
all.studies<-getCancerStudies(mycgds)

#'getDiseaseSampleMapping creats a unified mapping of all samples
#'to various cell lines and disease profiles so that when
#'all are joined we can compare one to another
getDiseaseSampleMapping<-function(dis=''){
  
}

#various disease types in cbioporta.
broad.cancer.types=c('brca','cellline','lcll','desm','dlbc','esca','hnsc','luad','mbl','skcm','mm','nsclc','es','prad')
mskcc.cancer.types=c('acyc','acbc','blca','coadread','luad','mpnst','thyroid','prad','hnc','sarc','scco')
tcga.cancer.types<-c('laml','acc','blca','lgg','brca','cesc','chol','coadread','esca','gbm','hnsc','kich','kirc','kirp','lihc','luad','lusc','dlbc','lgggbm','meso','ov','nsclc','paad','thca','pcpg','prad','sarc','skcm','stad','tgct','thym','ucs','ucec','uvm')

##not all have counts
cell.line.tiss<-c('CENTRAL_NERVOUS_SYSTEM','BONE','PROSTATE','STOMACH','URINARY_TRACT','OVARY','HAEMATOPOIETIC_AND_LYMPHOID_TISSUE','KIDNEY','THYROID','SKIN','SOFT_TISSUE','SALIVARY_GLAND','LUNG','PLEURA','LIVER','ENDOMETRIUM','PANCREAS','BREAST','UPPER_AERODIGESTIVE_TRACT','LARGE_INTESTINE','AUTONOMIC_GANGLIA','OESOPHAGUS','BILIARY_TRACT','SMALL_INTESTINE')

getDisMutationData<-function(dis='',study='tcga'){
  #if disease is blank will get all diseases
  ind=grep(paste(tolower(dis),study,sep='_'),all.studies$cancer_study_id)
  if(length(ind)==0)
    return(NULL)
  mycancerstudy<-all.studies$cancer_study_id[ind]
  expr.list<-lapply(mycancerstudy,function(cs){
    print(cs)
    caseLists<-getCaseLists(mycgds,cs)
    allprofs<-getGeneticProfiles(mycgds,cs)[,1]
    profile=allprofs[grep('mutations',allprofs)]
    seqSamps=caseLists$case_list_id[grep('sequenced',caseLists$case_list_id)]
    gene.groups=split(all.genes, ceiling(seq_along(all.genes)/500))
    dat<-lapply(gene.groups,function(g) getProfileData(mycgds,g,profile,seqSamps))
    ddat<-matrix()
    for(i in which(sapply(dat,nrow)!=0)){
      ddat<-cbind(ddat,dat[[i]])
    }
    nans<-which(apply(ddat,2,function(x) all(is.nan(x)||is.na(x))))
    # nas<-which(apply(ddat,2,function(x) all(is.na(x))))
    ddat<-ddat[,-nans]
    ##now set to binary matrix
    dfdat<-apply(ddat,1,function(x){
      sapply(unlist(x),function(y) !is.na(y) && y!='NaN')
    })
    
    return(dfdat)
    
  })
  
  comm.genes<-rownames(expr.list[[1]])
  for(i in 2:length(expr.list))
    comm.genes<-intersect(comm.genes,rownames(expr.list[[i]]))
  full.dat<-do.call('cbind',lapply(expr.list,function(x) x[comm.genes,]))
  return(full.dat)
  
}

#'formats TCGA expression data into a single matrix, often combining
#'samples from multiple studies
getDisExpressionData<-function(dis='',study='tcga',getZscores=FALSE){
  #if disease is blank will get all diseases
  ind=grep(paste(tolower(dis),study,sep='_'),all.studies$cancer_study_id)
  if(length(ind)==0)
    return(NULL)
  
  mycancerstudy<-all.studies$cancer_study_id[ind]
  expr.list<-lapply(mycancerstudy,function(cs){
    print(cs)
    caseLists<-getCaseLists(mycgds,cs)
    allprofs<-getGeneticProfiles(mycgds,cs)[,1]
    rnaseqs<-allprofs[grep('rna_seq',allprofs)]
    zscores<-grep('Zscores',rnaseqs)
    profile=rnaseqs[zscores]
    if(!getZscores)
      profile=rnaseqs[-zscores]
    if(length(profile)>1)
      profile=profile[grep('v2',profile)]
    mrnaSamps=caseLists$case_list_id[grep('rna_seq',caseLists$case_list_id)]
    if(length(mrnaSamps)>1)
      mrnaSamps=mrnaSamps[grep('v2',mrnaSamps)]
    gene.groups=split(all.genes, ceiling(seq_along(all.genes)/500))
    dat<-lapply(gene.groups,function(g) getProfileData(mycgds,g,profile,mrnaSamps))
    ddat<-matrix()
    for(i in which(sapply(dat,nrow)!=0)){
      ddat<-cbind(ddat,dat[[i]])
    }
    nans<-which(apply(ddat,2,function(x) all(is.na(x))||mean(x,na.rm=T)==0))
    if(length(nans)>0)
      ddat<-ddat[,-nans]
    #ddat<-ddat[,-1]
    ddat<-data.frame(t(ddat))
    ddat
  })
  comm.genes<-rownames(expr.list[[1]])
  for(i in 2:length(expr.list))
    comm.genes<-intersect(comm.genes,rownames(expr.list[[i]]))
  
  ##now combine all samples by doing a cbind
  full.dat<-do.call('cbind',lapply(expr.list,function(x) x[comm.genes,]))
  return(full.dat)
  
}


#'get CCLE expressiond ata. can get z score or affy data, not sure which to do yet
getCcleExpressionData<-function(tiss='',getZscores=FALSE){
  
  mycancerstudy='cellline_ccle_broad'
  
  profile='cellline_ccle_broad_mrna_median_Zscore' ##eventually test out both
  if(!getZscores)
    profile='cellline_ccle_broad_mrna'
  
  caseLists<-getCaseLists(mycgds,mycancerstudy)
  
  ##get those samples with mRNA expression data
  mrnaSamps=caseLists$case_list_id[grep('mrna',caseLists$case_list_id)]
  
  #cbio seems to handle chunks of 500 or so
  gene.groups=split(all.genes, ceiling(seq_along(all.genes)/500))
  dat<-lapply(gene.groups,function(g) getProfileData(mycgds,g,profile,mrnaSamps))

  ddat<-matrix()
  for(i in which(sapply(dat,nrow)!=0)){
    ddat<-cbind(ddat,dat[[i]])
  }
  
  nans<-which(apply(ddat,2,function(x) all(is.nan(x))))
  ddat<-ddat[,-nans]
  ddat<-ddat[,-1]
  ddat<-data.frame(t(ddat))
  
  ##tissue here
  if(tiss!=''){
    cols<-grep(tiss,colnames(ddat))
    print('Selecting',length(cols),'cell lines for tissue',tiss)
  }else{
    cols<-1:ncol(ddat)
  }
  
  
  return(ddat[,cols])
  
}
#'get CCLE mutation dat
getCcleMutationData<-function(tiss=''){
  mycancerstudy='cellline_ccle_broad'
  profile="cellline_ccle_broad_mutations" ##think about adding CNA data
  caseLists<-getCaseLists(mycgds,mycancerstudy)
  mutSamps<-caseLists$case_list_id[grep("sequenced",caseLists[,1])]
  gene.groups=split(all.genes, ceiling(seq_along(all.genes)/500))
  dat<-lapply(gene.groups,function(g) getProfileData(mycgds,g,profile,mutSamps))
  
  ddat<-matrix()
  for(i in which(sapply(dat,nrow)!=0)){
    ddat<-cbind(ddat,dat[[i]])
  }
  nans<-which(apply(ddat,2,function(x) all(is.nan(x)||is.na(x))))
  # nas<-which(apply(ddat,2,function(x) all(is.na(x))))
  ddat<-ddat[,-nans]
  ##now set to binary matrix
  dfdat<-apply(ddat,1,function(x){
    sapply(unlist(x),function(y) !is.na(y) && y!='NaN')
  })
  
  ##tissue here
  if(tiss!=''){
    cols<-grep(tiss,colnames(dfdat))
    print('Selecting',length(cols),'cell lines for tissue',tiss)
  }else{
    cols<-1:ncol(dfdat)
  }
  
  return(dfdat[,cols])
  

}