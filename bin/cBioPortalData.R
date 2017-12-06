##create a script to parse cBIOportal data
#perhaps this will facilitate analysis across the different datasets...
library(cgdsr)
library(tidyverse)
library(data.table)

script.dir <- dirname(sys.frame(1)$ofile)
all.genes<<-unique(fread('../data/ucsc_kgXref_hg19_2015_10_29.csv')$geneSymbol)


#'getSamplesForDisease creates a unified mapping of all samples
#'to various cell lines and disease profiles so that when
#'all are joined we can compare one to another
getSamplesForDisease<-function(dis='',study='tcga'){
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  all.studies<-getCancerStudies(mycgds)
  ind=grep(paste(tolower(dis),paste(study,'$',sep=''),sep='_'),all.studies$cancer_study_id)
    if(length(ind)==0)
    return(c())
  mycancerstudy<-all.studies$cancer_study_id[ind]

  sampList<-lapply(mycancerstudy,function(cs){
      caseLists<-getCaseLists(mycgds,cs)
      samps<-unlist(strsplit(caseLists[match(paste(cs,'all',sep='_'),caseLists[,1]),5],split=' '))
                                        #      allprofs<-getGeneticProfiles(mycgds,cs)[,1]
      print(paste('Found',length(samps),'for',cs))
      return(samps)
  })
  all.samps<-unique(unlist(sampList))
  print(paste("Found",length(all.samps),'samples for',study,dis))
  return(all.samps)
}

#various disease types in cbioporta.
broad.cancer.types=c('brca','cellline','lcll','desm','dlbc','esca','hnsc','luad','mbl','skcm','mm','nsclc','es','prad')
mskcc.cancer.types=c('acyc','acbc','blca','coadread','luad','mpnst','thyroid','prad','hnc','sarc','scco')
tcga.cancer.types<<-c('laml','acc','blca','lgg','brca','cesc','chol','coadread','esca','gbm','hnsc','kich','kirc','kirp','lihc','luad','lusc','dlbc','ov','paad','thca','pcpg','prad','sarc','skcm','stad','tgct','thym','ucs','ucec','uvm')#meso has no sequence data

##not all have counts
cell.line.tiss<-c('CENTRAL_NERVOUS_SYSTEM','BONE','PROSTATE','STOMACH','URINARY_TRACT','OVARY','HAEMATOPOIETIC_AND_LYMPHOID_TISSUE','KIDNEY','THYROID','SKIN','SOFT_TISSUE','SALIVARY_GLAND','LUNG','PLEURA','LIVER','ENDOMETRIUM','PANCREAS','BREAST','UPPER_AERODIGESTIVE_TRACT','LARGE_INTESTINE','AUTONOMIC_GANGLIA','OESOPHAGUS','BILIARY_TRACT','SMALL_INTESTINE')

getDisMutationData<-function(dis='',study='tcga',genelist=all.genes){

  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  all.studies<-getCancerStudies(mycgds)

  if(tolower(dis)=='alltcga')
    dval=''
  else
    dval=dis
  #if disease is blank will get all diseases
  ind=grep(paste(tolower(dval),paste(study,'$',sep=''),sep='_'),all.studies$cancer_study_id)
  print(paste('found',length(ind),study,'samples for disease',dis))

  if(length(ind)==0)
    return(NULL)

  mycancerstudy<-all.studies$cancer_study_id[ind]

  if(study=='tcga')
      mycancerstudy=intersect(mycancerstudy,paste(tcga.cancer.types,study,sep='_'))
  else if(study=='mskcc')
      mycancerstudy=intersect(mycancerstudy,paste(mskcc.cancer.types,study,sep='_'))
  else if(study=='broad')
      mycancerstudy=intersect(mycancerstudy,paste(broad.cancer.types,study,sep='_'))

  expr.list<-lapply(mycancerstudy,function(cs){
    print(paste(cs,study,'Mutation data'))
    caseLists<-getCaseLists(mycgds,cs)
    allprofs<-getGeneticProfiles(mycgds,cs)[,1]
    profile=allprofs[grep('mutations',allprofs)]
    seqSamps=caseLists$case_list_id[grep('sequenced',caseLists$case_list_id)]
    if(length(genelist)>800)
      gene.groups=split(genelist, ceiling(seq_along(genelist)/400))
    else
      gene.groups=list(genelist)
    dat<-lapply(gene.groups,function(g) getProfileData(mycgds,g,profile,seqSamps))
    ddat<-matrix()
    if(length(dat)==1){
      ddat<-t(as.data.frame(dat))
    }else{
    for(i in which(sapply(dat,nrow)!=0)){
      ddat<-cbind(ddat,dat[[i]])
    }
    }
    
#    nans<-which(apply(ddat,2,function(x) all(is.nan(x)||is.na(x))))
    # nas<-which(apply(ddat,2,function(x) all(is.na(x))))
#    if(nrow(ddat)>1)
#        ddat<-ddat[,-nans]
    ##now set to binary matrix
    dfdat<-apply(ddat,1,function(x){
      sapply(unlist(x),function(y) !is.na(y) && y!='NaN')
    })

    return(dfdat)

  })
  if(length(expr.list)>1){
    comm.genes<-c()#rownames(expr.list[[1]])
    for(i in 1:length(expr.list))
        comm.genes<-union(comm.genes,rownames(expr.list[[i]]))

    full.dat<-do.call('cbind',lapply(expr.list,function(x){
        missing<-setdiff(comm.genes,rownames(x))
        print(length(missing))
        if(length(missing)>0)
            dat<-rbind(x[intersect(rownames(x),comm.genes),],
                       t(sapply(missing,function(y,x) rep(FALSE,ncol(x)),x)))
        else{
            dat<-x[intersect(rownames(x),comm.genes),]
        }
        colnames(dat)<-colnames(x)
        dat<-dat[comm.genes,]
        return(dat)
       }))


   }
  else{
    full.dat<-expr.list[[1]]
  }
  return(full.dat)

}

#'get TCGA clientd ata required for survival analysis
#'
getDisClinicalData<-function(dis='',study='tcga'){
  mycgds = CGDS("http://www.cbioportal.org/")
  all.studies<-getCancerStudies(mycgds)
  # Get available case lists (collection of samples) for a given cancer study

  if(tolower(dis)=='alltcga')
    dval=''
  else
    dval=dis
  
  ind=grep(paste(tolower(dval),paste(study,'$',sep=''),sep='_'),all.studies$cancer_study_id)
  print(paste('found',length(ind),study,'studies for disease',dis))
  if(length(ind)==0){
    return(NULL)
  }
  
  mycancerstudy<-all.studies$cancer_study_id[ind]
  
  clin.list<-lapply(mycancerstudy,function(cs){
    print(paste(cs,study,'Clinical data'))
    caseLists<-getCaseLists(mycgds,cs)%>%select(case_list_id,case_ids)
    clinDat<-getClinicalData(mycgds,caseLists[grep('all',caseLists$case_list_id),1])
    clinDat%>%select(OS_MONTHS,OS_STATUS,DFS_MONTHS,DFS_STATUS)%>%mutate(Sample=rownames(clinDat),Study=cs)
  })
  
  full.dat<-do.call('rbind',clin.list)
  full.dat

  
}

#'formats TCGA expression data into a single matrix, often combining
#'samples from multiple studies
getDisExpressionData<-function(dis='',study='tcga',getZscores=FALSE,genelist=all.genes){
  #if disease is blank will get all diseases
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  all.studies<-getCancerStudies(mycgds)

  if(tolower(dis)=='alltcga')
    dval=''
  else
    dval=dis

  ind=grep(paste(tolower(dval),paste(study,'$',sep=''),sep='_'),all.studies$cancer_study_id)
  print(paste('found',length(ind),study,'samples for disease',dis))
  if(length(ind)==0){
    return(NULL)
  }

  mycancerstudy<-all.studies$cancer_study_id[ind]
  expr.list<-lapply(mycancerstudy,function(cs){
    print(paste(cs,study,'Expression data'))
    caseLists<-getCaseLists(mycgds,cs)
    allprofs<-getGeneticProfiles(mycgds,cs)[,1]
    rnaseqs<-allprofs[grep('rna_seq',allprofs)]
    if(length(rnaseqs)==0)
      rnaseqs<-allprofs[grep('mrna',allprofs)]
    if(length(grep('merged',rnaseqs))>0)
      rnaseqs<-rnaseqs[-grep('merged',rnaseqs)]

    zscores<-grep('Zscores',rnaseqs)
    profile=rnaseqs[zscores]
    if(!getZscores)
      profile=rnaseqs[-zscores]
    if(length(profile)>1)
      profile=profile[grep('v2',profile)]
    mrnaSamps=caseLists$case_list_id[grep('rna_seq',caseLists$case_list_id)]
    if(length(mrnaSamps)>1)
      mrnaSamps=mrnaSamps[grep('v2',mrnaSamps)]
    else if(length(mrnaSamps)==0)
      mrnaSamps=caseLists$case_list_id[grep('mrna',caseLists$case_list_id)]
    if(length(genelist)>400)
        gene.groups=split(all.genes, ceiling(seq_along(all.genes)/400))
    else
      gene.groups=list(genelist)
    dat<-lapply(gene.groups,function(g) getProfileData(mycgds,g,profile,mrnaSamps))
    ddat<-matrix()
    if(length(dat)==1){
      ddat<-t(as.data.frame(dat))
    }else{
      for(i in which(sapply(dat,nrow)!=0)){
        ddat<-cbind(ddat,dat[[i]])
      }
    }
    nans<-which(apply(ddat,2,function(x) all(is.na(x))||mean(x,na.rm=T)==0))
    if(length(nans)>0 && nrow(ddat)>1)
      ddat<-ddat[,-nans]
    #ddat<-ddat[,-1]
    ddat<-data.frame(t(ddat))
    ddat
  })
  if(length(expr.list)>1){
    comm.genes<-rownames(expr.list[[1]])
    for(i in 2:length(expr.list))
      comm.genes<-intersect(comm.genes,rownames(expr.list[[i]]))
    ##now combine all samples by doing a cbind
    full.dat<-do.call('cbind',lapply(expr.list,function(x) x[comm.genes,]))
  }
  else{
    full.dat<-expr.list[[1]]
  }
  new.dat<-full.dat
  #print(dim(full.dat))
  #new.dat<-apply(full.dat,2,function(x) as.numeric)
  #rownames(new.dat)<-rownames(full.dat)
  return(new.dat)

}


#'get CCLE expressiond ata. can get z score or affy data, not sure which to do yet
getCcleExpressionData<-function(tiss='',getZscores=FALSE){
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  #all.studies<-getCancerStudies(mycgds)

  mycancerstudy='cellline_ccle_broad'


  if(tolower(tiss)=='allccle')
    tval=''
  else
    tval=tiss

  mprofile<-'cellline_ccle_broad_mrna_median_Zscores' ##eventually test out both
  if(!getZscores)
    mprofile<-'cellline_ccle_broad_mrna'
  profile<-mprofile

  caseLists<-getCaseLists(mycgds,mycancerstudy)
  print(paste('Collecting CCLE expression data for',tiss,'tissue'))

  ##get those samples with mRNA expression data
  mrnaSamps<<-caseLists$case_list_id[grep('mrna',caseLists$case_list_id)]

  #cbio seems to handle chunks of 500 or so
  gene.groups=split(all.genes, ceiling(seq_along(all.genes)/400))
  dat<-lapply(gene.groups,function(g) getProfileData(mycgds,g,profile,mrnaSamps))

  ddat<-matrix()
  for(i in which(sapply(dat,nrow)!=0)){
    ddat<-cbind(ddat,dat[[i]])
  }

  nans<-which(apply(ddat,2,function(x) all(is.nan(x))))
  if(length(nans)>0)
    ddat<-ddat[,-nans]
  ddat<-ddat[,-1]
  ddat<-data.frame(t(ddat))

  ##tissue here
  if(tval!=''){
    cols<-grep(tval,colnames(ddat))
    print(paste('Selecting',length(cols),'cell lines for tissue',tiss))
  }else{
    cols<-1:ncol(ddat)
  }


  return(ddat[,cols])

}
#'get CCLE mutation dat
getCcleMutationData<-function(tiss=''){
  mycancerstudy<-'cellline_ccle_broad'
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  #all.studies<-getCancerStudies(mycgds)

  if(tolower(tiss)=='allccle')
    tval=''
  else
    tval=tiss


  profile<-"cellline_ccle_broad_mutations" ##think about adding CNA data
  caseLists<-getCaseLists(mycgds,mycancerstudy)
  print('Got caselists')
  mutSamps<-caseLists$case_list_id[grep("sequenced",caseLists[,1])]
  print(paste('Collecting CCLE mutation data for',tiss,'tissue'))
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
  if(tval!=''){
    cols<-grep(tval,colnames(dfdat))
    print(paste('Selecting',length(cols),'cell lines for tissue',tiss))
  }else{
    cols<-1:ncol(dfdat)
  }
  return(dfdat[,cols])

}
