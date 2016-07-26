##test cross-gene prediction
source("../../bin/elasticNetPred.R")

if(!exists('kras.mut'))
  kras.mut<<-getMutDataForGene("KRAS",TRUE)
if(!exists('nras.mut'))
  nras.mut<<-getMutDataForGene("NRAS",TRUE)
if(!exists("braf.mut"))
  braf.mut<<-getMutDataForGene("BRAF",TRUE)
if(!exists("nf1.mut"))
  nf1.mut<<-getMutDataForGene("NF1",TRUE)

##patient identifiers from the expression data
expr.pats<-toPatientId(colnames(alldat))

##TODO: get better list of tumors by disease.
tbd=tumsByDis#[which(names(tumsByDis)!='PANCAN')]
dis.inds<-lapply(tbd,function(x) which(expr.pats%in%toPatientId(x)))

##first get model of kras/nras
mutdata<-kras.mut
mut.pats=toPatientId(as.character(mutdata$Tumor))

exprdata<-alldat
expr.pats<-toPatientId(colnames(exprdata)[-1])

mut.vec=rep('WT',length(expr.pats))
mut.vec[match(mut.pats,expr.pats)]<-'MUTANT'

#build model
fit=model.build(exprdata,mut.vec,pref='KRAS')

##now try to predict model on NF1
nf1.vec<-rep("WT",length(expr.pats))
nf1.vec[match(toPatientId(as.character(nf1.mut$Tumor)),expr.pats)]<-'MUTANT'
res=model.pred(fit,exprdata,nf1.vec,pref='KRAS_to_NF1',doPlot=T)
df=res$Response
coeffs=res$Coeff

genes<-exprdata[which(coeffs[,1]>0),1]
gn<-sapply(as.character(genes),function(x) unlist(strsplit(x,split='|',fixed=T))[1])
