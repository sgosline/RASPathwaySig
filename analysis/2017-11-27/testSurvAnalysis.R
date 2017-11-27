##compare surv analysis calls for mutations and expression signatures across cancers/and NF1/NF2 signatures

source("../../bin/sigSurvAnalysis.R")

gene='NF1'

for(dis in tcga.cancer.types){
  mut.sig<-survivalAnalysisByMutation(dis,c(gene))
  expr.sig<-survivalAnalysisByExpression(dis,c(gene))
  print(paste('Significance of',gene,'in',dis,'is',expr.sig$pval,'(expression) and',mut.sig$pval,'(mutation)'))
  
}

gene='CREBBP'
for(dis in tcga.cancer.types){
  mut.sig<-survivalAnalysisByMutation(dis,c(gene))
  expr.sig<-survivalAnalysisByExpression(dis,c(gene))
  print(paste('Significance of',gene,'in',dis,'is',expr.sig$pval,'(expression) and',mut.sig$pval,'(mutation)'))
  
}
gene='CDC27'
for(dis in tcga.cancer.types){
  mut.sig<-survivalAnalysisByMutation(dis,c(gene))
  expr.sig<-survivalAnalysisByExpression(dis,c(gene))
  print(paste('Significance of',gene,'in',dis,'is',expr.sig$pval,'(expression) and',mut.sig$pval,'(mutation)'))
  
}
gene=c("CREBBP",'CDC27','NF1')
for(dis in tcga.cancer.types){
  mut.sig<-survivalAnalysisByMutation(dis,c(gene))
  expr.sig<-survivalAnalysisByExpression(dis,c(gene))
  print(paste('Significance of',paste(gene,collapse='_'),'in',dis,'is',expr.sig$pval,'(expression) and',mut.sig$pval,'(mutation)'))
  
}