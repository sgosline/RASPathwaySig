#test survival analysis that separates based on expression and mutation

source("../../bin/sigSurvAnalysis.R",chdir=T)

gene='NF1'

for(dis in tcga.cancer.types){
  sapply(list(c("CDC27","CREBBP"),"CREBBP",'CDC27'),function(g){
    mut.sig<-survivalAnalysisByMutationAndExpression(dis,mutGene=c(g),exprGene=c(gene))
  })
#  print(paste('Significance of',gene,'in',dis,'is',expr.sig$pval,'(expression) and',mut.sig$pval,'(mutation)'))
  
}