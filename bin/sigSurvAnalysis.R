source("../../bin/cBioPortalData.R")
library(survival)
library(survminer)
library(OIsurv) # Aumatically loads KMsurv

survivalAnalysisByMutation<-function(dis,genelist){
  
  #first get mutation status for genelist
  mut.data<-getDisMutationData(dis,genelist=genelist)
  
  ##then get patient data
  clin.data<-getDisClinicalData(dis=dis)
  
  ##then bracket patient data
  mutated<-data.frame(Mutation=apply(mut.data,1, any)*1)
  mutated$Sample=rownames(mutated)

  df<-clin.data%>%select(OS_MONTHS,OS_STATUS,Sample)%>%left_join(mutated,by='Sample')%>%mutate(event=(OS_STATUS=="DECEASED")*1)
  df<-df%>%mutate(MutLabel=ifelse(Mutation,paste(paste(genelist,collapse=','),'mutated'),'wildtype'))
  fit<-survfit(Surv(df$OS_MONTHS,df$event,)~df$MutLabel)
  res<-coxph(Surv(df$OS_MONTHS,df$event,)~df$MutLabel)
  
  pval=summary(res)$logtest[3]
  ##then compute curve, p-value, plot if signif.
  fname=''
  if(pval<0.05){
    ggsurvplot(fit = fit,data=df,pval=TRUE)+ggtitle(paste(paste(genelist,collapse='_'),'mutants in',dis))
    fname=paste(paste(genelist,collapse='_'),'mutantStatusIn',dis,'KMPlot.png',sep='')
    ggsave(file=fname)
    
  }
  return(list(pval=pval,file=fname))
}


survivalAnalysisByExpression<-function(dis,genelist){
  #first get gene list
  expr.data<-getDisExpressionData(dis=dis,getZscores=TRUE,genelist=genelist)
  mean.exp<-apply(expr.data,1,mean,na.rm=T)
  
  ##bracket by mean quartile
  quarts<-quantile(mean.exp,c(0.25,0.5,0.75))
  
  expr.status<-rep("Low",length(mean.exp))
  expr.status[which(mean.exp>quarts[2])]<-'Mid'
  expr.status[which(mean.exp>quarts[3])]<-'High'
  
  geneExpr<-data.frame(Expr=expr.status,Sample=names(mean.exp))
  
  ##then get patient data
  clin.data<-getDisClinicalData(dis=dis)
  
  df<-clin.data%>%select(OS_MONTHS,OS_STATUS,Sample)%>%left_join(geneExpr,by='Sample')%>%mutate(event=(OS_STATUS=="DECEASED")*1)
  #df<-df%>%mutate(ExprLabel=ifelse(Mutation,paste(paste(genelist,collapse=','),'mutated'),'wildtype'))
  fit<-survfit(Surv(df$OS_MONTHS,df$event,)~df$Expr)
  res<-coxph(Surv(df$OS_MONTHS,df$event,)~df$Expr)
  
  pval=summary(res)$logtest[3]
  ##then compute curve, p-value, plot if signif.
  if(pval<0.05){
    ggsurvplot(fit = fit,data=df,pval=TRUE)+ggtitle(paste(paste(genelist,collapse='_'),'expression in',dis))
    fname=paste(paste(genelist,collapse='_'),'exprStatusIn',dis,'KMPlot.png',sep='')
    ggsave(file=fname)
    
  }
  return(list(pval=pval,file=fname))
  
  
}