# 2022-04-06 AndyP
# Wrapper script to plot vmPFC data
library(pracma)
totest <- ''
toprocess <- 'symmetry'
toalign <- 'feedback'
behavmodel <- ''
if (strcmp(behavmodel,'explore')){
  behavmodel <- paste0(behavmodel,'-')
}
if (strcmp(totest,'online')){
  totest <- paste0('-',totest)
}
source('~/vmPFC/plot_mixed_by_vmPFC.R')
source('~/vmPFC/plot_emmeans_vmPFC.R')
source('~/vmPFC/plot_emtrends_vmPFC.R')
for (i in 1:2){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-',behavmodel,toprocess,'-',toalign,totest,'-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  if (strcmp(toprocess,'symmetry')){
    toprocess <- 'symmetry_group'
  }
  if (i==1 & strcmp(behavmodel,'explore')){
    totest <- 'entropy'
  } else if (i==2 & strcmp(behavmodel,'explore')){
    totest <- 'value'
  }
  if (strcmp(totest,'-online')){
    totest <- 'online'
  }
  plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  #plot_emmeans_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  #plot_emtrends_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  if (strcmp(toprocess,'symmetry_group')){
    toprocess <- 'symmetry'
  }  
  if (i==1 & strcmp(behavmodel,'explore')){
    totest <- ''
  } else if (i==2 & strcmp(behavmodel,'explore')){
    totest <- ''
  }
  if (strcmp(totest,'online')){
    totest <- paste0('-',totest)
  }
}
