# 2022-04-06 AndyP
# Wrapper script to plot vmPFC data

totest <- 'dHtest'
toprocess <- 'network'
toalign <- 'feedback'
behavmodel <- 'compressed'
source('~/vmPFC/plot_mixed_by_vmPFC.R')
source('~/vmPFC/plot_emmeans_vmPFC.R')
#source('~/vmPFC/plot_emtrends_vmPFC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-',toprocess,'-',toalign,'-',totest,'-',i,'.Rdata')
  model_str <- '2022-12-13-vmPFC-network-feedback-dHtest-3.Rdata'
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  if (strcmp(toprocess,'symmetry')){
    toprocess <- 'symmetry_group'
  }
  plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  #plot_emmeans_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  #plot_emtrends_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}
