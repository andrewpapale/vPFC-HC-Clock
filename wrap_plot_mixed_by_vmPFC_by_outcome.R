source('~/vmPFC/plot_mixed_by_vmPFC_by_outcome.R')
for (i in 1:8){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-feedback-split-by-condition-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model-'
  toprocess <- 'network'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC_by_outcome(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}
