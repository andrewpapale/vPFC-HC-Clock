# 2022-04-06 AndyP
# Wrapper script to plot vmPFC data

source('~/vmPFC/plot_mixed_by_vmPFC.R')
for (i in 2:10){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-last-model-w_expl_code-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'network'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}
