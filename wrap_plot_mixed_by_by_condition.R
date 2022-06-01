# 2022-05-12 AndyP
# Wrapper script to plot vmPFC data

library(pracma)
library(tidyverse)


# vmPFC-by-condition
source('~/vmPFC/plot_mixed_by_vmPFC_HC_by_condition.R')
source('~/vmPFC/plot_emtrends_vmPFC_HC_by_condition.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('vmPFC-HC-network-clock-split-by-condition-',i,'.Rdata')
  #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model-'
  toprocess <- 'network-by-HC'
  #toprocess <- 'network-by-HC'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_mixed_by_vmPFC_HC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  plot_emtrends_vmPFC_HC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}

# # vmPFC-by-condition
# source('~/vmPFC/plot_mixed_by_vmPFC_HC_by_condition.R')
# for (i in 1){
#   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   model_str <- paste0('vmPFC-HC-symmetry-clock-split-by-condition-',i,'.Rdata')
#   #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
#   model_str <- Sys.glob(paste0('*',model_str))
#   load(model_str)
#   model_iter <- i
#   totest <- 'final-model'
#   toprocess <- 'symmetry-by-HC'
#   #toprocess <- 'network-by-HC'
#   toalign <- 'clock'
#   behavmodel <- 'compressed'
#   hc_LorR <- 'LR'
#   plot_mixed_by_vmPFC_HC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
# }

# vmPFC-by-condition
source('~/vmPFC/plot_mixed_by_vmPFC_HC_by_condition.R')
source('~/vmPFC/plot_emtrends_vmPFC_HC_by_condition.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('vmPFC-HC-network-feedback-split-by-condition-',i,'.Rdata')
  #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model-'
  toprocess <- 'network-by-HC'
  #toprocess <- 'network-by-HC'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_mixed_by_vmPFC_HC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  plot_emtrends_vmPFC_HC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}

# vmPFC-by-condition
# source('~/vmPFC/plot_mixed_by_vmPFC_HC_by_condition.R')
# source('~/vmPFC/plot_emtrends_vmPFC_HC_by_condition.R')
# for (i in 1){
#   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   model_str <- paste0('vmPFC-HC-symmetry-feedback-split-by-condition-',i,'.Rdata')
#   #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
#   model_str <- Sys.glob(paste0('*',model_str))
#   load(model_str)
#   model_iter <- i
#   totest <- 'final-model'
#   toprocess <- 'symmetry-by-HC'
#   #toprocess <- 'network-by-HC'
#   toalign <- 'feedback'
#   behavmodel <- 'compressed'
#   hc_LorR <- 'LR'
#   plot_mixed_by_vmPFC_HC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
# }

source('~/vmPFC/plot_mixed_by_vmPFC_by_condition.R')
source('~/vmPFC/plot_emmeans_vmPFC_by_condition.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-clock-split-by-condition-trial_x_entropy-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model-trial_by_entropy-'
  toprocess <- 'network'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  plot_emmeans_vmPFC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}

source('~/vmPFC/plot_mixed_by_vmPFC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-symmetry-clock-split-by-condition-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'symmetry_group'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}

source('~/vmPFC/plot_mixed_by_vmPFC_by_condition.R')
source('~/vmPFC/plot_emmeans_vmPFC_by_condition.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-feedback-split-by-condition-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model-'
  toprocess <- 'network'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  plot_emmeans_vmPFC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}

# source('~/vmPFC/plot_mixed_by_vmPFC.R')
# for (i in 1){
#   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   model_str <- paste0('-vmPFC-symmetry-feedback-split-by-condition-',i,'.Rdata')
#   model_str <- Sys.glob(paste0('*',model_str))
#   load(model_str)
#   model_iter <- i
#   totest <- 'final-model'
#   toprocess <- 'symmetry_group'
#   toalign <- 'feedback'
#   behavmodel <- 'compressed'
#   plot_mixed_by_vmPFC_by_condition(ddf,toalign,toprocess,totest,behavmodel,model_iter)
# }
