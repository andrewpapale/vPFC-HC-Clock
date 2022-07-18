# 2022-04-06 AndyP
# Wrapper script to plot vmPFC data

library(pracma)
library(tidyverse)

source('~/vmPFC/plot_mixed_by_vmPFC.R')
source('~/vmPFC/plot_emmeans_vmPFC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-feedback-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'network'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  plot_emmeans_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}

source('~/vmPFC/plot_mixed_by_vmPFC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-HC-feedback-network-outcome-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'outcome'
  toprocess <- 'network'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}

source('~/vmPFC/plot_mixed_by_vmPFC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-symmetry-feedback-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'symmetry_group'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}
source('~/vmPFC/plot_mixed_by_vmPFC.R')
source('~/vmPFC/plot_emmeans_vmPFC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-clock-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'network'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  #plot_emmeans_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}

source('~/vmPFC/plot_mixed_by_vmPFC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-symmetry-clock-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'symmetry_group'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}


source('~/vmPFC/plot_mixed_by_HC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-HC-axis-feedback-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'axis'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  plot_mixed_by_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}
source('~/vmPFC/plot_mixed_by_HC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-HC-axis-clock-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'axis'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  plot_mixed_by_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
}
source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
  model_str <- paste0('-vmPFC-HC-network-feedback-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  #toprocess <- 'symmetry-by-HC'
  toprocess <- 'network-by-HC'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}

source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-HC-symmetry-feedback-',i,'.Rdata')
  #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'symmetry-by-HC'
  #toprocess <- 'network-by-HC'
  toalign <- 'feedback'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}
source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
source('~/vmPFC/plot_emmeans_vmPFC_HC.R')
source('~/vmPFC/plot_emmeans_vmPFC_HC2.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
  model_str <- paste0('vmPFC-HC-network-clock-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'testing-HCwithin-scaled-'
  #toprocess <- 'symmetry-by-HC'
  toprocess <- 'network-by-HC'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  plot_emmeans_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  plot_emmeans_vmPFC_HC2(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}

source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-HC-symmetry-clock-',i,'.Rdata')
  #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'symmetry-by-HC'
  #toprocess <- 'network-by-HC'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}


source('~/vmPFC/plot_mixed_by_vmPFC_HC_entropy_split.R')
source('~/vmPFC/plot_emtrends_vmPFC_HC_entropy_split.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-HC-network-clock-Hsplit-',i,'.Rdata')
  #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'network-by-HC'
  #toprocess <- 'network-by-HC'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}

# vmPFC-HC-network-testing
source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-HC-network-clock-ranslopes-',i,'.Rdata')
  #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'final-model'
  toprocess <- 'network-by-HC'
  #toprocess <- 'network-by-HC'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}

