# 2022-03-31 AndyP
# mixed_by model building for AIC comparison

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26




# message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
# load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
# vmPFC <- fb_comb
# vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
# rm(fb_comb)
# 
# vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
# vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
# vmPFC <- vmPFC %>% mutate(side_vmpfc = as.factor(substr(region, nchar(region), nchar(region))))
# 
# source('~/vmPFC/get_trial_data_vmPFC.R')
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# 
# 
# df <- df %>% select(id,run,run_trial,v_max_wi,v_entropy_wi,v_entropy_wi_change,rewFunc,rt_csv_sc,outcome,trial_neg_inv_sc,total_earnings,rt_vmax_change,log10kld3,rt_vmax_lag_sc,v_entropy_wi_full,rt_vmax_full)
# Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# rm(vmPFC)
# Q$outcome <- relevel(as.factor(Q$outcome),ref="Omission")
# 
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# splits = c("evt_time","network")
# 
# # working up to this...
# # decode_formula_1 = formula(~ trial_neg_inv_sc +
# #                            rt_csv_sc + 
# #                            v_entropy_wi*HCwithin + 
# #                            v_entropy_wi_change*HCwithin + 
# #                            v_max_wi*HCwithin  + log10kld3 +
# #                            outcome*HCwithin +
# #                            HCbetween + (1|id/run) + (1|region))  
# 
# # H1 = AIC decreases with each step here
# decode_formula = formula(~ (1|id))
# decode_formula[[1]] = formula(~ (1|id))
# decode_formula[[2]] = formula(~ (1|id/run))
# decode_formula[[3]] = formula(~ (1|id) + (1|region))
# decode_formula[[4]] = formula(~ (1|id/run) + (1|region))
# decode_formula[[5]] = formula(~ v_max_wi + (1|id))
# decode_formula[[6]] = formula(~ v_max_wi + (1|id/run))
# decode_formula[[7]] = formula(~ v_max_wi + (1|id/run) + (1|region))
# decode_formula[[8]] = formula(~ v_max_wi + v_entropy_wi + (1|id))
# decode_formula[[9]] = formula(~ v_max_wi + v_entropy_wi + (1|id/run))
# decode_formula[[10]] = formula(~ v_max_wi + v_entropy_wi + (1|id/run) + (1|region))
# decode_formula[[11]] = formula(~ v_max_wi + v_entropy_wi + outcome + (1|id))
# decode_formula[[12]] = formula(~ v_max_wi + v_entropy_wi + outcome + (1|id/run))
# decode_formula[[13]] = formula(~ v_max_wi + v_entropy_wi + outcome + (1|id/run) + (1|region))
# decode_formula[[14]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + (1|id))
# decode_formula[[15]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + (1|id/run))
# decode_formula[[16]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + (1|id/run) + (1|region))
# decode_formula[[17]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + (1|id))
# decode_formula[[18]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + (1|id/run))
# decode_formula[[19]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + (1|id/run) + (1|region))
# decode_formula[[20]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + (1|id))
# decode_formula[[21]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + (1|id/run))
# decode_formula[[22]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + (1|id/run) + (1|region))
# decode_formula[[23]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + (1|id))
# decode_formula[[24]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + (1|id/run))
# decode_formula[[25]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + (1|id/run) + (1|region))
# decode_formula[[26]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + log10kld3 + (1|id))
# decode_formula[[27]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + log10kld3 + (1|id/run))
# decode_formula[[28]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + log10kld3 + (1|id/run) + (1|region))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# 
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# epoch_label = paste("Time relative to",'feedback', "[s]")
# pal = wes_palette("FantasticFox1", 3, type = "discrete")
# pdf(paste0('AIC-Model-Comparison','.pdf'),width=15,height=12)
# gg <- ggplot(D0, aes(x=evt_time,y=AIC,group=network,color=network)) + facet_wrap(~model) +
#   geom_line(size=1) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) +
#   xlab(epoch_label) + ylab('AIC')
# print(gg)
# dev.off()
# 
# # repeat process without nesting in run
# # H1 = AIC decreases with each step here
# rm(decode_formula)
# decode_formula = formula(~ (1|id))
# decode_formula[[1]] = formula(~ (1|id))
# decode_formula[[2]] = formula(~ (1|id) + (1|region))
# decode_formula[[3]] = formula(~ v_max_wi + (1|id))
# decode_formula[[4]] = formula(~ v_max_wi + (1|id) + (1|region))
# decode_formula[[5]] = formula(~ v_max_wi + v_entropy_wi + (1|id))
# decode_formula[[6]] = formula(~ v_max_wi + v_entropy_wi + (1|id) + (1|region))
# decode_formula[[7]] = formula(~ v_max_wi + v_entropy_wi + outcome + (1|id))
# decode_formula[[8]] = formula(~ v_max_wi + v_entropy_wi + outcome + (1|id) + (1|region))
# decode_formula[[9]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + (1|id))
# decode_formula[[10]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + (1|id) + (1|region))
# decode_formula[[11]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + (1|id))
# decode_formula[[12]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + (1|id) + (1|region))
# decode_formula[[13]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + (1|id))
# decode_formula[[14]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + (1|id) + (1|region))
# decode_formula[[15]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + (1|id))
# decode_formula[[16]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + (1|id) + (1|region))
# decode_formula[[17]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + (1|id/run) + (1|region))
# decode_formula[[18]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + log10kld3 + (1|id))
# decode_formula[[19]] = formula(~ v_max_wi + v_entropy_wi + outcome + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + log10kld3 + (1|id) + (1|region))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-no-nesting-win-run',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-no-nesting-win-run',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# 
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# epoch_label = paste("Time relative to",'feedback', "[s]")
# pal = wes_palette("FantasticFox1", 3, type = "discrete")
# pdf(paste0('AIC-Model-Comparison','.pdf'),width=15,height=12)
# gg <- ggplot(D0, aes(x=evt_time,y=AIC,group=network,color=network)) + facet_wrap(~model) +
#   geom_line(size=1) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) +
#   xlab(epoch_label) + ylab('AIC')
# print(gg)
# dev.off()
# 
# # no by region
# rm(decode_formula)
# decode_formula = formula(~ (1|id))
# decode_formula[[1]] = formula(~ (1|id))
# decode_formula[[2]] = formula(~ v_max_wi + (1|id))
# decode_formula[[3]] = formula(~ v_max_wi + v_entropy_wi + (1|id))
# decode_formula[[4]] = formula(~ v_max_wi + v_entropy_wi + (1|id))
# decode_formula[[5]] = formula(~ v_max_wi + v_entropy_wi + trial_neg_inv_sc + (1|id))
# decode_formula[[6]] = formula(~ v_max_wi + v_entropy_wi + trial_neg_inv_sc + rt_csv_sc + (1|id))
# decode_formula[[7]] = formula(~ v_max_wi + v_entropy_wi + trial_neg_inv_sc + rt_csv_sc + total_earnings + (1|id))
# decode_formula[[8]] = formula(~ v_max_wi + v_entropy_wi + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + (1|id))
# decode_formula[[9]] = formula(~ v_max_wi + v_entropy_wi + trial_neg_inv_sc + rt_csv_sc + total_earnings + v_entropy_wi_change + log10kld3 + (1|id))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-no-region-no-outcome',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-no-region-no-outcome',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# 
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# epoch_label = paste("Time relative to",'feedback', "[s]")
# pal = wes_palette("FantasticFox1", 3, type = "discrete")
# pdf(paste0('AIC-Model-Comparison','.pdf'),width=15,height=12)
# gg <- ggplot(D0, aes(x=evt_time,y=AIC,group=network,color=network)) + facet_wrap(~model) +
#   geom_line(size=1) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) +
#   xlab(epoch_label) + ylab('AIC')
# print(gg)
# dev.off()
# 
# # repeat, move entropy change to H1, remove total_earnings,rt_csv_sc, and trial_neg_inv_sc, we can try variants on trial
# rm(decode_formula)
# decode_formula = formula(~ (1|id))
# decode_formula[[1]] = formula(~ (1|id))
# decode_formula[[2]] = formula(~ v_max_wi + (1|id))
# decode_formula[[3]] = formula(~ v_max_wi + v_entropy_wi + (1|id))
# decode_formula[[4]] = formula(~ v_max_wi + v_entropy_wi + v_entropy_wi_change + (1|id))
# decode_formula[[5]] = formula(~ v_max_wi + v_entropy_wi + v_entropy_wi_change + log10kld3 + (1|id))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-round4-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-round4-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# 
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# epoch_label = paste("Time relative to",'feedback', "[s]")
# pal = wes_palette("FantasticFox1", 3, type = "discrete")
# pdf(paste0('AIC-Model-Comparison','.pdf'),width=15,height=12)
# gg <- ggplot(D0, aes(x=evt_time,y=AIC,group=network,color=network)) + facet_wrap(~model) +
#   geom_line(size=1) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) +
#   xlab(epoch_label) + ylab('AIC')
# print(gg)
# dev.off()
# 
# # test cor of all variables
# dq <- df %>% select(v_max_wi, v_entropy_wi,v_entropy_wi_change,log10kld3)
# cor(dq,method='spearman',use='complete.obs')
# 
# # try some new variables
# source('~/vmPFC/get_trial_data_vmPFC.R')
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          pe_max_sc = scale(pe_max),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# 
# rm(Q)
# 
# message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
# load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
# vmPFC <- fb_comb
# vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
# rm(fb_comb)
# 
# vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
# vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
# vmPFC <- vmPFC %>% mutate(side_vmpfc = as.factor(substr(region, nchar(region), nchar(region))))
# 
# df <- df %>% select(id,run,run_trial,v_max_wi,v_entropy,v_entropy_wi_change,rewFunc,rt_csv_sc,outcome,run_trial,rt_vmax_change,log10kld3,rt_vmax_lag_sc,v_entropy_wi_full,rt_vmax_full)
# Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# rm(vmPFC)
# 
# # repeat, move entropy change to H1, remove total_earnings,rt_csv_sc, and trial_neg_inv_sc, we can try variants on trial
# rm(decode_formula)
# decode_formula = formula(~ (1|id))
# decode_formula[[1]] = formula(~ (1|id))
# decode_formula[[2]] = formula(~ v_max_wi + (1|id))
# decode_formula[[3]] = formula(~ v_max_wi + v_entropy + (1|id))
# decode_formula[[4]] = formula(~ v_max_wi + v_entropy + v_entropy_wi_change + run_trial + (1|id))
# decode_formula[[5]] = formula(~ v_max_wi + v_entropy + v_entropy_wi_change + run_trial + rt_vmax_change + (1|id))
# decode_formula[[6]] = formula(~ v_max_wi + v_entropy + v_entropy_wi_change + run_trial + rt_vmax_change + rt_vmax_lag_sc + (1|id))
# decode_formula[[7]] = formula(~ v_max_wi + v_entropy + v_entropy_wi_change + run_trial + rt_vmax_change + rt_vmax_lag_sc + log10kld3 + (1|id))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-round5-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-round5-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# 
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# epoch_label = paste("Time relative to",'feedback', "[s]")
# pal = wes_palette("FantasticFox1", 3, type = "discrete")
# pdf(paste0('AIC-Model-Comparison','.pdf'),width=15,height=12)
# gg <- ggplot(D0, aes(x=evt_time,y=AIC,group=network,color=network)) + facet_wrap(~model) +
#   geom_line(size=1) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) +
#   xlab(epoch_label) + ylab('AIC')
# print(gg)
# dev.off()
# 
# # add some value related variables
# source('~/vmPFC/get_trial_data_vmPFC.R')
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(id,run,run_trial,v_max_wi,score_csv,v_chosen,ev)
# 
# 
# message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
# load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
# vmPFC <- fb_comb
# vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
# rm(fb_comb)
# 
# vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
# vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
# vmPFC <- vmPFC %>% mutate(side_vmpfc = as.factor(substr(region, nchar(region), nchar(region))))
# 
# Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# rm(vmPFC)
# 
# rm(decode_formula)
# decode_formula = formula(~ (1|id))
# decode_formula[[1]] = formula(~ (1|id))
# decode_formula[[2]] = formula(~ v_max_wi + (1|id))
# decode_formula[[3]] = formula(~ v_max_wi + v_chosen + (1|id))
# decode_formula[[4]] = formula(~ v_max_wi + v_chosen + ev + (1|id))
# decode_formula[[5]] = formula(~ v_max_wi + v_chosen + ev + score_csv + (1|id))
# decode_formula[[6]] = formula(~ (1|id) + (1|region))
# decode_formula[[7]] = formula(~ v_max_wi + (1|id) + (1|region))
# decode_formula[[8]] = formula(~ v_max_wi + v_chosen + (1|id) + (1|region))
# decode_formula[[9]] = formula(~ v_max_wi + v_chosen + ev + (1|id) + (1|region))
# decode_formula[[10]] = formula(~ v_max_wi + v_chosen + ev + score_csv + (1|id) + (1|region))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-value-sensitivity-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-value-sensitivity-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# 
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# epoch_label = paste("Time relative to",'feedback', "[s]")
# pal = wes_palette("FantasticFox1", 3, type = "discrete")
# pdf(paste0('AIC-Model-Comparison','.pdf'),width=15,height=12)
# gg <- ggplot(D0, aes(x=evt_time,y=AIC,group=network,color=network)) + facet_wrap(~model) +
#   geom_line(size=1) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) +
#   xlab(epoch_label) + ylab('AIC')
# print(gg)
# dev.off()
# 
# # add pe signal to model
# # add some value related variables
# rm(Q)
# source('~/vmPFC/get_trial_data_vmPFC.R')
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# 
# message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
# load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
# vmPFC <- fb_comb
# vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
# rm(fb_comb)
# 
# vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
# vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
# vmPFC <- vmPFC %>% mutate(side_vmpfc = as.factor(substr(region, nchar(region), nchar(region))))
# 
# df <- df %>% select(id,run,run_trial,v_max_wi,v_entropy, v_entropy_wi_change,run_trial,rt_vmax_change,log10kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag)
# Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# rm(vmPFC)
# 
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~ v_max_wi + v_entropy + v_entropy_wi_change + run_trial + rt_vmax_change + log10kld3 + (1|id))
# decode_formula[[2]] = formula(~ v_max_wi + v_entropy + v_entropy_wi_change + run_trial + rt_vmax_change + log10kld3 + abs_pe_max_sc + (1|id))
# decode_formula[[3]] = formula(~ v_max_wi + v_entropy + v_entropy_wi_change + run_trial + rt_vmax_change + log10kld3 + pe_max + (1|id))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-wPE-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-wPE-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# 
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# epoch_label = paste("Time relative to",'feedback', "[s]")
# pal = wes_palette("FantasticFox1", 3, type = "discrete")
# pdf(paste0('AIC-Model-Comparison','.pdf'),width=15,height=12)
# gg <- ggplot(D0, aes(x=evt_time,y=AIC,group=network,color=network)) + facet_wrap(~model) +
#   geom_line(size=1) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) +
#   xlab(epoch_label) + ylab('AIC')
# print(gg)
# dev.off()
# 
# # test corr
# dq <- df %>% select(!id & !run & !run_trial)
# cor(dq,method='spearman',use='complete.obs')
# 
# # test rewFunc
# 
# rm(Q)
# 
# source('~/vmPFC/get_trial_data_vmPFC.R')
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# 
# message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
# load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
# vmPFC <- fb_comb
# vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
# rm(fb_comb)
# 
# vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
# vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
# vmPFC <- vmPFC %>% mutate(side_vmpfc = as.factor(substr(region, nchar(region), nchar(region))))
# 
# df <- df %>% select(id,run,run_trial,rewFunc,v_max_wi,v_entropy, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,log10kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag)
# Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# rm(vmPFC)
# 
# Q1 <- Q %>% filter(rewFunc == 'IEV')
# 
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~  (1|id))
# decode_formula[[2]] = formula(~  v_max_wi + (1|id))
# decode_formula[[3]] = formula(~  v_max_wi + v_entropy + (1|id))
# decode_formula[[4]] = formula(~  v_max_wi + v_entropy + v_entropy_wi_change + (1|id))
# decode_formula[[5]] = formula(~  v_max_wi + v_entropy + v_entropy_wi_change + run_trial + (1|id))
# decode_formula[[6]] = formula(~  v_max_wi + v_entropy + v_entropy_wi_change + run_trial + rt_vmax_change + (1|id))
# decode_formula[[7]] = formula(~  v_max_wi + v_entropy + v_entropy_wi_change + run_trial + rt_vmax_change + log10kld3 + (1|id))
# decode_formula[[8]] = formula(~  v_max_wi + v_entropy + v_entropy_wi_change + run_trial + rt_vmax_change + log10kld3 + pe_max + (1|id))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q1, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-wPE-IEV',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-wPE-IEV',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# write.csv(D0,'Model-Comparison-7-IEV.csv')
# 
# 
# Q1 <- Q %>% filter(rewFunc == 'DEV')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q1, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-wPE-DEV',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-wPE-DEV',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# write.csv(D0,'Model-Comparison-7-DEV.csv')
# 
# Q1 <- Q %>% filter(rewFunc == 'CEV')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q1, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-wPE-CEV',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-wPE-CEV',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# write.csv(D0,'Model-Comparison-7-CEV.csv')
# 
# Q1 <- Q %>% filter(rewFunc == 'CEVR')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q1, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-wPE-CEVR',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-wPE-CEVR',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# write.csv(D0,'Model-Comparison-7-CEVR.csv')
# 
# 
# ###############################
# #  test with trial_neg_inv_sc #
# ###############################
# 
# Q1 <- Q %>% filter(rewFunc == 'IEV')
# 
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~  (1|id))
# decode_formula[[2]] = formula(~  v_max_wi + (1|id))
# decode_formula[[3]] = formula(~  v_max_wi + v_entropy_wi_change + (1|id))
# decode_formula[[4]] = formula(~  v_max_wi + v_entropy_wi_change + trial_neg_inv_sc + (1|id))
# decode_formula[[5]] = formula(~  v_max_wi + v_entropy_wi_change + trial_neg_inv_sc + rt_vmax_change + (1|id))
# decode_formula[[6]] = formula(~  v_max_wi + v_entropy_wi_change + trial_neg_inv_sc + rt_vmax_change + log10kld3 + (1|id))
# decode_formula[[7]] = formula(~  v_max_wi + v_entropy_wi_change + trial_neg_inv_sc + rt_vmax_change + log10kld3 + pe_max + (1|id))
# decode_formula[[8]] = formula(~  v_max_wi + v_entropy_wi_change + trial_neg_inv_sc + rt_vmax_change + log10kld3 + pe_max + v_entropy + (1|id))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 4:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q1, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-wPE-IEV-2-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-wPE-IEV-2-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# write.csv(D0,'Model-Comparison-7-IEV-2.csv')
# 
# 
# Q1 <- Q %>% filter(rewFunc == 'DEV')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q1, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-wPE-DEV-2-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-wPE-DEV-2-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# write.csv(D0,'Model-Comparison-7-DEV-2.csv')
# 
# Q1 <- Q %>% filter(rewFunc == 'CEV')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q1, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-wPE-CEV-2-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-wPE-CEV-2-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# write.csv(D0,'Model-Comparison-7-CEV-2.csv')
# 
# Q1 <- Q %>% filter(rewFunc == 'CEVR')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q1, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-wPE-CEVR-2-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-wPE-CEVR-2-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# write.csv(D0,'Model-Comparison-7-CEVR-2.csv')
# 
# 
# 
# 
# 
# # rewFunc alone has a big impact, but once you introduce other variables it's impact diminishes
# 
# 
# # test choose sooner/later/stay
# source('~/vmPFC/get_trial_data_vmPFC.R')
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# df <- df %>% group_by(id,run) %>% mutate(expl_code =(case_when(
#   rt_csv - rt_lag > 1 ~ 'Longer',
#   rt_csv - rt_lag < -1 ~ 'Shorter',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ 'Same'
# )))
# df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
#   rt_csv_sc <= -1 ~ '-1',
#   rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
#   rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
#   rt_csv_sc > 1 ~ '1'
# )))
# df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 10 ~ 'Early',
#   run_trial > 10 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_code,rt_csv_sc,v_entropy, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,log10kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag)
# Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# Q$expl_code <- relevel(as.factor(Q$expl_code),ref='Same')
# Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
# Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~  (1|id))
# decode_formula[[2]] = formula(~  v_max_wi + (1|id))
# decode_formula[[3]] = formula(~  v_max_wi + v_entropy_wi_change + (1|id))
# decode_formula[[4]] = formula(~  v_max_wi + v_entropy_wi_change + trial_bin + (1|id))
# decode_formula[[5]] = formula(~  v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + (1|id))
# decode_formula[[6]] = formula(~  v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + log10kld3 + (1|id))
# decode_formula[[7]] = formula(~  v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + log10kld3 + pe_max + (1|id))
# decode_formula[[8]] = formula(~  v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + log10kld3 + pe_max + v_entropy + (1|id))
# decode_formula[[9]] = formula(~  v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + log10kld3 + pe_max + v_entropy + rt_bin + (1|id))
# decode_formula[[10]] = formula(~  v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + log10kld3 + pe_max + v_entropy + rt_bin + expl_code + (1|id))
# splits <- c('evt_time','network')
# 
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-w_expl_code-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-w_expl_code-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_BIC <- paste0('BIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,BIC) %>% rename(!! new_BIC := BIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_BIC <- paste0('BIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,BIC) %>% rename(!! new_BIC := BIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('BIC'),names_to='model',names_prefix="BIC-",values_to="BIC")
# write.csv(D0,'last-model-w_expl_code.csv')
# 
# rm(Q)
# Q <- inner_join(vmPFC,df,by=c('id','run','run_trial'))
# # test age & sex
# demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
# demo <- demo %>% rename(id=lunaid)
# demo <- demo %>% select(!adult & !scandate)
# Q <- inner_join(Q,demo,by=c('id'))
# Q <- Q %>% filter(Q$kld3 < mean(Q$kld3,na.rm=TRUE)+2.5*sd(Q$kld3,na.rm=TRUE))
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~  age + (1|id))
# decode_formula[[2]] = formula(~  age + female + (1|id))
# decode_formula[[3]] = formula(~  age + female + v_max_wi + (1|id))
# decode_formula[[4]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + (1|id))
# decode_formula[[5]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + (1|id))
# decode_formula[[6]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + (1|id))
# decode_formula[[7]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + (1|id))
# decode_formula[[8]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + pe_max + (1|id))
# decode_formula[[9]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + pe_max + v_entropy + (1|id))
# decode_formula[[10]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + pe_max + v_entropy + rt_bin + (1|id))
# decode_formula[[11]] = formula(~ age + female +  v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + pe_max + v_entropy + rt_bin + expl_code + (1|id))
# splits <- c('evt_time','network')
# 
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-last-model-w_expl_code&demo-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-last-model-w_expl_code&demo-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_BIC <- paste0('BIC','-',i)
#     D0 <- D0 %>% select(evt_time,network,BIC) %>% rename(!! new_BIC := BIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_BIC <- paste0('BIC','-',i)
#     D1 <- D1 %>% select(evt_time,network,BIC) %>% rename(!! new_BIC := BIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('BIC'),names_to='model',names_prefix="BIC-",values_to="BIC")
# write.csv(D0,'last-model-w_expl_code&demo.csv')
# 
# 
# 
# #################
# ### HC model  ###
# #################
# rm(Q)
# setwd('~/vmPFC')
# message('adding HC signals to models...')
# load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
# hc <- fb_comb
# hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
# rm(fb_comb)
# 
# hc <- hc %>% select(id,run,run_trial,decon_mean,evt_time,bin_num,side)
# hc <- hc %>% mutate(
#   HC_region = case_when(
#     bin_num <= 8 ~ 'AH',
#     bin_num >8 ~ 'PH'
#   ))
# hc <- hc %>% group_by(id,run,run_trial,evt_time,bin_num) %>% summarize(decon_mean1 = mean(decon_mean,na.rm=TRUE)) 
# hc <- hc %>% rename(decon_mean=decon_mean1)
# # hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>%
# #   summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
# # hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
# # hc <- hc %>% select(!decon1)
# 
# splits = c("evt_time","bin_num")
# 
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# df <- df %>% group_by(id,run) %>% mutate(expl_code =(case_when(
#   rt_csv - rt_lag > 1 ~ 'Longer',
#   rt_csv - rt_lag < -1 ~ 'Shorter',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ 'Same'
# )))
# df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
#   rt_csv_sc <= -1 ~ '-1',
#   rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
#   rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
#   rt_csv_sc > 1 ~ '1'
# )))
# df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 10 ~ 'Early',
#   run_trial > 10 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_code,rt_csv_sc,v_entropy, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag)
# Q <- merge(df, hc, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# Q$expl_code <- relevel(as.factor(Q$expl_code),ref='Same')
# Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
# Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# # test age & sex
# demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
# demo <- demo %>% rename(id=lunaid)
# demo <- demo %>% select(!adult & !scandate)
# Q <- inner_join(Q,demo,by=c('id'))
# Q <- Q %>% filter(Q$kld3 < mean(Q$kld3,na.rm=TRUE)+2.5*sd(Q$kld3,na.rm=TRUE))
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~  age + (1|id))
# decode_formula[[2]] = formula(~  age + female + (1|id))
# decode_formula[[3]] = formula(~  age + female + v_max_wi + (1|id))
# decode_formula[[4]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + (1|id))
# decode_formula[[5]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + (1|id))
# decode_formula[[6]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + (1|id))
# decode_formula[[7]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + (1|id))
# decode_formula[[8]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + pe_max + (1|id))
# decode_formula[[9]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + pe_max + v_entropy + (1|id))
# decode_formula[[10]] = formula(~  age + female + v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + pe_max + v_entropy + rt_bin + (1|id))
# decode_formula[[11]] = formula(~ age + female +  v_max_wi + v_entropy_wi_change + trial_bin + rt_vmax_change + pe_max + v_entropy + rt_bin + expl_code + (1|id))
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1:length(decode_formula)){
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-HC-last-model-w_expl_code&demo-',i,'.Rdata'))
# }
# 
# # plot
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-HC-last-model-w_expl_code&demo-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_BIC <- paste0('BIC','-',i)
#     D0 <- D0 %>% select(evt_time,bin_num,BIC) %>% rename(!! new_BIC := BIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_BIC <- paste0('BIC','-',i)
#     D1 <- D1 %>% select(evt_time,bin_num,BIC) %>% rename(!! new_BIC := BIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','bin_num'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('BIC'),names_to='model',names_prefix="BIC-",values_to="BIC")
# write.csv(D0,'HC-last-model-w_expl_code&demo.csv')
# 
# 
# 
# 
# ####################
# ### vmPFC - HC   ###
# ####################
# rm(Q)
# message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
# load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
# vmPFC <- fb_comb
# vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
# rm(fb_comb)
# vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
# vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
# load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
# hc <- fb_comb
# hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
# rm(fb_comb)
# hc <- hc %>% mutate(
#   HC_region = case_when(
#     bin_num <= 8 ~ 'AH',
#     bin_num >8 ~ 'PH'
#   ),
# )
# hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
# hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
# 
# Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
# Q <- Q %>% select(!decon1)
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          v_entropy_sc = scale(v_entropy),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
#   rt_csv - rt_lag > 1 ~ 'Longer',
#   rt_csv - rt_lag < -1 ~ '0',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
#   rt_csv - rt_lag > 1 ~ '0',
#   rt_csv - rt_lag < -1 ~ 'Shorter',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
#   rt_csv_sc <= -1 ~ '-1',
#   rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
#   rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
#   rt_csv_sc > 1 ~ '1'
# )))
# df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 10 ~ 'Early',
#   run_trial > 10 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_sc, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag)
# Q <- merge(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
# Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
# Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
# Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# # test age & sex
# demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
# demo <- demo %>% rename(id=lunaid)
# demo <- demo %>% select(!adult & !scandate)
# Q <- inner_join(Q,demo,by=c('id'))
# Q$female <- relevel(as.factor(Q$female),ref='0')
# Q$age <- scale(Q$age)
#   
#   
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
#                                 v_entropy_wi_change + trial_neg_inv_sc + 
#                                 pe_max*HCwithin + v_entropy_sc*HCwithin + rt_bin + 
#                                 expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
#                                 v_entropy_wi_change:expl_longer +
#                                 v_entropy_wi_change:expl_shorter +
#                                 HCbetween +
#                                 (1|id))
# 
# splits = c('evt_time','network','HC_region')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1){
#   setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                   tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-vmPFC-HC-full-network-',i,'.Rdata'))
# }
# 
# splits = c('evt_time','symmetry_group','HC_region')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1){
#   setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                   tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-vmPFC-HC-full-symmetry-',i,'.Rdata'))
# }
# 
# if (exists('D0')){ # start fresh
#   rm(D0)
# }
# for (i in 1:length(decode_formula)){
#   curr_model <- paste0('-vmPFC-HC-basic-',i,'.Rdata')
#   curr_model <- Sys.glob(paste0('*',curr_model))
#   load(curr_model)
#   if (i==1){
#     D0 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D0 <- D0 %>% select(evt_time,HC_region,network,AIC) %>% rename(!! new_AIC := AIC)
#   } else {
#     D1 <- ddf$fit_df
#     new_AIC <- paste0('AIC','-',i)
#     D1 <- D1 %>% select(evt_time,HC_region,network,AIC) %>% rename(!! new_AIC := AIC)
#     D0 <- inner_join(D0,D1,by=c('evt_time','HC_region','network'))
#     rm(D1)
#   }
#   rm(ddf)
#   disp(i)
# }
# D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC-",values_to="AIC")
# write.csv(D0,'vmPFC-HC-basic.csv')
# 
# rm(Q)
# message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
# load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
# vmPFC <- fb_comb
# vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
# rm(fb_comb)
# vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
# vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
# 
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          v_entropy_sc = scale(v_entropy),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
#   rt_csv - rt_lag > 1 ~ 'Longer',
#   rt_csv - rt_lag < -1 ~ '0',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
#   rt_csv - rt_lag > 1 ~ '0',
#   rt_csv - rt_lag < -1 ~ 'Shorter',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
#   rt_csv_sc <= -1 ~ '-1',
#   rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
#   rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
#   rt_csv_sc > 1 ~ '1'
# )))
# df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 10 ~ 'Early',
#   run_trial > 10 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_sc, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag)
# Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
# Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
# Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
# Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# # test age & sex
# demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
# demo <- demo %>% rename(id=lunaid)
# demo <- demo %>% select(!adult & !scandate)
# Q <- inner_join(Q,demo,by=c('id'))
# Q$female <- relevel(as.factor(Q$female),ref='0')
# Q$age <- scale(Q$age)
# 
# 
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~ age + female +  v_max_wi + 
#                                 v_entropy_wi_change + trial_neg_inv_sc + 
#                                 pe_max + v_entropy_sc + rt_bin + 
#                                 expl_shorter + expl_longer +  # binary expl_code incr / decr separate variables
#                                 v_entropy_wi_change:expl_longer +
#                                 v_entropy_wi_change:expl_shorter +
#                                 (1|id))
# 
# splits = c('evt_time','network')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1){
#   setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                   tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-vmPFC-network-',i,'.Rdata'))
# }
# 
# 
# 
# 
# rm(Q)
# setwd('~/vmPFC')
# message('adding HC signals to models...')
# load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
# hc <- clock_comb
# hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
# rm(clock_comb)
# 
# hc <- hc %>% select(id,run,run_trial,decon_mean,evt_time,bin_num,side)
# hc <- hc %>% mutate(
#   HC_region = case_when(
#     bin_num <= 8 ~ 'AH',
#     bin_num >8 ~ 'PH'
#   ))
# hc <- hc %>% group_by(id,run,run_trial,evt_time,bin_num) %>% summarize(decon_mean1 = mean(decon_mean,na.rm=TRUE)) 
# hc <- hc %>% rename(decon_mean=decon_mean1)
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          pe_max_sc = scale(pe_max),
#          pe_max_lag_sc = scale(lag(pe_max)),
#          rt_vmax_sc = scale(rt_vmax),
#          v_entropy_sc = scale(v_entropy),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
#   rt_csv - rt_lag > 1 ~ 'Longer',
#   rt_csv - rt_lag < -1 ~ '0',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
#   rt_csv - rt_lag > 1 ~ '0',
#   rt_csv - rt_lag < -1 ~ 'Shorter',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
#   rt_csv_sc <= -1 ~ '-1',
#   rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
#   rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
#   rt_csv_sc > 1 ~ '1'
# )))
# df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 10 ~ 'Early',
#   run_trial > 10 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_sc, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max_sc,pe_max_lag_sc)
# Q <- merge(df, hc, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
# Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
# Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
# Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# # test age & sex
# demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
# demo <- demo %>% rename(id=lunaid)
# demo <- demo %>% select(!adult & !scandate)
# Q <- inner_join(Q,demo,by=c('id'))
# Q$female <- relevel(as.factor(Q$female),ref='0')
# Q$age <- scale(Q$age)
# 
# 
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~ age + female +  v_max_wi + 
#                                 v_entropy_wi_change + trial_neg_inv_sc + 
#                                 pe_max_lag_sc + v_entropy_sc + rt_bin + 
#                                 expl_shorter + expl_longer +  # binary expl_code incr / decr separate variables
#                                 v_entropy_wi_change:expl_longer +
#                                 v_entropy_wi_change:expl_shorter +
#                                 (1|id))
# splits = c('evt_time','bin_num')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1){
#   setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                   tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-HC-axis-clock-',i,'.Rdata'))
# }
# 
# 
# rm(Q)
# message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
# load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
# vmPFC <- clock_comb
# vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
# rm(clock_comb)
# vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
# vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
# 
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          v_entropy_sc = scale(v_entropy),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
#   rt_csv - rt_lag > 1 ~ 'Longer',
#   rt_csv - rt_lag < -1 ~ '0',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
#   rt_csv - rt_lag > 1 ~ '0',
#   rt_csv - rt_lag < -1 ~ 'Shorter',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
#   rt_csv_sc <= -1 ~ '-1',
#   rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
#   rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
#   rt_csv_sc > 1 ~ '1'
# )))
# df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 10 ~ 'Early',
#   run_trial > 10 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_sc, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag)
# Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
# Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
# Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
# Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# # test age & sex
# demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
# demo <- demo %>% rename(id=lunaid)
# demo <- demo %>% select(!adult & !scandate)
# Q <- inner_join(Q,demo,by=c('id'))
# Q$female <- relevel(as.factor(Q$female),ref='0')
# Q$age <- scale(Q$age)
# 
# 
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~ age + female +  v_max_wi + 
#                                 v_entropy_wi_change + trial_neg_inv_sc + 
#                                 pe_max + v_entropy_sc + rt_bin + 
#                                 expl_shorter + expl_longer +  # binary expl_code incr / decr separate variables
#                                 v_entropy_wi_change:expl_longer +
#                                 v_entropy_wi_change:expl_shorter +
#                                 (1|id))
# 
# splits = c('evt_time','network')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1){
#   setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                   tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-vmPFC-network-clock-',i,'.Rdata'))
# }
# 
# 
# ####################
# ### vmPFC - HC   ###
# ####################
# rm(Q)
# message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
# load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
# vmPFC <- clock_comb
# vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
# rm(clock_comb)
# vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
# vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
# load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
# hc <- clock_comb
# hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
# rm(clock_comb)
# hc <- hc %>% mutate(
#   HC_region = case_when(
#     bin_num <= 8 ~ 'AH',
#     bin_num >8 ~ 'PH'
#   ),
# )
# hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
# hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
# 
# Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
# Q <- Q %>% select(!decon1)
# df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
# df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
#                     v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
#                     rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
#   group_by(id, run) %>% 
#   mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
#   mutate(v_chosen_sc = scale(v_chosen),
#          abs_pe_max_sc = scale(abs(pe_max)),
#          abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
#          rt_vmax_sc = scale(rt_vmax),
#          v_entropy_sc = scale(v_entropy),
#          rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
#            kld3 ==0 ~ NA_real_,
#            kld3 >0 ~ log10(kld3)
#          )) %>% mutate(log10kld3_lag = case_when(
#            kld3_lag==0 ~NA_real_,
#            kld3_lag>0 ~ log10(kld3_lag)
#          ))
# df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
#   rt_csv - rt_lag > 1 ~ 'Longer',
#   rt_csv - rt_lag < -1 ~ '0',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
#   rt_csv - rt_lag > 1 ~ '0',
#   rt_csv - rt_lag < -1 ~ 'Shorter',
#   rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
# )))
# df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
#   rt_csv_sc <= -1 ~ '-1',
#   rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
#   rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
#   rt_csv_sc > 1 ~ '1'
# )))
# df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 10 ~ 'Early',
#   run_trial > 10 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_sc, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag)
# Q <- merge(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
# Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
# Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
# Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
# Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# # test age & sex
# demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
# demo <- demo %>% rename(id=lunaid)
# demo <- demo %>% select(!adult & !scandate)
# Q <- inner_join(Q,demo,by=c('id'))
# Q$female <- relevel(as.factor(Q$female),ref='0')
# Q$age <- scale(Q$age)
# 
# 
# rm(decode_formula)
# decode_formula <- formula(~ (1|id))
# decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
#                                 v_entropy_wi_change + trial_neg_inv_sc + 
#                                 pe_max*HCwithin + v_entropy_sc*HCwithin + rt_bin + 
#                                 expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
#                                 v_entropy_wi_change:expl_longer +
#                                 v_entropy_wi_change:expl_shorter +
#                                 HCbetween +
#                                 (1|id))
# 
# splits = c('evt_time','network','HC_region')
# source("~/fmri.pipeline/R/mixed_by.R")
# for (i in 1){
#   setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   df0 <- decode_formula[[i]]
#   print(df0)
#   ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                   tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddf,file=paste0(curr_date,'-vmPFC-HC-full-network-clock-',i,'.Rdata'))
# }

###################################
#####      vmPFC -feedback    #####
###################################

rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- fb_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(fb_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
         ))
df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)


rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ age + female +  v_max_wi + 
                                v_entropy_wi_change + trial_neg_inv_sc + 
                                pe_max_sc + v_entropy_wi + rt_bin + 
                                expl_shorter + expl_longer +  # binary expl_code incr / decr separate variables
                                v_entropy_wi_change:expl_longer +
                                v_entropy_wi_change:expl_shorter +
                                (1|id))

splits = c('evt_time','symmetry_group')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-symmetry-feedback-',i,'.Rdata'))
}
splits = c('evt_time','network')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qPE = quantile(df$pe_max_sc,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emmeans_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', 
                                specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(qH))),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', 
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1',
                                specs=c('v_max_wi'), at=list(v_max_wi=c(qV))),
                    PE_HC = list(outcome='vmPFC_decon',model_name='model1',
                                 specs=c("pe_max_sc"), at = list(pe_max_sc=c(qPE))),
                    dH_HC = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=c('v_entropy_wi_change'), at=list(v_entropy_wi_change=c(qdH)))
                                )
                  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-network-feedback-',i,'.Rdata'))
}

###################################
#####      vmPFC - clock      #####
###################################

rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
         ))
df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)


rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ v_max_wi + 
                                trial_neg_inv_sc + 
                                v_entropy_wi + rt_bin + 
                                expl_shorter + expl_longer +  # binary expl_code incr / decr separate variables
                                (1|id))

splits = c('evt_time','symmetry_group')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-symmetry-clock-',i,'.Rdata'))
}

rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ age + female +  v_max_wi + 
                                v_entropy_wi_change_lag + trial_neg_inv_sc + 
                                v_entropy_wi + rt_bin + 
                                expl_shorter + expl_longer +  # binary expl_code incr / decr separate variables
                                v_entropy_wi_change_lag:expl_longer +
                                v_entropy_wi_change_lag:expl_shorter +
                                (1|id))

splits = c('evt_time','network')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qPE = quantile(df$pe_max_sc,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emmeans_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', 
                                specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(qH))),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', 
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1',
                                specs=c('v_max_wi'), at=list(v_max_wi=c(qV))),
                    dH_HC = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=c('v_entropy_wi_change_lag'), at=list(v_entropy_wi_change_lag=c(qdH)))
                  ))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-network-clock-',i,'.Rdata'))
}

#################################
#####   HC - feedback      ######
#################################
rm(Q)
setwd('~/vmPFC')
message('adding HC signals to models...')
load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
hc <- fb_comb
hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
rm(fb_comb)

hc <- hc %>% select(id,run,run_trial,decon_mean,evt_time,bin_num,side)
hc <- hc %>% mutate(
  HC_region = case_when(
    bin_num <= 8 ~ 'AH',
    bin_num >8 ~ 'PH'
  ))
hc <- hc %>% group_by(id,run,run_trial,evt_time,bin_num) %>% summarize(decon_mean1 = mean(decon_mean,na.rm=TRUE)) 
hc <- hc %>% rename(decon_mean=decon_mean1)
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
         ))
df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, hc, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)


rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ age + female +  v_max_wi + 
                                v_entropy_wi_change + trial_neg_inv_sc + 
                                pe_max_sc + v_entropy_wi + rt_bin + 
                                expl_shorter + expl_longer +  # binary expl_code incr / decr separate variables
                                v_entropy_wi_change:expl_longer +
                                v_entropy_wi_change:expl_shorter +
                                (1|id))
splits = c('evt_time','bin_num')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-HC-axis-feedback-',i,'.Rdata'))
}
###################################
#####      HC- clock          #####
###################################
rm(Q)
setwd('~/vmPFC')
message('adding HC signals to models...')
load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
hc <- clock_comb
hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)

hc <- hc %>% select(id,run,run_trial,decon_mean,evt_time,bin_num,side)
hc <- hc %>% mutate(
  HC_region = case_when(
    bin_num <= 8 ~ 'AH',
    bin_num >8 ~ 'PH'
  ))
hc <- hc %>% group_by(id,run,run_trial,evt_time,bin_num) %>% summarize(decon_mean1 = mean(decon_mean,na.rm=TRUE)) 
hc <- hc %>% rename(decon_mean=decon_mean1)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
         ))
df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, hc, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)


rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ age + female +  v_max_wi + 
                                v_entropy_wi_change_lag + trial_neg_inv_sc + 
                                v_entropy_wi + rt_bin + 
                                expl_shorter + expl_longer +  # binary expl_code incr / decr separate variables
                                v_entropy_wi_change_lag:expl_longer +
                                v_entropy_wi_change_lag:expl_shorter +
                                (1|id))
splits = c('evt_time','bin_num')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
                  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-HC-axis-clock-',i,'.Rdata'))
}
###################################
##### vmPFC - HC -feedback    #####
###################################
rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- fb_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(fb_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
hc <- fb_comb
hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
rm(fb_comb)
hc <- hc %>% mutate(
  HC_region = case_when(
    bin_num <= 8 ~ 'AH',
    bin_num >8 ~ 'PH'
  ),
)
hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()

Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
Q <- Q %>% select(!decon1)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(trial,v_max,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
         ))
df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max_sc,pe_max_lag_sc)
Q <- merge(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)



rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                v_entropy_wi_change*HCwithin + trial_neg_inv_sc*HCwithin + 
                                pe_max_sc*HCwithin + v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (1|id))

splits = c('evt_time','network','HC_region')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qPE = quantile(df$pe_max_sc,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emtrends_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                         specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(qH))),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                         specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=c('v_max_wi'), at=list(v_max_wi=c(qV))),
                    PE_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                 specs=c('pe_max_sc'), at=list(pe_max_sc=c(qPE))),
                    dH_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                 specs=c('v_entropy_wi_change'), at=list(v_entropy_wi_change=c(qdH)))
                    )
                  )
  
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-feedback-',i,'.Rdata'))
}

splits = c('evt_time','symmetry_group','HC_region')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qPE = quantile(df$pe_max_sc,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <-  mean(df$v_entropy_wi_change,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  qdH[2] <-  mean(df$v_entropy_wi_change,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-feedback-',i,'.Rdata'))
}
###################################
##### vmPFC - HC - clock      #####
###################################
rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
hc <- clock_comb
hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)
hc <- hc %>% mutate(
  HC_region = case_when(
    bin_num <= 8 ~ 'AH',
    bin_num >8 ~ 'PH'
  ),
)
hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()

Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
Q <- Q %>% select(!decon1)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
         ))
df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)


dP <- read_csv(file=file.path('~/vmPFC/','peaks.csv'),col_names = FALSE)
dP <- dP %>% rename(peaks=X1,id=X2,run=X3,run_trial=X4)
Q <- inner_join(Q,dP,by=c('id','run','run_trial'))
Q$peaks <- as.numeric(peaks)
Q <- Q %>% mutate(peaks2 = case_when(
  peaks==0 ~ NA_real_,
  peaks==1 ~ peaks,
  peaks==2 ~ peaks,
  peaks==3 ~ peaks,
  peaks==4 ~ peaks,
  peaks==5 ~ NA_real_,
  TRUE ~ as.numeric(peaks)
))
Q <- Q %>% select(!peaks)
Q$peaks2 <- relevel(as.factor(Q$peaks2),ref=1)

rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (1|id))


qHC <- Q %>% filter(atlas_value==55) %>% group_by(evt_time,HC_region) %>% 
  summarize(HC_p2SD = mean(HCwithin,na.rm=TRUE)+2*sd(HCwithin,na.rm=TRUE),
            HC_m2SD = mean(HCwithin,na.rm=TRUE)-2*sd(HCwithin,na.rm=TRUE),
            HC_p1SD = mean(HCwithin,na.rm=TRUE)+1*sd(HCwithin,na.rm=TRUE),
            HC_m1SD = mean(HCwithin,na.rm=TRUE)-1*sd(HCwithin,na.rm=TRUE),
            HC_p05SD = mean(HCwithin,na.rm=TRUE)+0.5*sd(HCwithin,na.rm=TRUE),
            HC_m05SD = mean(HCwithin,na.rm=TRUE)-0.5*sd(HCwithin,na.rm=TRUE))

qHC <- qHC %>% filter(HC_region=='PH')
qHC1 <- NULL
qHC1[1] <- mean(qHC$HC_m2SD,na.rm=TRUE)
qHC1[2] <- mean(qHC$HC_p2SD,na.rm=TRUE)
qHC1[3] <- mean(qHC$HC_m1SD,na.rm=TRUE)
qHC1[4] <- mean(qHC$HC_p1SD,na.rm=TRUE)
qHC1[5] <- mean(qHC$HC_m05SD,na.rm=TRUE)
qHC1[6] <- mean(qHC$HC_p05SD,na.rm=TRUE)
splits = c('evt_time','network','HC_region')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emtrends_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs='v_entropy_wi', at = list(v_entropy_wi=qH)),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=c('v_max_wi'), at=list(v_max_wi=c(qV)))
                  ), 
                  emmeans_spec = list(
                    H_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ v_entropy_wi * HCwithin), at=list(HCwithin=qHC1,v_entropy_wi=qH)),
                    T_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ trial_neg_inv_sc * HCwithin), at=list(HCwithin=qHC1,trial_neg_inv_sc=qT)),
                    V_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ v_max_wi * HCwithin), at=list(HCwithin=qHC1,v_max_wi=qV))
                  )
                  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-',i,'.Rdata'))
}
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emtrends_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(qH))),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=c('v_max_wi'), at=list(v_max_wi=c(qV)))
                  )
                  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-clock-',i,'.Rdata'))
}



##############################
### testing ...   ############
##############################
rm(Q)
Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
Q <- Q %>% select(!decon1)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(probability,magnitude,ev,v_chosen,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_chosen_sc = scale(v_chosen),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
         ))
df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,run_trial,v_chosen_sc,trial_neg_inv_sc,v_entropy_wi,v_entropy_wi_change_lag,expl_longer,expl_shorter,rt_bin,v_max_wi)
Q <- merge(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)


dP <- read_csv(file=file.path('~/vmPFC/','peaks.csv'),col_names = FALSE)
dP <- dP %>% rename(peaks=X1,id=X2,run=X3,run_trial=X4)
Q <- inner_join(Q,dP,by=c('id','run','run_trial'))
Q$peaks <- as.numeric(peaks)
Q <- Q %>% mutate(peaks2 = case_when(
  peaks==0 ~ NA_real_,
  peaks==1 ~ peaks,
  peaks==2 ~ peaks,
  peaks==3 ~ peaks,
  peaks==4 ~ peaks,
  peaks==5 ~ NA_real_,
  TRUE ~ as.numeric(peaks)
))
Q <- Q %>% select(!peaks)
Q$peaks <- relevel(as.factor(Q$peaks),ref=1)
Q <- Q %>% filter(!is.na(peaks))
Q$HCbetween <- scale(Q$HCbetween)

rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                rt_bin + v_entropy_wi*HCwithin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (1|id))
decode_formula[[2]] = formula(~ v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                peaks*HCwithin + rt_bin + v_entropy_wi*HCwithin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (1|id))
decode_formula[[3]] = formula(~ v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                peaks*HCwithin + rt_bin + v_entropy_wi*HCwithin + peaks:HCwithin:v_entropy_wi + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (1|id))

qHC <- Q %>% filter(atlas_value==55) %>% group_by(evt_time,HC_region) %>% 
  summarize(HC_p2SD = mean(HCwithin,na.rm=TRUE)+2*sd(HCwithin,na.rm=TRUE),
            HC_m2SD = mean(HCwithin,na.rm=TRUE)-2*sd(HCwithin,na.rm=TRUE),
            HC_p1SD = mean(HCwithin,na.rm=TRUE)+1*sd(HCwithin,na.rm=TRUE),
            HC_m1SD = mean(HCwithin,na.rm=TRUE)-1*sd(HCwithin,na.rm=TRUE),
            HC_p05SD = mean(HCwithin,na.rm=TRUE)+0.5*sd(HCwithin,na.rm=TRUE),
            HC_m05SD = mean(HCwithin,na.rm=TRUE)-0.5*sd(HCwithin,na.rm=TRUE))

qHC <- qHC %>% filter(HC_region=='PH')
qHC1 <- NULL
qHC1[1] <- mean(qHC$HC_m2SD,na.rm=TRUE)
qHC1[2] <- mean(qHC$HC_p2SD,na.rm=TRUE)
qHC1[3] <- mean(qHC$HC_m1SD,na.rm=TRUE)
qHC1[4] <- mean(qHC$HC_p1SD,na.rm=TRUE)
qHC1[5] <- mean(qHC$HC_m05SD,na.rm=TRUE)
qHC1[6] <- mean(qHC$HC_p05SD,na.rm=TRUE)

splits = c('evt_time','network','HC_region')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 3){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  #qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  #qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH <- c(1,2,3,4)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emtrends_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=formula(~peaks:HCwithin)),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=(~trial_neg_inv_sc:HCwithin), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=(~v_max_wi:HCwithin), at=list(v_max_wi=c(qV)))
                  ),
                  emmeans_spec = list(
                    H_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ peaks * HCwithin), at=list(HCwithin=qHC1)),
                    T_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ trial_neg_inv_sc * HCwithin), at=list(HCwithin=qHC1,trial_neg_inv_sc=qT)),
                    V_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ v_max_wi * HCwithin), at=list(HCwithin=qHC1,v_max_wi=qV))
                  )
  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-testing-wem',i,'.Rdata'))
}

