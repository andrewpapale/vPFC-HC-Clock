# 2022-03-31 AndyP
# mixed_by model building for AIC comparison

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
doTesting = FALSE
do_vPFC_fb = TRUE
do_vPFC_clock = TRUE
do_symmetry = TRUE
do_network = TRUE
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26
source("~/fmri.pipeline/R/mixed_by.R")



###################################
#####      vmPFC -feedback    #####
###################################
if (do_vPFC_fb){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
  vmPFC <- fb_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(fb_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  source('~/vmPFC/get_trial_data_vmPFC.R')
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
  df <- df %>% 
    group_by(id, run) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           score_sc = scale(score_csv),
           iti_sc = scale(iti_ideal),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           rt_vmax_sc = scale(rt_vmax),
           ev_sc = scale(ev),
           v_entropy_sc = scale(v_entropy),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
           abs_rt_vmax_change = abs(rt_vmax_change),
           rt_vmax_change_bin = case_when(
             abs_rt_vmax_change < 4/24 ~ 'No Change',
             abs_rt_vmax_change >= 4/24 ~ 'Change'
           ),
           v_entropy_wi_change_bin = case_when(
             v_entropy_wi_change < -0.5 ~ 'Decrease',
             v_entropy_wi_change > 0.5  ~ 'Increase',
             v_entropy_wi_change >= -0.5 & v_entropy_wi_change <= 0.5 ~ 'No Change'
           ),
           rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld4 = case_when(
             kld4 ==0 ~ NA_real_,
             kld4 >0 ~ log10(kld4)
           )) %>% mutate(log10kld4_lag = case_when(
             kld4_lag==0 ~NA_real_,
             kld4_lag>0 ~ log10(kld4_lag)
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
    run_trial <= 15 ~ 'Early',
    run_trial > 15 & run_trial < 30 ~ 'Middle',
    run_trial >=30 ~ 'Late',
  )))
  #df <- df %>% filter(!is.na(rt_vmax_change_bin) | !is.na(v_entropy_wi_change_bin))
  df <- df %>% select(id,run,run_trial,rt_vmax_change_sc,v_entropy_wi_change,trial_bin,iti_ideal,iti_prev,rt_csv,trial_neg_inv_sc,rt_csv_sc,rewFunc,v_entropy_sc,outcome,v_max_wi,score_sc,rt_bin,iti_sc,ev_sc,expl_longer,expl_shorter)
  Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$vmPFC_decon[Q$evt_time > Q$iti_ideal] = NA;
  #Q$vmPFC_decon[Q$evt_time < -(Q$rt_csv)] = NA;
  Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
  Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
  Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
  Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
  #Q$v_entropy_wi_change_bin <- relevel(as.factor(Q$v_entropy_wi_change_bin),ref='No Change')
  #Q$rt_vmax_change_bin <- relevel(as.factor(Q$rt_vmax_change_bin),ref='No Change')
  
  # test age & sex
  demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
  demo <- demo %>% rename(id=lunaid)
  demo <- demo %>% select(!adult & !scandate)
  Q <- inner_join(Q,demo,by=c('id'))
  Q$female <- relevel(as.factor(Q$female),ref='0')
  Q$age <- scale(Q$age)
  
  rm(decode_formula)
  decode_formula <- formula(~ (1|id))
  decode_formula[[1]] = formula(~ age + female + v_entropy_sc + v_entropy_wi_change + outcome + trial_neg_inv_sc + v_max_wi + rt_csv_sc + iti_sc + rt_vmax_change_sc + (1|id/run))
  decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_entropy_wi_change + outcome + trial_neg_inv_sc + v_max_wi + rt_csv_sc + iti_sc + rt_vmax_change_sc + (1 + v_entropy_sc |id/run))
  #decode_formula[[2]] = formula(~ v_entropy_sc*last_outcome + v_entropy_sc*trial_neg_inv_sc + v_max_wi*last_outcome + v_entropy_wi_change + rt_csv_sc + iti_sc + (1|id/run))
  #decode_formula[[2]] = formula(~ v_entropy_sc*last_outcome + v_entropy_sc*trial_neg_inv_sc + v_max_wi*last_outcome + v_entropy_wi_change + rt_csv_sc + iti_sc +  (1+v_max_wi + v_entropy_sc |id/run))
  qT <- c(-0.7,0.43)
  if (do_symmetry){
    splits = c('evt_time','symmetry_group')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', 
                                 specs=c("v_entropy_sc"), at = list(v_entropy_sc=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                        V = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5))),
                        dH = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=c("v_entropy_wi_change"), at=list(v_entropy_wi_change=c(-0.5,0.5)))
                      ),
                      emtrends_spec = list(
                        Tr = list(outcome='vmPFC_decon',model_name='model1', var = 'trial_neg_inv_sc',
                                  specs = formula(~trial_neg_inv_sc),at=list(trial_neg_inv_sc = qT)),
                        V = list(outcome='vmPFC_decon',model_name='model1', var = 'v_max_wi',
                                 specs = formula(~v_max_wi),at=list(v_max_wi=c(-1.5,1.5))),
                        H = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_sc',
                                 specs = formula(~v_entropy_sc),at=list(v_entropy_sc=c(-1.5,1.5))),
                        dH = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi_change',
                                  specs=formula(~v_entropy_wi_change), at=list(v_entropy_wi_change=c(-0.5,0.5)))
                      )
      )
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-symmetry-feedback-',i,'.Rdata'))
    }
  }
  if (do_network){
    splits = c('evt_time','network')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', 
                                 specs=c("v_entropy_sc"), at = list(v_entropy_sc=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                        V = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5))),
                        dH = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=c("v_entropy_wi_change"), at=list(v_entropy_wi_change=c(-0.5,0.5)))
                      ),
                      emtrends_spec = list(
                        Tr = list(outcome='vmPFC_decon',model_name='model1', var = 'trial_neg_inv_sc',
                                  specs = formula(~trial_neg_inv_sc),at=list(trial_neg_inv_sc = qT)),
                        V = list(outcome='vmPFC_decon',model_name='model1', var = 'v_max_wi',
                                 specs = formula(~v_max_wi),at=list(v_max_wi=c(-1.5,1.5))),
                        H = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_sc',
                                 specs = formula(~v_entropy_sc),at=list(v_entropy_sc=c(-1.5,1.5))),
                        dH = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi_change',
                                  specs=formula(~v_entropy_wi_change), at=list(v_entropy_wi_change=c(-0.5,0.5)))
                      )
      )
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-network-feedback-',i,'.Rdata'))
    }
  }
}
###################################
#####      vmPFC - clock      #####
###################################
if (do_vPFC_clock){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
  vmPFC <- clock_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(clock_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  source('~/vmPFC/get_trial_data_vmPFC.R')
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
  df <- df %>% 
    group_by(id, run) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           iti_lag_sc = scale(iti_prev),
           score_sc = scale(score_csv),
           score_lag_sc = scale(lag(score_csv)),
           iti_sc = scale(iti_ideal),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           rt_vmax_sc = scale(rt_vmax),
           ev_sc = scale(ev),
           ev_lag_sc = scale(lag(ev)),
           v_entropy_sc = scale(v_entropy),
           v_entropy_lag_sc = scale(lag(v_entropy)),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
           abs_rt_vmax_change = abs(rt_vmax_change),
           rt_vmax_change_bin = case_when(
             abs_rt_vmax_change < 4/24 ~ 'No Change',
             abs_rt_vmax_change >= 4/24 ~ 'Change'
           ),
           v_entropy_wi_change_lag_bin = case_when(
             v_entropy_wi_change_lag < -0.5 ~ 'Decrease',
             v_entropy_wi_change_lag > 0.5  ~ 'Increase',
             v_entropy_wi_change_lag >= -0.5 & v_entropy_wi_change_lag <= 0.5 ~ 'No Change'
           ),
           rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld4 = case_when(
             kld4 ==0 ~ NA_real_,
             kld4 >0 ~ log10(kld4)
           )) %>% mutate(log10kld4_lag = case_when(
             kld4_lag==0 ~NA_real_,
             kld4_lag>0 ~ log10(kld4_lag)
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
    run_trial <= 15 ~ 'Early',
    run_trial > 15 & run_trial < 30 ~ 'Middle',
    run_trial >=30 ~ 'Late',
  )))
  #df <- df %>% filter(!is.na(rt_vmax_change_bin) | !is.na(v_entropy_wi_change_lag_bin))
  df <- df %>% select(id,run,run_trial,rt_vmax_change_sc,v_entropy_wi_change_lag,iti_ideal, iti_prev, rt_csv, trial_bin,rewFunc,v_entropy_sc,expl_longer,rt_csv_sc, trial_neg_inv_sc,expl_shorter,rt_bin,trial_bin,last_outcome,v_max_wi,v_entropy_wi_change_lag,score_lag_sc,iti_sc,iti_lag_sc,ev_lag_sc)
  Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  #Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
  Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
  Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
  Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
  Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
  #Q$v_entropy_wi_change_lag_bin <- relevel(as.factor(Q$v_entropy_wi_change_lag_bin),ref='No Change')
  #Q$rt_vmax_change_bin <- relevel(as.factor(Q$rt_vmax_change_bin),ref='No Change')
  # test age & sex
  demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
  demo <- demo %>% rename(id=lunaid)
  demo <- demo %>% select(!adult & !scandate)
  Q <- inner_join(Q,demo,by=c('id'))
  Q$female <- relevel(as.factor(Q$female),ref='0')
  Q$age <- scale(Q$age)
  
  rm(decode_formula)
  decode_formula <- formula(~ (1|id))
  decode_formula[[1]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change_lag + trial_neg_inv_sc + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + rt_vmax_change_sc + (1|id/run))
  decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change_lag + trial_neg_inv_sc + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + rt_vmax_change_sc + (1 + v_entropy_sc |id/run))
  #decode_formula[[2]] = formula(~ v_entropy_lag_sc*last_outcome + v_entropy_lag_sc*trial_neg_inv_sc + v_max_wi_lag*last_outcome + v_entropy_wi_change_lag + rt_csv_sc + iti_lag_sc + (1|id/run))
  #decode_formula[[2]] = formula(~ v_entropy_lag_sc*last_outcome + v_entropy_lag_sc*trial_neg_inv_sc + v_max_wi_lag*last_outcome + v_entropy_wi_change_lag + rt_csv_sc + iti_lag_sc +  (1+v_max_wi_lag + v_entropy_lag_sc | id/run))
  
  qT <- c(-0.7,0.43)
  if (do_symmetry){
    splits = c('evt_time','symmetry_group')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', 
                                 specs=c("v_entropy_sc"), at = list(v_entropy_sc=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                        V = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5))),
                        dH = list(outcome='vmPFC_decon', model_name='model1',
                                  specs=c('v_entropy_wi_change_lag'),at = list(v_entropy_wi_change_lag=c(-0.5,0.5)))
                      ),
                      emtrends_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', var='v_entropy_sc',
                                 specs=formula(~v_entropy_sc), at = list(v_entropy_sc=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon', model_name='model1', var = 'trial_neg_inv_sc',
                                  specs=formula(~trial_neg_inv_sc), at = list(trial_neg_inv_sc=qT)),
                        V = list(outcome='vmPFC_decon', model_name='model1', var = 'v_max_wi',
                                 specs=formula(~v_max_wi), at=list(v_max_wi=c(-1.5,1.5))),
                        dH = list(outcome='vmPFC_decon', model_name='model1', var = 'v_entropy_wi_change_lag',
                                  specs=formula(~v_entropy_wi_change_lag),at = list(v_entropy_wi_change_lag=c(-0.5,0.5)))
                      )
      )
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-symmetry-clock-',i,'.Rdata'))
    }
  }
  if (do_network){
    splits = c('evt_time','network')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', 
                                 specs=c("v_entropy_sc"), at = list(v_entropy_sc=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                        V = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5))),
                        dH = list(outcome='vmPFC_decon', model_name='model1',
                                  specs=c('v_entropy_wi_change_lag'),at = list(v_entropy_wi_change_lag=c(-0.5,0.5)))
                      ),
                      emtrends_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', var='v_entropy_sc',
                                 specs=formula(~v_entropy_sc), at = list(v_entropy_sc=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon', model_name='model1', var = 'trial_neg_inv_sc',
                                  specs=formula(~trial_neg_inv_sc), at = list(trial_neg_inv_sc=qT)),
                        V = list(outcome='vmPFC_decon', model_name='model1', var = 'v_max_wi',
                                 specs=formula(~v_max_wi), at=list(v_max_wi=c(-1.5,1.5))),
                        dH = list(outcome='vmPFC_decon', model_name='model1', var = 'v_entropy_wi_change_lag',
                                  specs=formula(~v_entropy_wi_change_lag),at = list(v_entropy_wi_change_lag=c(-0.5,0.5)))
                      )
      )
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-network-clock-',i,'.Rdata'))
    }
  }
}