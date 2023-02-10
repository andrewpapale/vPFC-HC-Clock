# 2022-03-31 AndyP
# mixed_by model building for AIC comparison

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
doTesting = FALSE
do_vPFC_fb = FALSE
do_vPFC_clock = FALSE
do_HC_fb = FALSE
do_HC_clock = FALSE
do_HC2vPFC_fb = FALSE
do_HC2vPFC_clock = TRUE
do_anat_fb = FALSE
do_anat_clock = FALSE
do_symmetry = FALSE
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
    group_by(id) %>% 
    mutate(v_entropy_sc_r = scale(v_entropy)) %>%
    group_by(id, run) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           score_sc = scale(score_csv),
           iti_sc = scale(iti_ideal),
           iti_prev_sc = scale(iti_prev),
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
  df <- df %>% select(id,run,run_trial,last_outcome,trial_bin,rt_bin,v_entropy_wi,v_chosen_sc,rt_vmax_change_sc,v_entropy_wi_change,trial_bin,iti_ideal,iti_prev,iti_prev_sc,rt_csv,trial_neg_inv_sc,rt_csv_sc,rewFunc,v_entropy_sc,outcome,v_max_wi,score_sc,rt_bin,iti_sc,ev_sc,expl_longer,expl_shorter)
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
  decode_formula[[1]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + last_outcome + outcome +  (1|id/run))
  #decode_formula[[2]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + last_outcome + outcome +  (1 + v_entropy_wi |id/run))
  #decode_formula[[3]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + last_outcome + outcome +  (1 + v_entropy_wi | id) + (1 | id/run))
  #decode_formula[[4]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + last_outcome + outcome +  (1 + v_entropy_wi | run) + (1| id/run))
  # 
  decode_formula[[2]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + last_outcome + outcome + (1|id/run))
  #decode_formula[[6]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + last_outcome + outcome + (1 + v_max_wi |id/run))
  #decode_formula[[7]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + last_outcome + outcome + (1 + v_max_wi | id) + (1 | id/run))
  #decode_formula[[8]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + last_outcome + outcome + (1 + v_max_wi | run) + (1| id/run))
  qT <- c(-0.7,0.43)
  if (do_symmetry){
    splits = c('evt_time','symmetry_group')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      if (i < 2){
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon', model_name='model1',
                                  specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT))
                      )
      )
      } else{
        ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                        emmeans_spec = list(
                          Tr = list(outcome='vmPFC_decon', model_name='model1',
                                    specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                          V = list(outcome='vmPFC_decon', model_name='model1',
                                   specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5)))
                        )
        )       
      }
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
      if (i < 2){
        ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
        )
      } else{
        ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
        )       
      }
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
  df <- df %>% select(id,run,run_trial,rt_vmax_change_sc,v_entropy_wi,outcome,v_entropy_wi_change_lag,iti_ideal, iti_prev, rt_csv, trial_bin,rewFunc,v_entropy_sc,expl_longer,rt_csv_sc, trial_neg_inv_sc,expl_shorter,rt_bin,trial_bin,last_outcome,v_max_wi,v_entropy_wi_change_lag,score_lag_sc,iti_sc,iti_lag_sc,ev_lag_sc)
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
  decode_formula[[1]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + last_outcome + outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1|id/run))
  #decode_formula[[2]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + last_outcome + outcome+ rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_entropy_wi |id/run))
  #decode_formula[[3]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + last_outcome + outcome+ rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_entropy_wi | id) + (1 | id/run))
  #decode_formula[[4]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + last_outcome + outcome+ rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_entropy_wi | run) + (1 | id/run))
  decode_formula[[2]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + last_outcome + outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 |id/run))
  #decode_formula[[5]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + last_outcome + outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_max_wi |id/run))
  #decode_formula[[6]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + last_outcome + outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_max_wi | id) + (1 | id/run))
  #decode_formula[[7]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + last_outcome + outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_max_wi | run) + (1 | id/run))

  qT <- c(-0.7,0.43)
  if (do_symmetry){
    splits = c('evt_time','symmetry_group')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      if (i < 2){
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', 
                                 specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT))
                      )
                      )
      } else {
        ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                        emmeans_spec = list(
                          Tr = list(outcome='vmPFC_decon', model_name='model1', 
                                    specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                          V = list(outcome='vmPFC_decon', model_name='model1',
                                   specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5)))
                        )
                        )
      }
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
      if (i < 2){
        ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                        emmeans_spec = list(
                          H = list(outcome='vmPFC_decon', model_name='model1', 
                                   specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(-1.5,1.5))),
                          Tr = list(outcome='vmPFC_decon', model_name='model1', 
                                    specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT))
                        )
        )
      } else {
        ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                        emmeans_spec = list(
                          Tr = list(outcome='vmPFC_decon', model_name='model1', 
                                    specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                          V = list(outcome='vmPFC_decon', model_name='model1',
                                   specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5)))
                        )
        )        
      }
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-network-clock-',i,'.Rdata'))
    }
  }
}
#################################
#####   HC - feedback      ######
#################################
if (do_HC_fb){
  rm(Q)
  setwd('~/vmPFC')
  message('adding HC signals to models...')
  load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
  hc <- fb_comb
  hc <- hc %>% filter(evt_time > -4 & evt_time < 5)
  rm(fb_comb)
  
  hc <- hc %>% select(id,run,run_trial,decon_mean,evt_time,bin_num,side)
  hc <- hc %>% mutate(
    HC_region = case_when(
      bin_num <= 8 ~ 'AH',
      bin_num >8 ~ 'PH'
    ))
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon_mean1 = mean(decon_mean,na.rm=TRUE)) 
  hc <- hc %>% rename(decon_mean=decon_mean1)
  
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
  df <- df %>% select(id,run,run_trial,v_entropy_wi_change,rt_vmax_change_sc,trial_bin,iti_ideal,iti_prev, rt_csv,rt_csv_sc, trial_neg_inv_sc,rewFunc,v_entropy_sc,outcome,v_max_wi,v_entropy_wi_change,score_sc,rt_bin,iti_sc,ev_sc,expl_longer,expl_shorter)
  Q <- merge(df, hc, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$decon_mean[Q$evt_time > Q$iti_ideal] = NA;
  #Q$decon_mean[Q$evt_time < -(Q$rt_csv)] = NA;
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
  decode_formula[[1]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change + outcome + rt_csv_sc + iti_sc + (1|id/run))
  decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change + outcome + rt_csv_sc + iti_sc + (1 + v_entropy_sc |id/run))
  decode_formula[[3]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change + outcome + rt_csv_sc + iti_sc + (1 + v_entropy_sc |id) + (1|run))
  decode_formula[[4]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change + outcome + rt_csv_sc + iti_sc + (1 + v_entropy_sc |run) + (1|id))
  decode_formula[[5]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change + outcome + rt_csv_sc + iti_sc + (1 + v_max_wi |id/run))
  decode_formula[[6]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change + outcome + rt_csv_sc + iti_sc + (1 + v_max_wi |id) + (1|run))
  decode_formula[[7]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change + outcome + rt_csv_sc + iti_sc + (1 + v_max_wi |run) + (1|id))
  #decode_formula[[2]] = formula(~ v_entropy_sc*last_outcome + v_entropy_sc*trial_neg_inv_sc + v_max_wi*last_outcome + rt_csv_sc + iti_sc + (1|id/run))
  #decode_formula[[2]] = formula(~ v_entropy_sc*last_outcome + v_entropy_sc*trial_neg_inv_sc + v_max_wi*last_outcome + rt_csv_sc + iti_sc +  (1+v_max_wi + v_entropy_sc |id/run))
  
  qT <- c(-0.7,0.43)
  splits = c('evt_time','HC_region')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                    emmeans_spec = list(
                      H = list(outcome='decon_mean', model_name='model1',
                               specs=c("v_entropy_sc"), at = list(v_entropy_sc=c(-1.5,1.5))),
                      Tr = list(outcome='vmPFC_decon', model_name='model1',
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                      V = list(outcome='decon_mean', model_name='model1',
                               specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5))),
                      dH = list(outcome='decon_mean', model_name='model1',
                                specs=c('v_entropy_wi_change'),at=list(v_entropy_wi_change=c(-0.5,0.5))),
                      RT = list(outcome='decon_mean',model_name='model1',
                                specs=c('rt_vmax_change_sc'),at=list(rt_vmax_change_sc=c(-0.17,0.17))),
                      HxTr = list(outcome='decon_mean',model_name='model1',
                                  specs=formula(~v_entropy_sc:trial_neg_inv_sc),at=list(v_entropy_sc=c(-1.5,1.5),trial_neg_inv_sc=qT))
                    ),
                    emtrends_spec = list(
                      H = list(outcome='decon_mean', model_name='model1',var ='v_entropy_sc',
                               specs=formula(~v_entropy_sc), at = list(v_entropy_sc=c(-1.5,1.5))),
                      Tr = list(outcome='vmPFC_decon', model_name='model1', var='trial_neg_inv_sc',
                                specs=formula(~trial_neg_inv_sc), at = list(trial_neg_inv_sc=qT)),
                      V = list(outcome='decon_mean', model_name='model1',var='v_max_wi',
                               specs=formula(~v_max_wi), at=list(v_max_wi=c(-1.5,1.5))),
                      dH = list(outcome='decon_mean', model_name='model1', var='v_entropy_wi_change',
                                specs=formula(~v_entropy_wi_change),at=list(v_entropy_wi_change=c(-0.5,0.5))),
                      RT = list(outcome='decon_mean',model_name='model1', var = 'rt_vmax_change_sc',
                                specs=formula(~rt_vmax_change_sc),at=list(rt_vmax_change_sc=c(-0.17,0.17))),
                      HxTr = list(outcome='decon_mean',model_name='model1', var = 'v_entropy_sc',
                                  specs=formula(~v_entropy_sc:trial_neg_inv_sc),at=list(v_entropy_sc=c(-1.5,1.5),trial_neg_inv_sc=qT))
                    )
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-HC-axis-feedback-',i,'.Rdata'))
  }
}
###################################
#####      HC- clock          #####
###################################
if (do_HC_clock){
  rm(Q)
  setwd('~/vmPFC')
  message('adding HC signals to models...')
  load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
  hc <- clock_comb
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
  rm(clock_comb)
  
  hc <- hc %>% select(id,run,run_trial,decon_mean,evt_time,bin_num,side)
  hc <- hc %>% mutate(
    HC_region = case_when(
      bin_num <= 8 ~ 'AH',
      bin_num >8 ~ 'PH'
    ))
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon_mean1 = mean(decon_mean,na.rm=TRUE)) 
  hc <- hc %>% rename(decon_mean=decon_mean1)
  
  source('~/vmPFC/get_trial_data_vmPFC.R')
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
  df <- df %>% 
    group_by(id,run) %>%
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
  df <- df %>% select(id,run,run_trial,v_entropy_wi_change_lag,rt_vmax_change_sc,iti_ideal,iti_prev,rt_csv,trial_bin,rewFunc,v_entropy_sc,expl_longer,expl_shorter,rt_bin,trial_bin,last_outcome,v_max_wi,v_entropy_wi_change_lag,score_lag_sc,iti_lag_sc,iti_sc,ev_lag_sc,rt_csv_sc,trial_neg_inv_sc)
  Q <- merge(df, hc, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$decon_mean[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  #Q$decon_mean[Q$evt_time < -(Q$iti_prev)] = NA;
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
  decode_formula <- NULL
  decode_formula[[1]] = formula(~ age + female + v_entropy_wi + last_outcome + outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1|id/run))
  decode_formula[[2]] = formula(~ age + female + v_max_wi + last_outcome + outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1|id/run))
  
  # decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_max_wi + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_entropy_sc |id/run))
  # decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_max_wi + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_entropy_sc |id) + (1|run))
  # decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_max_wi + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_entropy_sc |run) + (1|id))
  # decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_max_wi + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_max_wi |id/run))
  # decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_max_wi + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_max_wi |id) + (1|run))
  # decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_max_wi + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + (1 + v_max_wi |run) + (1|id))
  splits = c('evt_time','HC_region')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                    emmeans_spec = list(
                      H = list(outcome='decon_mean', model_name='model1',
                               specs=c("v_entropy_wi"), at = list(v_entropy_sc=c(-1.5,1.5))),
                      Tr = list(outcome='vmPFC_decon', model_name='model1',
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                      V = list(outcome='decon_mean', model_name='model1',
                               specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5))),
                      dH = list(outcome='decon_mean', model_name='model1',
                                specs=c('v_entropy_wi_change_lag'),at=list(v_entropy_wi_change_lag=c(-0.5,0.5))),
                      RT = list(outcome='decon_mean',model_name='model1',
                                specs=c('rt_vmax_change_sc'),at=list(rt_vmax_change_sc=c(-0.17,0.17))),
                      HxTr = list(outcome='decon_mean',model_name='model1',
                                  specs=formula(~v_entropy_sc:trial_neg_inv_sc),at=list(v_entropy_sc=c(-1.5,1.5),trial_neg_inv_sc=qT))
                    ),
                    emtrends_spec = list(
                      H = list(outcome='decon_mean', model_name='model1',var ='v_entropy_sc',
                               specs=formula(~v_entropy_sc), at = list(v_entropy_sc=c(-1.5,1.5))),
                      Tr = list(outcome='vmPFC_decon', model_name='model1', var='trial_neg_inv_sc',
                                specs=formula(~trial_neg_inv_sc), at = list(trial_neg_inv_sc=qT)),
                      V = list(outcome='decon_mean', model_name='model1',var='v_max_wi',
                               specs=formula(~v_max_wi), at=list(v_max_wi=c(-1.5,1.5))),
                      dH = list(outcome='decon_mean', model_name='model1', var='v_entropy_wi_change_lag',
                                specs=formula(~v_entropy_wi_change_lag),at=list(v_entropy_wi_change_lag=c(-0.5,0.5))),
                      RT = list(outcome='decon_mean',model_name='model1', var = 'rt_vmax_change_sc',
                                specs=formula(~rt_vmax_change_sc),at=list(rt_vmax_change_sc=c(-0.17,0.17))),
                      HxTr = list(outcome='decon_mean',model_name='model1', var = 'v_entropy_sc',
                                  specs=formula(~v_entropy_sc:trial_neg_inv_sc),at=list(v_entropy_sc=c(-1.5,1.5),trial_neg_inv_sc=qT))
                    )
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-HC-axis-clock-',i,'.Rdata'))
  }
}
####################################
### vPFC - HC feedback - anatomy ###
####################################
if (do_anat_fb){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
  vmPFC <- fb_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(fb_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
  hc <- fb_comb
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
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
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
  df <- df %>% 
    group_by(id,run) %>%
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
  
  df <- df %>% select(id,run,run_trial,trial_neg_inv_sc,iti_ideal,iti_prev,rt_csv,rt_csv_sc,iti_sc,rt_csv_sc, rewFunc)
  
  Q <- inner_join(Q,df,by=c('id','run','run_trial'))
  Q$vmPFC_decon[Q$evt_time > Q$iti_ideal] = NA;
  #Q$vmPFC_decon[Q$evt_time < -(Q$rt_csv)] = NA;
  Q$decon_mean[Q$evt_time > Q$iti_ideal] = NA;
  #Q$decon_mean[Q$evt_time < -(Q$rt_csv)] = NA;  
 
  # test age & sex
  demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
  demo <- demo %>% rename(id=lunaid)
  demo <- demo %>% select(!adult & !scandate)
  Q <- inner_join(Q,demo,by=c('id'))
  Q$female <- relevel(as.factor(Q$female),ref='0')
  Q$age <- scale(Q$age)
  
  Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
  
  decode_formula <- formula(~ (1|id))
  decode_formula[[1]] <- formula(~ age + female + HCwithin*trial_neg_inv_sc + rt_csv_sc + iti_sc + HCbetween*HCwithin + (1|id/run))
  decode_formula[[2]] <- formula(~ age + female + HCwithin*trial_neg_inv_sc + rt_csv_sc + iti_sc + HCbetween*HCwithin + (1 + HCwithin |id/run))
  qT <- c(-0.7,0.43)
  
  if (do_network){
    
    splits = c('evt_time','network','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        HC = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=formula(~HCwithin), at = list(HCwithin=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5))),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',
                                    specs=formula(~HCbetween:HCwithin), at=list(HCwithin=c(-1.5,1.5),HCbetween=c(-1.5,1.5)))
                      )
      )
      
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-anatomy-',i,'.Rdata'))
    }
  }
  if (do_symmetry){
    
    splits = c('evt_time','symmetry_group','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        HC = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=formula(~HCwithin), at = list(HCwithin=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5))),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',
                                    specs=formula(~HCbetween:HCwithin), at=list(HCwithin=c(-1.5,1.5),HCbetween=c(-1.5,1.5)))
                      )
      )
      
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-clock-anatomy-',i,'.Rdata'))
    }
  }
}
  

####################################
### vPFC - HC clock - anatomy ###
####################################
if (do_anat_clock){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
  vmPFC <- clock_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(clock_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
  hc <- clock_comb
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
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
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
  df <- df %>% 
    group_by(id,run) %>%
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

  df <- df %>% select(id,run,run_trial,trial_neg_inv_sc,iti_ideal,iti_prev,rt_csv,rt_csv_sc,iti_sc,rt_csv_sc,rewFunc,iti_lag_sc)
  
  Q <- inner_join(Q,df,by=c('id','run','run_trial'))
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
  Q$decon_mean[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  Q$decon_mean[Q$evt_time < -(Q$iti_prev)] = NA;  

  # test age & sex
  demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
  demo <- demo %>% rename(id=lunaid)
  demo <- demo %>% select(!adult & !scandate)
  Q <- inner_join(Q,demo,by=c('id'))
  Q$female <- relevel(as.factor(Q$female),ref='0')
  Q$age <- scale(Q$age)
  
  Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
  
  decode_formula <- NULL
  decode_formula[[1]] <- formula(~ age + female + HCwithin*trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_lag_sc + HCbetween*HCwithin + (1|id/run))
  decode_formula[[2]] <- formula(~ age + female + HCwithin*trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_lag_sc + HCbetween*HCwithin + (1 + HCwithin | id/run))
  decode_formula[[3]] <- formula(~ age + female + HCwithin*trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_lag_sc + HCbetween*HCwithin + (1 + HCwithin | id) + (1|run))
  decode_formula[[4]] <- formula(~ age + female + HCwithin*trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_lag_sc + HCbetween*HCwithin + (1 + HCwithin | run) + (1|id))
  qT <- c(-0.7,0.43)
  
  if (do_network){
    
    splits = c('evt_time','network','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        HC = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=formula(~HCwithin), at = list(HCwithin=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5))),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',
                                    specs=formula(~HCbetween:HCwithin), at=list(HCwithin=c(-1.5,1.5),HCbetween=c(-1.5,1.5)))
                      )
      )
      
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-feedback-anatomy-',i,'.Rdata'))
    }
  }
  if (do_symmetry){
    
    splits = c('evt_time','symmetry_group','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        HC = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=formula(~HCwithin), at = list(HCwithin=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5))),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',
                                    specs=formula(~HCbetween:HCwithin), at=list(HCwithin=c(-1.5,1.5),HCbetween=c(-1.5,1.5)))
                      )
      )
      
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-feedback-anatomy-',i,'.Rdata'))
    }
  }
}

###################################
##### vmPFC - HC -feedback    #####
###################################
if (do_HC2vPFC_fb){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
  vmPFC <- fb_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(fb_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
  hc <- fb_comb
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
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
  df <- df %>% select(id,run,run_trial,v_entropy_wi_change,rt_vmax_change_sc,trial_bin,iti_ideal,iti_prev,rt_csv,rewFunc,trial_neg_inv_sc,rt_csv_sc,v_entropy_sc,outcome,v_max_wi,trial_neg_inv_sc,rt_csv_sc,v_entropy_wi_change,score_sc,rt_bin,iti_sc,ev_sc,expl_longer,expl_shorter)
  Q <- merge(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$vmPFC_decon[Q$evt_time > Q$iti_ideal] = NA;
  #Q$vmPFC_decon[Q$evt_time < -(Q$rt_csv)] = NA;
  Q$decon_mean[Q$evt_time > Q$iti_ideal] = NA;
  #Q$decon_mean[Q$evt_time < -(Q$rt_csv)] = NA;
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
  
  Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
  
  rm(decode_formula)
  decode_formula <- formula(~ (1|id))
  decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_sc*HCwithin + trial_neg_inv_sc*HCwithin +v_max_wi*HCwithin  + v_entropy_wi_change*HCwithin  + rt_csv_sc*HCwithin  + iti_sc*HCwithin + outcome*HCwithin + rt_vmax_change_sc*HCwithin + HCbetween*HCwithin + (1|id/run))
  decode_formula[[2]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_sc*HCwithin + trial_neg_inv_sc*HCwithin + v_max_wi*HCwithin  + v_entropy_wi_change*HCwithin   + rt_csv_sc*HCwithin  + iti_sc*HCwithin + outcome*HCwithin + rt_vmax_change_sc*HCwithin + HCbetween*HCwithin + (1 + HCwithin |id/run))
  #decode_formula[[2]] <- NULL
  
  qT <- c(-0.7,0.43)

  
  if (do_network){
    
    splits = c('evt_time','network','HC_region')
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
                                 specs=formula(~v_entropy_sc:HCwithin), at = list(v_entropy_sc=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        V = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=formula(~v_max_wi:HCwithin), at=list(v_max_wi=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5)))
                      ),
                      emtrends_spec = list(
                        H_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_entropy_sc:HCwithin), at=list(HCwithin=c(-1.5,1.5),v_entropy_sc=c(-1.5,1.5))),
                        T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~trial_neg_inv_sc:HCwithin),at=list(HCwithin=c(-1.5,1.5),trial_neg_inv_sc=qT)),
                        V_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_max_wi:HCwithin),at=list(HCwithin=c(-1.5,1.5),v_max_wi=c(-1.5,1,5))),
                        A_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~age:HCwithin),at=list(HCwithin=c(-1.5,1.5),age=c(-1.5,1.5))),
                        O_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~outcome:HCwithin),at=list(HCwithin=c(-1.5,1.5))),
                        dH_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                     specs=formula(~v_entropy_wi_change:HCwithin),at=list(HCwithin=c(-1.5,1.5),v_entropy_wi_change=c(-0.5,0.5))),
                        dRT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                      specs=formula(~rt_vmax_change_sc:HCwithin),at=list(HCwithin=c(-1.5,1.5),rt_vmax_change_sc=c(-0.17,0.17))),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~HCbetween:HCwithin),at=list(HCwithin=c(-1.5,1,5),HCbetween=c(-1.5,1,5))),
                        RT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                     specs=formula(~rt_csv_sc:HCwithin),at=list(rt_csv_sc=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))
                      )
      )
      
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-feedback-',i,'.Rdata'))
    }
  }
  if (do_symmetry){
    
    splits = c('evt_time','symmetry_group','HC_region')
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
                                 specs=formula(~v_entropy_sc:HCwithin), at = list(v_entropy_sc=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        V = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=formula(~v_max_wi:HCwithin), at=list(v_max_wi=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5)))
                      ),
                      emtrends_spec = list(
                        H_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_entropy_sc:HCwithin), at=list(HCwithin=c(-1.5,1.5),v_entropy_sc=c(-1.5,1.5))),
                        T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~trial_neg_inv_sc:HCwithin),at=list(HCwithin=c(-1.5,1.5),trial_neg_inv_sc=qT)),
                        V_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_max_wi:HCwithin),at=list(HCwithin=c(-1.5,1.5),v_max_wi=c(-1.5,1,5))),
                        A_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~age:HCwithin),at=list(HCwithin=c(-1.5,1.5),age=c(-1.5,1.5))),
                        O_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~outcome:HCwithin),at=list(HCwithin=c(-1.5,1.5))),
                        dH_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                     specs=formula(~v_entropy_wi_change:HCwithin),at=list(HCwithin=c(-1.5,1.5),v_entropy_wi_change=c(-0.5,0.5))),
                        dRT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                      specs=formula(~rt_vmax_change_sc:HCwithin),at=list(HCwithin=c(-1.5,1.5),rt_vmax_change_sc=c(-0.17,0.17))),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~HCbetween:HCwithin),at=list(HCwithin=c(-1.5,1,5),HCbetween=c(-1.5,1,5))),
                        RT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                     specs=formula(~rt_csv_sc:HCwithin),at=list(rt_csv_sc=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))
                      )
      )
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-feedback-',i,'.Rdata'))
    }
  }
}
###################################
##### vmPFC - HC - clock      #####
###################################
if (do_HC2vPFC_clock){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
  vmPFC <- clock_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(clock_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
  hc <- clock_comb
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
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
             abs_rt_vmax_change >= 4/24 ~ 'Change',
           ),
           v_entropy_wi_change_lag_bin = case_when(
             v_entropy_wi_change_lag < -0.5 ~ 'Decrease',
             v_entropy_wi_change_lag > 0.5  ~ 'Increase',
             v_entropy_wi_change_lag >= -0.5 & v_entropy_wi_change_lag <= 0.5 ~ 'No Change',
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
  df <- df %>% select(id,run,run_trial,iti_ideal,iti_prev,rt_csv,abs_pe_max_lag_sc,v_entropy_wi,rt_vmax_change_sc,trial_bin,rewFunc,trial_neg_inv_sc,rt_csv_sc,v_entropy_sc,expl_longer,expl_shorter,rt_bin,trial_bin,last_outcome,v_max_wi,v_entropy_wi_change_lag,score_lag_sc,iti_sc,iti_lag_sc,ev_lag_sc)
  Q <- merge(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  #Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
  Q$decon_mean[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  #Q$decon_mean[Q$evt_time < -(Q$iti_prev)] = NA;
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
  
  Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
  
  rm(decode_formula)
  decode_formula <- NULL
  decode_formula[[1]] = formula(~ age + female + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + rt_csv_sc*HCwithin  + iti_lag_sc*HCwithin + iti_sc + last_outcome*HCwithin + outcome + HCwithin*HCbetween + (1|id/run))
  decode_formula[[2]] = formula(~ age + female + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin  + rt_csv_sc*HCwithin  + iti_lag_sc*HCwithin + iti_sc + last_outcome*HCwithin + outcome + HCwithin*HCbetween + (1| id/run))
  decode_formula[[3]] = formula(~ age + female + abs_pe_max_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin  + rt_csv_sc*HCwithin  + iti_lag_sc*HCwithin + iti_sc + HCwithin*HCbetween + (1| id/run))
  qT <- c(-0.7,0.43)
  
  if (do_network){
    
    splits = c('evt_time','network','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      if (i ==1){
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', 
                                 specs=formula(~v_entropy_sc:HCwithin), at = list(v_entropy_sc=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        V = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=formula(~v_max_wi:HCwithin), at=list(v_max_wi=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5)))
                      ),
                      emtrends_spec = list(
                        H_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_entropy_sc:HCwithin), at=list(HCwithin=c(-1.5,1.5),v_entropy_sc=c(-1.5,1.5))),
                        T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~trial_neg_inv_sc:HCwithin),at=list(trial_neg_inv_sc=qT)),
                        V_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_max_wi:HCwithin),at=list(v_max_wi=c(-1.5,1,5))),
                        O_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~last_outcome:HCwithin)),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~HCbetween:HCwithin),at=list(HCbetween=c(-1.5,1,5))),
                        RT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                     specs=formula(~rt_csv_sc:HCwithin),at=list(rt_csv_sc=c(-1.5,1.5)))
                      )
      )
      } else if (i==2){
        
      } else if (i==3){
        ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                        emmeans_spec = list(
                          PE = list(outcome='vmPFC_decon',model_name='model1',
                                    specs=formula(~abs_pe_max_lag_sc:HCwithin), at=list(abs_pe_max_lag_sc=c(-2,-1,0,1,2),HCwithin=c(-2,-1,0,1,2)))
                        ),
                        emtrends_spec = list(
                          PE_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                       specs=formula(~abs_pe_max_lag_sc:HCwithin),at=list(abs_pe_max_lag_sc=c(-2,-1,0,1,2)))
                        )
                        )
      }
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-',i,'.Rdata'))
    }
  }
  
  if (do_symmetry){
    
    splits = c('evt_time','symmetry_group','HC_region')
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', 
                                 specs=formula(~v_entropy_sc:HCwithin), at = list(v_entropy_sc=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        V = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=formula(~v_max_wi:HCwithin), at=list(v_max_wi=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5)))
                      ),
                      emtrends_spec = list(
                        H_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_entropy_sc:HCwithin), at=list(v_entropy_sc=c(-1.5,1.5))),
                        T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~trial_neg_inv_sc:HCwithin),at=list(trial_neg_inv_sc=qT)),
                        V_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_max_wi:HCwithin),at=list(v_max_wi=c(-1.5,1,5))),
                        A_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~age:HCwithin),at=list(age=c(-1.5,1.5))),
                        O_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~last_outcome:HCwithin)),
                        dH_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                     specs=formula(~v_entropy_wi_change_lag:HCwithin)),
                        dRT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                      specs=formula(~rt_vmax_change_sc:HCwithin),at=list(HCwithin=c(-1.5,1.5))),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~HCbetween:HCwithin),at=list(HCwithin=c(-1.5,1,5),HCbetween=c(-1.5,1,5))),
                        RT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                     specs=formula(~rt_csv_sc:HCwithin),at=list(rt_csv_sc=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))
                      )
      )
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-clock-',i,'.Rdata'))
    }
  }
}

##############################
### testing ...   ############
##############################
if (doTesting){
  rm(Q)
  Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
  Q <- Q %>% select(!decon1)
  source('~/vmPFC/get_trial_data_vmPFC.R')
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
  df <- df %>%
    group_by(id, run) %>% 
    mutate(iti_prev = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           rt_vmax_sc = scale(rt_vmax),
           v_entropy_sc = scale(v_entropy),
           v_chosen_sc = scale(v_chosen),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
           abs_rt_vmax_change = abs(rt_vmax_change),
           rt_vmax_change_bin = case_when(
             abs_rt_vmax_change < 4/24 ~ 'No Change',
             abs_rt_vmax_change >= 4/24 ~ 'Change'
           ),
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
  df <- df %>% select(id,run,run_trial,rt_vmax_change_bin,v_chosen_sc,trial_neg_inv_sc,v_entropy_wi,v_entropy_wi_change_lag,expl_longer,expl_shorter,rt_bin,v_max_wi)
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
  
  
  decode_formula[[1]] = formula( ~ v_max_wi*HCwithin +
                                   trial_bin*HCwithin +
                                   rt_bin + v_entropy_wi*HCwithin +
                                   HCbetween + 
                                   (1|id))
  decode_formula[[2]] = formula( ~ v_max_wi*HCwithin +
                                   trial_bin*HCwithin +
                                   rt_bin + v_entropy_wi*HCwithin +
                                   expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                   HCbetween + 
                                   (1|id))
  decode_formula[[3]] = formula( ~ v_max_wi*HCwithin +
                                   trial_bin*HCwithin +
                                   rt_bin + v_entropy_wi*HCwithin +
                                   expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                   trial_bin:v_entropy_wi +
                                   HCbetween + 
                                   (1|id))
  decode_formula[[4]] = formula( ~ v_max_wi*HCwithin +
                                   trial_bin*HCwithin +
                                   rt_bin + v_entropy_wi*HCwithin +
                                   expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                   trial_bin:v_entropy_wi:HCwithin +
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
  for (i in 1:3){
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
                                  specs=formula(~v_entropy_wi:HCwithin), at=list(v_entropy_wi=qH)),
                      T_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                  specs=(~trial_bin:HCwithin)),
                      V_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                  specs=(~v_max_wi:HCwithin), at=list(v_max_wi=c(qV)))
                    ),
                    emmeans_spec = list(
                      #H_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ peaks * HCwithin), at=list(HCwithin=qHC1)),
                      T_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ trial_bin * HCwithin), at=list(HCwithin=qHC1)),
                      V_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ v_max_wi * HCwithin), at=list(HCwithin=qHC1,v_max_wi=qV))
                    )
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-testing-trxH',i,'.Rdata'))
  }
  
  
  
  #################################
  ### vPFC - Feedback - Outcome ###
  #################################
  
  Q$rewFunc <- factor(Q$rewFunc)
  Q$rewFunc <- relevel(Q$rewFunc,ref='CEV')
  
  decode_formula <- formula(~ (1|id))
  decode_formula[[1]] <- formula(~ age + female + rt_bin + trial_bin*outcome + expl_shorter + expl_longer + (1|id))
  
  splits = c('evt_time','network')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
    
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-vmPFC-HC-feedback-network-outcome-',i,'.Rdata'))
  }
}