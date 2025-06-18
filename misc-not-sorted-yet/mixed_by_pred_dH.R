# mixed_by_pred_dH
# 2022-12-09 AndyP
# model to predict entropy change


# 2022-03-31 AndyP
# mixed_by model building for AIC comparison

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26
toalign <- 'feedback'
do_rand_slopes = FALSE
do_rt_pred_fmri = TRUE
plot_rt_pred_fmri = TRUE
do_rt_pred_meg = TRUE
plot_rt_pred_meg = TRUE
#### clock ####

if (do_rand_slopes){
  load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
  vmPFC <- fb_comb
  vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
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
  df <- df %>% select(id,run,run_trial,rt_vmax_change_sc,v_entropy_wi_change,v_entropy_wi_change_lag, trial_bin,iti_ideal,iti_prev,rt_csv,trial_neg_inv_sc,rt_csv_sc,rewFunc,v_entropy_sc,outcome,v_max_wi,score_sc,rt_bin,iti_sc,ev_sc,expl_longer,expl_shorter)
  Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$vmPFC_decon[Q$evt_time > Q$iti_ideal] = NA;
  Q <- Q %>% mutate(online = case_when(
    Q$evt_time <= -Q$rt_csv-Q$iti_prev ~ 'online_pre',
    Q$evt_time < -Q$rt_csv & Q$evt_time <= 0 & Q$evt_time > -Q$rt_csv-Q$iti_prev ~ 'offline_pre',
    Q$evt_time >= -Q$rt_csv & Q$evt_time <= 0 ~ 'online',
    Q$evt_time > 0 ~ 'offline_post'
  ))
  Q <- Q %>% filter(!is.na(online))
  Q$online <- relevel(as.factor(Q$online),ref='offline_pre')
  #Q$vmPFC_decon[Q$evt_time < -(Q$rt_csv)] = NA;
  Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
  Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
  Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
  Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
  demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
  demo <- demo %>% rename(id=lunaid)
  demo <- demo %>% select(!adult & !scandate)
  Q <- inner_join(Q,demo,by=c('id'))
  Q$female <- relevel(as.factor(Q$female),ref='0')
  Q$age <- scale(Q$age)

  decode_formula = formula(~ age + female + v_entropy_sc + v_entropy_wi_change + v_max_wi + outcome + trial_neg_inv_sc + rt_csv_sc + iti_sc + (1 + evt_time | id/run))
  
  splits = {'network'}
  source("~/fmri.pipeline/R/mixed_by.R")
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,return_models=TRUE,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-vmPFC-network-',toalign,'-ranslopes-dH-',i,'.Rdata'))
}

if (do_rt_pred_fmri){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vmPFC-network-',toalign,'-ranslopes-dH-',1,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    
    ### Does random slope predict rt_vmax or rt_swing?
    rm(Q2)
  qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='evt_time' & group=='id')
  qdf <- qdf %>% group_by(network) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
  qdf <- qdf %>% rename(id=level)
  qdf$id <- as.character(qdf$id)
  qdf <- qdf %>% select(!outcome)
  source('~/vmPFC/get_trial_data_vmPFC.R')
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
  df <- df %>% select(rt_vmax_lag,rt_next,v_max_wi_lag,ev,iti_ideal,score_csv,v_max,last_outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                      v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                      rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
    group_by(id, run) %>% 
    mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           score_sc = scale(score_csv),
           iti_sc = scale(iti_ideal),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           ev_sc = scale(ev),
           rt_vmax_sc = scale(rt_vmax),
           v_entropy_sc = scale(v_entropy),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
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
  Q2 <- merge(qdf,df,by=c('id'))
  Q2 <- Q2 %>% rename(subj_level_rand_slope=estimate)
  rm(decode_formula)
  
  decode_formula = formula(~ (subj_level_rand_slope + trial_neg_inv_sc + iti_sc + rt_csv_sc + rt_lag_sc + v_entropy_sc)^2 + (1|id/run))
  
  qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
  qRTV <- quantile(df$rt_vmax_lag_sc,c(0.1,0.9),na.rm=TRUE)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
  #qH <- c(1,2,3,4)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qRS <- quantile(Q2$subj_level_rand_slope, c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
  splits = c('evt_time','network')
  source('~/fmri.pipeline/R/mixed_by.R')
  print(i)
  ddq <- mixed_by(Q2, outcomes = "v_entropy_wi_change", rhs_model_formulae = decode_formula, split_on = splits,return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed"),conf.int=TRUE)
  )
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-pred-dH-',i,'.Rdata'))
}

if (plot_rt_pred_fmri){
  source('~/vmPFC/plot_subject_level_random_slopes.R')
  #source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
  #source('~/vmPFC/plot_emmeans_subject_level_random_slopes.R')
  for (i in 1){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
    model_str <- paste0('-vmPFC-network-',toalign,'-ranslopes-dH-',1,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    model_iter <- i
    totest <- paste0('dH-',as.character(i))
    #toprocess <- 'symmetry'
    toprocess <- 'network'
    behavmodel <- 'compressed'
    hc_LorR <- 'LR'
    plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    #plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    #plot_emmeans_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  }
}

# replication

### Does random slope predict rt_vmax or rt_swing?
if (do_rt_pred_meg){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
  model_str <- paste0('-vmPFC-network-',toalign,'-ranslopes-',3,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  
  ### Does random slope predict rt_vmax or rt_swing?
  qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='evt_time' & group=='id')
  
  qdf <- qdf %>% group_by(network) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
  qdf <- qdf %>% rename(id=level)
  qdf$id <- as.character(qdf$id)
  qdf <- qdf %>% select(!outcome)
  source('~/vmPFC/get_trial_data_vmPFC.R')
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_meg')
  df <- df %>% select(rt_vmax_lag,rt_next,v_max_wi_lag,ev,score_csv,v_max,last_outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                      v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,
                      rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
    group_by(id, run) %>% 
    mutate(rt_sec = rt_csv/1000) %>% ungroup() %>%
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           score_sc = scale(score_csv),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           ev_sc = scale(ev),
           rt_vmax_sc = scale(rt_vmax),
           v_entropy_sc = scale(v_entropy),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
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
  Q2 <- merge(qdf,df,by=c('id'))
  Q2 <- Q2 %>% rename(subj_level_rand_slope=estimate)
  rm(decode_formula)
  
  decode_formula = formula(~ (subj_level_rand_slope + trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_entropy_sc)^2 + (1|id/run))
  
  qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
  qRTV <- quantile(df$rt_vmax_lag_sc,c(0.1,0.9),na.rm=TRUE)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
  #qH <- c(1,2,3,4)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qRS <- quantile(Q2$subj_level_rand_slope, c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
  splits = c('evt_time','network')
  source('~/fmri.pipeline/R/mixed_by.R')
  print(i)
  ddq <- mixed_by(Q2, outcomes = "v_entropy_wi_change", rhs_model_formulae = decode_formula, split_on = splits,return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed"),conf.int=TRUE)
  )
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-pred-dH-replication',i,'.Rdata'))
}
if (plot_rt_pred_meg){
  source('~/vmPFC/plot_subject_level_random_slopes.R')
  #source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
  #source('~/vmPFC/plot_emmeans_subject_level_random_slopes.R')
  for (i in 1){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
    model_str <- paste0('-vmPFC-network-ranslopes-',toalign,'-pred-dH-replication',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    model_iter <- i
    totest <- paste0('dH-replication-',as.character(i))
    #toprocess <- 'symmetry-by'
    toprocess <- 'network'
    behavmodel <- 'compressed'
    hc_LorR <- 'LR'
    plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    #plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    #plot_emmeans_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  }
}
