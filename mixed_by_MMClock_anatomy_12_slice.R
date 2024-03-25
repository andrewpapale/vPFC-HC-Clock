library(pracma)
library(tidyverse)

do_1to6 = TRUE
do_7to12= FALSE
ncores = 26;



if (do_1to6){
  rm(Q)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMclock_clock_Aug2023.Rdata')
  vmPFC <- clock_comb
  rm(clock_comb)
  vmPFC <- vmPFC %>% group_by(id,run,trial,atlas_value) %>% arrange(evt_time) %>% mutate(vmPFC_lag1 = dplyr::lag(decon_mean,1,order_by=evt_time),
                                                                                      vmPFC_lag2 = dplyr::lag(decon_mean,2,order_by=evt_time)) %>% ungroup()
  
  vmPFC <- vmPFC %>% select(id,run,trial,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  gc()
  # load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock.Rdata')
  # hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  # hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
  # hc <- hc %>% mutate(block = case_when(trial <= 40 ~ 1, 
  #                                     trial > 40 & trial <= 80 ~ 2,
  #                                     trial > 80 & trial <=120 ~ 3, 
  #                                     trial > 120 & trial <=160 ~ 4,
  #                                     trial > 160 & trial <=200 ~ 5,
  #                                     trial > 200 & trial <=240 ~ 6))
  # hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2  hc <- hc %>% rename(decon_mean=decon_mean1)
  #hc <- read_csv('/Volumes/Users/Bea/StriatumHippThalamus/clock_aligned_striatum_hipp_thalamus.csv.gz')
  #hc <- hc %>% mutate(run1 = as.integer(str_sub(run,4,4))) %>% select(!run) %>% rename(run=run1)
  #hc <- hc %>% filter(atlas_value %in% c(223,224,225,226,227,228,229,230))
  #hc #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_fb.Rdata')
  #hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  # hc <- hc %>% mutate(HC_region = case_when(atlas_value==223 ~ 'PH',
  #                                           atlas_value==224 ~ 'PH',
  #                                           atlas_value==225 ~ 'AH',
  #                                           atlas_value==226 ~ 'AH',
  #                                           atlas_value==227 ~ 'PH',
  #                                           atlas_value==228 ~ 'PH',
  #                                           atlas_value==229 ~ 'AH',
  #                                           atlas_value==230 ~ 'AH'))
  # hc <- hc %>% mutate(run_trial = case_when(trial <= 40 ~ trial,
  #                                           trial > 40 & trial <= 80 ~ trial - 40,
  #                                           trial > 80 & trial <=120 ~ trial - 80,
  #                                           trial > 120 & trial <= 160  ~ trial - 120,
  #                                           trial > 160 & trial <=200 ~ trial - 160,
  #                                           trial > 200 & trial <=240 ~ trial - 200))
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
  hc <- hc %>% select(!atlas_value0)
  #hc <- hc %>% group_by(id,run,trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  #hc <- hc %>% rename(decon_mean=decon1)
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon_mean),HCbetween=mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  hc <- hc %>% group_by(id,run,trial,atlas_value) %>% arrange(evt_time) %>% mutate(HC_lag1 = lag(HCwithin,1),
                                                                                   HC_lag2 = lag(HCwithin,2),
                                                                                   HC_lag3 = lag(HCwithin,3)) %>% ungroup()
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
  hc <- hc %>% filter(atlas_value %in% c(1,2,3,4,5,6))
  gc()
  hc <- hc %>% rename(bin = atlas_value)
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory="~/clock_analysis",dataset='mmclock_fmri')
  df <- df %>%
    group_by(id) %>% 
    mutate(v_entropy_sc_r = scale(v_entropy)) %>% ungroup() %>%
    group_by(id, run) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           score_sc = scale(score_csv),
           iti_sc = scale(iti_ideal),
           iti_lag_sc = lag(iti_sc),
           v_chosen_sc = scale(v_chosen),
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
    run_trial <= 13 ~ 'Early',
    run_trial > 13 & run_trial < 26 ~ 'Middle',
    run_trial >=26 ~ 'Late',
  )))
  df <- df %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                       total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))
  
  df <- df %>% select(iti_ideal,rt_lag_sc,iti_lag_sc,iti_prev,iti_sc,last_outcome,rt_csv,rt_csv_sc,id,run,trial)
  df$id <- as.numeric(df$id)
  Q <- full_join(vmPFC,df,by=c('id','run','trial'))
  rm(vmPFC)
  gc()
  hc$id <- as.numeric(hc$id)
  Q <- inner_join(Q,hc,by=c('id','run','trial','evt_time'))
  rm(hc)
  gc()
  Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCwithin[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_decon[Q$evt_time < -Q$iti_prev] = NA
  Q$HC_lag1[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HC_lag1[Q$evt_time < -Q$iti_prev] = NA
  Q$HC_lag2[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HC_lag2[Q$evt_time < -Q$iti_prev] = NA
  Q$HC_lag3[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HC_lag3[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_lag1[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_lag1[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_lag2[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_lag2[Q$evt_time < -Q$iti_prev] = NA
  Q <- Q %>% arrange(id,run,trial,evt_time)
  Q <- Q %>% filter(evt_time > -5 & evt_time < 5)
  
  
  rm(decode_formula)
  decode_formula <- NULL
  #decode_formula[[1]] = formula(~age * HCwithin + gender * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + v_max_wi * HCwithin + v_entropy_wi_change_lag * HCwithin + rt_csv_sc * HCwithin + iti_lag_sc * HCwithin + iti_sc * HCwithin + last_outcome * HCwithin + outcome*HCwithin + rt_vmax_change_sc * HCwithin +  HCbetween + (1 | id/run)) 
  #decode_formula[[1]] = formula(~age*HCwithin + v_entropy_wi*HCwithin + (1|id/run))
  #decode_formula[[2]] = formula(~age*HCwithin + v_max_wi*HCwithin + (1|id/run))
  #decode_formula[[3]] = formula(~age*HCwithin + gender*HCwithin + v_entropy_wi*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1|id/run))
  decode_formula <- NULL
  decode_formula[[1]] <- formula(~ HCwithin + HCbetween + (1|id/run))
  # decode_formula[[2]] <- formula(~ HCwithin + HC_lag1 + HCbetween + (1 | id/run))
  decode_formula[[2]] <- formula(~ HC_lag1 + HCbetween + (1 | id/run))
  decode_formula[[3]] <- formula(~ vmPFC_lag1 + HCwithin + (1|id/run))
  decode_formula[[4]] <- formula(~ HCwithin + (1|id/run))
  # decode_formula[[4]] <- formula(~ HCwithin + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[5]] <- formula(~ HCwithin + HC_lag1 + HC_lag2 +  HCbetween + (1 | id/run))
  # decode_formula[[6]] <- formula(~ HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[7]] <- formula(~ HCwithin + HC_lag1 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[8]] <- formula(~ HCwithin + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[9]] <- formula(~ HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[10]] <- formula(~ HC_lag1 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[11]] <- formula(~ vmPFC_lag1 + HCwithin + HCbetween + (1|id/run))
  # decode_formula[[12]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag1 + HCbetween + (1 | id/run))
  # decode_formula[[13]] <- formula(~ vmPFC_lag1 + HC_lag1 + HCbetween + (1 | id/run))
  # decode_formula[[14]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[15]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag1 + HC_lag2 +  HCbetween + (1 | id/run))
  # decode_formula[[16]] <- formula(~ vmPFC_lag1 + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[17]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag1 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[18]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[19]] <- formula(~ vmPFC_lag1 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[20]] <- formula(~ vmPFC_lag1 + HC_lag1 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[21]] <- formula(~ vmPFC_lag2 + HC_lag1 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[22]] <- formula(~ vmPFC_lag2 + HCwithin + HCbetween + (1|id/run))
  # decode_formula[[23]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag1 + HCbetween + (1 | id/run))
  # decode_formula[[24]] <- formula(~ vmPFC_lag2 + HC_lag1 + HCbetween + (1 | id/run))
  # decode_formula[[25]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[26]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag1 + HC_lag2 +  HCbetween + (1 | id/run))
  # decode_formula[[27]] <- formula(~ vmPFC_lag2 + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[28]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag1 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[29]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[30]] <- formula(~ vmPFC_lag2 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[31]] <- formula(~ vmPFC_lag2 + HC_lag1 + HC_lag3 + HCbetween + (1 | id/run))# 
  
  
  qT <- c(-0.8,0.46)
  splits = c('evt_time','network','bin')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 3:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-MMClock-vPFC-HC-network-clock-anatomy-12slice-1to6-',i,'.Rdata'))
  }
  rm(Q)
  gc()  
}

if (do_7to12){
  rm(Q)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMclock_clock_Aug2023.Rdata')
  vmPFC <- clock_comb
  rm(clock_comb)
  vmPFC <- vmPFC %>% group_by(id,run,trial,atlas_value) %>% arrange(evt_time) %>% mutate(vmPFC_lag1 = dplyr::lag(decon_mean,1,order_by=evt_time),
                                                                                         vmPFC_lag2 = dplyr::lag(decon_mean,2,order_by=evt_time)) %>% ungroup()
  
  vmPFC <- vmPFC %>% select(id,run,trial,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  gc()
  # load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock.Rdata')
  # hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  # hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
  # hc <- hc %>% mutate(block = case_when(trial <= 40 ~ 1, 
  #                                     trial > 40 & trial <= 80 ~ 2,
  #                                     trial > 80 & trial <=120 ~ 3, 
  #                                     trial > 120 & trial <=160 ~ 4,
  #                                     trial > 160 & trial <=200 ~ 5,
  #                                     trial > 200 & trial <=240 ~ 6))
  # hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2  hc <- hc %>% rename(decon_mean=decon_mean1)
  #hc <- read_csv('/Volumes/Users/Bea/StriatumHippThalamus/clock_aligned_striatum_hipp_thalamus.csv.gz')
  #hc <- hc %>% mutate(run1 = as.integer(str_sub(run,4,4))) %>% select(!run) %>% rename(run=run1)
  #hc <- hc %>% filter(atlas_value %in% c(223,224,225,226,227,228,229,230))
  #hc #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_fb.Rdata')
  #hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  # hc <- hc %>% mutate(HC_region = case_when(atlas_value==223 ~ 'PH',
  #                                           atlas_value==224 ~ 'PH',
  #                                           atlas_value==225 ~ 'AH',
  #                                           atlas_value==226 ~ 'AH',
  #                                           atlas_value==227 ~ 'PH',
  #                                           atlas_value==228 ~ 'PH',
  #                                           atlas_value==229 ~ 'AH',
  #                                           atlas_value==230 ~ 'AH'))
  # hc <- hc %>% mutate(run_trial = case_when(trial <= 40 ~ trial,
  #                                           trial > 40 & trial <= 80 ~ trial - 40,
  #                                           trial > 80 & trial <=120 ~ trial - 80,
  #                                           trial > 120 & trial <= 160  ~ trial - 120,
  #                                           trial > 160 & trial <=200 ~ trial - 160,
  #                                           trial > 200 & trial <=240 ~ trial - 200))
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
  hc <- hc %>% select(!atlas_value0)
  #hc <- hc %>% group_by(id,run,trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  #hc <- hc %>% rename(decon_mean=decon1)
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon_mean),HCbetween=mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  hc <- hc %>% group_by(id,run,trial,atlas_value) %>% arrange(evt_time) %>% mutate(HC_lag1 = lag(HCwithin,1),
                                                                                   HC_lag2 = lag(HCwithin,2),
                                                                                   HC_lag3 = lag(HCwithin,3)) %>% ungroup()
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
  hc <- hc %>% filter(atlas_value %in% c(7,8,9,10,11,12))
  gc()
  hc <- hc %>% rename(bin = atlas_value)
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory="~/clock_analysis",dataset='mmclock_fmri')
  df <- df %>%
    group_by(id) %>% 
    mutate(v_entropy_sc_r = scale(v_entropy)) %>% ungroup() %>%
    group_by(id, run) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           score_sc = scale(score_csv),
           iti_sc = scale(iti_ideal),
           iti_lag_sc = lag(iti_sc),
           v_chosen_sc = scale(v_chosen),
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
    run_trial <= 13 ~ 'Early',
    run_trial > 13 & run_trial < 26 ~ 'Middle',
    run_trial >=26 ~ 'Late',
  )))
  df <- df %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                       total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))
  
  df <- df %>% select(iti_ideal,rt_lag_sc,iti_lag_sc,iti_prev,iti_sc,last_outcome,rt_csv,rt_csv_sc,id,run,trial)
  df$id <- as.numeric(df$id)
  Q <- full_join(vmPFC,df,by=c('id','run','trial'))
  rm(vmPFC)
  gc()
  hc$id <- as.numeric(hc$id)
  Q <- inner_join(Q,hc,by=c('id','run','trial','evt_time'))
  rm(hc)
  gc()
  Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCwithin[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_decon[Q$evt_time < -Q$iti_prev] = NA
  Q$HC_lag1[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HC_lag1[Q$evt_time < -Q$iti_prev] = NA
  Q$HC_lag2[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HC_lag2[Q$evt_time < -Q$iti_prev] = NA
  Q$HC_lag3[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HC_lag3[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_lag1[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_lag1[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_lag2[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_lag2[Q$evt_time < -Q$iti_prev] = NA
  Q <- Q %>% arrange(id,run,trial,evt_time)
  Q <- Q %>% filter(evt_time > -5 & evt_time < 5)
  
  
  rm(decode_formula)
  decode_formula <- NULL
  #decode_formula[[1]] = formula(~age * HCwithin + gender * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + v_max_wi * HCwithin + v_entropy_wi_change_lag * HCwithin + rt_csv_sc * HCwithin + iti_lag_sc * HCwithin + iti_sc * HCwithin + last_outcome * HCwithin + outcome*HCwithin + rt_vmax_change_sc * HCwithin +  HCbetween + (1 | id/run)) 
  #decode_formula[[1]] = formula(~age*HCwithin + v_entropy_wi*HCwithin + (1|id/run))
  #decode_formula[[2]] = formula(~age*HCwithin + v_max_wi*HCwithin + (1|id/run))
  #decode_formula[[3]] = formula(~age*HCwithin + gender*HCwithin + v_entropy_wi*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1|id/run))
  decode_formula <- NULL
  decode_formula[[1]] <- formula(~ HCwithin + HCbetween + (1|id/run))
  # decode_formula[[2]] <- formula(~ HCwithin + HC_lag1 + HCbetween + (1 | id/run))
  decode_formula[[2]] <- formula(~ HC_lag1 + HCbetween + (1 | id/run))
  decode_formula[[3]] <- formula(~ vmPFC_lag1 + HCwithin + (1|id/run))
  decode_formula[[4]] <- formula(~ HCwithin + (1|id/run))
  # decode_formula[[4]] <- formula(~ HCwithin + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[5]] <- formula(~ HCwithin + HC_lag1 + HC_lag2 +  HCbetween + (1 | id/run))
  # decode_formula[[6]] <- formula(~ HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[7]] <- formula(~ HCwithin + HC_lag1 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[8]] <- formula(~ HCwithin + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[9]] <- formula(~ HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[10]] <- formula(~ HC_lag1 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[11]] <- formula(~ vmPFC_lag1 + HCwithin + HCbetween + (1|id/run))
  # decode_formula[[12]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag1 + HCbetween + (1 | id/run))
  # decode_formula[[13]] <- formula(~ vmPFC_lag1 + HC_lag1 + HCbetween + (1 | id/run))
  # decode_formula[[14]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[15]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag1 + HC_lag2 +  HCbetween + (1 | id/run))
  # decode_formula[[16]] <- formula(~ vmPFC_lag1 + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[17]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag1 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[18]] <- formula(~ vmPFC_lag1 + HCwithin + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[19]] <- formula(~ vmPFC_lag1 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[20]] <- formula(~ vmPFC_lag1 + HC_lag1 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[21]] <- formula(~ vmPFC_lag2 + HC_lag1 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[22]] <- formula(~ vmPFC_lag2 + HCwithin + HCbetween + (1|id/run))
  # decode_formula[[23]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag1 + HCbetween + (1 | id/run))
  # decode_formula[[24]] <- formula(~ vmPFC_lag2 + HC_lag1 + HCbetween + (1 | id/run))
  # decode_formula[[25]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[26]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag1 + HC_lag2 +  HCbetween + (1 | id/run))
  # decode_formula[[27]] <- formula(~ vmPFC_lag2 + HC_lag2 + HCbetween + (1 | id/run))
  # decode_formula[[28]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag1 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[29]] <- formula(~ vmPFC_lag2 + HCwithin + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[30]] <- formula(~ vmPFC_lag2 + HC_lag2 + HC_lag3 + HCbetween + (1 | id/run))
  # decode_formula[[31]] <- formula(~ vmPFC_lag2 + HC_lag1 + HC_lag3 + HCbetween + (1 | id/run))# 
  
  
  qT <- c(-0.8,0.46)
  splits = c('evt_time','network','bin')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 3:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-MMClock-vPFC-HC-network-clock-anatomy-12slice-7to12-',i,'.Rdata'))
  }
  rm(Q)
  gc()  
}

