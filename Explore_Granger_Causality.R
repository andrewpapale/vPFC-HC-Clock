
library(stringr)

ncores = 26
do_network = TRUE
do_symmetry = FALSE
do_vPFC_HC = TRUE


if (do_vPFC_HC){
  
  #source('/Volumes/Users/Andrew/MEDuSA_data_Explore/get_trial_data_explore.R')
  
  load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')
  
  md <- md %>% group_by(id,run,trial,atlas_value) %>% arrange(evt_time) %>% mutate(vmPFC_lag1 = dplyr::lag(decon_mean,1,order_by=evt_time),
                                                                                   vmPFC_lag2 = dplyr::lag(decon_mean,2, order_by=evt_time),
                                                                                   vmPFC_lag3 = dplyr::lag(decon_mean,3, order_by=evt_time)) %>% ungroup()
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
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock.Rdata')
  #hc <- read_csv('/Volumes/Users/Bea/StriatumHippThalamus/clock_aligned_striatum_hipp_thalamus.csv.gz')
  #hc <- hc %>% mutate(run1 = as.integer(str_sub(run,4,4))) %>% select(!run) %>% rename(run=run1)
  #hc <- hc %>% filter(atlas_value %in% c(223,224,225,226,227,228,229,230))
  #hc #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_fb.Rdata')
  #hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  hc <- hc %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                        trial > 40 & trial <= 80 ~ 2,
                                        trial > 80 & trial <=120 ~ 3, 
                                        trial > 120 & trial <=160 ~ 4,
                                        trial > 160 & trial <=200 ~ 5,
                                        trial > 200 & trial <=240 ~ 6))
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
  hc <- hc %>% select(!atlas_value)
  hc <- hc %>% group_by(id,run,trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  hc <- hc %>% rename(decon_mean=decon1)
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon_mean),HCbetween=mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  hc <- hc %>% group_by(id,run,trial,HC_region) %>% arrange(evt_time) %>% mutate(HC_lag1 = dplyr::lag(HCwithin,1,order_by=evt_time),
                                                                                 HC_lag2 = lag(HCwithin,2, order_by=evt_time),
                                                                                 HC_lag3 = lag(HCwithin,3, order_by=evt_time)) %>% ungroup()
  hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory='/Volumes/Users/Andrew/MEDuSA_data_Explore',dataset='explore')
  df <- df %>%
    group_by(id,run_number) %>%
    arrange(id, run_number, trial) %>%
    mutate(condition_trial = run_trial-floor(run_trial/40.5)*40,
           condition_trial_neg_inv = -1000 / condition_trial,
           condition_trial_neg_inv_sc = as.vector(scale(condition_trial_neg_inv)),
    ) %>% ungroup()
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
  
  df <- df %>% select(total_earnings_split,condition_trial_neg_inv_sc,iti_ideal,rt_lag_sc,iti_lag_sc,iti_prev,iti_sc,v_entropy_wi_change_lag,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
  df$id <- as.character(df$id)
  Q <- full_join(md,df,by=c('id','run','trial'))
  Q <- Q %>% rename(vmPFC_decon = decon_mean) %>% select(!decon_median & !decon_sd)
  rm(md)
  hc$id <- as.character(hc$id)
  Q <- Q %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                      trial > 40 & trial <= 80 ~ 2,
                                      trial > 80 & trial <=120 ~ 3, 
                                      trial > 120 & trial <=160 ~ 4,
                                      trial > 160 & trial <=200 ~ 5,
                                      trial > 200 & trial <=240 ~ 6))
  Q <- Q %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                           trial > 40 & trial <= 80 ~ trial-40,
                                           trial > 80 & trial <=120 ~ trial-80, 
                                           trial > 120 & trial <=160 ~ trial-120,
                                           trial > 160 & trial <=200 ~ trial-160,
                                           trial > 200 & trial <=240 ~ trial-200))
  Q <- Q %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                    run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                    run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))
  Q <- inner_join(Q,hc,by=c('id','run','trial','evt_time'))
  rm(hc)
  Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCwithin[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_decon[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_lag1[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_lag1[Q$evt_time < -Q$iti_prev] = NA
  Q$HC_lag1[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HC_lag1[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_lag2[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_lag2[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_lag3[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_lag3[Q$evt_time < -Q$iti_prev] = NA
  Q$HC_lag2[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HC_lag2[Q$evt_time < -Q$iti_prev] = NA
  Q$HC_lag3[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HC_lag3[Q$evt_time < -Q$iti_prev] = NA
  Q <- Q %>% mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM'
  )) %>% mutate(symmetry_group=case_when(
    atlas_value==67 | atlas_value==171 ~ 6,
    atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 1,
    atlas_value==89 | atlas_value==194 ~ 7,
    atlas_value==88 | atlas_value==192 ~ 5,
    atlas_value==84 | atlas_value==191 ~ 4,
    atlas_value==86 | atlas_value==161 ~ 2,
    atlas_value==55 | atlas_value==160 ~ 8,
    atlas_value==56 | atlas_value==159 ~ 3
  ))
  Q <- Q %>% arrange(id,run,trial,evt_time)
  
  Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
  Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
  
  demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
  demo <- demo %>% select(registration_redcapid,age,gender,registration_group,wtar,education_yrs)
  demo <- demo %>% rename(id=registration_redcapid,group=registration_group)
  demo$gender <- relevel(as.factor(demo$gender),ref='M')
  demo$age <- scale(demo$age)
  demo$wtar <- scale(demo$wtar)
  demo$education_yrs <- scale(demo$education_yrs)
  
  Q <- merge(demo,Q,by='id')
  #Q <- Q %>% filter(total_earnings_split=='richer')
  Q$group <- relevel(factor(Q$group),ref='HC')
  #Q <- Q %>% filter(group!='ATT')
  Q <- Q %>% filter(group=='HC')
  #Q <- Q %>% filter(!is.na(rewFunc))
  #Q <- Q %>% filter(trial > 10)
  
  rm(decode_formula)
  decode_formula <- NULL
  #decode_formula[[1]] = formula(~age * HCwithin + female * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + v_max_wi * HCwithin + v_entropy_wi_change_lag * HCwithin + rt_csv_sc * HCwithin + iti_lag_sc * HCwithin + iti_sc * HCwithin + last_outcome * HCwithin + rt_vmax_change_sc * HCwithin + HCwithin * HCbetween + (1 | id/run)) 
  #decode_formula[[3]] = formula(~ age + female + abs_pe_max_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin  + rt_csv_sc*HCwithin  + iti_lag_sc*HCwithin + iti_sc + HCwithin*HCbetween + (1| id/run))
  #decode_formula[[1]] = formula(~ v_max_wi*HCwithin + HCbetween + (1|id/run))
  #decode_formula[[2]] = formula(~ v_entropy_wi*HCwithin + HCbetween + (1|id/run))
  #decode_formula[[1]] = formula(~v_max_wi + (1|id/run))
  #decode_formula[[1]] = formula(~v_max_wi*HCwithin + v_entropy_wi*HCwithin + HCbetween + (1|id/run))
  #decode_formula[[1]] = formula(~age * HCwithin + female * HCwithin + v_entropy_wi * HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run))    
  #decode_formula[[2]] = formula(~HC_lag1 + age * HCwithin + female * HCwithin + v_entropy_wi * HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run))    
  decode_formula[[1]] = formula(~ HCwithin + vmPFC_lag1 + (1 | id/run))    
  decode_formula[[2]] = formula(~ HCwithin + vmPFC_lag2 + (1 | id/run))    
  decode_formula[[3]] = formula(~ HCwithin + vmPFC_lag3 + (1 | id/run))
  decode_formula[[4]] = formula(~ HCwithin + vmPFC_lag1*v_entropy_wi + (1 | id/run))    
  decode_formula[[5]] = formula(~ HCwithin + vmPFC_lag2*v_entropy_wi + (1 | id/run))    
  decode_formula[[6]] = formula(~ HCwithin + vmPFC_lag3*v_entropy_wi + (1 | id/run))
  decode_formula[[7]] = formula(~ vmPFC_lag1 + (1|id/run))
  decode_formula[[8]] = formula(~ vmPFC_lag2 + (1|id/run))
  decode_formula[[9]] = formula(~ vmPFC_lag3 + (1|id/run))
  decode_formula[[10]] = formula(~ vmPFC_lag1*v_entropy_wi + (1|id/run))
  decode_formula[[11]] = formula(~ vmPFC_lag2*v_entropy_wi + (1|id/run))
  decode_formula[[12]] = formula(~ vmPFC_lag3*v_entropy_wi + (1|id/run))
  #decode_formula[[4]] = formula(~HC_lag1 + HC_lag2 + age * HCwithin + female * HCwithin + v_entropy_wi * HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run))    
  qT <- c(-0.7,0.43)
  
  Q <- Q %>% filter(evt_time > -4 & evt_time < 4)
  splits = c('evt_time','network','HC_region')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits, 
                    calculate = c("parameter_estimates_reml","residuals","fit_statistics"),
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-Explore-Granger-Causality-pred-vmPFC-network-clock-',i,'.Rdata'))
  }
  
  
  rm(decode_formula)
  decode_formula <- NULL
  decode_formula[[1]] = formula(~ vmPFC_decon + HC_lag1 + (1 | id/run))    
  decode_formula[[2]] = formula(~ vmPFC_decon + HC_lag2 + (1 | id/run))    
  decode_formula[[3]] = formula(~ vmPFC_decon + HC_lag3 + (1 | id/run))
  decode_formula[[4]] = formula(~ vmPFC_decon + HC_lag1*v_entropy_wi + (1 | id/run))    
  decode_formula[[5]] = formula(~ vmPFC_decon + HC_lag2*v_entropy_wi + (1 | id/run))    
  decode_formula[[6]] = formula(~ vmPFC_decon + HC_lag3*v_entropy_wi + (1 | id/run))
  decode_formula[[7]] = formula(~ HC_lag1 + (1 | id/run))
  decode_formula[[8]] = formula(~ HC_lag2 + (1 | id/run))
  decode_formula[[9]] = formula(~ HC_lag3 + (1 | id/run))
  decode_formula[[10]] = formula(~ HC_lag1*v_entropy_wi + (1 | id/run))
  decode_formula[[11]] = formula(~ HC_lag2*v_entropy_wi + (1 | id/run))
  decode_formula[[12]] = formula(~ HC_lag3*v_entropy_wi + (1 | id/run))
  qT <- c(-0.7,0.43)
  
  Q <- Q %>% filter(evt_time > -4 & evt_time < 4)
  splits = c('evt_time','network','HC_region')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "HCwithin", rhs_model_formulae = df0 , split_on = splits, 
                    calculate = c("parameter_estimates_reml","residuals","fit_statistics"),
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-Explore-Granger-Causality-pred-HC-network-clock-',i,'.Rdata'))
  }
  
  
}
