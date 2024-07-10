# 2023-06-16 AndyP
# subject level random slopes of HC*entropy and HC*value for Explore dataset

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
ncores <- 26
toalign <- 'clock'
do_rand_slopes = TRUE
do_rt_pred_fmri = TRUE
simple_model = FALSE
trial_mod = TRUE
do_vif = FALSE

#### clock ####

if (do_rand_slopes){
  
  load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')
  md <- md %>% filter(evt_time > -4 & evt_time < 4)
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
  #hc <- hc %>% select(!decon_median & !decon_sd)
  #hc #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_fb.Rdata')
  #hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
  # hc <- hc %>% mutate(block = case_when(trial <= 40 ~ 1, 
  #                                       trial > 40 & trial <= 80 ~ 2,
  #                                       trial > 80 & trial <=120 ~ 3, 
  #                                       trial > 120 & trial <=160 ~ 4,
  #                                       trial > 160 & trial <=200 ~ 5,
  #                                       trial > 200 & trial <=240 ~ 6))
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
           v_max_wi_lag = lag(v_max_wi),
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
  Q <- inner_join(Q,hc,by=c('id','run','trial','evt_time'))
  rm(hc)
  Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCwithin[Q$evt_time < -Q$iti_prev] = NA
  Q$vmPFC_decon[Q$evt_time < -Q$iti_prev] = NA
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
  Q <- Q %>% filter(evt_time > -4 & evt_time < 4)
  
  
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
  Q <- Q %>% filter(!is.na(rewFunc))
  
  Q <- Q %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                           trial > 40 & trial <= 80 ~ trial-40,
                                           trial > 80 & trial <=120 ~ trial-80, 
                                           trial > 120 & trial <=160 ~ trial-120,
                                           trial > 160 & trial <=200 ~ trial-160,
                                           trial > 200 & trial <=240 ~ trial-200))
  Q <- Q %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                    run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                    run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv))
  )
  
  rm(decode_formula)
  decode_formula <- formula(~ (1|id))
  decode_formula[[1]] = formula(~ age + gender + v_entropy_wi + run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + last_outcome + HCbetween + (1 + HCwithin*v_entropy_wi |id) + (1|run))
  decode_formula[[2]] = formula(~ age + gender + run_trial0_neg_inv_sc + v_max_wi + rt_lag_sc  + iti_lag_sc + last_outcome + HCbetween + (1 + HCwithin*v_max_wi  |id) + (1|run))
  decode_formula[[3]] = formula(~ age + gender + run_trial0_neg_inv_sc + v_max_wi + v_entropy_wi + rt_lag_sc  + iti_lag_sc + last_outcome + HCbetween + (1 + HCwithin  |id) + (1|run))
  
  splits = c('evt_time','network','HC_region')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,return_models=TRUE,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-',toalign,'-Explore-ranslopes-HConly-trial_mod-trial1-10included-nofixedeffect-',i,'.Rdata'))
  }
}

if (do_rt_pred_fmri){
  for (i in 1:3){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vmPFC-HC-network-',toalign,'-Explore-ranslopes-HConly-trial_mod-trial1-10included-nofixedeffect-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    
    ### Does random slope predict rt_vmax or rt_swing?
    rm(Q2)
    if (strcmp(toalign,'clock')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin:v_entropy_wi')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin:v_max_wi')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')
      }
    } else if (strcmp(toalign,'feedback')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin:v_entropy_wi')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin:v_max_wi')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')
      }
    }
    #qdf <- qdf %>% group_by(network,HC_region) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
    qdf <- qdf %>% rename(id=level)
    qdf <- qdf %>% group_by(id,network,HC_region) %>% summarize(estimate = mean(estimate,na.rm=TRUE)) %>% ungroup()
    qdf <- qdf %>% group_by(network,HC_region) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
    qdf$id <- as.character(qdf$id)
    #qdf <- qdf %>% select(!outcome)
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
             v_max_wi_lag = lag(v_max_wi),
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
    
    df <- df %>% select(total_earnings_split,v_max_wi_lag,rt_vmax_lag_sc,condition_trial_neg_inv_sc,iti_ideal,rt_lag_sc,iti_lag_sc,iti_prev,iti_sc,v_entropy_wi_change_lag,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
    df$id <- as.character(df$id)
    Q2 <- inner_join(qdf,df,by='id')
    
    Q2$trial_bin <- factor(Q2$trial_bin)
    Q2$trial_bin <- relevel(Q2$trial_bin, ref='Middle')
    Q2 <- Q2 %>% rename(subj_level_rand_slope=estimate)
    Q2 <- Q2 %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                             trial > 40 & trial <= 80 ~ trial-40,
                                             trial > 80 & trial <=120 ~ trial-80, 
                                             trial > 120 & trial <=160 ~ trial-120,
                                             trial > 160 & trial <=200 ~ trial-160,
                                             trial > 200 & trial <=240 ~ trial-200))
    Q2 <- Q2 %>% mutate(run_trial0_c = run_trial-floor(run_trial/40.5)*40,
                        run_trial0_neg_inv = -1000 / run_trial0_c,
                        run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv))
    )
    #~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
    #  rt_lag:last_outcome:fmri_beta +
    #  rt_vmax_lag*trial_neg_inv*fmri_beta +
    #  (1 | id/run)
    decode_formula <- NULL
    #decode_formula[[1]] = formula(~ rt_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 | id/run))
    #decode_formula[[2]] = formula(~ rt_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 + rt_lag_sc | id/run))
    #decode_formula[[3]] = formula(~ rt_vmax_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 | id/run))  
    #decode_formula[[4]] = formula(~ rt_vmax_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 + rt_vmax_lag_sc | id/run))   
    if (!simple_model && !trial_mod){
      decode_formula[[1]] <- formula(~(condition_trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * condition_trial_neg_inv_sc * subj_level_rand_slope + (1 | id/run))
      decode_formula[[2]] <- formula(~(condition_trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * condition_trial_neg_inv_sc * subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc | id/run))
    } else if (simple_model){
      decode_formula[[1]] <- formula(~condition_trial_neg_inv_sc + rt_lag_sc*subj_level_rand_slope + last_outcome*subj_level_rand_slope + rt_vmax_lag_sc*subj_level_rand_slope + (1|id/run))
      decode_formula[[2]] <- formula(~condition_trial_neg_inv_sc + rt_lag_sc*subj_level_rand_slope + last_outcome*subj_level_rand_slope + rt_vmax_lag_sc*subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc |id/run))
    } else if (trial_mod){
      decode_formula[[1]] <- formula(~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope + (1 | id/run))
      #decode_formula[[2]] <- formula(~(run_trial0_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * run_trial0_neg_inv_sc * subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc | id/run))
    }
    qH <- NULL
    qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
    qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    #qH <- c(1,2,3,4)
    qT <- quantile(df$condition_trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
    splits = c('HC_region','network')
    source('~/fmri.pipeline/R/mixed_by.R')
    print(i)
    for (j in 1:length(decode_formula)){
      if (!simple_model && !trial_mod){
        ddq <- mixed_by(Q2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                        emmeans_spec = list(
                          RT = list(outcome='rt_csv_sc', model_name='model1', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2)))        
                          
                        ),
                        emtrends_spec = list(
                          RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                          
                        )
        )
        setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
        curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
        if (j==1){
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-Explore-ranslopes-',toalign,'-pred-int-trial1-10included-',i,'.Rdata'))
        } else {
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-Explore-ranslopes-',toalign,'-pred-slo-trial1-10included-',i,'.Rdata')) 
        }
      } else if (simple_model){
        ddq <- mixed_by(Q2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                        emmeans_spec = list(
                          RT = list(outcome='rt_csv_sc', model_name='model1', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2)))
                        ),
                        emtrends_spec = list(
                          RT = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                        )
        )
        setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
        curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
        if (j==1){
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-Explore-ranslopes-',toalign,'-pred-int-trial1-10included-simple-',i,'.Rdata'))
        } else {
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-Explore-ranslopes-',toalign,'-pred-slo-trial1-10included-simple-',i,'.Rdata')) 
        }  
      } else if (trial_mod){
        ddq <- mixed_by(Q2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                        emmeans_spec = list(
                          RT = list(outcome='rt_csv_sc', model_name='model1', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2)))#,
                          # TrxVmax = list(outcome='rt_csv',model_name='model1',
                          #                specs=formula(~rt_vmax_lag_sc:run_trial0_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2),run_trial0_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                          
                        ),
                        emtrends_spec = list(
                          RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2)))#,
                          # TrxVmax = list(outcome='rt_csv',model_name='model1', var = 'rt_vmax_lag_sc',
                          #                specs=formula(~rt_vmax_lag_sc:run_trial0_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),run_trial0_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                          # TrxVmax1 = list(outcome='rt_csv',model_name='model1', var = 'run_trial0_neg_inv_sc',
                          #                 specs=formula(~rt_vmax_lag_sc:run_trial0_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                          # TrxVmax2 = list(outcome='rt_csv',model_name='model1', var = 'subj_level_rand_slope',
                          #                 specs=formula(~rt_vmax_lag_sc:run_trial0_neg_inv_sc:subj_level_rand_slope), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),run_trial0_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                        )
        )
        setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
        curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
        if (j==1){
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-Explore-ranslopes-',toalign,'-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-',i,'.Rdata'))
        } else {
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-Explore-ranslopes-',toalign,'-pred-slo-trial_mod-trial1-10included-notimesplit-nofixedeffect-',i,'.Rdata')) 
        }
        
        if (do_vif){
          nE = unique(Q2$evt_time)
          nHC = unique(Q2$HC_region)
          nN = unique(Q2$network)
          vdf <- NULL
          network <- NULL
          HC_region <- NULL
          evt_time <- NULL
          for (iE in 1:length(nE)){
            for (iHC in 1:length(nHC)){
              for (iN in 1:length(nN)){
                  temp <- Q2 %>% filter(evt_time == nE[iE])
                  temp <- temp %>% filter(HC_region == nHC[iHC])
                  temp <- temp %>% filter(network == nN[iN])
                  m1 <- lmerTest::lmer(data=temp, rt_csv ~ (rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope + (1 | id/run))
                  v0 <- car::vif(m1)
                  vdf <- rbind(vdf,v0)
                  network <- rbind(network,nN[iN])
                  HC_region <- rbind(HC_region,nHC[iHC])
                  evt_time <- rbind(evt_time,nE[iE])
              }
            }
          }
          vdf1 <- data.frame(vdf,network=network,HC_regio=HC_region,evt_time=evt_time)
        }
        
        
      }
    }
  }
}

