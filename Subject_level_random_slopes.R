############
# 2022-05-23 AndyP
# Subject level random slopes

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26
toalign <- 'clock'
do_rand_slopes = FALSE
simple_model = FALSE
do_rt_pred_fmri = TRUE
plot_rt_pred_fmri = FALSE
do_rt_pred_meg = TRUE
plot_rt_pred_meg = FALSE
do_entropy_plot = FALSE
do_value_plot = FALSE
remove_1to10 = FALSE
do_vif = FALSE
#### clock ####

if (do_rand_slopes){
  
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  
  if (strcmp(toalign,'feedback')){
    load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
    vmPFC <- fb_comb
    vmPFC <- vmPFC %>% filter(evt_time > -4 & evt_time < 4)
    rm(fb_comb)
    vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
    vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
    load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
    hc <- fb_comb
    hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
    rm(fb_comb)
  } else if (strcmp(toalign,'clock')){
    load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMclock_clock_Aug2023.Rdata')
    vmPFC <- clock_comb
    vmPFC <- vmPFC %>% filter(evt_time > -4 & evt_time < 4)
    rm(clock_comb)
    vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
    vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
    
    load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
    #hc <- clock_comb
    hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
    #rm(clock_comb)
    
    hc <- hc %>% select(id,run,run_trial,decon_mean,evt_time,bin_num,side, HC_region)
  }
  
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
  gc()
  Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
  rm(hc,vmPFC)
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
  if (strcmp(toalign,'feedback')){
    df <- df %>% select(outcome,ev,iti_ideal,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                        v_max_wi, rt_lag_sc,v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                        rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
      group_by(id, run) %>% 
      mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
      mutate(v_chosen_sc = scale(v_chosen),
             abs_pe_max_sc = scale(abs(pe_max)),
             score_sc = scale(score_csv),
             iti_sc = scale(iti_ideal),
             iti_lag_sc = scale(iti_lag),
             pe_max_sc = scale(pe_max),
             pe_max_lag_sc = scale(lag(pe_max)),
             abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
             rt_vmax_sc = scale(rt_vmax),
             ev_sc = scale(ev),
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
    df <- df %>% select(id,run,trial,run_trial,trial_bin,outcome,rewFunc,trial_neg_inv_sc,rt_csv_sc,v_entropy_sc,last_outcome,v_max_wi,trial_neg_inv_sc,rt_csv_sc,v_entropy_wi_change,score_sc,rt_bin,iti_sc,iti_lag_sc,rt_lag_sc,ev_sc,expl_longer,expl_shorter)
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
    
    Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
    
    if (remove_1to10){
      Q <- Q %>% filter(trial > 10)
    }
    
    rm(decode_formula)
    decode_formula <- formula(~ (1|id))
    decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + last_outcome*HCwithin + HCbetween + (1 + HCwithin*v_entropy_wi |id) + (1|run))
    decode_formula[[2]] = formula(~ age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + v_max_wi*HCwithin  + rt_lag_sc*HCwithin + last_outcome*HCwithin + HCbetween + (1 + HCwithin*v_max_wi |id) + (1|run))
    
  } else if (strcmp(toalign,'clock')){
    df <- df %>% select(ev,iti_ideal,iti_prev,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                        v_max_wi, iti_prev,rt_lag_sc,v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                        rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
      group_by(id, run) %>% 
      mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>%
      mutate(v_chosen_sc = scale(v_chosen),
             abs_pe_max_sc = scale(abs(pe_max)),
             iti_lag_sc = scale(iti_lag),
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
    df <- df %>% select(id,run,run_trial,rt_lag_sc,iti_lag,trial_neg_inv_sc,rt_csv_sc,v_entropy_wi,expl_longer,expl_shorter,rt_bin,trial_bin,last_outcome,v_max_wi,v_entropy_wi_change_lag,score_lag_sc,iti_lag_sc,ev_lag_sc)
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
    Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
    Q$vmPFC_decon[Q$evt_time < -(Q$iti_lag)] = NA;
    Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_lag] = NA;
    Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_lag] = NA;
    Q$HCwithin[Q$evt_time < -(Q$iti_lag)] = NA;
    Q$HCbetween[Q$evt_time < -(Q$iti_lag)] = NA;
    Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
    Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
    
    rm(decode_formula)
    decode_formula <- formula(~ (1|id))
    #decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + iti_lag_sc*HCwithin + last_outcome*HCwithin + HCbetween + (1 + HCwithin*v_entropy_wi |id) + (1|run))
    #decode_formula[[2]] = formula(~ age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + v_max_wi*HCwithin + rt_lag_sc*HCwithin  + iti_lag_sc*HCwithin + last_outcome*HCwithin + HCbetween + (1 + HCwithin*v_max_wi  |id) + (1|run))
    #decode_formula[[3]] = formula(~ age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + v_max_wi*HCwithin + v_entropy_wi*HCwithin + rt_lag_sc*HCwithin  + iti_lag_sc*HCwithin + last_outcome*HCwithin + HCbetween + (1 + HCwithin  |id) + (1|run))
  
    decode_formula[[1]] = formula(~ age + female + trial_neg_inv_sc + v_entropy_wi + rt_lag_sc + iti_lag_sc + last_outcome + HCbetween + (1 + HCwithin*v_entropy_wi |id) + (1|run))
    decode_formula[[2]] = formula(~ age + female + trial_neg_inv_sc + v_max_wi + rt_lag_sc  + iti_lag_sc + last_outcome + HCbetween + (1 + HCwithin*v_max_wi  |id) + (1|run))
    decode_formula[[3]] = formula(~ age + female + trial_neg_inv_sc + v_max_wi + v_entropy_wi + rt_lag_sc  + iti_lag_sc + last_outcome + HCbetween + (1 + HCwithin  |id) + (1|run))
    
  }
  
  
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
    save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-',toalign,'-ranslopes-nofixedeffect-',i,'.Rdata'))
  }
}

if (do_rt_pred_fmri){
  for (i in 1:3){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vmPFC-HC-network-',toalign,'-ranslopes-nofixedeffect-',i,'.Rdata')
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
        qdf <- ddf$coef_df_reml %:% filter(effect=='ran_vals' & term=='HCwithin')
      } 
    }
    #qdf <- qdf %>% group_by(network,HC_region) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
    qdf <- qdf %>% rename(id=level)
    qdf <- qdf %>% group_by(id,network,HC_region) %>% summarize(estimate = mean(estimate,na.rm=TRUE)) %>% ungroup()
    qdf <- qdf %>% group_by(network,HC_region) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
    qdf$id <- as.character(qdf$id)
    #qdf <- qdf %>% select(!outcome)
    source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
    df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
    df <- df %>% select(rt_vmax_lag,iti_prev,ev,iti_ideal,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                        v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                        rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
      group_by(id, run) %>% 
      mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% 
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
             v_max_wi_lag = lag(v_max_wi),
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
    
    df$id <- as.character(df$id)
    Q2 <- inner_join(qdf,df,by='id')
    
    Q2$trial_bin <- factor(Q2$trial_bin)
    Q2$trial_bin <- relevel(Q2$trial_bin, ref='Middle')
    Q2 <- Q2 %>% rename(subj_level_rand_slope=estimate)
    #~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
    #  rt_lag:last_outcome:fmri_beta +
    #  rt_vmax_lag*trial_neg_inv*fmri_beta +
    #  (1 | id/run)
    decode_formula <- NULL
    #decode_formula[[1]] = formula(~ rt_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 | id/run))
    #decode_formula[[2]] = formula(~ rt_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 + rt_lag_sc | id/run))
    #decode_formula[[3]] = formula(~ rt_vmax_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 | id/run))  
    #decode_formula[[4]] = formula(~ rt_vmax_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 + rt_vmax_lag_sc | id/run))   
    
    if (!simple_model){
      decode_formula[[1]] <- formula(~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope * trial_neg_inv_sc + (1 | id/run))
      #decode_formula[[2]] <- formula(~(trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * trial_neg_inv_sc * subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc | id/run))
    } else if (simple_model){
      decode_formula[[1]] <- formula(~trial_neg_inv_sc + last_outcome*subj_level_rand_slope + rt_lag_sc*subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope + (1 | id/run))
      #decode_formula[[2]] <- formula(~trial_neg_inv_sc + last_outcome*subj_level_rand_slope + rt_lag_sc*subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc | id/run))
    }
    qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
    qRTV <- quantile(df$rt_vmax_lag,c(0.1,0.9),na.rm=TRUE)
    qH <- NULL
    qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
    qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    #qH <- c(1,2,3,4)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
    qRS <- quantile(Q2$estimate, c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    splits = c('HC_region','network')
    source('~/fmri.pipeline/R/mixed_by.R')
    print(i)
    for (j in 1:length(decode_formula)){
      if (!simple_model){
        ddq <- mixed_by(Q2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                        emmeans_spec = list(
                          RT = list(outcome='rt_csv_sc', model_name='model1', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),        
                          TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                                          specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                          
                          ),
                        emtrends_spec = list(
                          RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
                                          specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                          TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
                                           specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                          TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'subj_level_rand_slope',
                                           specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                          )
        )
        setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
        curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
        if (j==1){
          # AH
          zq <- Q2 %>% filter(network=='DMN' & HC_region =='AH')
          m_dmn_ah <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_dmn_ah <- anova(m_dmn_ah)
          zq <- Q2 %>% filter(network=='CTR' & HC_region =='AH')
          m_ctr_ah <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_ctr_ah <- anova(m_ctr_ah)
          zq <- Q2 %>% filter(network=='LIM' & HC_region =='AH')
          m_lim_ah <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_lim_ah <- anova(m_lim_ah)
          # PH
          zq <- Q2 %>% filter(network=='DMN' & HC_region =='PH')
          m_dmn_ph <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_dmn_ph <- anova(m_dmn_ph)
          zq <- Q2 %>% filter(network=='CTR' & HC_region =='PH')
          m_ctr_ph <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_ctr_ph <- anova(m_ctr_ph)
          zq <- Q2 %>% filter(network=='LIM' & HC_region =='PH')
          m_lim_ph <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_lim_ph <- anova(m_lim_ph)
          save(am_dmn_ah,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-dmn_ah-anova-',i,'.Rdata'))
          save(am_ctr_ah,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-ctr_ah-anova-',i,'.Rdata'))
          save(am_lim_ah,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-lim_ah-anova-',i,'.Rdata'))
          save(am_dmn_ph,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-dmn_ph-anova-',i,'.Rdata'))
          save(am_ctr_ph,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-ctr_ph-anova-',i,'.Rdata'))
          save(am_lim_ph,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-lim_ph-anova-',i,'.Rdata'))
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-',i,'.Rdata'))
        } else {
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-slo-notimesplit-nofixedeffect-',i,'.Rdata')) 
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
                m1 <- lmerTest::lmer(data=temp, rt_csv ~(trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * trial_neg_inv_sc * subj_level_rand_slope + (1 | id/run))
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
                          RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                        )
        )
        setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
        curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
        if (j==1){
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-int-simple-',i,'.Rdata'))
        } else {
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-slo-simple-',i,'.Rdata')) 
        }       
      }
    }
  }
}

if (plot_rt_pred_fmri){
  source('~/vmPFC/plot_subject_level_random_slopes.R')
  source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
  for (i in 3){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
    model_str <- paste0('-vmPFC-HC-network-ranslopes-',toalign,'-pred-rt_csv_sc-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    model_iter <- i
    totest <- paste0('rt_csv_sc-',as.character(i))
    #toprocess <- 'symmetry-by-HC'
    toprocess <- 'network-by-HC'
    behavmodel <- 'compressed'
    hc_LorR <- 'LR'
    plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    #plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  }
}

# replication

### Does random slope predict rt_vmax or rt_swing?
if (do_rt_pred_meg) {
  
  for (i in 1:3){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vmPFC-HC-network-',toalign,'-ranslopes-nofixedeffect-',i,'.Rdata')
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
      }  else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')
      } 
    }
    #qdf <- qdf %>% group_by(network,HC_region) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
    qdf <- qdf %>% rename(id=level)
    qdf$id <- as.character(qdf$id)
    qdf <- qdf %>% group_by(id,network,HC_region) %>% summarize(estimate = mean(estimate,na.rm=TRUE)) %>% ungroup()
    qdf <- qdf %>% group_by(network,HC_region) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
    #qdf <- qdf %>% select(!outcome)
    source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
    df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_meg')
    df <- df %>% select(rt_vmax_lag,iti_prev,ev,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                        v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,
                        rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
      group_by(id, run) %>% 
      mutate(rt_sec = rt_csv/1000) %>% 
      mutate(v_chosen_sc = scale(v_chosen),
             abs_pe_max_sc = scale(abs(pe_max)),
             score_sc = scale(score_csv),
             pe_max_sc = scale(pe_max),
             pe_max_lag_sc = scale(lag(pe_max)),
             abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
             rt_vmax_sc = scale(rt_vmax),
             ev_sc = scale(ev),
             v_entropy_sc = scale(v_entropy),
             v_max_wi_lag = lag(v_max_wi),
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
    
    df$id <- as.character(df$id)
    Q2 <- inner_join(qdf,df,by='id')
    
    Q2$trial_bin <- factor(Q2$trial_bin)
    Q2$trial_bin <- relevel(Q2$trial_bin, ref='Middle')
    Q2 <- Q2 %>% rename(subj_level_rand_slope=estimate)
    #~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
    #  rt_lag:last_outcome:fmri_beta +
    #  rt_vmax_lag*trial_neg_inv*fmri_beta +
    #  (1 | id/run)
    decode_formula <- NULL
    if (!simple_model){
      decode_formula[[1]] <- formula(~( rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope * trial_neg_inv_sc + (1 | id/run))
      #decode_formula[[2]] <- formula(~(trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * trial_neg_inv_sc * subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc | id/run))
    } else if (simple_model){
      decode_formula[[1]] <- formula(~trial_neg_inv_sc + last_outcome*subj_level_rand_slope + rt_lag_sc*subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope + (1 | id/run))
      decode_formula[[2]] <- formula(~trial_neg_inv_sc + last_outcome*subj_level_rand_slope + rt_lag_sc*subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc | id/run))
    }
    qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
    qRTV <- quantile(df$rt_vmax_lag,c(0.1,0.9),na.rm=TRUE)
    qH <- NULL
    qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
    qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    #qH <- c(1,2,3,4)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
    qRS <- quantile(Q2$estimate, c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    splits = c('HC_region','network')
    source('~/fmri.pipeline/R/mixed_by.R')
    print(i)
    for (j in 1:length(decode_formula)){
      if (!simple_model){
        ddq <- mixed_by(Q2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                        emmeans_spec = list(
                          RT = list(outcome='rt_csv_sc', model_name='model1', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),        
                          TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                                          specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                          
                        ),
                        emtrends_spec = list(
                          RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
                                          specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                          TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
                                           specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                          TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'subj_level_rand_slope',
                                           specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                        )
        )
        setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
        curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
        if (j==1){
          # AH
          zq <- Q2 %>% filter(network=='DMN' & HC_region =='AH')
          m_dmn_ah <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_dmn_ah <- anova(m_dmn_ah)
          zq <- Q2 %>% filter(network=='CTR' & HC_region =='AH')
          m_ctr_ah <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_ctr_ah <- anova(m_ctr_ah)
          zq <- Q2 %>% filter(network=='LIM' & HC_region =='AH')
          m_lim_ah <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_lim_ah <- anova(m_lim_ah)
          # PH
          zq <- Q2 %>% filter(network=='DMN' & HC_region =='PH')
          m_dmn_ph <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_dmn_ph <- anova(m_dmn_ph)
          zq <- Q2 %>% filter(network=='CTR' & HC_region =='PH')
          m_ctr_ph <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_ctr_ph <- anova(m_ctr_ph)
          zq <- Q2 %>% filter(network=='LIM' & HC_region =='PH')
          m_lim_ph <- lmer(rt_csv_sc ~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope*trial_neg_inv_sc + (1 | id/run), data=zq)
          am_lim_ph <- anova(m_lim_ph)
          save(am_dmn_ah,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-dmn_ah-anova-',i,'.Rdata'))
          save(am_ctr_ah,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-ctr_ah-anova-',i,'.Rdata'))
          save(am_lim_ah,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-lim_ah-anova-',i,'.Rdata'))
          save(am_dmn_ph,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-dmn_ph-anova-',i,'.Rdata'))
          save(am_ctr_ph,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-ctr_ph-anova-',i,'.Rdata'))
          save(am_lim_ph,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-lim_ph-anova-',i,'.Rdata'))
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-',i,'.Rdata'))
        } else {
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-slo-notimesplit-nofixedeffect-',i,'.Rdata')) 
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
                          RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                        )
        )
        setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
        curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
        if (j==1){
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-int-simple-',i,'.Rdata'))
        } else {
          save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-replication-pred-slo-simple-',i,'.Rdata')) 
        }       
      }
    }
  }
}

if (plot_rt_pred_meg){
  source('~/vmPFC/plot_subject_level_random_slopes.R')
  source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
  for (i in 1:2){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
    model_str <- paste0('-vmPFC-HC-network-ranslopes-',toalign,'-pred-rt_csv_sc-replication-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    model_iter <- i
    totest <- paste0('rt_csv_sc-replication-',as.character(i))
    #toprocess <- 'symmetry-by-HC'
    toprocess <- 'network-by-HC'
    behavmodel <- 'compressed'
    hc_LorR <- 'LR'
    plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  }
}

if (do_entropy_plot){
  
  #########################
  ###  Entropy - rt_lag ###
  #########################
  
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-pred-rt_csv_sc-rt_lag-1.Rdata')
  ddf1 <- ddq$coef_df_reml
  ddf1 <- ddf1 %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope')
  ddf1 <- ddf1 %>% mutate(dataset='fMRI')
  emt1 <- ddq$emtrends_list$RT
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-replication-pred-rt_csv_sc-rt_lag-1.Rdata')
  ddf2 <- ddq$coef_df_reml
  ddf2 <- ddf2 %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope')
  ddf2 <- ddf2 %>% mutate(dataset='Replication')
  emt2 <- ddq$emtrends_list$RT
  emt1 <- emt1 %>% mutate(dataset='fMRI')
  emt2 <- emt2 %>% mutate(dataset='Replication')
  emt <- rbind(emt1,emt2)
  ddf <- rbind(ddf1,ddf2)
  emt <- emt %>% filter(subj_level_rand_slope==-2 | subj_level_rand_slope==2)
  ddq <- merge(emt,ddf,by=c('evt_time','network','HC_region','dataset'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddqA <- ddq %>% filter(HC_region=='AH')
  ddqP <- ddq %>% filter(HC_region=='PH')
  pdf('Entropy-rt_lag-AH.pdf',height=3.5,width=9)
  gg1 <- ggplot(ddqA,aes(x=evt_time,y=rt_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error.x,ymax=rt_lag_sc.trend+std.error.x))
  print(gg1)
  dev.off()
  pdf('Entropy-rt_lag-PH.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqP,aes(x=evt_time,y=rt_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error.x,ymax=rt_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  
}


if (do_value_plot){
  
  #######################
  ###  Value - rt_lag ###
  #######################
  
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-pred-rt_csv_sc-rt_lag-1.Rdata')
  ddf1 <- ddq$coef_df_reml
  ddf1 <- ddf1 %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope')
  ddf1 <- ddf1 %>% mutate(dataset='fMRI')
  emt1 <- ddq$emtrends_list$RT
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-replication-pred-rt_csv_sc-rt_lag-1.Rdata')
  ddf2 <- ddq$coef_df_reml
  ddf2 <- ddf2 %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope')
  ddf2 <- ddf2 %>% mutate(dataset='Replication')
  emt2 <- ddq$emtrends_list$RT
  emt1 <- emt1 %>% mutate(dataset='fMRI')
  emt2 <- emt2 %>% mutate(dataset='Replication')
  emt <- rbind(emt1,emt2)
  ddf <- rbind(ddf1,ddf2)
  emt <- emt %>% filter(subj_level_rand_slope==-2 | subj_level_rand_slope==2)
  ddq <- merge(emt,ddf,by=c('evt_time','network','HC_region','dataset'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddqA <- ddq %>% filter(HC_region=='AH')
  ddqP <- ddq %>% filter(HC_region=='PH')
  pdf('Value-rt_lag-AH.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqA,aes(x=evt_time,y=rt_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error.x,ymax=rt_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  pdf('Value-rt_lag-PH.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqP,aes(x=evt_time,y=rt_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error.x,ymax=rt_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  
}


if (do_entropy_plot){
  
  ##########################
  ###  Entropy - rt_vmax ###
  ##########################
  
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-pred-rt_csv_sc-rt_vmax_lag-1.Rdata')
  ddf1 <- ddq$coef_df_reml
  ddf1 <- ddf1 %>% filter(effect=='fixed' & term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddf1 <- ddf1 %>% mutate(dataset='fMRI')
  emt1 <- ddq$emtrends_list$Vmax
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-replication-pred-rt_csv_sc-rt_vmax_lag-1.Rdata')
  ddf2 <- ddq$coef_df_reml
  ddf2 <- ddf2 %>% filter(effect=='fixed' & term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddf2 <- ddf2 %>% mutate(dataset='Replication')
  emt2 <- ddq$emtrends_list$Vmax
  emt1 <- emt1 %>% mutate(dataset='fMRI')
  emt2 <- emt2 %>% mutate(dataset='Replication')
  emt <- rbind(emt1,emt2)
  ddf <- rbind(ddf1,ddf2)
  emt <- emt %>% filter(subj_level_rand_slope==-2 | subj_level_rand_slope==2)
  ddq <- merge(emt,ddf,by=c('evt_time','network','HC_region','dataset'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddqA <- ddq %>% filter(HC_region=='AH')
  ddqP <- ddq %>% filter(HC_region=='PH')
  pdf('Entropy-rt_vmax-AH.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqA,aes(x=evt_time,y=rt_vmax_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_vmax_lag_sc.trend-std.error.x,ymax=rt_vmax_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  pdf('Entropy-rt_vmax-PH.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqP,aes(x=evt_time,y=rt_vmax_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_vmax_lag_sc.trend-std.error.x,ymax=rt_vmax_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  
}


if (do_value_plot){
  
  ########################
  ###  Value - rt_vmax ###
  ########################
  
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-pred-rt_csv_sc-rt_vmax_lag-2.Rdata')
  ddf1 <- ddq$coef_df_reml
  ddf1 <- ddf1 %>% filter(effect=='fixed' & term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddf1 <- ddf1 %>% mutate(dataset='fMRI')
  emt1 <- ddq$emtrends_list$Vmax
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-replication-pred-rt_csv_sc-rt_vmax_lag-2.Rdata')
  ddf2 <- ddq$coef_df_reml
  ddf2 <- ddf2 %>% filter(effect=='fixed' & term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddf2 <- ddf2 %>% mutate(dataset='Replication')
  emt2 <- ddq$emtrends_list$Vmax
  emt1 <- emt1 %>% mutate(dataset='fMRI')
  emt2 <- emt2 %>% mutate(dataset='Replication')
  emt <- rbind(emt1,emt2)
  ddf <- rbind(ddf1,ddf2)
  emt <- emt %>% filter(subj_level_rand_slope==-2 | subj_level_rand_slope==2)
  ddq <- merge(emt,ddf,by=c('evt_time','network','HC_region','dataset'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddqA <- ddq %>% filter(HC_region=='AH')
  ddqP <- ddq %>% filter(HC_region=='PH')
  pdf('Value-rt_vmax-AH.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqA,aes(x=evt_time,y=rt_vmax_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_vmax_lag_sc.trend-std.error.x,ymax=rt_vmax_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  pdf('Value-rt_vmax-PH.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqP,aes(x=evt_time,y=rt_vmax_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_vmax_lag_sc.trend-std.error.x,ymax=rt_vmax_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  
}

if (do_entropy_plot){
  
  #########################
  ###  Entropy - rt_lag ###
  #########################
  
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-ranslope-pred-rt_csv_sc-rt_lag-1.Rdata')
  ddf1 <- ddq$coef_df_reml
  ddf1 <- ddf1 %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope')
  ddf1 <- ddf1 %>% mutate(dataset='fMRI')
  emt1 <- ddq$emtrends_list$RT
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-replication-ranslope-pred-rt_csv_sc-rt_lag-1.Rdata')
  ddf2 <- ddq$coef_df_reml
  ddf2 <- ddf2 %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope')
  ddf2 <- ddf2 %>% mutate(dataset='Replication')
  emt2 <- ddq$emtrends_list$RT
  emt1 <- emt1 %>% mutate(dataset='fMRI')
  emt2 <- emt2 %>% mutate(dataset='Replication')
  emt <- rbind(emt1,emt2)
  ddf <- rbind(ddf1,ddf2)
  emt <- emt %>% filter(subj_level_rand_slope==-2 | subj_level_rand_slope==2)
  ddq <- merge(emt,ddf,by=c('evt_time','network','HC_region','dataset'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddqA <- ddq %>% filter(HC_region=='AH')
  ddqP <- ddq %>% filter(HC_region=='PH')
  pdf('Entropy-rt_lag-AH-ranslope.pdf',height=3.5,width=9)
  gg1 <- ggplot(ddqA,aes(x=evt_time,y=rt_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error.x,ymax=rt_lag_sc.trend+std.error.x))
  print(gg1)
  dev.off()
  pdf('Entropy-rt_lag-PH-ranslope.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqP,aes(x=evt_time,y=rt_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error.x,ymax=rt_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  
}


if (do_value_plot){
  
  #######################
  ###  Value - rt_lag ###
  #######################
  
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-ranslope-pred-rt_csv_sc-rt_lag-1.Rdata')
  ddf1 <- ddq$coef_df_reml
  ddf1 <- ddf1 %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope')
  ddf1 <- ddf1 %>% mutate(dataset='fMRI')
  emt1 <- ddq$emtrends_list$RT
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-replication-ranslope-pred-rt_csv_sc-rt_lag-1.Rdata')
  ddf2 <- ddq$coef_df_reml
  ddf2 <- ddf2 %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope')
  ddf2 <- ddf2 %>% mutate(dataset='Replication')
  emt2 <- ddq$emtrends_list$RT
  emt1 <- emt1 %>% mutate(dataset='fMRI')
  emt2 <- emt2 %>% mutate(dataset='Replication')
  emt <- rbind(emt1,emt2)
  ddf <- rbind(ddf1,ddf2)
  emt <- emt %>% filter(subj_level_rand_slope==-2 | subj_level_rand_slope==2)
  ddq <- merge(emt,ddf,by=c('evt_time','network','HC_region','dataset'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddqA <- ddq %>% filter(HC_region=='AH')
  ddqP <- ddq %>% filter(HC_region=='PH')
  pdf('Value-rt_lag-AH-ranslope.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqA,aes(x=evt_time,y=rt_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error.x,ymax=rt_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  pdf('Value-rt_lag-PH-ranslope.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqP,aes(x=evt_time,y=rt_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error.x,ymax=rt_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  
}


if (do_entropy_plot){
  
  ##########################
  ###  Entropy - rt_vmax ###
  ##########################
  
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-ranslope-pred-rt_csv_sc-rt_vmax_lag-1.Rdata')
  ddf1 <- ddq$coef_df_reml
  ddf1 <- ddf1 %>% filter(effect=='fixed' & term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddf1 <- ddf1 %>% mutate(dataset='fMRI')
  emt1 <- ddq$emtrends_list$Vmax
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-replication-ranslope-pred-rt_csv_sc-rt_vmax_lag-1.Rdata')
  ddf2 <- ddq$coef_df_reml
  ddf2 <- ddf2 %>% filter(effect=='fixed' & term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddf2 <- ddf2 %>% mutate(dataset='Replication')
  emt2 <- ddq$emtrends_list$Vmax
  emt1 <- emt1 %>% mutate(dataset='fMRI')
  emt2 <- emt2 %>% mutate(dataset='Replication')
  emt <- rbind(emt1,emt2)
  ddf <- rbind(ddf1,ddf2)
  emt <- emt %>% filter(subj_level_rand_slope==-2 | subj_level_rand_slope==2)
  ddq <- merge(emt,ddf,by=c('evt_time','network','HC_region','dataset'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddqA <- ddq %>% filter(HC_region=='AH')
  ddqP <- ddq %>% filter(HC_region=='PH')
  pdf('Entropy-rt_vmax-AH-ranslope.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqA,aes(x=evt_time,y=rt_vmax_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_vmax_lag_sc.trend-std.error.x,ymax=rt_vmax_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  pdf('Entropy-rt_vmax-PH-ranslope.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqP,aes(x=evt_time,y=rt_vmax_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_vmax_lag_sc.trend-std.error.x,ymax=rt_vmax_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  
}


if (do_value_plot){
  
  ########################
  ###  Value - rt_vmax ###
  ########################
  
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-ranslope-pred-rt_csv_sc-rt_vmax_lag-2.Rdata')
  ddf1 <- ddq$coef_df_reml
  ddf1 <- ddf1 %>% filter(effect=='fixed' & term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddf1 <- ddf1 %>% mutate(dataset='fMRI')
  emt1 <- ddq$emtrends_list$Vmax
  load('2023-01-10-vmPFC-HC-network-ranslopes-clock-replication-ranslope-pred-rt_csv_sc-rt_vmax_lag-2.Rdata')
  ddf2 <- ddq$coef_df_reml
  ddf2 <- ddf2 %>% filter(effect=='fixed' & term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddf2 <- ddf2 %>% mutate(dataset='Replication')
  emt2 <- ddq$emtrends_list$Vmax
  emt1 <- emt1 %>% mutate(dataset='fMRI')
  emt2 <- emt2 %>% mutate(dataset='Replication')
  emt <- rbind(emt1,emt2)
  ddf <- rbind(ddf1,ddf2)
  emt <- emt %>% filter(subj_level_rand_slope==-2 | subj_level_rand_slope==2)
  ddq <- merge(emt,ddf,by=c('evt_time','network','HC_region','dataset'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddqA <- ddq %>% filter(HC_region=='AH')
  ddqP <- ddq %>% filter(HC_region=='PH')
  pdf('Value-rt_vmax-AH-ranslope.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqA,aes(x=evt_time,y=rt_vmax_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_vmax_lag_sc.trend-std.error.x,ymax=rt_vmax_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  pdf('Value-rt_vmax-PH-ranslope.pdf',height=3.5,width=9)
  gg1<-ggplot(ddqP,aes(x=evt_time,y=rt_vmax_lag_sc.trend)) + 
    facet_grid(network~dataset) +
    geom_line(aes(group=subj_level_rand_slope,color=subj_level_rand_slope)) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr,color=subj_level_rand_slope)) + 
    geom_errorbar(aes(ymin=rt_vmax_lag_sc.trend-std.error.x,ymax=rt_vmax_lag_sc.trend+std.error.x))
  print(gg1); dev.off()
  
}


