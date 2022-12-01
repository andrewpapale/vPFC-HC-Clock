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
do_rt_pred_fmri = FALSE
plot_rt_pred_fmri = TRUE
do_rt_pred_meg = FALSE
plot_rt_pred_meg = TRUE
#### clock ####

if (do_rand_slopes){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  
  if (strcmp(toalign,'feedback')){
    load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
    vmPFC <- fb_comb
    vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
    rm(fb_comb)
    vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
    vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
    #load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
    #hc <- fb_comb
    #hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
    #rm(fb_comb)
  } else if (strcmp(toalign,'clock')){
    load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
    vmPFC <- clock_comb
    vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
    rm(clock_comb)
    vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
    vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
    #load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
    #hc <- clock_comb
    #hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
    #rm(clock_comb)
  }
  #hc <- hc %>% mutate(
  #  HC_region = case_when(
  #    bin_num <= 8 ~ 'AH',
  #    bin_num >8 ~ 'PH'
  #  ),
  #)
  #hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
  #hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
  
  #Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
  Q <- vmPFC
  rm(vmPFC)
  
  source('~/vmPFC/get_trial_data_vmPFC.R')
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
  
  if (strcmp(toalign,'feedback')){
    df <- df %>% select(v_max_wi_lag,ev,iti_ideal,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
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
    df <- df %>% select(id,run,run_trial,trial_bin,rt_vmax_change_sc,iti_lag_sc,rewFunc,trial_neg_inv_sc,rt_csv_sc,v_entropy_sc,last_outcome,v_max_wi,trial_neg_inv_sc,rt_csv_sc,v_entropy_wi_change,score_sc,rt_bin,iti_sc,ev_sc,expl_longer,expl_shorter)
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
    
    #Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
    
    rm(decode_formula)
    decode_formula <- formula(~ (1|id))
    #decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_sc*HCwithin + trial_neg_inv_sc*HCwithin + v_max_wi*HCwithin  + v_entropy_wi_change*HCwithin   + rt_csv_sc  + iti_sc + last_outcome*HCwithin + HCbetween + (1 + HCwithin*v_entropy_sc |id/run))
    #decode_formula[[2]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_sc*HCwithin + trial_neg_inv_sc*HCwithin + v_max_wi*HCwithin  + v_entropy_wi_change*HCwithin   + rt_csv_sc  + iti_sc + last_outcome*HCwithin + HCbetween + (1 + HCwithin*v_max_wi |id/run))
    decode_formula[[1]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change_lag + trial_neg_inv_sc + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + rt_vmax_change_sc + (1 + v_entropy_sc |id/run))
    decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change_lag + trial_neg_inv_sc + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + rt_vmax_change_sc + (1 + v_max_wi |id/run))
    decode_formula[[3]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change_lag + trial_neg_inv_sc + last_outcome + rt_csv_sc + iti_sc + iti_lag_sc + rt_vmax_change_sc + (1 + v_entropy_wi_change_lag |id/run))
    
  } else if (strcmp(toalign,'clock')){
    df <- df %>% select(v_max_wi_lag,ev,iti_ideal,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                        v_max_wi, v_entropy_wi, v_entropy_wi_change,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                        rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
      group_by(id, run) %>% 
      mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
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
    df <- df %>% select(id,run,run_trial,trial_bin,iti_sc,v_max_wi,v_entropy_wi_change_lag,rt_vmax_change_sc,rewFunc,trial_neg_inv_sc,rt_csv_sc,v_entropy_sc,expl_longer,expl_shorter,rt_bin,trial_bin,outcome,score_lag_sc,iti_lag_sc,ev_lag_sc)
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
    
    #Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
    
    rm(decode_formula)
    decode_formula <- formula(~ (1|id))
    #decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin +  v_max_wi_lag*HCwithin  + v_entropy_wi_change_lag*HCwithin  + rt_csv_sc  + iti_lag_sc + last_outcome*HCwithin + HCwithin + HCbetween + (1 + HCwithin*v_entropy_lag_sc |id/run))
    #decode_formula[[2]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + v_max_wi_lag*HCwithin  + v_entropy_wi_change_lag*HCwithin + rt_csv_sc  + iti_lag_sc + last_outcome*HCwithin + HCwithin + HCbetween + (1 + HCwithin*v_max_wi_lag  |id/run))
    decode_formula[[1]] = formula(~ age + female + v_entropy_sc + v_entropy_wi_change_lag + outcome + trial_neg_inv_sc + v_max_wi + rt_csv_sc + iti_sc + rt_vmax_change_sc + (1 + v_entropy_sc |id/run))
    decode_formula[[2]] = formula(~ age + female + v_entropy_sc + v_entropy_wi_change_lag + outcome + trial_neg_inv_sc + v_max_wi + rt_csv_sc + iti_sc + rt_vmax_change_sc + (1 + v_max_wi |id/run))
    decode_formula[[3]] = formula(~ age + female + v_entropy_sc + v_entropy_wi_change_lag + outcome + trial_neg_inv_sc + v_max_wi + rt_csv_sc + iti_sc + rt_vmax_change_sc + (1 + v_entropy_wi_change_lag |id/run))
  }
  
  
  splits = c('evt_time','network')
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
    save(ddf,file=paste0(curr_date,'-vmPFC-network-',toalign,'-ranslopes-',i,'.Rdata'))
  }
}

if (do_rt_pred_fmri){
  for (i in 1:3){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vmPFC-network-',toalign,'-ranslopes-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    
    ### Does random slope predict rt_vmax or rt_swing?
    rm(Q2)
    if (strcmp(toalign,'clock')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_sc' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change_lag' & group=='id')
      }
    } else if (strcmp(toalign,'feedback')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_sc' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change' & group=='id')
      } 
    }
    qdf <- qdf %>% group_by(network) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
    qdf <- qdf %>% rename(id=level)
    qdf$id <- as.character(qdf$id)
    qdf <- qdf %>% select(!outcome)
    source('~/vmPFC/get_trial_data_vmPFC.R')
    df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
    df <- df %>% select(rt_vmax_lag,v_max_wi_lag,ev,iti_ideal,score_csv,v_max,last_outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
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
    
    #~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
    #  rt_lag:last_outcome:fmri_beta +
    #  rt_vmax_lag*trial_neg_inv*fmri_beta +
    #  (1 | id/run)
    decode_formula = formula(~ (trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_sc + subj_level_rand_slope + last_outcome)^2 +
                               rt_lag_sc*subj_level_rand_slope*last_outcome + rt_vmax_lag_sc*subj_level_rand_slope*trial_neg_inv_sc + (1|id/run))
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
    ddq <- mixed_by(Q2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula, split_on = splits,return_models=TRUE,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed"),conf.int=TRUE),
                    emmeans_spec = list(
                      LO = list(outcome='rt_csv_sc',model_name='model1',
                                specs=formula(~rt_lag_sc:subj_level_rand_slope:last_outcome),at=list(subj_level_rand_slope=qRS)),
                      TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                                     specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope:trial_neg_inv_sc),at=list(trial_neg_inv_sc=qT,subj_level_rand_slope=qRS)),
                      RT = list(outcome='rt_csv_sc', model_name='model1', 
                                specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=qRS)),
                      Vmax = list(outcome='rt_csv_sc', model_name='model1',
                                  specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at=list(subj_level_rand_slope=qRS))
                    ),
                    emtrends_spec = list(
                      LO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                specs=formula(~rt_lag_sc:subj_level_rand_slope:last_outcome),at=list(subj_level_rand_slope=qRS)),
                      TrxVmax = list(outcome='rt_csv_sc',model_name='model1',var='rt_vmax_lag_sc',
                                     specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope:trial_neg_inv_sc),at=list(trial_neg_inv_sc=qT,subj_level_rand_slope=qRS)),
                      RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                                specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=qRS)),
                      Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc',
                                  specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at=list(subj_level_rand_slope=qRS))
                    )
    )
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-',i,'.Rdata'))
  }
}
if (plot_rt_pred_fmri){
  source('~/vmPFC/plot_subject_level_random_slopes.R')
  source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
  for (i in 1:3){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
    model_str <- paste0('-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    model_iter <- i
    totest <- paste0('rt_csv_sc-',as.character(i))
    #toprocess <- 'symmetry'
    toprocess <- 'network'
    behavmodel <- 'compressed'
    hc_LorR <- 'LR'
    plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  }
}
  # #test rt_swing
  # splits = c('evt_time','HC_region','network')
  # ddq <- mixed_by(Q2, outcomes = "rt_swing", rhs_model_formulae = decode_formula[[1]] , split_on = splits,
  #                 padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
  #                 tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
  #                 emtrends_spec = list(
  #                   V_HC = list(outcome='rt_swing', model_name='model1', var='estimate', 
  #                               specs='v_max_wi', at = list(v_max_wi=qV)),
  #                   T_HC = list(outcome='rt_swing', model_name='model1', var='estimate', 
  #                               specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
  #                   O_HC = list(outcome='rt_swing', model_name='model1', var='estimate',
  #                               specs=c('last_outcome'))
  #                 )
  # )
  # setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  # curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  # save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-pred-rt_swing-',i,'.Rdata'))
  # 
  # 
  # source('~/vmPFC/plot_subject_level_random_slopes.R')
  # source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
  # for (i in 1){
  #   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  #   #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
  #   model_str <- paste0('-vmPFC-HC-network-ranslopes-pred-rt_csv_sc-1','.Rdata')
  #   model_str <- Sys.glob(paste0('*',model_str))
  #   load(model_str)
  #   model_iter <- i
  #   totest <- 'rt_csv-HCwithin'
  #   #toprocess <- 'symmetry-by-HC'
  #   toprocess <- 'network-by-HC'
  #   toalign <- 'clock'
  #   behavmodel <- 'compressed'
  #   hc_LorR <- 'LR'
  #   plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  #   plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  # }
  
  
  # replication
  
  ### Does random slope predict rt_vmax or rt_swing?
if (do_rt_pred_meg){
  for (i in 1:3){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vmPFC-network-',toalign,'-ranslopes-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    
    ### Does random slope predict rt_vmax or rt_swing?
    rm(Q2)
    if (strcmp(toalign,'clock')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_sc' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change_lag' & group=='id')
      }
    } else if (strcmp(toalign,'feedback')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_sc' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change' & group=='id')
      } 
    }
    qdf <- qdf %>% group_by(network) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
    qdf <- qdf %>% rename(id=level)
    qdf$id <- as.character(qdf$id)
    qdf <- qdf %>% select(!outcome)
    
    source('~/vmPFC/get_trial_data_vmPFC.R')
    df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_meg')
    df <- df %>% select(rt_vmax_lag,v_max_wi_lag,ev,score_csv,v_max,last_outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
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
    Q2$trial_bin <- factor(Q2$trial_bin)
    Q2$trial_bin <- relevel(Q2$trial_bin, ref='Middle')
    Q2 <- Q2 %>% rename(subj_level_rand_slope=estimate)
    
    decode_formula = formula(~ (trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_sc + subj_level_rand_slope + last_outcome)^2 +
                               rt_lag_sc*subj_level_rand_slope*last_outcome + rt_vmax_lag_sc*subj_level_rand_slope*trial_neg_inv_sc + (1|id/run))
    qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
    qRTV <- quantile(df$rt_vmax_lag,c(0.1,0.9),na.rm=TRUE)
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
    ddq <- mixed_by(Q2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula, split_on = splits,return_models=TRUE,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed"),conf.int=TRUE),
                    emmeans_spec = list(
                      LO = list(outcome='rt_csv_sc',model_name='model1',
                                specs=formula(~rt_lag_sc:subj_level_rand_slope:last_outcome),at=list(subj_level_rand_slope=qRS)),
                      TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                                     specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope:trial_neg_inv_sc),at=list(trial_neg_inv_sc=qT,subj_level_rand_slope=qRS)),
                      RT = list(outcome='rt_csv_sc', model_name='model1', 
                                specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=qRS)),
                      Vmax = list(outcome='rt_csv_sc', model_name='model1',
                                  specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at=list(subj_level_rand_slope=qRS))
                    ),
                    emtrends_spec = list(
                      LO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                specs=formula(~rt_lag_sc:subj_level_rand_slope:last_outcome),at=list(subj_level_rand_slope=qRS)),
                      TrxVmax = list(outcome='rt_csv_sc',model_name='model1',var='rt_vmax_lag_sc',
                                     specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope:trial_neg_inv_sc),at=list(trial_neg_inv_sc=qT,subj_level_rand_slope=qRS)),
                      RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                                specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=qRS)),
                      Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc',
                                  specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at=list(subj_level_rand_slope=qRS))
                    )
    )
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-replication-',i,'.Rdata'))
  }
}
if (plot_rt_pred_meg){
  source('~/vmPFC/plot_subject_level_random_slopes.R')
  source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
  for (i in 1:3){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
    model_str <- paste0('-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-replication-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    model_iter <- i
    totest <- paste0('rt_csv_sc-replication-',as.character(i))
    #toprocess <- 'symmetry-by'
    toprocess <- 'network-by'
    behavmodel <- 'compressed'
    hc_LorR <- 'LR'
    plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  }
}
# #test rt_swing
# splits = c('evt_time','HC_region','network')
# ddq <- mixed_by(Q2, outcomes = "rt_swing", rhs_model_formulae = decode_formula[[1]] , split_on = splits,
#                 padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                 tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
#                 emtrends_spec = list(
#                   V_HC = list(outcome='rt_swing', model_name='model1', var='estimate', 
#                               specs='v_max_wi', at = list(v_max_wi=qV)),
#                   T_HC = list(outcome='rt_swing', model_name='model1', var='estimate', 
#                               specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
#                   O_HC = list(outcome='rt_swing', model_name='model1', var='estimate',
#                               specs=c('last_outcome'))
#                 )
# )
# setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
# curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
# save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-pred-rt_swing-replication-',i,'.Rdata'))
# 
# 
# source('~/vmPFC/plot_subject_level_random_slopes.R')
# source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
# for (i in 1){
#   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
#   model_str <- paste0('-vmPFC-HC-network-ranslopes-pred-rt_swing-replication-',i,'.Rdata')
#   model_str <- Sys.glob(paste0('*',model_str))
#   load(model_str)
#   model_iter <- i
#   totest <- 'rt_swing-replication'
#   #toprocess <- 'symmetry-by-HC'
#   toprocess <- 'network-by-HC'
#   toalign <- 'clock'
#   behavmodel <- 'compressed'
#   hc_LorR <- 'LR'
#   plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
#   plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
# }
# 
# rm(decode_formula)
# decode_formula <- formula(~ (1|id/run))
# decode_formula[[1]] = formula(~ estimate + (1|id))
# decode_formula[[2]] = formula(~ estimate +  + (1|id))
# decode_formula[[3]] = formula(~ estimate +  + last_outcome + (1|id))
# decode_formula[[4]] = formula(~ estimate +  + last_outcome:estimate + (1|id))
# decode_formula[[5]] = formula(~ estimate + v_entropy_wi + (1|id))
# 
# Q2 <- Q2 %>% mutate(Swing_Category = case_when(rt_swing <= 0.333 ~ 'Exploitation', 
#                                                rt_swing > 0.333  ~ 'Exploration' 
#                                                ))
# Q2$Swing_Category <- as.factor(Q2$Swing_Category)
# Q2$Swing_Category <- relevel(Q2$Swing_Category,ref='Exploitation')
# source('~/vmPFC/mixed_by.R')
# for (i in 1:length(decode_formula)){
#   print(i)
#   ddq <- mixed_by(Q2, outcomes = "Swing_Category", rhs_model_formulae = decode_formula[[i]], 
#                   split_on = splits,return_models=TRUE,padjust_by="term",padjust_method="fdr",
#                   refit_on_nonconvergence = 0, calculate=c('parameter_estimates_ml','fit_statistics'),ncores=ncores,
#                   tidy_args = list(effects=c("coefficients")), doclmm=TRUE)
#   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-pred-Swing_Category-',i,'.Rdata'))
# }


