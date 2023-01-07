############
# 2022-05-23 AndyP
# Subject level random slopes

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(ggnewscale)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26
toalign <- 'clock'
do_rand_slopes = FALSE
do_rt_pred_fmri = FALSE
plot_rt_pred_fmri = FALSE
do_rt_pred_meg = FALSE
plot_rt_pred_meg = FALSE
do_entropy_plot = TRUE
do_value_plot = TRUE

#### clock ####

if (do_rand_slopes){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  
  if (strcmp(toalign,'feedback')){
    load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
    vmPFC <- fb_comb
    vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
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
    vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
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
    df <- df %>% select(rewFunc,v_max_wi_lag,ev,iti_ideal,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                        v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                        rt_vmax_lag_sc,iti_prev,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
      group_by(id, run) %>% 
      mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>%
      mutate(v_chosen_sc = scale(v_chosen),
             abs_pe_max_sc = scale(abs(pe_max)),
             score_sc = scale(score_csv),
             iti_sc = scale(iti_ideal),
             iti_lag_sc = scale(iti_prev),
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
    df <- df %>% select(id,run,run_trial,trial_bin,iti_prev_sc,v_entropy_wi_change,rt_vmax_change_sc,iti_lag_sc,rewFunc,trial_neg_inv_sc,rt_csv_sc,v_entropy_sc,outcome,v_max_wi,trial_neg_inv_sc,rt_csv_sc,v_entropy_wi_change,score_sc,rt_bin,iti_sc,ev_sc,expl_longer,expl_shorter)
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
    decode_formula[[1]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + outcome + last_outcome + rt_csv_sc + iti_sc + iti_prev_sc + (1 + v_entropy_wi | id) + (1|id/run))
    decode_formula[[2]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + outcome + last_outcome + rt_csv_sc + iti_sc + iti_prev_sc + (1 + v_max_wi |id) + (1 |  id/run))
    decode_formula[[3]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + outcome + last_outcome + rt_csv_sc + iti_sc + iti_prev_sc + (1 + v_entropy_wi | rewFunc) + (1|id/rewFunc))
    decode_formula[[4]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + outcome + last_outcome + rt_csv_sc + iti_sc + iti_prev_sc + (1 + v_max_wi | rewFunc) + (1 |  id/rewFunc))
    
    #decode_formula[[3]] = formula(~ age + female + v_entropy_sc + v_max_wi + v_entropy_wi_change + trial_neg_inv_sc + outcome + rt_csv_sc + iti_sc + iti_lag_sc + rt_vmax_change_sc + (1 + v_entropy_wi_change |id/run))
    
  } else if (strcmp(toalign,'clock')){
    df <- df %>% select(rewFunc,v_max_wi_lag,ev,iti_ideal,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                        v_max_wi, v_entropy_wi,iti_prev, v_entropy_wi_change,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                        rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
      group_by(id, run) %>% 
      mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
      mutate(v_chosen_sc = scale(v_chosen),
             abs_pe_max_sc = scale(abs(pe_max)),
             iti_prev_sc = scale(iti_prev),
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
    df <- df %>% select(id,run,run_trial,trial_bin,iti_sc,iti_prev_sc,outcome,v_max_wi,v_entropy_wi,rt_vmax_change_sc,rewFunc,trial_neg_inv_sc,rt_csv_sc,v_entropy_sc,expl_longer,expl_shorter,rt_bin,trial_bin,last_outcome,score_lag_sc,ev_lag_sc)
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
    
    Q$vmPFC_decon[Q$evt_time > Q$iti_ideal] = NA
    #Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
    
    rm(decode_formula)
    decode_formula <- formula(~ (1|id))
    #decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin +  v_max_wi_lag*HCwithin  + v_entropy_wi_change_lag*HCwithin  + rt_csv_sc  + iti_lag_sc + last_outcome*HCwithin + HCwithin + HCbetween + (1 + HCwithin*v_entropy_lag_sc |id/run))
    #decode_formula[[2]] = formula(~ age*HCwithin + female*HCwithin + v_entropy_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + v_max_wi_lag*HCwithin  + v_entropy_wi_change_lag*HCwithin + rt_csv_sc  + iti_lag_sc + last_outcome*HCwithin + HCwithin + HCbetween + (1 + HCwithin*v_max_wi_lag  |id/run))
    decode_formula[[1]] = formula(~ age + female + v_entropy_wi + last_outcome + outcome + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + (1 + v_entropy_wi |id) + (1|run))
    decode_formula[[2]] = formula(~ age + female + v_max_wi + last_outcome + outcome + trial_neg_inv_sc + rt_csv_sc + iti_sc + iti_prev_sc + (1 + v_max_wi |id) + (1|run))
    decode_formula[[3]] = formula(~ age + female + v_entropy_wi + trial_neg_inv_sc + outcome + last_outcome + rt_csv_sc + iti_sc + iti_prev_sc + (1 + v_entropy_wi | rewFunc) + (1|id/rewFunc))
    decode_formula[[4]] = formula(~ age + female + v_max_wi + trial_neg_inv_sc + outcome + last_outcome + rt_csv_sc + iti_sc + iti_prev_sc + (1 + v_max_wi | rewFunc) + (1 |  id/rewFunc))
    
    #decode_formula[[3]] = formula(~ age + female + v_entropy_sc + v_entropy_wi_change_lag + last_outcome + trial_neg_inv_sc + v_max_wi + rt_csv_sc + iti_sc + rt_vmax_change_sc + (1 + v_entropy_wi_change_lag |id/run))
  }
  
  
  splits = c('evt_time','network')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 3:length(decode_formula)){
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
  for (i in 1:2){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vmPFC-network-',toalign,'-ranslopes-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    
    ### Does random slope predict rt_vmax or rt_swing?
    rm(Q2)
    if (strcmp(toalign,'clock')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi' & group=='rewFunc')
        #qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change_lag' & group=='id')
      } else if (i==4){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='rewFunc')
      }
    } else if (strcmp(toalign,'feedback')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change' & group=='id')
      } 
    }
    qdf <- qdf %>% group_by(network) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate & !rhs) %>% rename(estimate=estimate1)
    qdf <- qdf %>% rename(id=level)
    qdf$id <- as.character(qdf$id)
    qdf <- qdf %>% select(!outcome)
    source('~/vmPFC/get_trial_data_vmPFC.R')
    df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
    df <- df %>% select(rewFunc,rt_vmax_lag,v_max_wi_lag,ev,score_csv,v_max,last_outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                        v_max_wi, iti_ideal, iti_prev,v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,
                        rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
      group_by(id, run) %>% 
      mutate(rt_sec = rt_csv/1000) %>% 
      mutate(v_chosen_sc = scale(v_chosen),
             abs_pe_max_sc = scale(abs(pe_max)),
             iti_sc = scale(iti_ideal),
             iti_prev_sc = scale(iti_prev),
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
             )) %>% ungroup()
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
    
    if (strcmp(toalign,'clock')){
      decode_formula <- NULL
      decode_formula[[1]] = formula(~ rt_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 | id/run))
      decode_formula[[2]] = formula(~ rt_vmax_lag_sc*subj_level_rand_slope + iti_sc + iti_prev_sc + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 | id/run))
    } else if (strcmp(toalign,'feedback')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + last_outcome + outcome + iti_sc + iti_prev_sc +
                                 rt_csv_sc*subj_level_rand_slope + rt_vmax_sc*subj_level_rand_slope*trial_neg_inv_sc + (1|id/run))     
    }
    qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
    qRTV <- quantile(df$rt_vmax_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    qH <- NULL
    qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
    qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    #qH <- c(1,2,3,4)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    qRS <- quantile(Q2$subj_level_rand_slope, c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    
    #write.csv(Q2,file=paste0('vPFC_network_ranslopes_pred_rt_csv','_',i,'.csv'))
    
    splits = c('evt_time','network')
    source('~/fmri.pipeline/R/mixed_by.R')
    print(i)
    for (j in 1:length(decode_formula)){
      print(j)
      if (j==1){
        ddq <- mixed_by(Q2, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                        emmeans_spec = list(
                          RT = list(outcome='rt_csv', model_name='model1', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2)))
                        ),
                        emtrends_spec = list(
                          RT = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                        )
        )
      } else if (j==2){
             ddq <- mixed_by(Q2, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                               padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                               tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                               emmeans_spec = list(
                                 Vmax = list(outcome='rt_csv', model_name='model1', 
                                           specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2)))
                               ),
                               emtrends_spec = list(
                                 Vmax = list(outcome='rt_csv', model_name='model1', var='rt_vmax_lag_sc', 
                                           specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                               )
      )
      }
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      if (j==1){
        save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-rt_lag-',i,'.Rdata'))
      } else if (j==2){
        save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-rt_vmax_lag-',i,'.Rdata'))
      }
  }
  }
}
if (plot_rt_pred_fmri){
  source('~/vmPFC/plot_subject_level_random_slopes.R')
  source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
  source('~/vmPFC/plot_emmeans_subject_level_random_slopes.R')
  for (i in 1:2){
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
    if (strcmp(toprocess,'network-by-HC')){
      if (strcmp(toalign,"clock")){
        setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction')
      } else if (strcmp(toalign,"feedback")){
        setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC_interaction')
      }
    } else if (strcmp(toprocess,'network')){
      if (strcmp(toalign,"clock")){
        setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock')
      } else if (strcmp(toalign,"feedback")){
        setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback')
      }
    }
    #plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    plot_emmeans_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  }
}


# replication

### Does random slope predict rt_vmax or rt_swing?
if (do_rt_pred_meg){
  for (i in 1:2){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vmPFC-network-',toalign,'-ranslopes-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    
    ### Does random slope predict rt_vmax or rt_swing?
    rm(Q2)
    if (strcmp(toalign,'clock')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change_lag' & group=='id')
      }
    } else if (strcmp(toalign,'feedback')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change' & group=='id')
      } 
    }
    qdf <- qdf %>% group_by(network) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate & !rhs) %>% rename(estimate=estimate1)
    qdf <- qdf %>% rename(id=level)
    qdf$id <- as.character(qdf$id)
    qdf <- qdf %>% select(!outcome)
    
    source('~/vmPFC/get_trial_data_vmPFC.R')
    df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_meg')
    df <- df %>% select(rewFunc,rt_vmax_lag,v_max_wi_lag,ev,score_csv,v_max,last_outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
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
      run_trial <= 21 ~ 'Early',
      run_trial > 21 & run_trial < 42 ~ 'Middle',
      run_trial >=42 ~ 'Late',
    )))
    
    Q2 <- merge(qdf,df,by=c('id'))
    Q2$trial_bin <- factor(Q2$trial_bin)
    Q2$trial_bin <- relevel(Q2$trial_bin, ref='Middle')
    Q2 <- Q2 %>% rename(subj_level_rand_slope=estimate)
    
    if (strcmp(toalign,'clock')){
      decode_formula <- NULL
      decode_formula[[1]] = formula(~ rt_lag_sc*subj_level_rand_slope + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 | id/run))
      decode_formula[[2]] = formula(~ rt_vmax_lag_sc*subj_level_rand_slope + last_outcome + outcome + v_entropy_wi + v_max_wi +  (1 | id/run))
    } else if (strcmp(toalign,'feedback')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + subj_level_rand_slope + last_outcome + outcome +
                                 rt_csv_sc*subj_level_rand_slope*outcome + rt_vmax_sc*subj_level_rand_slope*trial_neg_inv_sc + (1|id/run))      
    }
    qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
    qRTV <- quantile(df$rt_vmax_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    qH <- NULL
    qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    #qH <- c(1,2,3,4)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
    qRS <- quantile(Q2$subj_level_rand_slope, c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    
    #write.csv(Q2,file=paste0('vPFC_network_ranslopes_replication_pred_rt_csv','_',i,'.csv'))
    
    splits = c('evt_time','network')
    source('~/fmri.pipeline/R/mixed_by.R')
    print(i)
    for (j in 1:length(decode_formula)){
      print(j)
      if (j==1){
        ddq <- mixed_by(Q2, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                        emmeans_spec = list(
                          RT = list(outcome='rt_csv', model_name='model1', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2)))
                        ),
                        emtrends_spec = list(
                          RT = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                        )
        )
      } else if (j==2){
        ddq <- mixed_by(Q2, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                               padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                               tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                               emmeans_spec = list(
                                 Vmax = list(outcome='rt_csv', model_name='model1', 
                                           specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2)))
                               ),
                               emtrends_spec = list(
                                 Vmax = list(outcome='rt_csv', model_name='model1', var='rt_vmax_lag_sc', 
                                           specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                               )
      )
      }
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      if (j==1){
        save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-replication-',toalign,'-pred-rt_csv_sc-rt_lag-',i,'.Rdata'))
      } else if (j==2){
        save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-replication-',toalign,'-pred-rt_csv_sc-rt_vmax_lag-',i,'.Rdata'))
      }
    }
  }
}
  
if (plot_rt_pred_meg){
  source('~/vmPFC/plot_subject_level_random_slopes.R')
  source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
  source('~/vmPFC/plot_emmeans_subject_level_random_slopes.R')
  for (i in 1:2){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
    model_str <- paste0('-vmPFC-network-ranslopes-replication-',toalign,'-pred-rt_csv_sc-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    model_iter <- i
    totest <- paste0('rt_csv_sc-replication-',as.character(i))
    #toprocess <- 'symmetry-by'
    toprocess <- 'network'
    behavmodel <- 'compressed'
    hc_LorR <- 'LR'
    if (strcmp(toprocess,'network-by-HC')){
      if (strcmp(toalign,"clock")){
        setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction')
      } else if (strcmp(toalign,"feedback")){
        setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC_interaction')
      }
    } else if (strcmp(toprocess,'network')){
      if (strcmp(toalign,"clock")){
        setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock')
      } else if (strcmp(toalign,"feedback")){
        setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback')
      }
    }
    plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    plot_emmeans_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  }
}

if (do_entropy_plot){
  # plot nice figure
  library(wesanderson)
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1 = palette()
  pal1[1] <- pal[1]
  pal1[2] <- pal[2]
  pal1[3] <- pal[3]
  i <- 1 # entropy
  toalign <- 'clock' #
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-ranslopes-replication-',toalign,'-pred-rt_csv_sc-rt_vmax_lag-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_r <- ddq
  rm(ddq)
  model_str <- paste0('-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-rt_vmax_lag-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_f <- ddq
  rm(ddq)
  
  ddq_r_emm <- ddq_r$emmeans_list$Vmax
  ddq_r_emt <- ddq_r$emtrends_list$Vmax
  ddq_f_emm <- ddq_f$emmeans_list$Vmax
  ddq_f_emt <- ddq_f$emtrends_list$Vmax
  ddq_r <- ddq_r$coef_df_reml
  ddq_f <- ddq_f$coef_df_reml
  
  ddq_r <- ddq_r %>% filter(effect=='fixed') 
  ddq_f <- ddq_f %>% filter(effect=='fixed')
  
  ddq_r_emm <- ddq_r_emm %>% mutate(dataset = 'Replication') %>% rename(rt_vmax_lag_sc_r = rt_vmax_lag_sc, slrs_r=subj_level_rand_slope) 
  uS <- sort(unique(ddq_r_emm$slrs_r))
  uRTV <- sort(unique(ddq_r_emm$rt_vmax_lag_sc_r))
  uT = sort(unique(ddq_r_emm$trial_neg_inv_sc))
  ddq_r_emm <- ddq_r_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_r==uS[1] ~ '-2 std',
    slrs_r==uS[2] ~ '-1 std',
    slrs_r==uS[3] ~ 'mean',
    slrs_r==uS[4] ~ '+1 std',
    slrs_r==uS[5] ~ '+2 std'
  ), RT_vmax_bin = case_when(
    rt_vmax_lag_sc_r==uRTV[1] ~ '-2 std',
    rt_vmax_lag_sc_r==uRTV[2] ~ '-1 std',
    rt_vmax_lag_sc_r==uRTV[3] ~ 'mean',
    rt_vmax_lag_sc_r==uRTV[4] ~ '+1 std',
    rt_vmax_lag_sc_r==uRTV[5] ~ '+2 std'
  ))
  #ddq_r_emt <- ddq_r_emt %>% mutate(dataset_r = 'replication')%>% rename(rt_vmax_lag_sc_r.trend = rt_vmax_lag_sc.trend, slrs_r=subj_level_rand_slope) %>% select(network,dataset_r,rt_vmax_lag_sc_r.trend,slrs_r)
  ddq_f_emm <- ddq_f_emm %>% mutate(dataset = 'fMRI') %>% rename(rt_vmax_lag_sc_f = rt_vmax_lag_sc, slrs_f=subj_level_rand_slope)
  uS <- sort(unique(ddq_f_emm$slrs_f))
  uRTV <- sort(unique(ddq_f_emm$rt_vmax_lag_sc_f))
  ddq_f_emm <- ddq_f_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_f==uS[1] ~ '-2 std',
    slrs_f==uS[2] ~ '-1 std',
    slrs_f==uS[3] ~ 'mean',
    slrs_f==uS[4] ~ '+1 std',
    slrs_f==uS[5] ~ '+2 std'
  ), RT_vmax_bin = case_when(
    rt_vmax_lag_sc_f==uRTV[1] ~ '-2 std',
    rt_vmax_lag_sc_f==uRTV[2] ~ '-1 std',
    rt_vmax_lag_sc_f==uRTV[3] ~ 'mean',
    rt_vmax_lag_sc_f==uRTV[4] ~ '+1 std',
    rt_vmax_lag_sc_f==uRTV[5] ~ '+2 std'
  ))
  #ddq_f_emt <- ddq_f_emt %>% mutate(dataset_f = 'fMRI')%>% rename(rt_vmax_lag_sc_f.trend = rt_vmax_lag_sc.trend, slrs_f=subj_level_rand_slope) %>% select(vPFC_Entropy_random_slope,network,dataset_f,rt_vmax_lag_sc_f.trend,slrs_f,evt_time)
  ddq_r <- ddq_r %>% mutate(dataset = 'Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddq <- ddq%>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN' | network1=='CTR')
  
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Entropy-rt_vmax-convergence-DC-main.pdf',width=10,height=6)
  gg1 <- ggplot(ddq,aes(x=network1,y=estimate,color=network1)) +
    geom_violin() + facet_wrap(~dataset) +
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr)) +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    ylab('Convergence on Best RT (AU)') + xlab('Entropy Response') +
    scale_color_manual(values = pal1,labels=c('DMN','CTR')) + geom_hline(yintercept=0) +
    guides(color=guide_legend(title='Network'))
  #gg1 <- ggplot(ddq, aes(x=network,y=estimate, color=network)) + geom_bar(aes(alpha=p_level_fdr),stat='identity') + geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) + facet_grid(~dataset)
  print(gg1)
  dev.off()
  
  ddq_r_emm <- ddq_r_emm %>% select(!slrs_r & !rt_vmax_lag_sc_r)
  ddq_f_emm <- ddq_f_emm %>% select(!slrs_f & !rt_vmax_lag_sc_f)
  ddf <- rbind(ddq_r_emm,ddq_f_emm)
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std')) #& RT_vmax_bin=='mean')
  ddf <- ddf %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR', network=='L'~'LIM')) %>% select(!std.error)
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std')) #& RT_vmax_bin=='mean')
  ddf <- ddf %>% filter(network1=='DMN' | network1=='CTR')
  ddf <- inner_join(ddf,ddq,by=c('evt_time','dataset','network1'))

  ddf$RT_vmax_bin <- factor(ddf$RT_vmax_bin,levels=c('-2 std','-1 std','mean','+1 std','+2 std'))
  
   setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
   pdf('Entropy-rt_vmax-convergence-DC-nodiff.pdf',width=10,height=12)
   gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
     geom_violin() +
     facet_grid(RT_vmax_bin~dataset) + ylab('Convergence on Best RT (AU)') + xlab('Entropy Response') +
     theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
     guides(fill = guide_legend(byrow=TRUE)) +
     ylab('Convergence on Best RT (AU)') + xlab('Entropy Response') +
     scale_color_manual(values = pal1,labels=c('DMN','CTR')) +
     guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
  
  rm(ddf)
  rm(ddq)
  
  i <- 1 # entropy
  toalign <- 'clock' #
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-ranslopes-replication-',toalign,'-pred-rt_csv_sc-rt_lag-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_r <- ddq
  rm(ddq)
  model_str <- paste0('-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-rt_lag-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_f <- ddq
  rm(ddq)
  
  ddq_r_emm <- ddq_r$emmeans_list$RT
  ddq_r_emt <- ddq_r$emtrends_list$RT
  ddq_f_emm <- ddq_f$emmeans_list$RT
  ddq_f_emt <- ddq_f$emtrends_list$RT
  ddq_r <- ddq_r$coef_df_reml
  ddq_f <- ddq_f$coef_df_reml
  
  ddq_r <- ddq_r %>% filter(effect=='fixed') 
  ddq_f <- ddq_f %>% filter(effect=='fixed') 
  
  ddq_r_emm <- ddq_r_emm %>% mutate(dataset = 'Replication') %>% rename(rt_lag_sc_r = rt_lag_sc, slrs_r=subj_level_rand_slope)
  uS <- sort(unique(ddq_r_emm$slrs_r))
  uRTV <- sort(unique(ddq_r_emm$rt_lag_sc_r))
  uT = sort(unique(ddq_r_emm$trial_neg_inv_sc))
  ddq_r_emm <- ddq_r_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_r==uS[1] ~ '-2 std',
    slrs_r==uS[2] ~ '-1 std',
    slrs_r==uS[3] ~ 'mean',
    slrs_r==uS[4] ~ '+1 std',
    slrs_r==uS[5] ~ '+2 std'
  ), RT_lag_bin = case_when(
    rt_lag_sc_r==uRTV[1] ~ '-2 std',
    rt_lag_sc_r==uRTV[2] ~ '-1 std',
    rt_lag_sc_r==uRTV[3] ~ 'mean',
    rt_lag_sc_r==uRTV[4] ~ '+1 std',
    rt_lag_sc_r==uRTV[5] ~ '+2 std'
  ))
  #ddq_r_emt <- ddq_r_emt %>% mutate(dataset_r = 'replication')%>% rename(rt_vmax_lag_sc_r.trend = rt_vmax_lag_sc.trend, slrs_r=subj_level_rand_slope) %>% select(network,dataset_r,rt_vmax_lag_sc_r.trend,slrs_r)
  ddq_f_emm <- ddq_f_emm %>% mutate(dataset = 'fMRI')%>% rename(rt_lag_sc_f = rt_lag_sc, slrs_f=subj_level_rand_slope)
  uS <- sort(unique(ddq_f_emm$slrs_f))
  uRTV <- sort(unique(ddq_f_emm$rt_lag_sc_f))
  ddq_f_emm <- ddq_f_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_f==uS[1] ~ '-2 std',
    slrs_f==uS[2] ~ '-1 std',
    slrs_f==uS[3] ~ 'mean',
    slrs_f==uS[4] ~ '+1 std',
    slrs_f==uS[5] ~ '+2 std'
  ), RT_lag_bin = case_when(
    rt_lag_sc_f==uRTV[1] ~ '-2 std',
    rt_lag_sc_f==uRTV[2] ~ '-1 std',
    rt_lag_sc_f==uRTV[3] ~ 'mean',
    rt_lag_sc_f==uRTV[4] ~ '+1 std',
    rt_lag_sc_f==uRTV[5] ~ '+2 std'
  ))
  #ddq_f_emt <- ddq_f_emt %>% mutate(dataset_f = 'fMRI')%>% rename(rt_vmax_lag_sc_f.trend = rt_vmax_lag_sc.trend, slrs_f=subj_level_rand_slope) %>% select(vPFC_Entropy_random_slope,network,dataset_f,rt_vmax_lag_sc_f.trend,slrs_f,evt_time)
  ddq_r <- ddq_r %>% mutate(dataset = 'Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
  ddq <- ddq %>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN' | network1=='CTR')
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Entropy-rt_csv-divergence-DC-main.pdf',width=10,height=6)
  gg1 <- ggplot(ddq,aes(x=network1,y=estimate,color=network1)) +
    geom_violin() + facet_wrap(~dataset) +
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr)) +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) + scale_y_reverse() +
    ylab('RT swings (AU)') + xlab('Entropy Response') +
    scale_color_manual(values = pal1,labels=c('DMN')) + geom_hline(yintercept=0) +
    guides(color=guide_legend(title='Network'))
  #gg1 <- ggplot(ddq, aes(x=network,y=estimate, color=network)) + geom_bar(aes(alpha=p_level_fdr),stat='identity') + geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) + facet_grid(~dataset)
  print(gg1)
  dev.off()
  
  
  
  ddq_r_emm <- ddq_r_emm %>% select(!slrs_r & !rt_lag_sc_r)
  ddq_f_emm <- ddq_f_emm %>% select(!slrs_f & !rt_lag_sc_f)
  
  ddf <- rbind(ddq_r_emm,ddq_f_emm)
  #ddf <- ddf %>% filter(`vPFC Entropy response`=='10th %ile' | `vPFC Entropy response`=='90th %ile')
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std'))
  ddf <- ddf %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR',network=='L'~'LIM')) %>% select(!std.error)
  ddf <- ddf %>% filter(network=='D' | network=='C')
  ddf <- inner_join(ddf,ddq,by=c('evt_time','dataset','network1'))
  #ddf <- ddf %>% filter(RT_lag_bin=='mean')
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Entropy-rt_csv-divergence-D-nodiff.pdf',width=15,height=15)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin() + 
    facet_grid(RT_lag_bin~dataset) + ylab('RT Swings (AU)') +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    xlab('Entropy Response') + scale_y_reverse() +
    scale_color_manual(values = pal1,labels=c('DMN')) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
  
  ddq_r_emt <- ddq_r_emt %>% mutate(dataset = 'Replication') %>% rename(rt_lag_sc_r = rt_lag_sc.trend, slrs_r=subj_level_rand_slope)
  uS <- sort(unique(ddq_r_emt$slrs_r))
  uRTV <- sort(unique(ddq_r_emt$rt_lag_sc_r))
  ddq_r_emt <- ddq_r_emt %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_r==uS[1] ~ '-2 std',
    slrs_r==uS[2] ~ '-1 std',
    slrs_r==uS[3] ~ 'mean',
    slrs_r==uS[4] ~ '+1 std',
    slrs_r==uS[5] ~ '+2 std'
  ))
  
  ddq_f_emt <- ddq_f_emt %>% mutate(dataset = 'fMRI')%>% rename(rt_lag_sc_f = rt_lag_sc.trend, slrs_f=subj_level_rand_slope) 
  uS <- sort(unique(ddq_f_emt$slrs_f))
  uRTV <- sort(unique(ddq_f_emt$rt_lag_sc_f))
  ddq_f_emt <- ddq_f_emt %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_f==uS[1] ~ '-2 std',
    slrs_f==uS[2] ~ '-1 std',
    slrs_f==uS[3] ~ 'mean',
    slrs_f==uS[4] ~ '+1 std',
    slrs_f==uS[5] ~ '+2 std'
  ))
  #ddq_f_emt <- ddq_f_emt %>% mutate(dataset_f = 'fMRI')%>% rename(rt_vmax_lag_sc_f.trend = rt_vmax_lag_sc.trend, slrs_f=subj_level_rand_slope) %>% select(vPFC_Entropy_random_slope,network,dataset_f,rt_vmax_lag_sc_f.trend,slrs_f,evt_time)
  ddq_r <- ddq_r %>% mutate(dataset = 'Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
  ddq <- ddq%>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN')
  
  
  
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Entropy-rt_csv-divergence-D-nodiff-ss.pdf',width=15,height=15)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin(position=position_dodge(width=1)) +
    facet_grid(RT_lag_bin~dataset) + ylab('RT Swings (AU)') + xlab('Entropy Response') +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    xlab('Entropy Response') + scale_y_reverse() +
    scale_color_manual(values = pal1,labels=c('DMN')) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
}

if (do_value_plot){
  
  # plot nice figure
  library(wesanderson)
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1 = palette()
  pal1[1] <- pal[1]
  pal1[2] <- pal[2]
  pal1[3] <- pal[3]
  i <- 2 # value
  toalign <- 'clock' #
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-ranslopes-replication-',toalign,'-pred-rt_csv_sc-rt_vmax_lag-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_r <- ddq
  rm(ddq)
  model_str <- paste0('-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-rt_vmax_lag-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_f <- ddq
  rm(ddq)
  
  ddq_r_emm <- ddq_r$emmeans_list$Vmax
  ddq_r_emt <- ddq_r$emtrends_list$Vmax
  ddq_f_emm <- ddq_f$emmeans_list$Vmax
  ddq_f_emt <- ddq_f$emtrends_list$Vmax
  ddq_r <- ddq_r$coef_df_reml
  ddq_f <- ddq_f$coef_df_reml
  
  ddq_r <- ddq_r %>% filter(effect=='fixed') 
  ddq_f <- ddq_f %>% filter(effect=='fixed') 
  
  ddq_r_emm <- ddq_r_emm %>% mutate(dataset = 'Replication') %>% rename(rt_vmax_lag_sc_r = rt_vmax_lag_sc, slrs_r=subj_level_rand_slope)
  uS <- sort(unique(ddq_r_emm$slrs_r))
  uRTV <- sort(unique(ddq_r_emm$rt_vmax_lag_sc_r))
  uT = sort(unique(ddq_r_emm$trial_neg_inv_sc))
  ddq_r_emm <- ddq_r_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_r==uS[1] ~ '-2 std',
    slrs_r==uS[2] ~ '-1 std',
    slrs_r==uS[3] ~ 'mean',
    slrs_r==uS[4] ~ '+1 std',
    slrs_r==uS[5] ~ '+2 std'
  ), RT_vmax_bin = case_when(
    rt_vmax_lag_sc_r==uRTV[1] ~ '-2 std',
    rt_vmax_lag_sc_r==uRTV[2] ~ '-1 std',
    rt_vmax_lag_sc_r==uRTV[3] ~ 'mean',
    rt_vmax_lag_sc_r==uRTV[4] ~ '+1 std',
    rt_vmax_lag_sc_r==uRTV[5] ~ '+2 std'
  ))
  #ddq_r_emt <- ddq_r_emt %>% mutate(dataset_r = 'replication')%>% rename(rt_vmax_lag_sc_r.trend = rt_vmax_lag_sc.trend, slrs_r=subj_level_rand_slope) %>% select(network,dataset_r,rt_vmax_lag_sc_r.trend,slrs_r)
  ddq_f_emm <- ddq_f_emm %>% mutate(dataset = 'fMRI')%>% rename(rt_vmax_lag_sc_f = rt_vmax_lag_sc, slrs_f=subj_level_rand_slope)
  uS <- sort(unique(ddq_f_emm$slrs_f))
  uRTV <- sort(unique(ddq_f_emm$rt_vmax_lag_sc_f))
  ddq_f_emm <- ddq_f_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_f==uS[1] ~ '-2 std',
    slrs_f==uS[2] ~ '-1 std',
    slrs_f==uS[3] ~ 'mean',
    slrs_f==uS[4] ~ '+1 std',
    slrs_f==uS[5] ~ '+2 std'
  ), RT_vmax_bin = case_when(
    rt_vmax_lag_sc_f==uRTV[1] ~ '-2 std',
    rt_vmax_lag_sc_f==uRTV[2] ~ '-1 std',
    rt_vmax_lag_sc_f==uRTV[3] ~ 'mean',
    rt_vmax_lag_sc_f==uRTV[4] ~ '+1 std',
    rt_vmax_lag_sc_f==uRTV[5] ~ '+2 std'
  ))
  #ddq_f_emt <- ddq_f_emt %>% mutate(dataset_f = 'fMRI')%>% rename(rt_vmax_lag_sc_f.trend = rt_vmax_lag_sc.trend, slrs_f=subj_level_rand_slope) %>% select(vPFC_Entropy_random_slope,network,dataset_f,rt_vmax_lag_sc_f.trend,slrs_f,evt_time)
  ddq_r <- ddq_r %>% mutate(dataset = 'Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='rt_vmax_lag_sc:subj_level_rand_slope')
  ddq <- ddq%>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN' | network1=='CTR')
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Value-rt_vmax-convergence-D-main.pdf',width=10,height=6)
  gg1 <- ggplot(ddq,aes(x=network1,y=estimate,color=network1)) +
    geom_violin() + facet_wrap(~dataset) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr)) +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    ylab('Convergence on Best RT (AU)') + xlab('Value Response') +
    scale_color_manual(values = pal1,labels=c('DMN','CTR')) + geom_hline(yintercept=0) +
    guides(color=guide_legend(title='Network'))
  #gg1 <- ggplot(ddq, aes(x=network,y=estimate, color=network)) + geom_bar(aes(alpha=p_level_fdr),stat='identity') + geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) + facet_grid(~dataset)
  print(gg1)
  dev.off()
  
  ddq_r_emm <- ddq_r_emm %>% select(!slrs_r & !rt_vmax_lag_sc_r)
  ddq_f_emm <- ddq_f_emm %>% select(!slrs_f & !rt_vmax_lag_sc_f)
  ddf <- rbind(ddq_r_emm,ddq_f_emm)
  ddf <- ddf %>% filter((network=='D') & (vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std'))# & RT_vmax_bin=='mean')
  ddf <- ddf %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR', network=='L'~'LIM')) %>% select(!std.error)
  
  
  #ddf <- ddf %>% filter(`vPFC Entropy response`=='10th %ile' | `vPFC Entropy response`=='90th %ile')
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std') & RT_vmax_bin=='mean')
  ddf <- ddf %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR', network=='L'~'LIM'))
  ddf <- ddf %>% filter(network1=='DMN' | network1=='CTR')
  ddf <- inner_join(ddf,ddq,by=c('evt_time','dataset','network1'))
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Value-rt_vmax-convergence-DC-nodiff.pdf',width=10,height=6)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin() + 
    facet_grid(RT_vmax_bin~dataset) + ylab('Convergence on Best RT (AU)') +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    ylab('Convergence on Best RT (AU)') + xlab('Value Response') +
    scale_color_manual(values = pal1,labels=c('DMN','CTR')) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
  
  rm(ddf)
  rm(ddq)
  
  i <- 2 # value
  toalign <- 'clock' #
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-ranslopes-replication-',toalign,'-pred-rt_csv_sc-rt_lag-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_r <- ddq
  rm(ddq)
  model_str <- paste0('-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-rt_lag-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_f <- ddq
  rm(ddq)
  
  ddq_r_emm <- ddq_r$emmeans_list$RT
  ddq_r_emt <- ddq_r$emtrends_list$RT
  ddq_f_emm <- ddq_f$emmeans_list$RT
  ddq_f_emt <- ddq_f$emtrends_list$RT
  ddq_r <- ddq_r$coef_df_reml
  ddq_f <- ddq_f$coef_df_reml
  
  ddq_r <- ddq_r %>% filter(effect=='fixed') 
  ddq_f <- ddq_f %>% filter(effect=='fixed') 
  
  ddq_r_emm <- ddq_r_emm %>% mutate(dataset = 'Replication') %>% rename(rt_lag_sc_r = rt_lag_sc, slrs_r=subj_level_rand_slope)
  uS <- sort(unique(ddq_r_emm$slrs_r))
  uRTV <- sort(unique(ddq_r_emm$rt_lag_sc_r))
  uT = sort(unique(ddq_r_emm$trial_neg_inv_sc))
  ddq_r_emm <- ddq_r_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_r==uS[1] ~ '-2 std',
    slrs_r==uS[2] ~ '-1 std',
    slrs_r==uS[3] ~ 'mean',
    slrs_r==uS[4] ~ '+1 std',
    slrs_r==uS[5] ~ '+2 std'
  ), RT_lag_bin = case_when(
    rt_lag_sc_r==uRTV[1] ~ '-2 std',
    rt_lag_sc_r==uRTV[2] ~ '-1 std',
    rt_lag_sc_r==uRTV[3] ~ 'mean',
    rt_lag_sc_r==uRTV[4] ~ '+1 std',
    rt_lag_sc_r==uRTV[5] ~ '+2 std'
  ))
  #ddq_r_emt <- ddq_r_emt %>% mutate(dataset_r = 'replication')%>% rename(rt_vmax_lag_sc_r.trend = rt_vmax_lag_sc.trend, slrs_r=subj_level_rand_slope) %>% select(network,dataset_r,rt_vmax_lag_sc_r.trend,slrs_r)
  ddq_f_emm <- ddq_f_emm %>% mutate(dataset = 'fMRI')%>% rename(rt_lag_sc_f = rt_lag_sc, slrs_f=subj_level_rand_slope)
  uS <- sort(unique(ddq_f_emm$slrs_f))
  uRTV <- sort(unique(ddq_f_emm$rt_lag_sc_f))
  ddq_f_emm <- ddq_f_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_f==uS[1] ~ '-2 std',
    slrs_f==uS[2] ~ '-1 std',
    slrs_f==uS[3] ~ 'mean',
    slrs_f==uS[4] ~ '+1 std',
    slrs_f==uS[5] ~ '+2 std'
  ), RT_lag_bin = case_when(
    rt_lag_sc_f==uRTV[1] ~ '-2 std',
    rt_lag_sc_f==uRTV[2] ~ '-1 std',
    rt_lag_sc_f==uRTV[3] ~ 'mean',
    rt_lag_sc_f==uRTV[4] ~ '+1 std',
    rt_lag_sc_f==uRTV[5] ~ '+2 std'
  ))
  #ddq_f_emt <- ddq_f_emt %>% mutate(dataset_f = 'fMRI')%>% rename(rt_vmax_lag_sc_f.trend = rt_vmax_lag_sc.trend, slrs_f=subj_level_rand_slope) %>% select(vPFC_Entropy_random_slope,network,dataset_f,rt_vmax_lag_sc_f.trend,slrs_f,evt_time)
  ddq_r <- ddq_r %>% mutate(dataset = 'Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
  ddq <- ddq%>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN' | network1=='CTR')
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Value-rt_csv-divergence-D-main.pdf',width=10,height=6)
  gg1 <- ggplot(ddq,aes(x=network1,y=estimate,color=network1)) +
    geom_violin() + facet_wrap(~dataset) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr)) +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) + scale_y_reverse() +
    ylab('RT swings (AU)') + xlab('Value Response') +
    scale_color_manual(values = pal1,labels=c('DMN')) + geom_hline(yintercept=0) +
    guides(color=guide_legend(title='Network'))
  #gg1 <- ggplot(ddq, aes(x=network,y=estimate, color=network)) + geom_bar(aes(alpha=p_level_fdr),stat='identity') + geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) + facet_grid(~dataset)
  print(gg1)
  dev.off()
  
  
  
  ddq_r_emm <- ddq_r_emm %>% select(!slrs_r & !rt_lag_sc_r)
  ddq_f_emm <- ddq_f_emm %>% select(!slrs_f & !rt_lag_sc_f)
  
  ddf <- rbind(ddq_r_emm,ddq_f_emm)
  #ddf <- ddf %>% filter(`vPFC Entropy response`=='10th %ile' | `vPFC Entropy response`=='90th %ile')
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std'))
  ddf <- ddf %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR',network=='L'~'LIM')) %>% select(!std.error)
  ddf <- ddf %>% filter(network=='D')
  ddf <- inner_join(ddf,ddq,by=c('evt_time','dataset','network1'))
  #ddf <- ddf %>% filter(RT_lag_bin=='mean')
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Value-rt_csv-divergence-D-nodiff.pdf',width=15,height=15)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin(position=position_dodge(width=1)) +
    facet_grid(RT_lag_bin~dataset) + ylab('RT Swings (AU)') + xlab('Value Response') +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    xlab('Value Response') + scale_y_reverse() +
    scale_color_manual(values = pal1,labels=c('DMN')) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
  
  ddq_r_emt <- ddq_r_emt %>% mutate(dataset = 'Replication') %>% rename(rt_lag_sc_r = rt_lag_sc.trend, slrs_r=subj_level_rand_slope) 
  uS <- sort(unique(ddq_r_emt$slrs_r))
  uRTV <- sort(unique(ddq_r_emt$rt_lag_sc_r))
  ddq_r_emt <- ddq_r_emt %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_r==uS[1] ~ '-2 std',
    slrs_r==uS[2] ~ '-1 std',
    slrs_r==uS[3] ~ 'mean',
    slrs_r==uS[4] ~ '+1 std',
    slrs_r==uS[5] ~ '+2 std'
  ))
  
  ddq_f_emt <- ddq_f_emt %>% mutate(dataset = 'fMRI')%>% rename(rt_lag_sc_f = rt_lag_sc.trend, slrs_f=subj_level_rand_slope) 
  uS <- sort(unique(ddq_f_emt$slrs_f))
  uRTV <- sort(unique(ddq_f_emt$rt_lag_sc_f))
  ddq_f_emt <- ddq_f_emt %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_f==uS[1] ~ '-2 std',
    slrs_f==uS[2] ~ '-1 std',
    slrs_f==uS[3] ~ 'mean',
    slrs_f==uS[4] ~ '+1 std',
    slrs_f==uS[5] ~ '+2 std'
  ))
  #ddq_f_emt <- ddq_f_emt %>% mutate(dataset_f = 'fMRI')%>% rename(rt_vmax_lag_sc_f.trend = rt_vmax_lag_sc.trend, slrs_f=subj_level_rand_slope) %>% select(vPFC_Entropy_random_slope,network,dataset_f,rt_vmax_lag_sc_f.trend,slrs_f,evt_time)
  ddq_r <- ddq_r %>% mutate(dataset = 'Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
  ddq <- ddq%>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN')
  
  
  
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Value-rt_csv-divergence-D-nodiff-ss.pdf',width=15,height=15)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin(position=position_dodge(width=1)) + geom_point(position=position_dodge(width=1),aes(alpha=p_level_fdr,size=p_level_fdr)) +
    facet_grid(RT_lag_bin~dataset) + ylab('RT Swings (AU)') + xlab('Value Response') +
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr)) +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    xlab('Value Response') + scale_y_reverse() +
    scale_color_manual(values = pal1,labels=c('DMN')) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
}
