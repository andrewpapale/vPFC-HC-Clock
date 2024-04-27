############
# 2022-05-23 AndyP
# Subject level random slopes

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(ggnewscale)
library(quest)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26
toalign <- 'clock'
do_rand_slopes = TRUE
do_rt_pred_fmri = TRUE
simple_model = FALSE
trial_mod_model = TRUE

#### clock ####

if (do_rand_slopes){
  load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')  
  
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
           iti_lag_sc = scale(iti_prev),
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
  df <- df %>% select(iti_ideal,condition_trial_neg_inv_sc,rt_lag_sc,iti_prev,iti_sc,v_entropy_wi_change,iti_lag_sc,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
  df$id <- as.character(df$id)
  Q <- full_join(md,df,by=c('id','run','trial'))
  
  Q$decon_mean[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$decon_mean[Q$evt_time < -Q$iti_prev] = NA
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
  Q <- Q %>% rename(vmPFC_decon = decon_mean)
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
  Q <- Q %>% filter(group!='ATT')
  Q$group <- relevel(factor(Q$group),ref='HC')
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
  Q <- Q %>% mutate(run_trial0_c = run_trial-floor(run_trial/40.5)*40,
                    run_trial0_neg_inv = -1000 / run_trial0_c,
                    run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv))
  )
  Q <- Q %>% filter(group=='HC')
  Q <- Q %>% filter(!is.na(rewFunc))
  decode_formula <- NULL
  decode_formula[[1]] = formula(~ age + last_outcome + run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + v_entropy_wi + (1 + v_entropy_wi|id) + (1|run))
  decode_formula[[2]] = formula(~ age + last_outcome + run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + v_max_wi + (1 + v_max_wi |id) + (1|run))
  #decode_formula[[3]] = formula(~ age + last_outcome + rt_lag_sc + iti_lag_sc + v_max_wi + v_entropy_wi + (1 + v_entropy_wi + v_max_wi |id) + (1|run))
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
    save(ddf,file=paste0(curr_date,'-vPFC-network-',toalign,'-Explore-ranslopes-HConly-trial_mod-trial1-10included-',i,'.Rdata'))
  }
}

if (do_rt_pred_fmri){
  for (i in 1:2){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vPFC-network-',toalign,'-Explore-ranslopes-HConly-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    
    ### Does random slope predict rt_vmax or rt_swing?
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
    #qdf <- qdf %>% group_by(network) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate & !rhs) %>% rename(estimate=estimate1)
    qdf <- qdf %>% rename(id=level)
    qdf$id <- as.integer(qdf$id)
    qdf <- qdf %>% select(!outcome)
    qdf <- qdf %>% group_by(network) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate = estimate1)
    #qdf %>% group_by(evt_time) %>% mutate(estimate = quest::winsor(estimate,z.min=2,z.max=2,to.na=TRUE)) %>% ungroup()
    source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
    df <- get_trial_data(repo_directory='/Volumes/Users/Andrew/MEDuSA_data_Explore',dataset='explore')
    df <- df %>%
      group_by(id,run_number) %>%
      arrange(id, run_number, trial) %>%
      mutate(condition_trial = run_trial-floor(run_trial/40.5)*40,
             condition_trial_neg_inv = -1000 / condition_trial,
             condition_trial_neg_inv_sc = scale(as.vector(condition_trial_neg_inv)),
      ) %>% ungroup()
    df <- df %>%
      group_by(id) %>% 
      mutate(v_entropy_sc_r = scale(v_entropy)) %>% ungroup() %>%
      group_by(id, run) %>% 
      mutate(v_chosen_sc = scale(v_chosen),
             abs_pe_max_sc = scale(abs(pe_max)),
             score_sc = scale(score_csv),
             iti_sc = scale(iti_ideal),
             iti_prev_sc = scale(iti_prev),
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
    df <- df %>% select(iti_ideal,condition_trial_neg_inv_sc,rt_vmax_lag_sc, rt_lag_sc,iti_prev,iti_sc,v_entropy_wi_change,iti_prev_sc,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
    df$id <- as.integer(df$id)
    df <- df %>% group_by(id,run) %>% mutate(v_max_wi_lag = lag(v_max_wi)) %>% ungroup()
    Q <- inner_join(qdf,df,by=c('id'))
    
    Q <- Q %>% arrange(id,run,trial,evt_time)
    Q <- Q %>% filter(evt_time > -4 & evt_time < 4)
    
    
    Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
    Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
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
                        run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv))
    )
    Q <- Q %>% rename(subj_level_rand_slope=estimate)
    Q$last_outcome <- relevel(factor(Q$last_outcome),ref='Omission')
    decode_formula <- NULL
    #decode_formula[[1]] = formula(~ rt_lag*subj_level_rand_slope + last_outcome + v_entropy_wi + v_max_wi + trial_neg_inv_sc +  (1+rt_lag | id/run))
    #decode_formula[[2]] = formula(~ rt_vmax_lag*subj_level_rand_slope + last_outcome + v_entropy_wi + v_max_wi + trial_neg_inv_sc +  (1 + rt_vmax_lag | id/run))
    if (!simple_model && !trial_mod_model){
      decode_formula[[1]] <- formula(~(condition_trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * condition_trial_neg_inv_sc * subj_level_rand_slope + (1 | id/run))
      decode_formula[[2]] <- formula(~(condition_trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * condition_trial_neg_inv_sc * subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc | id/run))
    } else if (simple_model){
      decode_formula[[1]] <- formula(~condition_trial_neg_inv_sc + rt_lag_sc*subj_level_rand_slope + last_outcome*subj_level_rand_slope + rt_vmax_lag_sc*subj_level_rand_slope + (1|id/run))
      decode_formula[[2]] <- formula(~condition_trial_neg_inv_sc + rt_lag_sc*subj_level_rand_slope + last_outcome*subj_level_rand_slope + rt_vmax_lag_sc*subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc |id/run))
    } else if (trial_mod_model){
      decode_formula[[1]] <- formula(~(run_trial0_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * run_trial0_neg_inv_sc * subj_level_rand_slope + (1 | id/run))
      decode_formula[[2]] <- formula(~(run_trial0_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * run_trial0_neg_inv_sc * subj_level_rand_slope + (1 + rt_vmax_lag_sc + rt_lag_sc | id/run))
    }
    splits = c('evt_time','network')
    source('~/fmri.pipeline/R/mixed_by.R')
    print(i)
    for (j in 1:length(decode_formula)){
      if (!simple_model && !trial_mod_model){
        ddq <- mixed_by(Q, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                        emmeans_spec = list(
                          RT = list(outcome='rt_csv', model_name='model1', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv', model_name='model1', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv',model_name='model1',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2)))        
                          
                        ),
                        emtrends_spec = list(
                          RT = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv', model_name='model1', var='rt_vmax_lag_sc', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          RTxO = list(outcome='rt_csv',model_name='model1',var='rt_lag_sc',
                                      specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                          
                        )
        )
        
        setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
        curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
        if (j==1){
          save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-Explore-pred-rt_csv_sc-int-HConly-trial1-10included-',i,'.Rdata'))
        } else {
          save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-Explore-pred-rt_csv_sc-slo-HConly-trial1-10included-',i,'.Rdata'))
        }
      } else if (simple_model){
        ddq <- mixed_by(Q, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                        emmeans_spec = list(
                          RT = list(outcome='rt_csv', model_name='model1', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv', model_name='model1', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2)))
                        ),
                        emtrends_spec = list(
                          RT = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
                                    specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                          Vmax = list(outcome='rt_csv', model_name='model1', var='rt_vmax_lag_sc', 
                                      specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2)))
                          )
        )
        
        setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
        curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
        if (j==1){
          save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-Explore-pred-rt_csv_sc-int-HConly-simple-trial1-10included-',i,'.Rdata'))
        } else {
          save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-Explore-pred-rt_csv_sc-slo-HConly-simple-trial1-10included-',i,'.Rdata'))
        }        
    } else if (trial_mod_model){
      ddq <- mixed_by(Q, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                      emmeans_spec = list(
                        RT = list(outcome='rt_csv', model_name='model1', 
                                  specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                        Vmax = list(outcome='rt_csv', model_name='model1', 
                                    specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                        RTxO = list(outcome='rt_csv',model_name='model1',
                                    specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                        TrxVmax = list(outcome='rt_csv',model_name='model1',
                                       specs=formula(~rt_vmax_lag_sc:run_trial0_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2),run_trial0_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                        
                      ),
                      emtrends_spec = list(
                        RT = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
                                  specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                        Vmax = list(outcome='rt_csv', model_name='model1', var='rt_vmax_lag_sc', 
                                    specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                        RTxO = list(outcome='rt_csv',model_name='model1',var='rt_lag_sc',
                                    specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                        TrxVmax = list(outcome='rt_csv',model_name='model1', var = 'rt_vmax_lag_sc',
                                       specs=formula(~rt_vmax_lag_sc:run_trial0_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),run_trial0_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                        TrxVmax1 = list(outcome='rt_csv',model_name='model1', var = 'run_trial0_neg_inv_sc',
                                        specs=formula(~rt_vmax_lag_sc:run_trial0_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                        TrxVmax2 = list(outcome='rt_csv',model_name='model1', var = 'subj_level_rand_slope',
                                        specs=formula(~rt_vmax_lag_sc:run_trial0_neg_inv_sc:subj_level_rand_slope), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),run_trial0_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                      )
      )
      
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      if (j==1){
        save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-',i,'.Rdata'))
      } else {
        save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-Explore-pred-rt_csv_sc-slo-HConly-trial_mod-trial1-10included-',i,'.Rdata'))
      }
    }
  }
}
}