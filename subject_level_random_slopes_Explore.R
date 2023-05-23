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
do_entropy_plot = TRUE
do_value_plot = TRUE

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
  Q$decon_mean[Q$evt_time < -Q$rt_csv] = NA
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
  #Q <- Q %>% filter(group=='HC')
  Q <- Q %>% filter(!is.na(rewFunc))
  decode_formula <- NULL
  decode_formula[[1]] = formula(~ age + wtar + last_outcome + rt_lag_sc + iti_lag_sc + v_entropy_wi + (1 + v_entropy_wi|id) + (1|run))
  decode_formula[[2]] = formula(~ age + wtar + last_outcome + rt_lag_sc + iti_lag_sc + v_max_wi + (1 + v_max_wi |id) + (1|run))
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
    save(ddf,file=paste0(curr_date,'-vPFC-network-',toalign,'-Explore-ranslopes-',i,'.Rdata'))
  }
}

if (do_rt_pred_fmri){
  for (i in 1:2){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vPFC-network-',toalign,'-Explore-ranslopes-',i,'.Rdata')
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
    Q <- full_join(qdf,df,by=c('id'))
    Q$decon_mean[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
    Q$decon_mean[Q$evt_time < -Q$rt_csv] = NA

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
    Q <- Q %>% select(!group)
    Q <- merge(demo,Q,by='id')
    Q <- Q %>% filter(group!='ATT')
    Q$group <- relevel(factor(Q$group),ref='HC')
    Q <- Q %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                        trial > 40 & trial <= 80 ~ 2,
                                        trial > 80 & trial <=120 ~ 3, 
                                        trial > 120 & trial <=160 ~ 4,
                                        trial > 160 & trial <=200 ~ 5,
                                        trial > 200 & trial <=240 ~ 6))
    #Q <- Q %>% filter(group=='HC')
    Q <- Q %>% filter(!is.na(rewFunc))
    Q <- Q %>% rename(subj_level_rand_slope=estimate)
    Q$last_outcome <- relevel(factor(Q$last_outcome),ref='Reward')
      decode_formula <- NULL
      #decode_formula[[1]] = formula(~ rt_lag*subj_level_rand_slope + last_outcome + v_entropy_wi + v_max_wi + trial_neg_inv_sc +  (1+rt_lag | id/run))
      #decode_formula[[2]] = formula(~ rt_vmax_lag*subj_level_rand_slope + last_outcome + v_entropy_wi + v_max_wi + trial_neg_inv_sc +  (1 + rt_vmax_lag | id/run))
      decode_formula[[1]] <- formula(~(condition_trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * trial_neg_inv_sc * subj_level_rand_slope + (1 | id/run))
  
    splits = c('evt_time','network')
    source('~/fmri.pipeline/R/mixed_by.R')
    print(i)
        ddq <- mixed_by(Q, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[1]], split_on = splits,return_models=TRUE,
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
      save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-replication-',toalign,'-Explore-pred-rt_csv_sc-both-rs-',i,'.Rdata'))
}
}

if (do_entropy_plot){
  # plot nice figure
  library(wesanderson)
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1 = palette()
  pal1[1] <- pal[2]
  pal1[2] <- pal[1]
  pal1[3] <- pal[3]
  i <- 1 # entropy
  toalign <- 'clock' #
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-ranslopes-replication-clock-Explore-pred-rt_csv_sc-both-rs-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_r <- ddq
  rm(ddq)
  model_str <- paste0('-vmPFC-network-ranslopes-replication-clock-Explore-pred-rt_csv_sc-both-rs-',i,'.Rdata')
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
  
  ddq_r_emm <- ddq_r_emm %>% mutate(dataset = 'non-fMRI Replication') %>% rename(rt_vmax_lag_r = rt_vmax_lag_sc, slrs_r=subj_level_rand_slope) 
  uS <- sort(unique(ddq_r_emm$slrs_r))
  uRTV <- sort(unique(ddq_r_emm$rt_vmax_lag_r))
  uT = sort(unique(ddq_r_emm$trial_neg_inv_sc))
  ddq_r_emm <- ddq_r_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_r==uS[1] ~ '-2 std',
    slrs_r==uS[2] ~ '-1 std',
    slrs_r==uS[3] ~ 'mean',
    slrs_r==uS[4] ~ '+1 std',
    slrs_r==uS[5] ~ '+2 std'
  ), RT_vmax_bin = case_when(
    rt_vmax_lag_r==uRTV[1] ~ '-2 std',
    rt_vmax_lag_r==uRTV[2] ~ '-1 std',
    rt_vmax_lag_r==uRTV[3] ~ 'mean',
    rt_vmax_lag_r==uRTV[4] ~ '+1 std',
    rt_vmax_lag_r==uRTV[5] ~ '+2 std'
  ))
  #ddq_r_emt <- ddq_r_emt %>% mutate(dataset_r = 'replication')%>% rename(rt_vmax_lag_sc_r.trend = rt_vmax_lag_sc.trend, slrs_r=subj_level_rand_slope) %>% select(network,dataset_r,rt_vmax_lag_sc_r.trend,slrs_r)
  ddq_f_emm <- ddq_f_emm %>% mutate(dataset = 'fMRI') %>% rename(rt_vmax_lag_f = rt_vmax_lag_sc, slrs_f=subj_level_rand_slope)
  uS <- sort(unique(ddq_f_emm$slrs_f))
  uRTV <- sort(unique(ddq_f_emm$rt_vmax_lag_f))
  ddq_f_emm <- ddq_f_emm %>% mutate(vPFC_Entropy_random_slope = case_when(
    slrs_f==uS[1] ~ '-2 std',
    slrs_f==uS[2] ~ '-1 std',
    slrs_f==uS[3] ~ 'mean',
    slrs_f==uS[4] ~ '+1 std',
    slrs_f==uS[5] ~ '+2 std'
  ), RT_vmax_bin = case_when(
    rt_vmax_lag_f==uRTV[1] ~ '-2 std',
    rt_vmax_lag_f==uRTV[2] ~ '-1 std',
    rt_vmax_lag_f==uRTV[3] ~ 'mean',
    rt_vmax_lag_f==uRTV[4] ~ '+1 std',
    rt_vmax_lag_f==uRTV[5] ~ '+2 std'
  ))
  #ddq_f_emt <- ddq_f_emt %>% mutate(dataset_f = 'fMRI')%>% rename(rt_vmax_lag_sc_f.trend = rt_vmax_lag_sc.trend, slrs_f=subj_level_rand_slope) %>% select(vPFC_Entropy_random_slope,network,dataset_f,rt_vmax_lag_sc_f.trend,slrs_f,evt_time)
  ddq_r <- ddq_r %>% mutate(dataset = 'non-fMRI Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
  ddq <- ddq%>% mutate(network1 = case_when(network=='DMN'~'DMN', network=='CTR'~'CTR',network=='LIM'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term <= .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term <= .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <=.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <=.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN' | network1=='CTR')
  
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Explore-Entropy-rt_vmax-convergence-DC-main.pdf',width=14,height=7)
  gg1 <- ggplot(ddq,aes(x=network1,y=estimate,color=network1)) +
    geom_violin() + facet_wrap(~dataset) +
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr)) +
    theme(axis.text.x = element_text(size=25), axis.text.y = element_text(size=25), axis.title = element_text(size=35), strip.text.x = element_text(size=35),strip.text.y = element_text(size=35),legend.title=element_text(size=35),legend.text=element_text(size=30),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    ylab('Convergence on Best RT (AU)') + xlab('') +
    scale_color_manual(values = pal1) + geom_hline(yintercept=0) +
    guides(color=guide_legend(title='Network'))
  #gg1 <- ggplot(ddq, aes(x=network,y=estimate, color=network)) + geom_bar(aes(alpha=p_level_fdr),stat='identity') + geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) + facet_grid(~dataset)
  print(gg1)
  dev.off()
  
  ddq_r_emm <- ddq_r_emm %>% select(!slrs_r & !rt_vmax_lag_r)
  ddq_f_emm <- ddq_f_emm %>% select(!slrs_f & !rt_vmax_lag_f)
  ddf <- rbind(ddq_r_emm,ddq_f_emm)
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std')) #& RT_vmax_bin=='mean')
  ddf <- ddf %>% mutate(network1 = case_when(network=='DMN'~'DMN',network=='CTR'~'CTR', network=='LIM'~'LIM')) %>% select(!std.error)
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std')) #& RT_vmax_bin=='mean')
  ddf <- ddf %>% filter(network1=='DMN' | network1=='CTR')
  ddf <- inner_join(ddf,ddq,by=c('evt_time','dataset','network1'))
  
  ddf$RT_vmax_bin <- factor(ddf$RT_vmax_bin,levels=c('-2 std','-1 std','mean','+1 std','+2 std'))
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Explore-Entropy-rt_vmax-convergence-DC-nodiff.pdf',width=14,height=12)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin() +
    facet_grid(RT_vmax_bin~dataset) + ylab('Convergence on Best RT (AU)') + xlab('Entropy Response') +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    ylab('Convergence on Best RT (AU)') + xlab('Entropy Response') +
    scale_color_manual(values = pal1) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
  
  rm(ddf)
  rm(ddq)
  
  i <- 1 # entropy
  toalign <- 'clock' #
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-ranslopes-replication-clock-Explore-pred-rt_csv_sc-both-rs-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_r <- ddq
  rm(ddq)
  model_str <- paste0('-vmPFC-network-ranslopes-replication-clock-Explore-pred-rt_csv_sc-both-rs-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_f <- ddq
  rm(ddq)
  
  ddq_r_emm <- ddq_r$emmeans_list$RTxO
  ddq_r_emt <- ddq_r$emtrends_list$RTxO
  ddq_f_emm <- ddq_f$emmeans_list$RTxO
  ddq_f_emt <- ddq_f$emtrends_list$RTxO
  ddq_r <- ddq_r$coef_df_reml
  ddq_f <- ddq_f$coef_df_reml
  
  ddq_r <- ddq_r %>% filter(effect=='fixed') 
  ddq_f <- ddq_f %>% filter(effect=='fixed') 
  
  ddq_r_emm <- ddq_r_emm %>% mutate(dataset = 'non-fMRI Replication') %>% rename(rt_lag_sc_r = rt_lag_sc, slrs_r=subj_level_rand_slope)
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
  ddq_r <- ddq_r %>% mutate(dataset = 'non-fMRI Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
  ddq <- ddq %>% mutate(network1 = case_when(network=='DMN'~'DMN', network=='CTR'~'CTR',network=='LIM'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term <= .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term <= .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <=.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <=.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN' | network1=='CTR')
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Explore-Explore-Entropy-rt_csv-divergence-DC-main.pdf',width=14,height=7)
  gg1 <- ggplot(ddq,aes(x=network1,y=estimate,color=network1)) +
    geom_violin() + facet_wrap(~dataset) +
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr)) +
    theme(axis.text.x = element_text(size=25), axis.text.y = element_text(size=25), axis.title = element_text(size=35), strip.text.x = element_text(size=35),strip.text.y = element_text(size=35),legend.title=element_text(size=35),legend.text=element_text(size=30),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) + scale_y_reverse() +
    ylab('RT swings (AU)') + xlab('') +
    scale_color_manual(values = pal1) + geom_hline(yintercept=0) +
    guides(color=guide_legend(title='Network'))
  #gg1 <- ggplot(ddq, aes(x=network,y=estimate, color=network)) + geom_bar(aes(alpha=p_level_fdr),stat='identity') + geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) + facet_grid(~dataset)
  print(gg1)
  dev.off()
  
  
  
  ddq_r_emm <- ddq_r_emm %>% select(!slrs_r & !rt_lag_sc_r)
  ddq_f_emm <- ddq_f_emm %>% select(!slrs_f & !rt_lag_sc_f)
  
  ddf <- rbind(ddq_r_emm,ddq_f_emm)
  #ddf <- ddf %>% filter(`vPFC Entropy response`=='10th %ile' | `vPFC Entropy response`=='90th %ile')
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std'))
  ddf <- ddf %>% mutate(network1 = case_when(network=='DMN'~'DMN',network=='CTR'~'CTR',network=='LIM'~'LIM')) %>% select(!std.error)
  ddf <- ddf %>% filter(network1=='DMN' | network1=='CTR')
  ddf <- inner_join(ddf,ddq,by=c('evt_time','dataset','network1'))
  ddf <- ddf %>% filter(last_outcome=='Reward')
  #ddf <- ddf %>% filter(RT_lag_bin=='mean')
  
  ddf$RT_lag_bin <- factor(ddf$RT_lag_bin,levels=c('-2 std','-1 std','mean','+1 std','+2 std'))
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Explore-Explore-Entropy-rt_csv-divergence-DC-nodiff.pdf',width=17,height=15)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin() + 
    facet_grid(RT_lag_bin~dataset) + ylab('RT Swings (AU)') +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    xlab('Entropy Response') + scale_y_reverse() +
    scale_color_manual(values = pal1) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
  
  ddq_r_emt <- ddq_r_emt %>% mutate(dataset = 'non-fMRI Replication') %>% rename(rt_lag_sc_r = rt_lag_sc.trend, slrs_r=subj_level_rand_slope)
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
  ddq_r <- ddq_r %>% mutate(dataset = 'non-fMRI Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
  ddq <- ddq%>% mutate(network1 = case_when(network=='DMN'~'DMN', network=='CTR'~'CTR',network=='LIM'~'LIM'))
  #ddq <- ddq %>% mutate(network1 = network)
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term <= .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term <= .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <=.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <=.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN')
  
  
  ddf$RT_lag_bin <- factor(ddf$RT_lag_bin,levels=c('-2 std','-1 std','mean','+1 std','+2 std'))
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Explore-Explore-Entropy-rt_csv-divergence-D-nodiff-ss.pdf',width=15,height=15)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin(position=position_dodge(width=1)) +
    facet_grid(RT_lag_bin~dataset) + ylab('RT Swings (AU)') + xlab('Entropy Response') +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    xlab('Entropy Response') + scale_y_reverse() +
    scale_color_manual(values = pal1) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
}

if (do_value_plot){
  
  # plot nice figure
  library(wesanderson)
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1 = palette()
  pal1[1] <- pal[2]
  pal1[2] <- pal[1]
  pal1[3] <- pal[3]
  i <- 2 # value
  toalign <- 'clock' #
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-ranslopes-replication-clock-Explore-pred-rt_csv_sc-both-rs-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_r <- ddq
  rm(ddq)
  model_str <- paste0('-vmPFC-network-ranslopes-replication-clock-Explore-pred-rt_csv_sc-both-rs-',i,'.Rdata')
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
  
  ddq_r_emm <- ddq_r_emm %>% mutate(dataset = 'non-fMRI Replication') %>% rename(rt_vmax_lag_sc_r = rt_vmax_lag_sc, slrs_r=subj_level_rand_slope)
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
  ddq_r <- ddq_r %>% mutate(dataset = 'non-fMRI Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
  ddq <- ddq%>% mutate(network1 = case_when(network=='DMN'~'DMN', network=='CTR'~'CTR',network=='LIM'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term <= .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term <= .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <=.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <=.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  ddq <- ddq %>% filter(network1=='DMN' | network1=='CTR' | network1=='LIM')
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Explore-Value-rt_vmax-convergence-DCL-main.pdf',width=14,height=7)
  gg1 <- ggplot(ddq,aes(x=network1,y=estimate,color=network1)) +
    geom_violin() + facet_wrap(~dataset) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr)) +
    theme(axis.text.x = element_text(size=25), axis.text.y = element_text(size=25), axis.title = element_text(size=35), strip.text.x = element_text(size=35),strip.text.y = element_text(size=35),legend.title=element_text(size=35),legend.text=element_text(size=30),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    ylab('Convergence on Best RT (AU)') + xlab('') +
    scale_color_manual(values = pal1) + geom_hline(yintercept=0) +
    guides(color=guide_legend(title='Network'))
  #gg1 <- ggplot(ddq, aes(x=network,y=estimate, color=network)) + geom_bar(aes(alpha=p_level_fdr),stat='identity') + geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) + facet_grid(~dataset)
  print(gg1)
  dev.off()
  
  ddq_r_emm <- ddq_r_emm %>% select(!slrs_r & !rt_vmax_lag_sc_r)
  ddq_f_emm <- ddq_f_emm %>% select(!slrs_f & !rt_vmax_lag_sc_f)
  ddf <- rbind(ddq_r_emm,ddq_f_emm)
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std'))# & RT_vmax_bin=='mean')
  ddf <- ddf %>% mutate(network1 = case_when(network=='DMN'~'DMN',network=='CTR'~'CTR', network=='LIM'~'LIM')) %>% select(!std.error)
  
  #ddf <- ddf %>% filter(`vPFC Entropy response`=='10th %ile' | `vPFC Entropy response`=='90th %ile')
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std'))
  ddf <- ddf %>% mutate(network1 = case_when(network=='DMN'~'DMN',network=='CTR'~'CTR', network=='LIM'~'LIM'))
  ddf <- ddf %>% filter(network1=='DMN' | network1=='CTR' | network1=='LIM')
  ddf <- inner_join(ddf,ddq,by=c('evt_time','dataset','network1'))
  ddf$RT_vmax_bin <- factor(ddf$RT_vmax_bin,levels=c('-2 std','-1 std','mean','+1 std','+2 std'))
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Explore-Value-rt_vmax-convergence-DCL-nodiff.pdf',width=16,height=12)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin() + 
    facet_grid(RT_vmax_bin~dataset) + ylab('Convergence on Best RT (AU)') +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    ylab('Convergence on Best RT (AU)') + xlab('Value Response') +
    scale_color_manual(values = pal1) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
  
  rm(ddf)
  rm(ddq)
  
  i <- 2 # value
  toalign <- 'clock' #
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  model_str <- paste0('-vmPFC-network-ranslopes-replication-clock-Explore-pred-rt_csv_sc-both-rs-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_r <- ddq
  rm(ddq)
  model_str <- paste0('-vmPFC-network-ranslopes-replication-clock-Explore-pred-rt_csv_sc-both-rs-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  ddq_f <- ddq
  rm(ddq)
  
  ddq_r_emm <- ddq_r$emmeans_list$RTxO
  ddq_r_emt <- ddq_r$emtrends_list$RTxO
  ddq_f_emm <- ddq_f$emmeans_list$RTxO
  ddq_f_emt <- ddq_f$emtrends_list$RTxO
  ddq_r <- ddq_r$coef_df_reml
  ddq_f <- ddq_f$coef_df_reml
  
  ddq_r <- ddq_r %>% filter(effect=='fixed') 
  ddq_f <- ddq_f %>% filter(effect=='fixed') 
  
  ddq_r_emm <- ddq_r_emm %>% mutate(dataset = 'non-fMRI Replication') %>% rename(rt_lag_sc_r = rt_lag_sc, slrs_r=subj_level_rand_slope)
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
  ddq_r <- ddq_r %>% mutate(dataset = 'non-fMRI Replication')
  ddq_f <- ddq_f %>% mutate(dataset = 'fMRI')
  ddq <- rbind(ddq_r,ddq_f)
  ddq <- ddq %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
  ddq <- ddq%>% mutate(network1 = case_when(network=='DMN'~'DMN', network=='CTR'~'CTR',network=='LIM'~'LIM'))
  ddq <- ddq  %>% 
    mutate(p_level_fdr = case_when(
      padj_fdr_term > .05 ~ 'NS',
      padj_fdr_term <= .05 & padj_fdr_term > .01 ~ 'p < 0.05',
      padj_fdr_term <= .01 & padj_fdr_term > .001 ~ 'p < 0.01',
      padj_fdr_term <=.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
      padj_fdr_term <=.0001 ~ 'p < 0.0001'
    ))
  ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
  
  ddq <- ddq %>% filter(network1=='DMN' | network1=='CTR' | network1=='LIM')
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Explore-Value-rt_csv-divergence-DCL-main.pdf',width=14,height=7)
  gg1 <- ggplot(ddq,aes(x=network1,y=estimate,color=network1)) +
    geom_violin() + facet_wrap(~dataset) + 
    geom_point(aes(alpha=p_level_fdr,size=p_level_fdr)) +
    theme(axis.text.x = element_text(size=25), axis.text.y = element_text(size=25), axis.title = element_text(size=35), strip.text.x = element_text(size=35),strip.text.y = element_text(size=35),legend.title=element_text(size=35),legend.text=element_text(size=30),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) + scale_y_reverse() +
    ylab('RT swings (AU)') + xlab('') +
    scale_color_manual(values = pal1) + geom_hline(yintercept=0) +
    guides(color=guide_legend(title='Network'))
  #gg1 <- ggplot(ddq, aes(x=network,y=estimate, color=network)) + geom_bar(aes(alpha=p_level_fdr),stat='identity') + geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) + facet_grid(~dataset)
  print(gg1)
  dev.off()
  
  
  
  ddq_r_emm <- ddq_r_emm %>% select(!slrs_r & !rt_lag_sc_r)
  ddq_f_emm <- ddq_f_emm %>% select(!slrs_f & !rt_lag_sc_f)
  
  ddf <- rbind(ddq_r_emm,ddq_f_emm)
  #ddf <- ddf %>% filter(`vPFC Entropy response`=='10th %ile' | `vPFC Entropy response`=='90th %ile')
  ddf <- ddf %>% filter((vPFC_Entropy_random_slope=='-2 std' | vPFC_Entropy_random_slope=='+2 std'))
  ddf <- ddf %>% mutate(network1 = case_when(network=='DMN'~'DMN', network=='CTR'~'CTR',network=='LIM'~'LIM')) %>% select(!std.error)
  ddf <- ddf %>% filter(network1=='DMN' | network1=='CTR' | network1=='LIM')
  ddf <- inner_join(ddf,ddq,by=c('evt_time','dataset','network1'))
  ddf$RT_lag_bin <- factor(ddf$RT_lag_bin,levels=c('-2 std','-1 std','mean','+1 std','+2 std'))
  ddf <- ddf %>% filter(last_outcome=='Reward')
  #ddf <- ddf %>% filter(RT_lag_bin=='mean')
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
  pdf('Explore-Value-rt_csv-divergence-DCL-nodiff.pdf',width=15,height=15)
  gg1 <- ggplot(ddf,aes(x=vPFC_Entropy_random_slope,y=estimate.x,color=network1)) +
    geom_violin(position=position_dodge(width=1)) +
    facet_grid(RT_lag_bin~dataset) + ylab('RT Swings (AU)') + xlab('Value Response') +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title = element_text(size=30), strip.text.x = element_text(size=30),strip.text.y = element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20),legend.spacing.y=unit(1.0,'cm')) +
    guides(fill = guide_legend(byrow=TRUE)) +
    xlab('Value Response') + scale_y_reverse() +
    scale_color_manual(values = pal1) +
    guides(color=guide_legend(title='Network'))
  print(gg1)
  dev.off()
  
  
  
}
