# 2024-09-10 AndyP
# vPFC -> HC & HC -> vPFC exploring age effects

library(stringr)

ncores = 26
do_vPFC = TRUE
do_network = TRUE
do_symmetry = FALSE
do_anat_clock = TRUE
do_anat_feedback = TRUE

####################################
### vPFC - HC clock - anatomy ###
####################################
if (do_anat_clock){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMclock_clock_Aug2023.Rdata')
  vmPFC <- clock_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(clock_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
  
  Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
  Q <- Q %>% select(!decon1)
  
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
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
  
  df <- df %>% select(id,run,trial,run_trial,trial_neg_inv_sc,iti_ideal,iti_prev,rt_csv,rt_csv_sc,iti_sc,rt_lag_sc,rt_csv_sc,rewFunc,iti_lag_sc)
  
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
  decode_formula[[1]] <- formula(~ HCwithin*age + HCbetween + (1|id/run))
  decode_formula[[2]] <- formula(~ HCwithin*age + HCwithin*last_outcome + HCbetween + (1 | id/run))
  decode_formula[[3]] <- formula(~ HCwithin*age + HCwithin*female + HCwithin*trial_neg_inv_sc + HCwithin*rt_lag_sc + HCwithin*iti_lag_sc + HCwithin*last_outcome + HCbetween + (1 | id/run))
  decode_formula[[4]] <- formula(~ HCwithin*female*age + HCwithin*trial_neg_inv_sc*age + HCwithin*rt_lag_sc*age + HCwithin*iti_lag_sc*age + HCwithin*last_outcome*age + HCbetween*age + (1 | id/run)) 
  #decode_formula[[4]] <- formula(~ HCwithin*age + HCbetween + (1 + HCwithin |id/run))
  #decode_formula[[5]] <- formula(~ HCwithin*age + HCwithin*last_outcome + HCbetween + (1 + HCwithin | id/run))
  #decode_formula[[6]] <- formula(~ HCwithin*age + HCwithin*female + HCwithin*trial_neg_inv_sc + HCwithin*rt_lag_sc + HCwithin*iti_lag_sc + HCwithin*last_outcome + HCbetween + (1 + HCwithin | id/run))
  qT <- c(-0.7,0.43)
  
  if (do_network){
    
    splits = c('evt_time','network','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        HC = list(outcome='vmPFC_decon', model_name='model1', 
                                  specs=formula(~HCwithin), at = list(HCwithin=c(-1.5,1.5))),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',
                                    specs=formula(~HCbetween:HCwithin), at=list(HCwithin=c(-1.5,1.5),HCbetween=c(-1.5,1.5)))
                      )
      )
      
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-age-',i,'.Rdata'))
    }
  }
  if (do_symmetry){
    
    splits = c('evt_time','symmetry_group','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
      )
      
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-clock-age-',i,'.Rdata'))
    }
  }
}

if (do_anat_feedback){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMClock_vmPFC_fb_Aug2023.csv.gz')
  vmPFC <- fb_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(fb_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_fb_Aug2023.Rdata')
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
  
  Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
  Q <- Q %>% select(!decon1)
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
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
  df <- df %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                       total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))
  
  #df <- df %>% filter(!is.na(rt_vmax_change_bin) | !is.na(v_entropy_wi_change_bin))
  df <- df %>% select(total_earnings_split,id,run,run_trial,trial_bin,iti_ideal,iti_prev,rt_csv,trial_neg_inv_sc,rt_csv_sc,v_entropy_wi,outcome,v_max_wi,rt_csv_sc,v_entropy_wi_change,score_sc,rt_bin,iti_sc,ev_sc,expl_longer,expl_shorter)
  Q <- merge(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$vmPFC_decon[Q$evt_time > Q$iti_ideal] = NA;
  Q$vmPFC_decon[Q$evt_time < -(Q$rt_csv)] = NA;
  Q$decon_mean[Q$evt_time > Q$iti_ideal] = NA;
  Q$decon_mean[Q$evt_time < -(Q$rt_csv)] = NA;
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
  
  decode_formula <- NULL
  decode_formula[[1]] <- formula(~ HCwithin*age + HCbetween + (1|id/run))
  decode_formula[[2]] <- formula(~ HCwithin*age + HCwithin*outcome + HCbetween + (1 | id/run))
  decode_formula[[3]] <- formula(~ HCwithin*age + HCwithin*female + HCwithin*trial_neg_inv_sc + HCwithin*rt_lag_sc + HCwithin*iti_sc + HCwithin*outcome + HCbetween + (1 | id/run))
  decode_formula[[4]] <- formula(~ HCwithin*female*age + HCwithin*trial_neg_inv_sc*age + HCwithin*rt_lag_sc*age + HCwithin*iti_sc*age + HCwithin*outcome*age + HCbetween*age + (1 | id/run)) 
  #decode_formula[[4]] <- formula(~ HCwithin*age + HCbetween + (1 + HCwithin |id/run))
  #decode_formula[[5]] <- formula(~ HCwithin*age + HCwithin*outcome + HCbetween + (1 + HCwithin | id/run))
  #decode_formula[[6]] <- formula(~ HCwithin*age + HCwithin*female + HCwithin*trial_neg_inv_sc + HCwithin*rt_lag_sc + HCwithin*iti_lag_sc + HCwithin*outcome + HCbetween + (1 + HCwithin | id/run))
  #decode_formula[[2]] <- NULL
  
  qT <- c(-0.7,0.43)
  
  
  if (do_network){
    
    splits = c('evt_time','network','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
      )
      
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-feedback-age-',i,'.Rdata'))
    }
  }
  if (do_symmetry){
    
    splits = c('evt_time','symmetry_group','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)
      )
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-feedback-age-',i,'.Rdata'))
    }
  }
}