
library(stringr)
library(fmri.pipeline)
library(tidyverse)
library(pracma)

ncores = 26
do_vPFC = TRUE
do_network = TRUE
do_symmetry = FALSE
do_HC = TRUE
do_vPFC_HC = TRUE
do_vPFC_HC_fb = FALSE
do_HC_anatomy = FALSE
if (do_vPFC){

#load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')  
md <- read_csv('/Volumes/bierka_root/datamesh/PROC/EXP/MEDuSA/2024-TRdiv2/clock_aligned_explore_vmPFC.csv.gz')
md <- md %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
md$id <- as.character(md$id)
    
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
df <- df %>% select(iti_ideal,condition_trial_neg_inv_sc,rt_lag_sc,iti_prev,iti_sc,v_entropy_wi_change,iti_prev_sc,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
#df$id <- as.character(df$id)
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

scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = as.vector(t(scale(difftime(scaninfo$scan_date,min(scaninfo$scan_date)))))) %>% select(!registration_redcapid)

demo <- inner_join(demo,scaninfo,by='id')

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
Q <- Q %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                  run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                  run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))
Q <- Q %>% filter(group=='HC')
Q <- Q %>% filter(!is.na(rewFunc))
Q$scan_which <- relevel(as.factor(Q$scan_which),ref='P1')
#Q <- Q %>% filter(trial > 10)

rm(decode_formula)
decode_formula <- NULL
#decode_formula[[1]] = formula(~ v_max_wi*group + (1 |id/run))
#decode_formula[[2]] = formula(~ age + gender + wtar + education_yrs + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc |id/run))
#decode_formula[[3]] = formula(~ age + gender + wtar + education_yrs + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | id) + (1 | run))
#decode_formula[[4]] = formula(~ age + gender + wtar + education_yrs + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | run) + (1| id))

#decode_formula[[2]] = formula(~ v_max_wi + trial_neg_inv_sc + rt_csv_sc + gender + wtar + education_yrs + age + iti_sc + iti_prev_sc + last_outcome*iti_prev_sc + outcome*rt_csv_sc + (1 |id/run))
#decode_formula[[6]] = formula(~ age + gender + wtar + education_yrs + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi |id/run))
#decode_formula[[7]] = formula(~ age + gender + wtar + education_yrs + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | id) + (1 | run))
#decode_formula[[8]] = formula(~ age + gender + wtar + education_yrs + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | run) + (1| id))
decode_formula[[1]] = formula(~age + v_entropy_wi + last_outcome + run_trial0_neg_inv_sc + rt_lag_sc + iti_prev_sc + (1|id/run))
decode_formula[[2]] = formula(~age + v_max_wi + last_outcome + run_trial0_neg_inv_sc + rt_lag_sc + iti_prev_sc + (1|id/run))
decode_formula[[3]] = formula(~age + v_max_wi + v_entropy_wi + last_outcome + run_trial0_neg_inv_sc + rt_lag_sc + iti_prev_sc + (1|id/run))
decode_formula[[4]] = formula(~v_max_wi + (1|id))
decode_formula[[5]] = formula(~v_entropy_wi+ (1|id))
# decode_formula[[1]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1|id/run))
# decode_formula[[2]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc |id/run))
# decode_formula[[3]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | id) + (1 | run))
# decode_formula[[4]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | run) + (1| id))
# 
# decode_formula[[5]] = formula(~ age + gender + v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1|id/run))
# decode_formula[[6]] = formula(~ age + gender + v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi |id/run))
# decode_formula[[7]] = formula(~ age + gender + v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | id) + (1 | run))
# decode_formula[[8]] = formula(~ age + gender + v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | run) + (1| id))

# decode_formula[[1]] = formula(~ v_entropy_sc + (1|id))
# decode_formula[[2]] = formula(~ v_max_wi + (1|id))

qT <- c(-0.8,0.46)
if (do_network){
  splits = c('evt_time','network')
  #source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-Explore-vmPFC-network-clock-HConly-trial1-10included-RTcorrected-scannersplit-',i,'.Rdata'))
  }
}
if (do_symmetry){
  splits = c('evt_time','symmetry_group')
  #source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-Explore-vmPFC-symmetry-clock-HConly-trial1-10included-RTcorrected-',i,'.Rdata'))
  }
}
}

if (do_HC){
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock_TRdiv2.Rdata')
  #hc <- read_csv('/Volumes/Users/Bea/StriatumHippThalamus/clock_aligned_striatum_hipp_thalamus.csv.gz')
  #hc <- hc %>% mutate(run1 = as.integer(str_sub(run,4,4))) %>% select(!run) %>% rename(run=run1)
  #hc <- hc %>% filter(atlas_value %in% c(223,224,225,226,227,228,229,230))
  #hc #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_fb.Rdata')
  #hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
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
  #hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2  hc <- hc %>% rename(decon_mean=decon_mean1)
  #hc <- hc %>% group_by(id,run,HC_region) %>% mutate(HCwithin = scale(decon_mean),HCbetween=mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  hc <- hc %>% rename(HCwithin = decon_mean)
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
  df <- df %>% select(total_earnings,condition_trial_neg_inv_sc,trial,iti_ideal,rt_lag_sc,iti_lag_sc,iti_prev,iti_sc,v_entropy_wi_change_lag,v_entropy_wi_change,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
  df$id <- as.character(df$id)
  df <- df %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                      trial > 40 & trial <= 80 ~ 2,
                                      trial > 80 & trial <=120 ~ 3, 
                                      trial > 120 & trial <=160 ~ 4,
                                      trial > 160 & trial <=200 ~ 5,
                                      trial > 200 & trial <=240 ~ 6))
  hc$id <- as.character(hc$id)
  Q <- inner_join(hc,df,by=c('id','run','trial'))
  rm(hc)
  Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCwithin[Q$evt_time < -Q$iti_prev] = NA

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
  
  scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
  scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = as.vector(t(scale(difftime(scaninfo$scan_date,min(scaninfo$scan_date)))))) %>% select(!registration_redcapid)
  
  demo <- inner_join(demo,scaninfo,by='id')
  
  Q <- merge(demo,Q,by='id')
  #Q <- Q %>% filter(group!='ATT')
  Q$group <- relevel(factor(Q$group),ref='HC')
  Q <- Q %>% filter(group=='HC')
  #Q <- Q %>% filter(trial > 10)
  Q <- Q %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                           trial > 40 & trial <= 80 ~ trial-40,
                                           trial > 80 & trial <=120 ~ trial-80, 
                                           trial > 120 & trial <=160 ~ trial-120,
                                           trial > 160 & trial <=200 ~ trial-160,
                                           trial > 200 & trial <=240 ~ trial-200))
  Q <- Q %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                    run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                    run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))
  Q$scan_which <- relevel(as.factor(Q$scan_which),ref='P1')
  rm(decode_formula)
  decode_formula <- NULL
  #decode_formula[[1]] = formula(~ v_entropy_wi + last_outcome + rt_csv_sc + iti_lag_sc + (1|id/run))
  #decode_formula[[2]] = formula(~ v_max_wi + last_outcome + rt_csv_sc +  iti_lag_sc + (1|id/run))
  #decode_formula[[1]] = formula(~ v_entropy_wi_change +  (1|id/run))
  #decode_formula[[1]] = formula(~ group + condition_trial_neg_inv_sc + v_max_wi + last_outcome + age + gender + iti_lag_sc + rt_lag_sc + (1|id/run)) 
  decode_formula[[1]] = formula(~age + run_trial0_neg_inv_sc + gender + v_max_wi + rt_lag_sc + iti_lag_sc + last_outcome  + (1|id/run))
  decode_formula[[2]] = formula(~age + run_trial0_neg_inv_sc + gender + v_entropy_wi + rt_lag_sc + iti_lag_sc + last_outcome + (1|id/run))
  decode_formula[[3]] = formula(~v_max_wi +  (1|id/run))
  decode_formula[[4]] = formula(~v_entropy_wi +  (1|id/run))
  qT <- c(-0.8,0.46)
  splits = c('evt_time','HC_region','scan_which')
  #source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "HCwithin", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))#,
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-Explore-HC-region-clock-HConly-trial1-10included-RTcorrected-scannersplit-',i,'.Rdata'))
  }
}

if (do_vPFC_HC){
  
  #source('/Volumes/Users/Andrew/MEDuSA_data_Explore/get_trial_data_explore.R')
  
  #load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')
  md <- read_csv('/Volumes/bierka_root/datamesh/PROC/EXP/MEDuSA/2024-TRdiv2/clock_aligned_explore_vmPFC.csv.gz')
  md <- md %>% filter(evt_time > -4 & evt_time < 4)
  md <- md %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
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
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock_TRdiv2.Rdata')
  #hc <- read_csv('/Volumes/Users/Bea/StriatumHippThalamus/clock_aligned_striatum_hipp_thalamus.csv.gz')
  #hc <- hc %>% mutate(run1 = as.integer(str_sub(run,4,4))) %>% select(!run) %>% rename(run=run1)
  #hc <- hc %>% filter(atlas_value %in% c(223,224,225,226,227,228,229,230))
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
  
  df <- df %>% select(total_earnings_split,condition_trial_neg_inv_sc,iti_ideal,rt_lag_sc,iti_lag_sc,iti_prev,iti_sc,v_entropy_wi_change,outcome,ev_sc,v_chosen_sc,last_outcome,outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
  df$id <- as.character(df$id)
  md$id <- as.character(md$id)
  Q <- full_join(md,df,by=c('id','run','trial'))
  Q <- Q %>% rename(vmPFC_decon = decon_mean) 
  rm(md)
  hc$id <- as.character(hc$id)
  Q <- Q %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                      trial > 40 & trial <= 80 ~ 2,
                                      trial > 80 & trial <=120 ~ 3, 
                                      trial > 120 & trial <=160 ~ 4,
                                      trial > 160 & trial <=200 ~ 5,
                                      trial > 200 & trial <=240 ~ 6))
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
  
  scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
  scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = as.vector(t(scale(difftime(scaninfo$scan_date,min(scaninfo$scan_date)))))) %>% select(!registration_redcapid)
  
  demo <- inner_join(demo,scaninfo,by='id')
  
  Q <- merge(demo,Q,by='id')
  #Q <- Q %>% filter(total_earnings_split=='richer')
  Q$group <- relevel(factor(Q$group),ref='HC')
  Q <- Q %>% filter(group!='ATT')
  Q <- Q %>% filter(group=='HC')
  Q <- Q %>% filter(!is.na(rewFunc))
  #Q <- Q %>% filter(trial > 10)
  Q <- Q %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                           trial > 40 & trial <= 80 ~ trial-40,
                                           trial > 80 & trial <=120 ~ trial-80, 
                                           trial > 120 & trial <=160 ~ trial-120,
                                           trial > 160 & trial <=200 ~ trial-160,
                                           trial > 200 & trial <=240 ~ trial-200))
  Q <- Q %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                    run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                    run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))
  Q$scan_which <- relevel(as.factor(Q$scan_which),ref='P1')
  rm(decode_formula)
  decode_formula <- NULL
  #decode_formula[[1]] = formula(~age * HCwithin + gender * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + v_max_wi * HCwithin + v_entropy_wi_change_lag * HCwithin + rt_csv_sc * HCwithin + iti_lag_sc * HCwithin + iti_sc * HCwithin + last_outcome * HCwithin + outcome*HCwithin + rt_vmax_change_sc * HCwithin +  HCbetween + (1 | id/run)) 
  #decode_formula[[1]] = formula(~age*HCwithin + v_entropy_wi*HCwithin + (1|id/run))
  #decode_formula[[2]] = formula(~age*HCwithin + v_max_wi*HCwithin + (1|id/run))
  #decode_formula[[3]] = formula(~age*HCwithin + gender*HCwithin + v_entropy_wi*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1|id/run))
  decode_formula[[1]] = formula(~trial_neg_inv_sc*HCwithin + age*HCwithin + run_trial0_neg_inv_sc*HCwithin + gender*HCwithin + v_max_wi*HCwithin + rt_lag_sc*HCwithin + HCbetween + last_outcome*HCwithin  + (1|id/run))
  decode_formula[[2]] = formula(~trial_neg_inv_sc*HCwithin + age*HCwithin + run_trial0_neg_inv_sc*HCwithin + gender*HCwithin + v_entropy_wi*HCwithin + rt_lag_sc*HCwithin + HCbetween + last_outcome*HCwithin + (1|id/run))
  decode_formula[[3]] = formula(~trial_neg_inv_sc*HCwithin + age*HCwithin + run_trial0_neg_inv_sc*HCwithin + gender*HCwithin + v_max_wi*HCwithin + v_entropy_wi*HCwithin + rt_lag_sc*HCwithin + HCbetween + last_outcome*HCwithin + (1|id/run))
  decode_formula[[4]] = formula(~HCwithin*v_entropy_wi + HCbetween  + (1|id/run))
  #decode_formula[[4]] = formula(~trial_neg_inv_sc*HCwithin + v_max_wi*HCwithin + HCbetween + (1|id/run))
  #decode_formula[[5]] = formula(~trial_neg_inv_sc*HCwithin + v_entropy_wi*HCwithin + HCbetween + (1|id/run))
  #decode_formula[[6]] = formula(~trial_neg_inv_sc*HCwithin + v_max_wi*HCwithin + outcome + HCbetween + (1|id/run))
  #decode_formula[[7]] = formula(~trial_neg_inv_sc*HCwithin + v_entropy_wi*HCwithin + outcome + HCbetween + (1|id/run))
  #decode_formula[[8]] = formula(~trial_neg_inv_sc*HCwithin + v_entropy_wi*HCwithin + v_entropy_wi_change + HCbetween + (1|id/run))
  #decode_formula[[2]] = formula(~ v_entropy_wi*HCwithin + HCbetween + (1|id/run))
  # decode_formula[[3]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_c  (1|id/run))
  # decode_formula[[4]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc +   (1|id/run))
  # 
  # decode_formula[[5]] = formula(~ (1|id/run))
  # decode_formula[[6]] = formula(~  (1|id/run))
  # decode_formula[[7]] = formula(~  (1|id/run))
  # decode_formula[[8]] = formula(~  (1|id/run))
  
  # decode_formula[[1]] = formula(~ v_entropy_sc + (1|id))
  # decode_formula[[2]] = formula(~ v_max_wi + (1|id))
  
  qT2 <- c(-2.62,-0.544,0.372, 0.477)
  qT1 <- c(-2.668, -0.12, 0.11, 0.258, 0.288, 0.308, 0.323, 0.348)
  splits = c('evt_time','network','HC_region','scan_which')
  #source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)#,
                      # emmeans_spec = list(
                      #   A = list(outcome='vmPFC_decon',model_name='model1',
                      #            specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
                      #   W = list(outcome='vmPFC_decon',model_name='model1',
                      #            specs=formula(~wtar),at=list(wtar=c(-1,-0.5,0,0.5,1))),
                      #   H = list(outcome='vmPFC_decon', model_name='model1',
                      #            specs=formula(~v_entropy_wi), at = list(v_entropy_wi=c(-2,-1,0,1,2))),
                      #   Y = list(outcome='vmPFC_decon',model_name='model1',
                      #            specs=formula(~education_yrs), at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
                      # )#,
                      # emtrends_spec = list(
                      # #   HxW = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
                      # #              specs = formula(~v_entropy_wi:wtar),at=list(wtar = c(-1,-0.5,0,0.5,1))),
                      #   T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                      #               specs=formula(~trial_neg_inv_sc:HCwithin),at=list(trial_neg_inv_sc=qT1)),
                      #   T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                      #               specs=formula(~run_trial0_neg_inv_sc:HCwithin),at=list(run_trial0_neg_inv_sc=qT2))
                      # #   H = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
                      # #            specs = formula(~v_entropy_wi),at=list(v_entropy_wi=c(-2,-1,0,1,2))),
                      # #   HxY = list(outcome='vmPFC_decon',model_name='model1',var='v_entropy_wi',
                      # #              specs=formula(~v_entropy_wi:education_yrs),at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
                      # )
      )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-Explore-vPFC-HC-network-clock-HConly-trial_and_run_trial0-RTcorrected-scannersplit-',i,'.Rdata'))
  }
  
  
}

if (do_vPFC_HC_fb){
  
  #source('/Volumes/Users/Andrew/MEDuSA_data_Explore/get_trial_data_explore.R')
  
  load('/Volumes/Users/Andrew/MEDuSA_data_Explore/fb-vPFC.Rdata')
  md <- md %>% filter(evt_time > -4 & evt_time < 4)
  md <- md %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
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
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_fb.Rdata')
  #hc <- read_csv('/Volumes/Users/Bea/StriatumHippThalamus/clock_aligned_striatum_hipp_thalamus.csv.gz')
  #hc <- hc %>% mutate(run1 = as.integer(str_sub(run,4,4))) %>% select(!run) %>% rename(run=run1)
  #hc <- hc %>% filter(atlas_value %in% c(223,224,225,226,227,228,229,230))
  #hc <- hc %>% select(!decon_median & !decon_sd)
  #hc #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_fb.Rdata')
  #hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
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
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup()
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
  
  df <- df %>% select(total_earnings_split,condition_trial_neg_inv_sc,iti_ideal,rt_csv_sc,iti_lag_sc,iti_prev,iti_sc,v_entropy_wi_change_lag,outcome,ev_sc,v_chosen_sc,outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
  df$id <- as.character(df$id)
  md$id <- as.character(md$id)
  Q <- full_join(md,df,by=c('id','run','trial'))
  Q <- Q %>% rename(vmPFC_decon = decon_mean) 
  rm(md)
  hc$id <- as.character(hc$id)
  Q <- Q %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                      trial > 40 & trial <= 80 ~ 2,
                                      trial > 80 & trial <=120 ~ 3, 
                                      trial > 120 & trial <=160 ~ 4,
                                      trial > 160 & trial <=200 ~ 5,
                                      trial > 200 & trial <=240 ~ 6))
  Q <- inner_join(Q,hc,by=c('id','run','run_trial','evt_time'))
  rm(hc)
  Q$HCwithin[Q$evt_time > Q$iti_ideal] = NA
  Q$vmPFC_decon[Q$evt_time > Q$iti_ideal] = NA
  Q$HCbetween[Q$evt_time > Q$iti_ideal] = NA
  Q$HCwithin[Q$evt_time < -Q$rt_csv] = NA
  Q$vmPFC_decon[Q$evt_time < -Q$rt_csv] = NA
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
  Q <- Q %>% filter(group!='ATT')
  Q <- Q %>% filter(group=='HC')
  Q <- Q %>% filter(!is.na(rewFunc))
  Q <- Q %>% filter(trial > 10)
  
  rm(decode_formula)
  decode_formula <- NULL
  #decode_formula[[1]] = formula(~age * HCwithin + gender * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + v_max_wi * HCwithin + v_entropy_wi_change_lag * HCwithin + rt_csv_sc * HCwithin + iti_lag_sc * HCwithin + iti_sc * HCwithin + last_outcome * HCwithin + outcome*HCwithin + rt_vmax_change_sc * HCwithin +  HCbetween + (1 | id/run)) 
  #decode_formula[[1]] = formula(~age*HCwithin + v_entropy_wi*HCwithin + (1|id/run))
  #decode_formula[[2]] = formula(~age*HCwithin + v_max_wi*HCwithin + (1|id/run))
  #decode_formula[[3]] = formula(~age*HCwithin + gender*HCwithin + v_entropy_wi*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1|id/run))
  decode_formula[[1]] = formula(~age*HCwithin + condition_trial_neg_inv_sc*HCwithin + gender*HCwithin + v_max_wi*HCwithin + rt_csv_sc*HCwithin + HCbetween + outcome*HCwithin + HCwithin*iti_sc + (1|id))
  decode_formula[[2]] = formula(~age*HCwithin + condition_trial_neg_inv_sc*HCwithin + gender*HCwithin + v_entropy_wi*HCwithin + rt_csv_sc*HCwithin + HCbetween + outcome*HCwithin + HCwithin*iti_sc + (1|id))
  decode_formula[[3]] = formula(~age*HCwithin + condition_trial_neg_inv_sc*HCwithin + gender*HCwithin + v_max_wi*HCwithin + v_entropy_wi*HCwithin + rt_csv_sc*HCwithin + HCbetween + outcome*HCwithin + HCwithin*iti_sc +(1|id))
  decode_formula[[4]] = formula(~v_max_wi*HCwithin + HCbetween + (1|id))
  #decode_formula[[2]] = formula(~ v_entropy_wi*HCwithin + HCbetween + (1|id/run))
  # decode_formula[[3]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_c  (1|id/run))
  # decode_formula[[4]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc +   (1|id/run))
  # 
  # decode_formula[[5]] = formula(~ (1|id/run))
  # decode_formula[[6]] = formula(~  (1|id/run))
  # decode_formula[[7]] = formula(~  (1|id/run))
  # decode_formula[[8]] = formula(~  (1|id/run))
  
  # decode_formula[[1]] = formula(~ v_entropy_sc + (1|id))
  # decode_formula[[2]] = formula(~ v_max_wi + (1|id))
  
  qT <- c(-0.8,0.46)
  splits = c('evt_time','network','HC_region')
  #source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))#,
    # emmeans_spec = list(
    #   A = list(outcome='vmPFC_decon',model_name='model1',
    #            specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
    #   W = list(outcome='vmPFC_decon',model_name='model1',
    #            specs=formula(~wtar),at=list(wtar=c(-1,-0.5,0,0.5,1))),
    #   H = list(outcome='vmPFC_decon', model_name='model1',
    #            specs=formula(~v_entropy_wi), at = list(v_entropy_wi=c(-2,-1,0,1,2))),
    #   Y = list(outcome='vmPFC_decon',model_name='model1',
    #            specs=formula(~education_yrs), at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
    # )#,
    # emtrends_spec = list(
    #   HxW = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
    #              specs = formula(~v_entropy_wi:wtar),at=list(wtar = c(-1,-0.5,0,0.5,1))),
    #   HxA = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
    #             specs = formula(~v_entropy_wi:age),at=list(age = c(-1,-0.5,0,0.5,1))),
    #   H = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
    #            specs = formula(~v_entropy_wi),at=list(v_entropy_wi=c(-2,-1,0,1,2))),
    #   HxY = list(outcome='vmPFC_decon',model_name='model1',var='v_entropy_wi',
    #              specs=formula(~v_entropy_wi:education_yrs),at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
    # )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-Explore-vPFC-HC-network-feedback-HConly-',i,'.Rdata'))
  }
  
  
}

if (do_HC_anatomy){
  
  #source('/Volumes/Users/Andrew/MEDuSA_data_Explore/get_trial_data_explore.R')
  
  #load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')
  md <- read_csv('/Volumes/bierka_root/datamesh/PROC/EXP/MEDuSA/2024-TRdiv2/clock_aligned_explore_vmPFC.csv.gz')
  md <- md %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
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
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock_TRdiv2.Rdata')
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
  hc <- hc %>% group_by(id,run,trial,HC_region) %>% arrange(evt_time) %>% mutate(HC_lag1 = lag(HCwithin,1),
                                                                                     HC_lag2 = lag(HCwithin,2),
                                                                                     HC_lag3 = lag(HCwithin,3),
                                                                                     HC_lag4 = lag(HCwithin,4)) %>% ungroup()
  hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
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
  md$id <- as.character(md$id)
  Q <- full_join(md,df,by=c('id','run','trial'))
  Q <- Q %>% rename(vmPFC_decon = decon_mean)
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
  Q <- Q %>% filter(group!='ATT')
  Q <- Q %>% filter(group=='HC')
  Q <- Q %>% filter(!is.na(rewFunc))
  #Q <- Q %>% filter(trial > 10)
  
  rm(decode_formula)
  decode_formula <- NULL
  #decode_formula[[1]] = formula(~age * HCwithin + gender * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + v_max_wi * HCwithin + v_entropy_wi_change_lag * HCwithin + rt_csv_sc * HCwithin + iti_lag_sc * HCwithin + iti_sc * HCwithin + last_outcome * HCwithin + outcome*HCwithin + rt_vmax_change_sc * HCwithin +  HCbetween + (1 | id/run)) 
  #decode_formula[[1]] = formula(~age*HCwithin + v_entropy_wi*HCwithin + (1|id/run))
  #decode_formula[[2]] = formula(~age*HCwithin + v_max_wi*HCwithin + (1|id/run))
  #decode_formula[[3]] = formula(~age*HCwithin + gender*HCwithin + v_entropy_wi*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1|id/run))
  decode_formula <- NULL
  decode_formula[[1]] <- formula(~ HCwithin + HCbetween + (1|id/run))
  decode_formula[[2]] <- formula(~ HCwithin + run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + HCbetween + (1 | id/run))
  decode_formula[[3]] <- formula(~ age*HCwithin + gender*HCwithin + HCwithin*run_trial0_neg_inv_sc + rt_lag_sc*HCwithin + iti_lag_sc*HCwithin + last_outcome*HCwithin + HCbetween + (1 | id/run))
  #decode_formula[[3]] <- formula(~ HC_lag1 + age + gender + HCwithin*run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + HCbetween + (1|id))
  #decode_formula[[4]] <- formula(~ HC_lag1 + age + gender + HCwithin*run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + HCbetween + (1 + HCwithin | id))
  #decode_formula[[5]] <- formula(~ HC_lag2 + age + gender + HCwithin*run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + HCbetween + (1|id))
  #decode_formula[[6]] <- formula(~ HC_lag2 + age + gender + HCwithin*run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + HCbetween + (1 + HCwithin | id))
  #decode_formula[[7]] <- formula(~ HC_lag1 + HC_lag2 + age + gender + HCwithin*run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + HCbetween + (1|id))
  #decode_formula[[8]] <- formula(~ HC_lag1 + HC_lag2 + age + gender + HCwithin*run_trial0_neg_inv_sc + rt_lag_sc + iti_lag_sc + HCbetween + (1 + HCwithin | id))
  #decode_formula[[2]] = formula(~ v_entropy_wi*HCwithin + HCbetween + (1|id/run))
  # decode_formula[[3]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_c  (1|id/run))
  # decode_formula[[4]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc +   (1|id/run))
  # 
  # decode_formula[[5]] = formula(~ (1|id/run))
  # decode_formula[[6]] = formula(~  (1|id/run))
  # decode_formula[[7]] = formula(~  (1|id/run))
  # decode_formula[[8]] = formula(~  (1|id/run))
  
  # decode_formula[[1]] = formula(~ v_entropy_sc + (1|id))
  # decode_formula[[2]] = formula(~ v_max_wi + (1|id))
  
  qT <- c(-0.8,0.46)
  splits = c('evt_time','network','HC_region')
  #source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))#,
    # emmeans_spec = list(
    #   A = list(outcome='vmPFC_decon',model_name='model1',
    #            specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
    #   W = list(outcome='vmPFC_decon',model_name='model1',
    #            specs=formula(~wtar),at=list(wtar=c(-1,-0.5,0,0.5,1))),
    #   H = list(outcome='vmPFC_decon', model_name='model1',
    #            specs=formula(~v_entropy_wi), at = list(v_entropy_wi=c(-2,-1,0,1,2))),
    #   Y = list(outcome='vmPFC_decon',model_name='model1',
    #            specs=formula(~education_yrs), at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
    # )#,
    # emtrends_spec = list(
    #   HxW = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
    #              specs = formula(~v_entropy_wi:wtar),at=list(wtar = c(-1,-0.5,0,0.5,1))),
    #   HxA = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
    #             specs = formula(~v_entropy_wi:age),at=list(age = c(-1,-0.5,0,0.5,1))),
    #   H = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
    #            specs = formula(~v_entropy_wi),at=list(v_entropy_wi=c(-2,-1,0,1,2))),
    #   HxY = list(outcome='vmPFC_decon',model_name='model1',var='v_entropy_wi',
    #              specs=formula(~v_entropy_wi:education_yrs),at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
    # )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-Explore-vPFC-HC-network-clock-HConly-anatomy-trial1-10included-',i,'.Rdata'))
  }
  
  
}