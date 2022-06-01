# 2022-05-12 AndyP
# split_by_condition
library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26
###################################
#####      vmPFC -feedback    #####
###################################

rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- fb_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(fb_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
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
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
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


rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ v_max_wi + trial_neg_inv_sc + pe_max_sc + v_entropy_wi + v_entropy_wi_change + (1|id))

splits = c('evt_time','network','rewFunc')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qPE = quantile(df$pe_max_sc,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emmeans_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', 
                                specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(qH))),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', 
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1',
                                specs=c('v_max_wi'), at=list(v_max_wi=c(qV))),
                    PE_HC = list(outcome='vmPFC_decon',model_name='model1',
                                 specs=c("pe_max_sc"), at = list(pe_max_sc=c(qPE))),
                    dH_HC = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=c('v_entropy_wi_change'), at=list(v_entropy_wi_change=c(qdH)))
                  )
  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-network-feedback-split-by-condition-simple-',i,'.Rdata'))
}

###################################
#####      vmPFC - clock      #####
###################################

rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
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
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
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


rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ v_max_wi + v_entropy_wi + trial_neg_inv_sc + v_entropy_wi_change_lag + (1|id))

splits = c('evt_time','network','rewFunc')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qPE = quantile(df$pe_max_sc,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emmeans_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', 
                                specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(qH))),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', 
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1',
                                specs=c('v_max_wi'), at=list(v_max_wi=c(qV))),
                    dH_HC = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=c('v_entropy_wi_change_lag'), at=list(v_entropy_wi_change_lag=c(qdH)))
                  ))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-network-clock-split-by-condition-simple-',i,'.Rdata'))
}


###################################
##### vmPFC - HC -feedback    #####
###################################
rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- fb_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(fb_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
hc <- fb_comb
hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
rm(fb_comb)
hc <- hc %>% mutate(
  HC_region = case_when(
    bin_num <= 8 ~ 'AH',
    bin_num >8 ~ 'PH'
  ),
)
hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()

Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
Q <- Q %>% select(!decon1)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(trial,v_max,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
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
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max_sc,pe_max_lag_sc)
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



rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin +  
                                pe_max_sc*HCwithin + v_entropy_wi*HCwithin + 
                                v_entropy_wi_change*HCwithin + 
                                HCbetween +
                                (1|id))

splits = c('evt_time','network','HC_region','rewFunc')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qPE = quantile(df$pe_max_sc,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emtrends_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(qH))),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=c('v_max_wi'), at=list(v_max_wi=c(qV))),
                    PE_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                 specs=c('pe_max_sc'), at=list(pe_max_sc=c(qPE))),
                    dH_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                 specs=c('v_entropy_wi_change'), at=list(v_entropy_wi_change=c(qdH)))
                  )
  )
  
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-feedback-split-by-condition-simple-',i,'.Rdata'))
}

###################################
##### vmPFC - HC - clock      #####
###################################
rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
hc <- clock_comb
hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)
hc <- hc %>% mutate(
  HC_region = case_when(
    bin_num <= 8 ~ 'AH',
    bin_num >8 ~ 'PH'
  ),
)
hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()

Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
Q <- Q %>% select(!decon1)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
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
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
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


rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + v_entropy_wi_change_lag*HCwithin +
                                HCbetween +
                                (1|id))

splits = c('evt_time','network','HC_region','rewFunc')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emtrends_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=c("v_entropy_wi"), at = list(v_entropy_wi=c(qH))),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=c('v_max_wi'), at=list(v_max_wi=c(qV))),
                    dH_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                 specs=c('v_entropy_wi_change_lag'), at=list(v_entropy_wi_change_lag=c(qdH)))
                  )
  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-split-by-condition-simple-',i,'.Rdata'))
}
