



library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(fmri.pipeline)


repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26

rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMclock_clock_Aug2023.Rdata')
#load('/Volumes/Users/Andrew/MMclock_clock.Rdata') % drop2
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% 
  group_by(id, run) %>% 
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
         rt_vmax_change_sc = scale(rt_vmax_change),
         rt_vmax_change_bin = case_when(
           abs_rt_vmax_change < 4/24 ~ 'No Change',
           abs_rt_vmax_change >= 4/24 ~ 'Change'
         ),
         v_entropy_wi_change_lag_bin = case_when(
           v_entropy_wi_change_lag < -0.5 ~ 'Decrease',
           v_entropy_wi_change_lag > 0.5  ~ 'Increase',
           v_entropy_wi_change_lag >= -0.5 & v_entropy_wi_change_lag <= 0.5 ~ 'No Change'
         )) %>% arrange(id, run, run_trial) %>% mutate(log10kld4 = case_when(
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
#df <- df %>% filter(!is.na(rt_vmax_change_bin) | !is.na(v_entropy_wi_change_lag_bin))
df <- df %>% select(id,run,trial,run_trial,rt_lag_sc,rt_vmax_change_sc,abs_pe_max_lag_sc,v_entropy_wi,outcome,v_entropy_wi_change_lag, rt_vmax_lag_sc,iti_ideal, iti_prev, rt_csv, trial_bin,rewFunc,v_entropy_sc,expl_longer,rt_csv_sc, trial_neg_inv_sc,expl_shorter,rt_bin,trial_bin,last_outcome,v_max_wi,v_entropy_wi_change_lag,score_lag_sc,iti_sc,iti_lag_sc,ev_lag_sc)
Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
#Q$v_entropy_wi_change_lag_bin <- relevel(as.factor(Q$v_entropy_wi_change_lag_bin),ref='No Change')
#Q$rt_vmax_change_bin <- relevel(as.factor(Q$rt_vmax_change_bin),ref='No Change')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')

Q <- Q %>% mutate(sex = case_when(female==0 ~ 'M',
                                  female==1 ~ 'F'))
Q$sex <- relevel(as.factor(Q$sex),ref='M')
Q$age <- scale(Q$age)

load('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/mmclock_fmri_peaks.Rdata')
mmclock_fmri_peaks <- mmclock_fmri_peaks %>% group_by(id,trial) %>% summarize(num_peaks = length(peaks.peak)) %>% ungroup()
Q <- inner_join(Q,mmclock_fmri_peaks,by=c('id','trial'))
Q$num_peaks <- Q$num_peaks > 1
Q$num_peaks <- as.factor(Q$num_peaks)
Q$num_peaks <- relevel(Q$num_peaks,ref='FALSE')


rm(decode_formula)
decode_formula <- NULL
decode_formula[[1]] = formula(~ num_peaks + (1|id/run))
decode_formula[[2]] = formula(~ age + sex + num_peaks + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc +  (1 |id/run))
#decode_formula[[3]] = formula(~ age + sex + v_entropy_wi*sex + v_max_wi*sex + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc + (1|id/run))
#decode_formula[[4]] = formula(~ age + sex + v_entropy_wi*age + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc + (1|id/run))
#decode_formula[[5]] = formula(~ age + sex + v_max_wi*age + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc +  (1 |id/run))
#decode_formula[[6]] = formula(~ age + sex + v_entropy_wi*age + v_max_wi*age + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc + (1|id/run))
#decode_formula[[1]] = formula(~ sex + v_entropy_wi_change_lag*age + rt_vmax_lag_sc*age + abs_pe_max_lag_sc*age + rt_vmax_change_sc*age +  + rt_lag_sc + iti_lag_sc + (1|id/run))
#decode_formula[[2]] = formula(~ age + v_entropy_wi_change_lag*sex + rt_vmax_lag_sc*sex + abs_pe_max_lag_sc*sex + rt_vmax_change_sc*sex +  + rt_lag_sc + iti_lag_sc + (1|id/run))

qT <- c(-0.7,0.43)

splits = c('evt_time','network')
#source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  print(decode_formula[[i]])
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = decode_formula[[i]] , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emmeans_spec = list(
                     P = list(outcome='vmPFC_decon',model_name='model1',
                              specs=formula(~num_peaks)))
  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-network-clock-numpeaks-',i,'.Rdata'))
}

rm(decode_formula)
decode_formula <- NULL
decode_formula[[1]] = formula(~ ev_lag_sc + score_lag_sc + v_max_wi + (1|id/run))
decode_formula[[2]] = formula(~ age + sex + ev_lag_sc + score_lag_sc + v_max_wi + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc +  (1 |id/run))
#decode_formula[[3]] = formula(~ age + sex + v_entropy_wi*sex + v_max_wi*sex + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc + (1|id/run))
#decode_formula[[4]] = formula(~ age + sex + v_entropy_wi*age + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc + (1|id/run))
#decode_formula[[5]] = formula(~ age + sex + v_max_wi*age + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc +  (1 |id/run))
#decode_formula[[6]] = formula(~ age + sex + v_entropy_wi*age + v_max_wi*age + trial_neg_inv_sc + last_outcome + rt_lag_sc + iti_lag_sc + (1|id/run))
#decode_formula[[1]] = formula(~ sex + v_entropy_wi_change_lag*age + rt_vmax_lag_sc*age + abs_pe_max_lag_sc*age + rt_vmax_change_sc*age +  + rt_lag_sc + iti_lag_sc + (1|id/run))
#decode_formula[[2]] = formula(~ age + v_entropy_wi_change_lag*sex + rt_vmax_lag_sc*sex + abs_pe_max_lag_sc*sex + rt_vmax_change_sc*sex +  + rt_lag_sc + iti_lag_sc + (1|id/run))

qT <- c(-0.7,0.43)

splits = c('evt_time','network')
#source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  print(decode_formula[[i]])
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = decode_formula[[i]] , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emmeans_spec = list(
                    Vmax = list(outcome='vmPFC_decon',model_name='model1',
                             specs=formula(~v_max_wi), at=list(v_max_wi=c(-1.5,1.5))))
  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-network-clock-value-test-ev-score-',i,'.Rdata'))
}

