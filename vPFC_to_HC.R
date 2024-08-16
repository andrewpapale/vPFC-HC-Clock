# 2024-08-16 AndyP
# Test vPFC --> HC anatomy

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26
source("~/fmri.pipeline/R/mixed_by.R")

rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMclock_clock_Aug2023.Rdata')
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
vmPFC <- vmPFC %>% group_by(id,run,run_trial,evt_time,network) %>% summarize(decon1 = mean(vmPFC_decon,na.rm=TRUE)) %>% ungroup()
vmPFC <- vmPFC %>% group_by(id,run) %>% mutate(vmPFC_within = scale(decon1), vmPFC_between = mean(decon1,na.rm=TRUE)) %>% ungroup()
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
hc <- hc %>% rename(HC_decon = decon_mean)
#hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
#hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()

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

#Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
Q <- Q %>% group_by(network,HC_region) %>% mutate(vmPFC_between1 = scale(vmPFC_between)) %>% select(!vmPFC_between) %>% rename(vmPFC_between = vmPFC_between1)

decode_formula <- NULL
decode_formula[[1]] <- formula(~ vmPFC_within + vmPFC_between + (1|id/run))
decode_formula[[2]] <- formula(~ age + female + vmPFC_within + vmPFC_between + (1 + vmPFC_within | id/run))
qT <- c(-0.7,0.43)

splits = c('evt_time','network','HC_region')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "HC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emmeans_spec = list(
                    HC = list(outcome='HC_decon', model_name='model1', 
                              specs=formula(~vmPFC_within), at = list(vmPFC_within=c(-1.5,1.5)))
                  )
  )
  
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vPFC-to-HC-network-clock-anatomy-',i,'.Rdata'))
}



rm(Q)
load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')
vmPFC <- md %>% filter(evt_time > -5 & evt_time < 5)
rm(md)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
vmPFC <- vmPFC %>% group_by(id,run,run_trial,evt_time,network) %>% summarize(decon1 = mean(vmPFC_decon,na.rm=TRUE)) %>% ungroup()
vmPFC <- vmPFC %>% group_by(id,run) %>% mutate(vmPFC_within = scale(decon1), vmPFC_between = mean(decon1,na.rm=TRUE)) %>% ungroup()
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock.Rdata')
hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
hc <- hc %>% rename(HC_decon = decon_mean)
#hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
#hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()

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

#Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
Q <- Q %>% group_by(network,HC_region) %>% mutate(vmPFC_between1 = scale(vmPFC_between)) %>% select(!vmPFC_between) %>% rename(vmPFC_between = vmPFC_between1)

decode_formula <- NULL
decode_formula[[1]] <- formula(~ vmPFC_within + vmPFC_between + (1|id/run))
decode_formula[[2]] <- formula(~ age + female + vmPFC_within + vmPFC_between + (1 + vmPFC_within | id/run))
qT <- c(-0.7,0.43)

splits = c('evt_time','network','HC_region')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "HC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emmeans_spec = list(
                    HC = list(outcome='HC_decon', model_name='model1', 
                              specs=formula(~vmPFC_within), at = list(vmPFC_within=c(-1.5,1.5)))
                  )
  )
  
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-Explore-vPFC-to-HC-network-clock-anatomy-',i,'.Rdata'))
}
