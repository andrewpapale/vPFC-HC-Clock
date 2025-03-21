# 2025-03-20 AndyP
# test random slope effects on RT~RT_lag correlations to compare to MPlus MSEM

library(tidyverse)
library(lmerTest)

toalign <- 'clock'
repo_directory <- '~/clock_analysis'


#################
### MMClock  ####
#################

# from Subject_level_random_slopes.R
# decode_formula[[3]] = formula(~ age + female + trial_neg_inv_sc + v_max_wi + v_entropy_wi + rt_lag_sc  + iti_lag_sc + last_outcome + (1 + HCwithin  |id) + (1|run))
setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
model_str <- paste0('-vmPFC-HC-network-',toalign,'-ranslopes-nofixedeffect-noHCbetween-',3,'.Rdata') #3'rd iteration is HCwithin
model_str <- Sys.glob(paste0('*',model_str))
load(model_str)

qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')
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

Q2 <- Q2 %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)

Q2 <- inner_join(Q2,demo,by=c('id'))
Q2$age <- scale(Q2$age)
Q2 <- Q2 %>% mutate(sex = case_when(female == 1 ~ 'F',
                                    female == 0 ~ 'M'))
Q2$sex <- relevel(as.factor(Q2$sex),ref='F')

mDMNsex_mmclock <- lmerTest::lmer(data=Q2 %>% filter(network=='DMN'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRsex_mmclock <- lmerTest::lmer(data=Q2 %>% filter(network=='CTR'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMsex_mmclock <- lmerTest::lmer(data=Q2 %>% filter(network=='LIM'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))


#################
### BSocial  ####
#################




qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')
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

Q2 <- Q2 %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)

Q2 <- inner_join(Q2,demo,by=c('id'))
Q2$age <- scale(Q2$age)
Q2 <- Q2 %>% mutate(sex = case_when(female == 1 ~ 'F',
                                    female == 0 ~ 'M'))
Q2$sex <- relevel(as.factor(Q2$sex),ref='F')

mDMNsex_mmclock <- lmerTest::lmer(data=Q2 %>% filter(network=='DMN'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRsex_mmclock <- lmerTest::lmer(data=Q2 %>% filter(network=='CTR'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMsex_mmclock <- lmerTest::lmer(data=Q2 %>% filter(network=='LIM'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))


