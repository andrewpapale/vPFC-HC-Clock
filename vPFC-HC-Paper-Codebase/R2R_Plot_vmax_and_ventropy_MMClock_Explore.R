# 2025-04-24 AndyP
# Plot v_entropy and v_max timecourses for R2R

library(tidyverse)
library(fmri.pipeline)

ncores <- 26
repo_directory <- "~/clock_analysis"

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


df <- df %>% select(!run_trial) %>% rename(run_trial=condition_trial)

df <- df %>% select(!run) %>% mutate(run = case_when(trial <= 40 ~ 1, 
                                      trial > 40 & trial <= 80 ~ 2,
                                      trial > 80 & trial <=120 ~ 3, 
                                      trial > 120 & trial <=160 ~ 4,
                                      trial > 160 & trial <=200 ~ 5,
                                      trial > 200 & trial <=240 ~ 6))

df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  # run_trial <= 13 ~ 'Early',
  # run_trial > 13 & run_trial < 26 ~ 'Middle',
  # run_trial >=26 ~ 'Late',
  run_trial <= 5 ~ '1-5',
  run_trial > 5 & run_trial <= 10 ~ '6-10',
  run_trial > 10 & run_trial <= 15 ~ '11-15',
  run_trial > 15 & run_trial <= 20 ~ '16-20',
  run_trial > 20 & run_trial <= 25 ~ '21-25',
  run_trial > 25 & run_trial <= 30 ~ '26-30',
  run_trial > 30 & run_trial <= 35 ~ '30-35',
  run_trial > 35 & run_trial <= 40 ~ '36-40',
  run_trial > 40 & run_trial <= 45 ~ '40-45',
  run_trial > 45 & run_trial <= 50 ~ '46-50',
  run_trial > 50 & run_trial <= 55 ~ '51-55',
  run_trial > 55 & run_trial <= 63 ~ '55-63'
)))

df1 <- df
rm(df)

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
         rt_vmax_change_bin = case_when(
           abs_rt_vmax_change < 4/24 ~ 'No Change',
           abs_rt_vmax_change >= 4/24 ~ 'Change'
         ),
         v_entropy_wi_change_lag_bin = case_when(
           v_entropy_wi_change_lag < -0.5 ~ 'Decrease',
           v_entropy_wi_change_lag > 0.5  ~ 'Increase',
           v_entropy_wi_change_lag >= -0.5 & v_entropy_wi_change_lag <= 0.5 ~ 'No Change'
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
  # run_trial <= 15 ~ 'Early',
  # run_trial > 15 & run_trial < 30 ~ 'Middle',
  # run_trial >=30 ~ 'Late',
  run_trial <= 5 ~ '1-5',
  run_trial > 5 & run_trial <= 10 ~ '6-10',
  run_trial > 10 & run_trial <= 15 ~ '11-15',
  run_trial > 15 & run_trial <= 20 ~ '16-20',
  run_trial > 20 & run_trial <= 25 ~ '21-25',
  run_trial > 25 & run_trial <= 30 ~ '26-30',
  run_trial > 30 & run_trial <= 35 ~ '30-35',
  run_trial > 35 & run_trial <= 40 ~ '36-40',
  run_trial > 40 & run_trial <= 45 ~ '40-45',
  run_trial > 45 & run_trial <= 50 ~ '46-50',
  run_trial > 50 & run_trial <= 55 ~ '51-55',
  run_trial > 55 & run_trial <= 63 ~ '55-63'
)))

df2 <- df
rm(df)

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_meg')
df <- df %>% 
  group_by(id, run) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         iti_lag_sc = scale(iti_prev),
         score_sc = scale(score_csv),
         score_lag_sc = scale(lag(score_csv)),
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
         v_entropy_wi_change_lag_bin = case_when(
           v_entropy_wi_change_lag < -0.5 ~ 'Decrease',
           v_entropy_wi_change_lag > 0.5  ~ 'Increase',
           v_entropy_wi_change_lag >= -0.5 & v_entropy_wi_change_lag <= 0.5 ~ 'No Change'
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
  # run_trial <= 15 ~ 'Early',
  # run_trial > 15 & run_trial < 30 ~ 'Middle',
  # run_trial >=30 ~ 'Late',
  run_trial <= 5 ~ '1-5',
  run_trial > 5 & run_trial <= 10 ~ '6-10',
  run_trial > 10 & run_trial <= 15 ~ '11-15',
  run_trial > 15 & run_trial <= 20 ~ '16-20',
  run_trial > 20 & run_trial <= 25 ~ '21-25',
  run_trial > 25 & run_trial <= 30 ~ '26-30',
  run_trial > 30 & run_trial <= 35 ~ '30-35',
  run_trial > 35 & run_trial <= 40 ~ '36-40',
  run_trial > 40 & run_trial <= 45 ~ '40-45',
  run_trial > 45 & run_trial <= 50 ~ '46-50',
  run_trial > 50 & run_trial <= 55 ~ '51-55',
  run_trial > 55 & run_trial <= 63 ~ '55-63'
)))

df3 <- df
rm(df)

df1 <- df1 %>% mutate(dataset = 'Experiment 2')
df2 <- df2 %>% mutate(dataset = 'Experiment 1 - fMRI')
df3 <- df3 %>% mutate(dataset = 'Experiment 1 - MEG')


demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
demo <- demo %>% select(registration_redcapid,age,gender,registration_group,wtar,education_yrs)
demo <- demo %>% rename(id=registration_redcapid,group=registration_group)
demo$gender <- relevel(as.factor(demo$gender),ref='M')
demo$age <- scale(demo$age)
demo$wtar <- scale(demo$wtar)
demo$education_yrs <- scale(demo$education_yrs)

demo$id <- as.character(demo$id)

df1 <- inner_join(df1, demo, by=c('id'))

demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$female <- relevel(as.factor(demo$female),ref='0')
demo$age <- scale(demo$age)

demo$id <- as.character(demo$id)
df2$id <- as.character(df2$id)

df2 <- inner_join(df2, demo, by=c('id'))
df3 <- inner_join(df3, demo, by=c('id'))

df1 <- df1 %>% filter(group == 'HC')

df1 <- df1 %>% select(id,run,trial,run_trial,trial_bin,rt_swing,rt_csv,rt_lag_sc,last_outcome,outcome,ev,magnitude,v_max_above_median,rewFunc,v_entropy,v_entropy_wi,v_max,v_max_wi,dataset)
df2 <- df2 %>% select(id,run,trial,run_trial,trial_bin,rt_swing,rt_csv,rt_lag_sc,last_outcome,outcome,ev,magnitude,v_max_above_median,rewFunc,v_entropy,v_entropy_wi,v_max,v_max_wi,dataset)
df3 <- df3 %>% select(id,run,trial,run_trial,trial_bin,rt_swing,rt_csv,rt_lag_sc,last_outcome,outcome,ev,magnitude,v_max_above_median,rewFunc,v_entropy,v_entropy_wi,v_max,v_max_wi,dataset)

df <- rbind(df1,df2,df3)
#rm(df1,df2,df3)
df0 <- df %>% group_by(rewFunc, dataset, run_trial) %>% summarize(mean_v_max = mean(v_max,na.rm=TRUE), sd_vmax = sd(v_max,na.rm=TRUE),mean_v_entropy = mean(v_entropy,na.rm=TRUE), sd_ventropy = sd(v_entropy,na.rm=TRUE), N = n()) %>% ungroup()

ggplot(df0 %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV'), aes(x=run_trial,y=mean_v_max,ymin = mean_v_max-sd_vmax/sqrt(N), ymax = mean_v_max+sd_vmax/sqrt(N), color=dataset,group=dataset)) + geom_errorbar() + facet_grid(~rewFunc)
ggplot(df0 %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV'), aes(x=run_trial,y=mean_v_entropy,ymin = mean_v_entropy-sd_ventropy/sqrt(N), ymax = mean_v_entropy+sd_ventropy/sqrt(N), color=dataset,group=dataset)) + geom_errorbar() + geom_line() + facet_grid(~rewFunc)


df0 <- df %>% group_by(rewFunc, dataset, run_trial) %>% summarize(mean_v_max = mean(v_max_wi,na.rm=TRUE), sd_vmax = sd(v_max_wi,na.rm=TRUE),mean_v_entropy = mean(v_entropy_wi,na.rm=TRUE), sd_ventropy = sd(v_entropy_wi,na.rm=TRUE), N = n()) %>% ungroup()

#ggplot(df0 %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV'), aes(x=run_trial,y=mean_v_max,ymin = mean_v_max-sd_vmax/sqrt(N), ymax = mean_v_max+sd_vmax/sqrt(N), color=dataset,group=dataset)) + geom_errorbar() + geom_line() + facet_grid(~rewFunc)
#ggplot(df0 %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV'), aes(x=run_trial,y=mean_v_entropy,ymin = mean_v_entropy-sd_ventropy/sqrt(N), ymax = mean_v_entropy+sd_ventropy/sqrt(N), color=dataset,group=dataset)) + geom_errorbar() + geom_line() + facet_grid(~rewFunc)

# # The relative scaled v_max versus v_entropy could be interesting, the variables are correlated
# #ggplot(df0 %>% filter(rewFunc=='IEV' | rewFunc=='DEV'), aes(x=run_trial,y=mean_v_max-mean_v_entropy,color=dataset,group=dataset)) + geom_line() + facet_wrap(~rewFunc)
#
#df0 <- df %>% group_by(rewFunc,dataset,run_trial,last_outcome) %>% summarize(mean_rt_swing = mean(rt_swing,na.rm=TRUE), sd_rt_swing = sd(rt_swing,na.rm=TRUE), N=n()) %>% ungroup()
#df0 <- df0 %>% filter(!is.na(last_outcome))
#
#ggplot(df0 %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV'), aes(x=run_trial,y=mean_rt_swing,ymin=mean_rt_swing-sd_rt_swing/sqrt(N),ymax=mean_rt_swing+sd_rt_swing/sqrt(N),color=dataset,group=dataset)) + geom_line() + geom_errorbar() + facet_grid(rewFunc~last_outcome)
#


df$trial_bin = factor(df$trial_bin,levels = c('1-5','6-10','11-15','16-20','21-25','26-30','30-35','36-40','40-45','46-50','51-55','55-63'))

df <- df %>%  group_by(id,run,dataset) %>% 
  mutate(ev_above_median = case_when(mean(ev,na.rm=TRUE) > median(ev,na.rm=TRUE) ~ 'TRUE',
                                    mean(ev,na.rm=TRUE) < median(ev,na.rm=TRUE) ~ 'FALSE')) %>% ungroup()

#df$trial_bin <- factor(df$trial_bin,levels = c('Early','Middle','Late'))

df <- df %>%  group_by(id,run,dataset) %>% 
  mutate(magnitude_above_median = case_when(mean(magnitude,na.rm=TRUE) > median(magnitude,na.rm=TRUE) ~ 'TRUE',
                                            mean(magnitude,na.rm=TRUE) < median(magnitude,na.rm=TRUE) ~ 'FALSE')) %>% ungroup()

dfev <- df %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV') %>% filter(last_outcome == 'Reward')

# df0 <- df0 %>% group_by(dataset,trial_bin,rewFunc,ev_above_median) %>% summarize(mean_rt_swing = mean(rt_swing,na.rm=TRUE), sd_rt_swing = sd(rt_swing,na.rm=TRUE),N=n()) %>% ungroup()
# 
# # this is a boxplot
# ggplot(df %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV' & outcome == 'Reward'), aes(x=trial_bin,y=rt_swing,color=dataset,group=interaction(dataset,trial_bin))) + geom_boxplot(notch=TRUE,outlier.shape=NA) + facet_grid(rewFunc~ev_above_median) + ylim(0,2) + ggtitle('Expected Value Median Split')
# 
# df0 <- df %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV') %>% filter(outcome == 'Reward') %>% group_by(dataset,run_trial,rewFunc,v_max_above_median) %>% summarize(mean_rt_swing = mean(rt_swing,na.rm=TRUE), sd_rt_swing = sd(rt_swing,na.rm=TRUE),N=n()) %>% ungroup()
# 
# # this is a line plot
# ggplot(df0, aes(x=run_trial,y=mean_rt_swing,ymin=mean_rt_swing-sd_rt_swing/sqrt(N),ymax=mean_rt_swing+sd_rt_swing/sqrt(N),color=dataset,group=dataset)) + geom_line() + geom_errorbar() + facet_grid(rewFunc~v_max_above_median) + ggtitle('Vmax Median Split')
# 
dfmag <- df %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV') %>% filter(last_outcome == 'Reward')
# 
# df0 <- df %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV') %>% filter(last_outcome == 'Reward') %>% group_by(dataset,run_trial,rewFunc,magnitude_above_median) %>% summarize(mean_rt_swing = mean(rt_swing,na.rm=TRUE), sd_rt_swing = sd(rt_swing,na.rm=TRUE),N=n()) %>% ungroup()
# dfmag <- dfmag %>% filter(!is.na(magnitude_above_median))
# 
# # this is a line plot
# ggplot(df0, aes(x=trial_bin,y=mean_rt_swing,ymin=mean_rt_swing-sd_rt_swing/sqrt(N),ymax=mean_rt_swing+sd_rt_swing/sqrt(N),color=dataset,group=dataset)) + geom_line() + geom_errorbar() + facet_grid(rewFunc~magnitude_above_median) + ggtitle('Magnitude Median Split')

dfvmax <- df %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV') %>% filter(last_outcome == 'Reward')

# decode_formula <- NULL
# decode_formula[[1]] <- formula(~ rt_lag_sc*v_entropy_wi*last_outcome + rt_lag_sc*v_max_wi*last_outcome + rt_lag_sc*trial_bin*last_outcome + (1 | id))
# 

decode_formula <- NULL
decode_formula[[1]] <- formula(~ rt_lag_sc + (1|id/run))
splits1 = c('dataset','trial_bin','v_max_above_median')
splits2 = c('dataset','trial_bin','ev_above_median')
splits3 = c('dataset','trial_bin','magnitude_above_median')
splits4 = c('dataset','trial_bin')
nsplits = 4;
# for (j in 1:nsplits){
#   if (j==1){
#     curr_splits <- splits1
#     df0 <- df %>% filter(!is.na(ev_above_median))
#   } else if (j==2){
#     curr_splits <- splits2
#     df0 <- df %>% filter(!is.na(magnitude_above_median))
#   } else if (j==3){
#     curr_splits <- splits3
#   }
#   ddq <- mixed_by(df0, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[1]], split_on = curr_splits,return_models=TRUE,
#                   padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                   tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
#                   emtrends_spec = list(
#                     ventropy = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
#                               specs=formula(~rt_lag_sc:v_entropy_wi:last_outcome), at=list(v_entropy_wi = c(-2,-1,0,1,2))),
#                     vmax = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
#                                 specs=formula(~rt_lag_sc:v_max_wi:last_outcome), at=list(v_max_wi = c(-2,-1,0,1,2))),
#                     trial = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
#                                 specs=formula(~rt_lag_sc:trial_bin:last_outcome))
#                     )
#   )
#   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
#   save(ddq,file=paste0(curr_date,'-vPFC-HC-R2R-rtpred-',j,'.Rdata'))
#   
# }
# 
# load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2025-04-26-vPFC-HC-R2R-rtpred-3.Rdata')
# 
# ventropy <- ddq$emtrends_list$ventropy
# vmax <- ddq$emtrends_list$vmax
# trial <- ddq$emtrends_list$trial
# 
# ventropy <- ventropy %>% filter(v_entropy_wi %in% c(-1,1))%>% mutate(entropy = case_when(v_entropy_wi == -1 ~ 'low entropy',
#                                                                                v_entropy_wi == +1 ~ 'high entropy'))
# vmax <- vmax %>% filter(v_max_wi %in% c(-1,1)) %>% mutate(vmax = case_when(v_max_wi == -1 ~ 'low vmax',
#                                                                     v_max_wi == +1 ~ 'high vmax'))
# trial$trial_bin <- factor(trial$trial_bin,levels=c('Early','Middle','Late'))
# 
# ggplot(ventropy, aes(x=entropy,y=rt_lag_sc.trend,ymin = rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color = dataset, group=dataset)) + geom_errorbar() + geom_point() + facet_grid(dataset~last_outcome,scales = 'free_y') + scale_y_reverse()
# 
# ggplot(vmax, aes(x=vmax,y=rt_lag_sc.trend,ymin = rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color = dataset, group=dataset)) + geom_errorbar() + geom_point() + facet_grid(dataset~last_outcome,scales = 'free_y') + scale_y_reverse()
# 
# ggplot(trial, aes(x=trial_bin,y=rt_lag_sc.trend,ymin = rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color = dataset, group=dataset)) + geom_errorbar() + geom_point() + facet_grid(dataset~last_outcome,scales = 'free_y') + scale_y_reverse()
# 
# 
# 
# load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2025-04-26-vPFC-HC-R2R-rtpred-1.Rdata')
# 
# ventropy <- ddq$emtrends_list$ventropy
# vmax <- ddq$emtrends_list$vmax
# trial <- ddq$emtrends_list$trial
# 
# ventropy <- ventropy %>% filter(v_entropy_wi %in% c(-1,1))%>% mutate(entropy = case_when(v_entropy_wi == -1 ~ 'low entropy',
#                                                                                          v_entropy_wi == +1 ~ 'high entropy'))
# vmax <- vmax %>% filter(v_max_wi %in% c(-1,1)) %>% mutate(vmax = case_when(v_max_wi == -1 ~ 'low vmax',
#                                                                            v_max_wi == +1 ~ 'high vmax'))
# trial$trial_bin <- factor(trial$trial_bin,levels=c('Early','Middle','Late'))
# 
# ggplot(ventropy, aes(x=interaction(entropy,ev_above_median),y=rt_lag_sc.trend,ymin = rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color = ev_above_median, group=ev_above_median)) + geom_errorbar() + geom_point() + facet_grid(dataset~last_outcome,scales = 'free_y') + scale_y_reverse() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggplot(vmax, aes(x=interaction(vmax,ev_above_median),y=rt_lag_sc.trend,ymin = rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color = ev_above_median, group=ev_above_median)) + geom_errorbar() + geom_point() + facet_grid(dataset~last_outcome,scales = 'free_y') + scale_y_reverse() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggplot(trial, aes(x=interaction(trial_bin,ev_above_median),y=rt_lag_sc.trend,ymin = rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color = ev_above_median, group=ev_above_median)) + geom_errorbar() + geom_point() + facet_grid(dataset~last_outcome,scales = 'free_y') + scale_y_reverse() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

for (j in 4:nsplits){
  if (j==1){
    curr_splits <- splits1
    df0 <- dfvmax %>% filter(!is.na(v_max_above_median))
  } else if (j==2){
    curr_splits <- splits2
    df0 <- dfev %>% filter(!is.na(ev_above_median))
  } else if (j==3){
    curr_splits <- splits3
    df0 <- dfmag %>% filter(!is.na(magnitude_above_median))
  } else if (j==4){
    curr_splits <- splits4
    df0 <- df %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV') %>% filter(last_outcome == 'Reward')
  }
  ddq <- mixed_by(df0, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[1]], split_on = curr_splits,return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE))
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-vPFC-HC-R2R-rtpred-Omission',5,'.Rdata'))
}
ddf <- ddq$coef_df_reml %>% filter(effect=='fixed' & term=='rt_lag_sc')
ggplot(ddf, aes(x=trial_bin,y=estimate,ymin = estimate-std.error,ymax=estimate+std.error,color=dataset,group=dataset)) + geom_errorbar() + geom_line() + scale_y_reverse() + geom_point() + ylab('<-- less -- RT swings -- more -->') + facet_wrap(~v_max_above_median)
