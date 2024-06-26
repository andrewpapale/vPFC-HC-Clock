# 2024-06-25 AndyP
# Investigating RT~RTlag correlations in MMClock and Explore

library(tidyverse)
library(lmerTest)
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

df <- df %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                     total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))


df <- df %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                    trial > 40 & trial <= 80 ~ 2,
                                    trial > 80 & trial <=120 ~ 3, 
                                    trial > 120 & trial <=160 ~ 4,
                                    trial > 160 & trial <=200 ~ 5,
                                    trial > 200 & trial <=240 ~ 6))
df <- df %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                         trial > 40 & trial <= 80 ~ trial-40,
                                         trial > 80 & trial <=120 ~ trial-80, 
                                         trial > 120 & trial <=160 ~ trial-120,
                                         trial > 160 & trial <=200 ~ trial-160,
                                         trial > 200 & trial <=240 ~ trial-200))
df <- df %>% mutate(run_trial_c = run_trial-floor(run_trial/40.5),
                  run_trial_neg_inv = -(1 / run_trial_c) * 100,
                  run_trial_neg_inv_sc = as.vector(scale(run_trial_neg_inv)))

df <- df %>% mutate(dataset = 'Explore')

df$rt_bin <- relevel(as.factor(df$rt_bin),ref='-0.5')
df$trial_bin <- relevel(as.factor(df$trial_bin),ref='Middle')
df$expl_longer <- relevel(as.factor(df$expl_longer),ref='0')
df$expl_shorter <- relevel(as.factor(df$expl_shorter),ref='0')
df$rt_bin <- relevel(as.factor(df$rt_bin),ref='-0.5')
df$trial_bin <- relevel(as.factor(df$trial_bin),ref='Middle')

demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
demo$gender <- relevel(as.factor(demo$gender),ref='M')
demo$age_sc <- scale(demo$age)
demo$wtar_sc <- scale(demo$wtar)
demo$education_yrs_sc <- scale(demo$education_yrs)

df <- merge(df,demo,by='id')

lmer_control$optimizer <- 'bobyqa'
lmer_control$optCtrl = list(maxfun=2E5)
m1 <- lmerTest::lmer(rt_csv_sc ~ gender*rt_lag_sc + age_sc * rt_lag_sc + v_entropy_wi * rt_lag_sc * group_bin +  v_max_wi * rt_lag_sc * group_bin + trial_neg_inv_sc * rt_lag_sc *group_bin + group_bin * last_outcome * rt_lag_sc + (1 + rt_lag_sc | id),data=df, REML=TRUE,control=lmer_control)

df <- df %>% filter(Group=='Controls')

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df1 <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')

df1 <- df1 %>%
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         iti_lag_sc = scale(iti_lag),
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
df1 <- df1 %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df1 <- df1 %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df1 <- df1 %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df1 <- df1 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 15 ~ 'Early',
  run_trial > 15 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))

df1 <- df1 %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                     total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))

df1 <- df1 %>% mutate(dataset = 'MMClock')

df1$rt_bin <- relevel(as.factor(df1$rt_bin),ref='-0.5')
df1$trial_bin <- relevel(as.factor(df1$trial_bin),ref='Middle')
df1$expl_longer <- relevel(as.factor(df1$expl_longer),ref='0')
df1$expl_shorter <- relevel(as.factor(df1$expl_shorter),ref='0')
df1$rt_bin <- relevel(as.factor(df1$rt_bin),ref='-0.5')
df1$trial_bin <- relevel(as.factor(df1$trial_bin),ref='Middle')


demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
df1 <- inner_join(df1,demo,by=c('id'))
df1$female <- relevel(as.factor(df1$female),ref='0')
df1 <- df1 %>% mutate(gender = case_when(female==1 ~ 'F',
                                       female==0 ~ 'M'))
df1$gender <- relevel(as.factor(df1$gender),ref='M')
df1$age_sc <- scale(df1$age)


common_cols = intersect(colnames(df),colnames(df1))
ddf <- rbind(subset(df,select = common_cols),subset(df1,select = common_cols))


ddf <- ddf %>% filter(rewFunc == 'IEV' | rewFunc=='DEV')

ddf$dataset <- relevel(as.factor(ddf$dataset),ref='MMClock')

ddf$age_sc_merged <- scale(ddf$age)

df <- df %>% mutate(group_bin = case_when(registration_lethality=='ll' ~ 'LL', registration_lethality=='hl' ~ 'HL', is.na(registration_lethality) ~ 'NON'))
df$group_bin = relevel(as.factor(df$group_bin),ref='NON')

# LL and HL affected by task complexity or information load (entropy term).  For LL this depends on interacting RTlag:group:Vmax
m1 <- lmerTest::lmer(rt_csv_sc ~ age_sc*rt_lag_sc + v_entropy_wi*rt_lag_sc*group_bin + v_max_wi*rt_lag_sc*group_bin  + trial_neg_inv_sc*rt_lag_sc + group_bin*last_outcome*rt_lag_sc + (1+rt_lag_sc | id),data=df)

m2 <- lmerTest::lmer(rt_csv_sc ~ gender + age_sc + v_entropy_wi * rt_lag_sc * group_bin *last_outcome +  v_max_wi + trial_neg_inv_sc + (1 + rt_lag_sc | id),data=df, REML=TRUE,control=lmer_control)


