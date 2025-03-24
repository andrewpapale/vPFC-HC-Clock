# 2025-03-20 AndyP
# test random slope effects on RT~RT_lag correlations to compare to MPlus MSEM

library(tidyverse)
library(lmerTest)

toalign <- 'clock'
repo_directory <- '~/clock_analysis'


#######################
### MMClock - Sex  ####
#######################

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

mDMNsex_mmclock_fmri_AH <- lmerTest::lmer(data=Q2 %>% filter(network=='DMN' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRsex_mmclock_fmri_AH <- lmerTest::lmer(data=Q2 %>% filter(network=='CTR' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMsex_mmclock_fmri_AH <- lmerTest::lmer(data=Q2 %>% filter(network=='LIM' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mDMNsex_mmclock_fmri_PH <- lmerTest::lmer(data=Q2 %>% filter(network=='DMN' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRsex_mmclock_fmri_PH <- lmerTest::lmer(data=Q2 %>% filter(network=='CTR' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMsex_mmclock_fmri_PH <- lmerTest::lmer(data=Q2 %>% filter(network=='LIM' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))

####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNsex_mmclock_fmri_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-AH fMRI')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mCTRsex_mmclock_fmri_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-AH fMRI')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMsex_mmclock_fmri_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-AH fMRI')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNsex_mmclock_fmri_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-PH fMRI')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mCTRsex_mmclock_fmri_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-PH fMRI')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMsex_mmclock_fmri_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-PH fMRI')


source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_meg')
df <- df %>% select(rt_vmax_lag,ev,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(rt_sec = rt_csv/1000) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         score_sc = scale(score_csv),
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
Q3 <- inner_join(qdf,df,by='id')

Q3$trial_bin <- factor(Q3$trial_bin)
Q3$trial_bin <- relevel(Q3$trial_bin, ref='Middle')
Q3 <- Q3 %>% rename(subj_level_rand_slope=estimate)

Q3 <- Q3 %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)

Q3 <- inner_join(Q3,demo,by=c('id'))
Q3$age <- scale(Q3$age)
Q3 <- Q3 %>% mutate(sex = case_when(female == 1 ~ 'F',
                                    female == 0 ~ 'M'))
Q3$sex <- relevel(as.factor(Q3$sex),ref='F')

mDMNsex_mmclock_meg_AH <- lmerTest::lmer(data=Q3 %>% filter(network=='DMN' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRsex_mmclock_meg_AH <- lmerTest::lmer(data=Q3 %>% filter(network=='CTR' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMsex_mmclock_meg_AH <- lmerTest::lmer(data=Q3 %>% filter(network=='LIM' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mDMNsex_mmclock_meg_PH <- lmerTest::lmer(data=Q3 %>% filter(network=='DMN' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRsex_mmclock_meg_PH <- lmerTest::lmer(data=Q3 %>% filter(network=='CTR' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMsex_mmclock_meg_PH <- lmerTest::lmer(data=Q3 %>% filter(network=='LIM' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))

####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNsex_mmclock_meg_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-AH MEG')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mCTRsex_mmclock_meg_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-AH MEG')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMsex_mmclock_meg_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-AH MEG')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNsex_mmclock_meg_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-PH MEG')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mCTRsex_mmclock_meg_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-PH MEG')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMsex_mmclock_meg_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-PH MEG')

#################
### BSocial  ####
#################


setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
model_str <- paste0('Bsocial-vPFC-HC-network-',toalign,'-ranslopes-nofixedeffect-',3,'.Rdata')
model_str <- Sys.glob(paste0('*',model_str))
load(model_str)

qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')
qdf <- qdf %>% rename(id=level)
qdf <- qdf %>% group_by(id,network,HC_region) %>% summarize(estimate = mean(estimate,na.rm=TRUE)) %>% ungroup()
qdf <- qdf %>% group_by(network,HC_region) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
qdf$id <- as.character(qdf$id)
#qdf <- qdf %>% select(!outcome)

df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
split_ksoc_bsoc <- df %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300])
bsoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240])
ksoc <- rbind(ksoc,221193,221611)
bsoc <- rbind(bsoc,221973,219757,220419,220691,221507,221842,440223)
df_bsoc <- df %>% filter(id %in% bsoc$id) %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                                                       trial > 40 & trial <= 80 ~ 2,
                                                                       trial > 80 & trial <=120 ~ 3, 
                                                                       trial > 120 & trial <=160 ~ 4,
                                                                       trial > 160 & trial <=200 ~ 5,
                                                                       trial > 200 & trial <=240 ~ 6))
df_bsoc <- df_bsoc %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                                     trial > 40 & trial <= 80 ~ trial-40,
                                                     trial > 80 & trial <=120 ~ trial-80, 
                                                     trial > 120 & trial <=160 ~ trial-120,
                                                     trial > 160 & trial <=200 ~ trial-160,
                                                     trial > 200 & trial <=240 ~ trial-200))
df_bsoc <- df_bsoc %>% mutate(protocol = 'bsocial',
                              run_trial0_c = run_trial0-floor(run_trial0/40.5),
                              run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                              run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))

df_ksoc <- df %>% filter(id %in% ksoc$id) %>% mutate(block = case_when(trial <= 50 ~ 1, 
                                                                       trial > 50 & trial <= 100 ~ 2,
                                                                       trial > 100 & trial <=150 ~ 3, 
                                                                       trial > 150 & trial <=200 ~ 4,
                                                                       trial > 200 & trial <=250 ~ 5,
                                                                       trial > 250 & trial <=300 ~ 6))
df_ksoc <- df_ksoc %>% mutate(run_trial0 = case_when(trial <= 50 ~ trial, 
                                                     trial > 50 & trial <= 100 ~ trial-50,
                                                     trial > 100 & trial <=150 ~ trial-100, 
                                                     trial > 150 & trial <=200 ~ trial-150,
                                                     trial > 200 & trial <=250 ~ trial-200,
                                                     trial > 250 & trial <=300 ~ trial-250))
df_ksoc <- df_ksoc %>% mutate(protocol = 'ksocial',
                              run_trial0_c = run_trial0-floor(run_trial0/50.5),
                              run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                              run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))

df <- rbind(df_bsoc,df_ksoc)
# select and scale variables of interest
df <- df %>% 
  group_by(id, scanner_run) %>% 
  mutate(id = as.character(id),
         v_chosen_sc = scale(v_chosen),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
         iti_lag_sc = scale(iti_prev),
         v_max_sc = scale(v_max),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         rt_swing_sc = scale(rt_swing)) %>% ungroup()

# select only vars of interest and merge into MRI data
df <- df %>% select(!run) %>% rename(run = scanner_run) %>% select(!run_trial) %>% rename(run_trial = run_trial0)
Qb <- inner_join(df, qdf, by = c("id"))


# add in age and sex variables
demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(repo_directory,'2025-02-27-Partial-demo-pull-KSOC.csv'))
demo$id <- as.character(demo$id)
demo1$id <- as.character(demo1$registration_redcapid)
demo <- demo %>% rename(sex=registration_birthsex,
                        gender=registration_gender,
                        group=registration_group) %>%
  select(id,group,age,sex,gender)
demo1 <- demo1 %>% rename(sex=registration_birthsex,
                          gender=registration_gender,
                          group=registration_group) %>%
  select(id,group,age,sex,gender)
demo2 <- rbind(demo,demo1)
Qb <- inner_join(Qb,demo2,by=c('id'))
Qb$female <- ifelse(Qb$sex==1,1,0)
Qb <- Qb %>% select(!sex)
Qb$age <- scale(Qb$age)
Qb <- Qb %>% filter(group=='HC')
#Qb$group <- relevel(factor(Qb$group),ref='HC')

Qb <- Qb %>% rename(subj_level_rand_slope=estimate)

Qb <- Qb %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

Qb$age <- scale(Qb$age)
Qb <- Qb %>% mutate(sex = case_when(female == 1 ~ 'F',
                                    female == 0 ~ 'M'))
Qb$sex <- relevel(as.factor(Qb$sex),ref='F')

Qb <- Qb %>% filter(run_trial < 41 & rt_csv > 0.2 & rt_csv < 4)

mDMNsex_bsocial_fmri_AH <- lmerTest::lmer(data=Qb %>% filter(network=='DMN' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRsex_bsocial_fmri_AH <- lmerTest::lmer(data=Qb %>% filter(network=='CTR' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMsex_bsocial_fmri_AH <- lmerTest::lmer(data=Qb %>% filter(network=='LIM' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mDMNsex_bsocial_fmri_PH <- lmerTest::lmer(data=Qb %>% filter(network=='DMN' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRsex_bsocial_fmri_PH <- lmerTest::lmer(data=Qb %>% filter(network=='CTR' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMsex_bsocial_fmri_PH <- lmerTest::lmer(data=Qb %>% filter(network=='LIM' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*sex*subj_level_rand_slope + rt_vmax_lag_sc*sex*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNsex_bsocial_fmri_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-AH Bsocial')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mCTRsex_bsocial_fmri_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-AH Bsocial')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMsex_bsocial_fmri_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-AH Bsocial')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNsex_bsocial_fmri_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-PH Bsocial')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMsex_bsocial_fmri_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-PH Bsocial')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMsex_bsocial_fmri_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "sex") ,at = list(subj_level_rand_slope = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = sex)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-PH Bsocial')



# 

#######################
### MMClock - Age  ####
#######################

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

mDMNage_mmclock_fmri_AH <- lmerTest::lmer(data=Q2 %>% filter(network=='DMN' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRage_mmclock_fmri_AH <- lmerTest::lmer(data=Q2 %>% filter(network=='CTR' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMage_mmclock_fmri_AH <- lmerTest::lmer(data=Q2 %>% filter(network=='LIM' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mDMNage_mmclock_fmri_PH <- lmerTest::lmer(data=Q2 %>% filter(network=='DMN' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRage_mmclock_fmri_PH <- lmerTest::lmer(data=Q2 %>% filter(network=='CTR' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMage_mmclock_fmri_PH <- lmerTest::lmer(data=Q2 %>% filter(network=='LIM' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))

####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNage_mmclock_fmri_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-AH fMRI')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mCTRage_mmclock_fmri_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-AH fMRI')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMage_mmclock_fmri_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-AH fMRI')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNage_mmclock_fmri_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-PH fMRI')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mCTRage_mmclock_fmri_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-PH fMRI')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMage_mmclock_fmri_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-PH fMRI')


source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_meg')
df <- df %>% select(rt_vmax_lag,ev,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(rt_sec = rt_csv/1000) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         score_sc = scale(score_csv),
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
Q3 <- inner_join(qdf,df,by='id')

Q3$trial_bin <- factor(Q3$trial_bin)
Q3$trial_bin <- relevel(Q3$trial_bin, ref='Middle')
Q3 <- Q3 %>% rename(subj_level_rand_slope=estimate)

Q3 <- Q3 %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)

Q3 <- inner_join(Q3,demo,by=c('id'))
Q3$age <- scale(Q3$age)
Q3 <- Q3 %>% mutate(sex = case_when(female == 1 ~ 'F',
                                    female == 0 ~ 'M'))
Q3$sex <- relevel(as.factor(Q3$sex),ref='F')

mDMNage_mmclock_meg_AH <- lmerTest::lmer(data=Q3 %>% filter(network=='DMN' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRage_mmclock_meg_AH <- lmerTest::lmer(data=Q3 %>% filter(network=='CTR' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMage_mmclock_meg_AH <- lmerTest::lmer(data=Q3 %>% filter(network=='LIM' & HC_region == 'AH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mDMNage_mmclock_meg_PH <- lmerTest::lmer(data=Q3 %>% filter(network=='DMN' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mCTRage_mmclock_meg_PH <- lmerTest::lmer(data=Q3 %>% filter(network=='CTR' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))
mLIMage_mmclock_meg_PH <- lmerTest::lmer(data=Q3 %>% filter(network=='LIM' & HC_region == 'PH'), rt_csv ~ rt_lag_sc*age*subj_level_rand_slope + rt_vmax_lag_sc*age*subj_level_rand_slope + rt_lag_sc*trial_neg_inv_sc + rt_lag_sc*reward_lag_rec + (1|id/run))

####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNage_mmclock_meg_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-AH MEG')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mCTRage_mmclock_meg_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-AH MEG')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMage_mmclock_meg_AH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-AH MEG')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mDMNage_mmclock_meg_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('DMN-PH MEG')



####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mCTRage_mmclock_meg_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('CTR-PH MEG')


####PLOTTING RT SWINGS
em3a <- as_tibble(emmeans::emtrends(mLIMage_mmclock_meg_PH, var = "rt_lag_sc", specs = c("subj_level_rand_slope", "age") ,at = list(subj_level_rand_slope = c(-1.5,1.5), age = c(-1.5,1.5))))
em3a$subj_level_rand_slope <- factor(em3a$subj_level_rand_slope, levels = c(-1.5, 1.5), labels = c("low", "high"))
em3a$age <- factor(em3a$age, levels = c(-1.5, 1.5), labels = c("younger", "older"))
##basic palette plot
ggplot(em3a, aes(x = subj_level_rand_slope, y = rt_lag_sc.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = age)) +
  geom_point(position = position_dodge(width = .6), size = 2.5) + 
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.4, size = 0.9) + 
  theme_bw(base_size = 12) +
  ylab("RT swings (AU)\n Small <---------> Large") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8.5, color = "grey10")) +  
  scale_x_discrete() +
  scale_y_reverse() +
  xlab('HC-vPFC Modulation') +
  ggtitle('LIM-PH MEG')
