library(tidyverse)
library(fmri.pipeline)
library(broom)

repo_directory_mmclock <- "~/clock_analysis"
ncores <- 26

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
dfmmc <- get_trial_data(repo_directory=repo_directory_mmclock,dataset='mmclock_fmri')
demo <- read.table(file=file.path(repo_directory_mmclock, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.double(demo$id)
dfmmc <- inner_join(dfmmc,demo,by=c('id'))
dfmmc <- dfmmc %>% mutate(dataset = 'Experiment 1 - fMRI')

stats <- read_csv('/Users/dnplserv/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats <- stats %>% select(!dataset)
dfmmc <- inner_join(dfmmc,stats,by='id')
dfmmc <- dfmmc %>% mutate(sex = case_when(female==1 ~ 'F',
                                          female==0 ~ 'M'))
dfmmc$sex <- relevel(as.factor(dfmmc$sex),ref='M')
dfmmc$age <- scale(dfmmc$age)

dfmmc <- dfmmc %>% select(dataset,id,trial,age,sex,run_trial,iti_onset,rewFunc,ev,score_csv,feedback_onset,rt_lag_sc,rt_csv,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)

dfmmc_meg <- get_trial_data(repo_directory=repo_directory_mmclock,dataset='mmclock_meg')
dfmmc_meg <- dfmmc_meg %>% mutate(dataset = 'Experiment 1 - MEG Replication')

stats <- read_csv('/Users/dnplserv/clock_analysis/meg/data/mmclock_meg_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats <- stats %>% select(!dataset)
stats <- stats %>% mutate(id2 = str_extract(id, "^[0-9]+")) %>% select(!id) %>% rename(id = id2)

dfmmc_meg <- inner_join(dfmmc_meg,stats,by=c('id'))

demo <- read.table(file=file.path(repo_directory_mmclock, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)
dfmmc_meg <- inner_join(dfmmc_meg,demo,by='id')

dfmmc_meg <- dfmmc_meg %>% mutate(sex = case_when(female==1 ~ 'F',
                                                  female==0 ~ 'M'))
dfmmc_meg$sex <- relevel(as.factor(dfmmc_meg$sex),ref='M')
dfmmc_meg$age <- scale(dfmmc_meg$age)
dfmmc_meg <- dfmmc_meg %>% mutate(iti_onset = 0, feedback_onset = 0)
dfmmc_meg <- dfmmc_meg %>% select(dataset,trial,run_trial,id,age,sex,iti_onset,ev,score_csv,rewFunc,feedback_onset,rt_csv,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)


dfmer <- rbind(dfmmc,dfmmc_meg)
dfmer <- dfmer %>% mutate(ddate = 0)

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
demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
demo <- demo %>% select(!id) %>% rename(id=registration_redcapid,group=registration_group)
demo$gender <- relevel(as.factor(demo$gender),ref='M')
demo$age <- scale(demo$age)

scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = scan_date) %>% select(!registration_redcapid)

demo <- inner_join(demo,scaninfo,by='id')

demo$id <- as.character(demo$id)
df_exp <- inner_join(demo,df,by=c('id'))
df_exp <- df_exp %>% mutate(sex = gender) %>% select(!gender)
df_exp$sex <- relevel(as.factor(df_exp$sex),ref='M')
df_exp$age <- scale(df_exp$age)
df_exp <- df_exp %>% filter(group == 'HC')
df_exp <- df_exp %>% mutate(dataset = 'Experiment 2')
df_exp <- df_exp %>% select(!run_trial)
df_exp <- df_exp %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                         trial > 40 & trial <= 80 ~ trial-40,
                                         trial > 80 & trial <=120 ~ trial-80, 
                                         trial > 120 & trial <=160 ~ trial-120,
                                         trial > 160 & trial <=200 ~ trial-160,
                                         trial > 200 & trial <=240 ~ trial-200))


df_exp <- df_exp %>% select(dataset,trial,id,age,sex,run_trial,iti_onset,score_csv,ev,rewFunc,ddate,rt_csv,feedback_onset,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)

df3 <- rbind(df_exp,dfmer)

df3 <- df3 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
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

df3$trial_bin <- factor(df3$trial_bin, levels = c('1-5','6-10','11-15','16-20','21-25','26-30','30-35','36-40','40-45','46-50','51-55','55-63'))

df0 <- df3 %>% group_by(id,dataset,rewFunc,trial_bin) %>% summarize(rt_csv = mean(rt_csv,na.rm=TRUE)) %>% ungroup()
df0 <- df0 %>% pivot_wider(names_from='rewFunc',values_from='rt_csv') %>% mutate(rt_diff_IEV = IEV - (CEV+CEVR)/2, rt_diff_DEV = DEV - (CEV+CEVR)/2, rt_unlearnable = (CEV+CEVR)/2 - mean((CEV+CEVR)/2,na.rm=TRUE)) %>% pivot_longer(cols = c('rt_diff_IEV','rt_diff_DEV','rt_unlearnable'))
df0 <- df0 %>% select(!DEV & !IEV & !CEV & !CEVR) %>% filter(dataset != 'Experiment 2')
ggplot(df0,aes(x=trial_bin,y=value)) + geom_smooth(aes(group=name,color=name)) + facet_wrap(~dataset) + theme(axis.text.x = element_text(angle = 45, vjust=0.5,hjust=0.5)) + ylab('RT difference')


df0 <- df3 %>% group_by(id,dataset,rewFunc) %>% summarize(score = sum(score_csv,na.rm=TRUE)/max(run_trial)) %>% ungroup()
df0 <- df0 %>% group_by(dataset) %>% pivot_wider(names_from='rewFunc',values_from='score') %>% mutate(score_diff_IEV = IEV - (CEV+CEVR)/2, score_diff_DEV = DEV - (CEV+CEVR)/2) %>% pivot_longer(cols = c('score_diff_IEV','score_diff_DEV')) %>% ungroup()
df0 <- df0 %>% select(!DEV & !IEV & !CEV & !CEVR) %>% filter(dataset != 'Experiment 2')
df0$name <- factor(df0$name,levels=c('score_diff_DEV','score_diff_IEV'))
ggplot(df0,aes(x=value,group=name,color=name)) + geom_histogram() + facet_grid(name~dataset) + theme(axis.text.x = element_text(angle = 45, vjust=0.5,hjust=0.5)) + scale_y_continuous(breaks = c(0,25,50,75,100,125,150))
  
bins <- seq(from = 0.5, to = 65, by = 5)
#bin labels for each interval (1-5, 6-10, ..., 36-40)
bin_labels <- paste(seq(1, 60, by = 5), seq(5, 60, by = 5), sep = "-")
                                                                                                                                                                       
#cut to create the bin variable as a factor
#set right = TRUE to include the upper limit and exclude the lower limit
df3$bin <- cut(df3$run_trial, breaks = bins, include.lowest = TRUE, right = TRUE, labels = bin_labels)
#make sure the factor is unordered
df3$trial_bin <- factor(df3$trial_bin, ordered = FALSE)
### RT ~ RT_lag * trial_bin, plot interaction.

#df3 <- df3 %>% filter(!is.na(bin))

df3$reward_lag_rec<-df3$last_outcome
df3 <- df3 %>% dplyr::mutate(reward_lag_rec = case_when(reward_lag_rec=="reward" ~ 0.5, reward_lag_rec=="omission" ~ -0.5))


decode_formula <- NULL
decode_formula[[1]] <- formula(~rt_lag_sc*trial_bin*last_outcome + (1 | id/run))

splits = c('dataset','rewFunc')
df0 = decode_formula[[1]]
setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
ddf <- mixed_by(df3, outcomes = "rt_csv", rhs_model_formulae = df0 , split_on = splits,
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
save(ddf,file=paste0(curr_date,'-All-datasets-rt-rt_lag_by_trial_bin',1,'.Rdata'))

#term_list <- c("rt_lag_sc:trial_bin6-10","rt_lag_sc:trial_bin11-15","rt_lag_sc:trial_bin16-20","rt_lag_sc:trial_bin21-25","rt_lag_sc:trial_bin26-30","rt_lag_sc:trial_bin30-35","rt_lag_sc:trial_bin36-40","rt_lag_sc:trial_bin40-45","rt_lag_sc:trial_bin46-50","rt_lag_sc:trial_bin51-55","rt_lag_sc:trial_bin55-63")
term_list <- c("rt_lag_sc:trial_bin6-10:last_outcomeReward","rt_lag_sc:trial_bin11-15:last_outcomeReward","rt_lag_sc:trial_bin16-20:last_outcomeReward","rt_lag_sc:trial_bin21-25:last_outcomeReward","rt_lag_sc:trial_bin26-30:last_outcomeReward","rt_lag_sc:trial_bin30-35:last_outcomeReward","rt_lag_sc:trial_bin36-40:last_outcomeReward","rt_lag_sc:trial_bin40-45:last_outcomeReward","rt_lag_sc:trial_bin46-50:last_outcomeReward","rt_lag_sc:trial_bin51-55:last_outcomeReward","rt_lag_sc:trial_bin55-63:last_outcomeReward")


ddq <- ddf$coef_df_reml %>% filter(effect=='fixed') %>% filter(term %in% term_list)
ddq$term_list <- factor(ddq$term, levels=c("rt_lag_sc:trial_bin6-10:last_outcomeReward","rt_lag_sc:trial_bin11-15:last_outcomeReward:last_outcomeReward","rt_lag_sc:trial_bin16-20:last_outcomeReward","rt_lag_sc:trial_bin21-25:last_outcomeReward","rt_lag_sc:trial_bin26-30:last_outcomeReward","rt_lag_sc:trial_bin30-35:last_outcomeReward","rt_lag_sc:trial_bin36-40:last_outcomeReward","rt_lag_sc:trial_bin40-45:last_outcomeReward","rt_lag_sc:trial_bin46-50:last_outcomeReward","rt_lag_sc:trial_bin51-55:last_outcomeReward","rt_lag_sc:trial_bin55-63:last_outcomeReward"))
ddq <- ddq %>% mutate(p_level = case_when(padj_fdr_term >= 0.05 ~ 'NS', padj_fdr_term < 0.05 & padj_fdr_term > 0.001 ~ 'p < 0.05 & p > 0.001', padj_fdr_term <= 0.001 ~ 'p < 0.001'))
ddq$p_level <- factor(ddq$p_level,levels = c('NS','p < 0.05 & p > 0.001','p < 0.001'))
ggplot(ddq,aes(x=term_list,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error,color=rewFunc,group=rewFunc)) + 
  geom_errorbar() + 
  geom_line() + 
  geom_point(aes(alpha=p_level,size=p_level)) + 
  facet_wrap(~dataset) + scale_y_reverse() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  ylab('<-- less -- Exploration -- more -->')



###########################################################################
##### Plot RT Swings, lmer of rt~rt_lag*entropy and rt_~rt_lag*vmax  ######
###########################################################################


repo_directory_mmclock <- "~/clock_analysis"

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
dfmmc <- get_trial_data(repo_directory=repo_directory_mmclock,dataset='mmclock_fmri')
demo <- read.table(file=file.path(repo_directory_mmclock, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.double(demo$id)
dfmmc <- inner_join(dfmmc,demo,by=c('id'))
dfmmc <- dfmmc %>% mutate(dataset = 'Experiment 1 - fMRI')

stats <- read_csv('/Users/dnplserv/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats <- stats %>% select(!dataset)
dfmmc <- inner_join(dfmmc,stats,by='id')
dfmmc <- dfmmc %>% mutate(sex = case_when(female==1 ~ 'F',
                                          female==0 ~ 'M'))
dfmmc$sex <- relevel(as.factor(dfmmc$sex),ref='M')
dfmmc$age <- scale(dfmmc$age)

dfmmc <- dfmmc %>% select(dataset,id,trial,age,sex,run_trial,iti_onset,rt_swing,rewFunc,ev,score_csv,feedback_onset,rt_lag_sc,rt_csv,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)

dfmmc_meg <- get_trial_data(repo_directory=repo_directory_mmclock,dataset='mmclock_meg')
dfmmc_meg <- dfmmc_meg %>% mutate(dataset = 'Experiment 1 - MEG Replication')

stats <- read_csv('/Users/dnplserv/clock_analysis/meg/data/mmclock_meg_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats <- stats %>% select(!dataset)
stats <- stats %>% mutate(id2 = str_extract(id, "^[0-9]+")) %>% select(!id) %>% rename(id = id2)

dfmmc_meg <- inner_join(dfmmc_meg,stats,by=c('id'))

demo <- read.table(file=file.path(repo_directory_mmclock, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)
dfmmc_meg <- inner_join(dfmmc_meg,demo,by='id')

dfmmc_meg <- dfmmc_meg %>% mutate(sex = case_when(female==1 ~ 'F',
                                                  female==0 ~ 'M'))
dfmmc_meg$sex <- relevel(as.factor(dfmmc_meg$sex),ref='M')
dfmmc_meg$age <- scale(dfmmc_meg$age)
dfmmc_meg <- dfmmc_meg %>% mutate(iti_onset = 0, feedback_onset = 0)
dfmmc_meg <- dfmmc_meg %>% select(dataset,trial,run_trial,id,age,rt_swing,sex,iti_onset,ev,score_csv,rewFunc,feedback_onset,rt_csv,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)


dfmer <- rbind(dfmmc,dfmmc_meg)
dfmer <- dfmer %>% mutate(ddate = 0)

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
demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
demo <- demo %>% select(!id) %>% rename(id=registration_redcapid,group=registration_group)
demo$gender <- relevel(as.factor(demo$gender),ref='M')
demo$age <- scale(demo$age)

scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = scan_date) %>% select(!registration_redcapid)

demo <- inner_join(demo,scaninfo,by='id')

demo$id <- as.character(demo$id)
df_exp <- inner_join(demo,df,by=c('id'))
df_exp <- df_exp %>% mutate(sex = gender) %>% select(!gender)
df_exp$sex <- relevel(as.factor(df_exp$sex),ref='M')
df_exp$age <- scale(df_exp$age)
df_exp <- df_exp %>% filter(group == 'HC')
df_exp <- df_exp %>% mutate(dataset = 'Experiment 2')
df_exp <- df_exp %>% select(!run_trial)
df_exp <- df_exp %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                                  trial > 40 & trial <= 80 ~ trial-40,
                                                  trial > 80 & trial <=120 ~ trial-80, 
                                                  trial > 120 & trial <=160 ~ trial-120,
                                                  trial > 160 & trial <=200 ~ trial-160,
                                                  trial > 200 & trial <=240 ~ trial-200))


df_exp <- df_exp %>% select(dataset,trial,id,age,sex,run_trial,rt_swing,iti_onset,score_csv,ev,rewFunc,ddate,rt_csv,feedback_onset,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)

df3 <- rbind(df_exp,dfmer)

df3 <- df3 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
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

df3$trial_bin <- factor(df3$trial_bin, levels = c('1-5','6-10','11-15','16-20','21-25','26-30','30-35','36-40','40-45','46-50','51-55','55-63'))


df3 <- df3 %>% group_by(dataset) %>% mutate(median_split = total_earnings >= median(total_earnings,na.rm=TRUE)) %>% ungroup()

ggplot(data = df3 %>% filter(rewFunc == 'IEV' | rewFunc=='DEV'), aes(x=run_trial,y=rt_swing,color=rewFunc,group=rewFunc)) + geom_smooth(method='loess',span=0.7) + facet_grid(last_outcome~dataset)


# df3 <- df3 %>% group_by(dataset,id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 13 ~ 'Early',
#   run_trial > 13 & run_trial < 26 ~ 'Middle',
#   run_trial >=26 ~ 'Late',
# )))

#df3$trial_bin <- relevel(as.factor(df3$trial_bin),ref='Middle')

df3$last_outcome <- relevel(as.factor(df3$last_outcome),ref='Reward')

decode_formula <- NULL
decode_formula[[1]] <- formula(~rt_lag_sc*v_entropy_wi*last_outcome + rt_lag_sc*v_max_wi*last_outcome + rt_lag_sc*run_trial + (1 | id/run))

splits = c('dataset')
df0 = decode_formula[[1]]
setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
ddf <- mixed_by(df3, outcomes = "rt_csv", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
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
                emtrends_spec = list(
                  HTrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                              specs = formula(~rt_lag_sc:v_entropy_wi:last_outcome),at=list(v_entropy_wi = c(-1.5,0,1.5))),
                  HVrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                              specs = formula(~rt_lag_sc:v_max_wi:last_outcome),at=list(v_max_wi = c(-1.5,0,1.5)))
                )
)
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf,file=paste0(curr_date,'-All-datasets-entropy-vmax-rtlag-emtrends-',1,'.Rdata'))


ddq <- ddf$emtrends_list$HTrt %>% filter(v_entropy_wi != 0)
ddq <- ddq %>% mutate(entropy = case_when(v_entropy_wi==-1.5 ~ 'low entropy',
                                          v_entropy_wi==1.5 ~ 'high entropy'))

ggplot(data=ddq %>% filter(last_outcome=='Reward') ,aes(x=entropy,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  facet_wrap(~dataset,scales='free_y') + geom_errorbar(width=0.5) + scale_y_reverse() + ylab('<-- less -- Exploration -- more -->')


ddq <- ddf$emtrends_list$HVrt %>% filter(v_max_wi != 0)
ddq <- ddq %>% mutate(value_maximum = case_when(v_max_wi==-1.5 ~ 'low value max',
                                          v_max_wi==1.5 ~ 'high value max'))

ggplot(data=ddq %>% filter(last_outcome=='Omission') ,aes(x=value_maximum,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  facet_wrap(~dataset,scales='free_y') + geom_errorbar(width=0.5) + scale_y_reverse() + ylab('<-- less -- Exploration -- more -->')
