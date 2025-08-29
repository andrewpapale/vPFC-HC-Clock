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

dfmmc <- dfmmc %>% select(dataset,id,trial,age,sex,run_trial,rt_vmax,iti_onset,rewFunc,ev,score_csv,feedback_onset,rt_lag_sc,rt_csv,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)

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
dfmmc_meg <- dfmmc_meg %>% select(dataset,trial,run_trial,id,age,sex,rt_vmax,iti_onset,ev,score_csv,rewFunc,feedback_onset,rt_csv,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)


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


df_exp <- df_exp %>% select(dataset,trial,id,age,sex,run_trial,iti_onset,rt_vmax,score_csv,ev,rewFunc,ddate,rt_csv,feedback_onset,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)

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


df1 <- df3 %>% group_by(id,dataset,rewFunc) %>% summarize(score = sum(score_csv,na.rm=TRUE)/max(run_trial)) %>% ungroup()

pal <- NULL
pal[2] <- '#009B72'
pal[1] <- '#D65B26'
pal[3] <- '#6F6AB0'
pal[4] <- '#E50886'

setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
pdf('Figure_S1A.pdf',height = 8, width = 14)
gg1 <- ggplot(df1, aes(x=score,fill=rewFunc)) + geom_histogram(bins=30,color="black") + 
  facet_grid(rewFunc~dataset) + 
  theme_minimal(base_size = 18) + theme(axis.text.x = element_text(size = 18,angle = 45, hjust = 0.75),
                          axis.text.y = element_text(size = 18),
                          axis.title.x = element_text(size = 18),
                          axis.title.y = element_text(size = 18),
                          strip.text.x = element_text(size = 18),
                          strip.text.y = element_text(size = 18),
                          legend.text = element_text(size = 18),
                          legend.spacing.y = unit(20.0, 'cm'),
                          legend.title = element_blank()) + 
  scale_fill_manual(values = pal) +
  xlab('Average Reward Per Trial')
print(gg1)
dev.off()

df0 <- df1 %>% group_by(dataset) %>% pivot_wider(names_from='rewFunc',values_from='score') %>% mutate(IEV_minus_unlearnable = IEV - (CEV+CEVR)/2, DEV_minus_unlearnable = DEV - (CEV+CEVR)/2) %>% pivot_longer(cols = c('IEV_minus_unlearnable','DEV_minus_unlearnable')) %>% ungroup()
df0 <- df0 %>% select(!DEV & !IEV & !CEV & !CEVR) %>% filter(dataset != 'Experiment 2')
df0 <- df0 %>% mutate(name1 = case_when(name == 'DEV_minus_unlearnable' ~ 'DEV-Unlearnable',
                                        name == 'IEV_minus_unlearnable' ~ 'IEV-Unlearnable'))

df0$name1 <- factor(df0$name1,levels=c('DEV-unlearnable','IEV-unlearnable'))
ggplot(df0,aes(x=value,group=name,color=name)) + geom_histogram() + facet_grid(name~dataset) + theme(axis.text.x = element_text(angle = 45, vjust=0.5,hjust=0.5)) + scale_y_continuous(breaks = c(0,25,50,75,100,125,150))

setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
pdf('Figure_S1B.pdf',height = 8, width = 14)
gg1 <- ggplot(df0, aes(x=value,group=name1,fill=name1)) + geom_histogram(bins=30,color="black") + 
  facet_grid(name1~dataset) + 
  theme_minimal(base_size = 18) + theme(axis.text.x = element_text(size = 18,angle = 45, hjust = 0.75),
                                        axis.text.y = element_text(size = 18),
                                        axis.title.x = element_text(size = 18),
                                        axis.title.y = element_text(size = 18),
                                        strip.text.x = element_text(size = 18),
                                        strip.text.y = element_text(size = 18),
                                        legend.text = element_text(size = 18),
                                        legend.spacing.y = unit(20.0, 'cm'),
                                        legend.title = element_blank()) + 
  scale_fill_manual(values = pal) +
  xlab('Score Difference')
print(gg1)
dev.off()

  
setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
pdf('Figure_S1C.pdf',height = 8, width = 14)
gg1 <- ggplot(df1 %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV'), aes(x=score,fill=rewFunc)) + geom_histogram(bins=30,color="black") + 
  facet_grid(rewFunc~dataset) + 
  theme_minimal(base_size = 18) + theme(axis.text.x = element_text(size = 18,angle = 45, hjust = 0.75),
                                        axis.text.y = element_text(size = 18),
                                        axis.title.x = element_text(size = 18),
                                        axis.title.y = element_text(size = 18),
                                        strip.text.x = element_text(size = 18),
                                        strip.text.y = element_text(size = 18),
                                        legend.text = element_text(size = 18),
                                        legend.spacing.y = unit(20.0, 'cm'),
                                        legend.title = element_blank()) + 
  scale_fill_manual(values = pal) +
  xlab('Average Reward Per Trial')
print(gg1)
dev.off()

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

dfmmc <- dfmmc %>% select(dataset,id,trial,age,sex,score_csv,rt_vmax,probability,run_trial,v_entropy_lag,iti_onset,rt_swing,rewFunc,ev,score_csv,feedback_onset,rt_lag_sc,rt_csv,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)

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
dfmmc_meg <- dfmmc_meg %>% select(dataset,trial,run_trial,rt_vmax, probability,score_csv,v_entropy_lag,id,age,rt_swing,sex,iti_onset,ev,score_csv,rewFunc,feedback_onset,rt_csv,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)


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


df_exp <- df_exp %>% select(dataset,trial,id,age,sex,rt_vmax,probability,run_trial,v_entropy_lag,score_csv,rt_swing,iti_onset,score_csv,ev,rewFunc,ddate,rt_csv,feedback_onset,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings,v_max,v_entropy,ev,magnitude,probability)

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

df3 <- df3 %>% group_by(dataset,id,run) %>% mutate(v_entropy_wi_lag = lag(v_entropy_wi), v_max_wi_lag = lag(v_max_wi)) %>% ungroup()

ggplot(data = df3 %>% filter(rewFunc == 'IEV' | rewFunc=='DEV'), aes(x=run_trial,y=rt_swing,color=rewFunc,group=rewFunc)) + geom_smooth(method='loess',span=0.7) + facet_grid(last_outcome~dataset)


# df3 <- df3 %>% group_by(dataset,id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 13 ~ 'Early',
#   run_trial > 13 & run_trial < 26 ~ 'Middle',
#   run_trial >=26 ~ 'Late',
# )))

#df3$trial_bin <- relevel(as.factor(df3$trial_bin),ref='Middle')

df3$last_outcome <- relevel(as.factor(df3$last_outcome),ref='Reward')

df3 <- df3 %>% group_by(dataset,id) %>% mutate(v_entropy_lag_sc = as.vector(scale(v_entropy_lag))) %>% ungroup()


decode_formula <- NULL
decode_formula[[1]] <- formula(~ v_entropy_lag_sc + v_max_wi_lag + run_trial + last_outcome  + (1 | id/run))

splits = c('dataset')
df0 = decode_formula[[1]]
setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
ddf <- mixed_by(df3, outcomes = "rt_swing", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                 emmeans_spec = list(
                #   A = list(outcome='vmPFC_decon',model_name='model1',
                #            specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
                   W = list(outcome='rt_swing',model_name='model1',
                            specs=formula(~v_max_wi_lag),at=list(v_max_wi_lag=c(-2,0,2))),
                   H = list(outcome='rt_swing', model_name='model1',
                            specs=formula(~v_entropy_lag_sc), at = list(v_entropy_lag_sc=c(-2,0,2)))
                #   Y = list(outcome='vmPFC_decon',model_name='model1',
                #            specs=formula(~education_yrs), at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
                )
                # emtrends_spec = list(
                #   HTrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:v_entropy_wi:last_outcome),at=list(v_entropy_wi = c(-1.5,0,1.5))),
                #   HVrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:v_max_wi:last_outcome),at=list(v_max_wi = c(-1.5,0,1.5)))
                # )
)
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf,file=paste0(curr_date,'-All-datasets-entropy-vmax-rtlag-emtrends-',1,'.Rdata'))


ddq <- ddf$emmeans_list$H %>% filter(v_entropy_lag_sc != 0)
ddq <- ddq %>% mutate(entropy = case_when(v_entropy_lag_sc==-2 ~ 'low entropy',
                                          v_entropy_lag_sc==2 ~ 'high entropy'))

ggplot(data=ddq ,aes(x=entropy,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->')


ddq <- ddf$emmeans_list$W %>% filter(v_max_wi_lag != 0)
ddq <- ddq %>% mutate(value_maximum = case_when(v_max_wi_lag==-2 ~ 'low value max',
                                                v_max_wi_lag==2 ~ 'high value max'))

ggplot(data=ddq ,aes(x=value_maximum,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->')

ggplot(df3 %>% filter(rewFunc == 'IEV' | rewFunc=='DEV'),aes(x=run_trial,y=rt_swing,color=rewFunc,group=rewFunc)) + 
  geom_smooth(method='loess',span = 0.5) + facet_wrap(~dataset)


################################
### low vs high outcomes #######
################################

df4 <- df3 %>% filter(score_csv != 0)
df4 <- df4 %>% group_by(dataset,id,run) %>% mutate(score_csv_wi_lag = lag(scale(score_csv))) %>% ungroup()

decode_formula <- NULL
decode_formula[[1]] <- formula(~ score_csv_wi_lag + v_entropy_lag_sc + v_max_wi_lag + run_trial + last_outcome + (1 | id/run))

splits = c('dataset')
df0 = decode_formula[[1]]
setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
ddf <- mixed_by(df4, outcomes = "rt_swing", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                emmeans_spec = list(
                  #   A = list(outcome='vmPFC_decon',model_name='model1',
                  #            specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
                  W = list(outcome='rt_swing',model_name='model1',
                           specs=formula(~v_max_wi_lag),at=list(v_max_wi_lag=c(-2,0,2))),
                  H = list(outcome='rt_swing', model_name='model1',
                           specs=formula(~v_entropy_lag_sc), at = list(v_entropy_lag_sc=c(-2,0,2))),
                  Y = list(outcome='rt_swing',model_name='model1',
                              specs=formula(~score_csv_wi_lag), at=list(score_csv_wi_lag=c(-2,0,2)))
                )
                # emtrends_spec = list(
                #   HTrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:v_entropy_wi),at=list(v_entropy_wi = c(-1.5,0,1.5))),
                #   HVrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:v_max_wi),at=list(v_max_wi = c(-1.5,0,1.5))),
                #   HSrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:score_csv_sc),at=list(score_csv_sc = c(-1.5,0,1.5)))
                # )
)
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf,file=paste0(curr_date,'-All-score_csv-vmax-rtlag-emtrends-',1,'.Rdata'))

ddq <- ddf$emmeans_list$H %>% filter(v_entropy_lag_sc != 0)
ddq <- ddq %>% mutate(entropy = case_when(v_entropy_lag_sc==-2 ~ 'low entropy',
                                          v_entropy_lag_sc==2 ~ 'high entropy'))

ggplot(data=ddq ,aes(x=entropy,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->')


ddq <- ddf$emmeans_list$W %>% filter(v_max_wi_lag != 0)
ddq <- ddq %>% mutate(value_maximum = case_when(v_max_wi_lag==-2 ~ 'low value max',
                                                v_max_wi_lag==2 ~ 'high value max'))

ggplot(data=ddq ,aes(x=value_maximum,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->')


ddq <- ddf$emmeans_list$Y %>% filter(score_csv_wi_lag != 0)
ddq <- ddq %>% mutate(score = case_when(score_csv_wi_lag==-2 ~ 'low score',
                                        score_csv_wi_lag==2 ~ 'high score'))

ggplot(data=ddq ,aes(x=score,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->')

ggplot(df4 %>% filter(rewFunc == 'IEV' | rewFunc=='DEV'),aes(x=run_trial,y=rt_swing,color=rewFunc,group=rewFunc)) + 
  geom_smooth(method='loess',span = 0.5) + facet_wrap(~dataset)



#df5 <- df4 %>% filter(rewFunc == 'IEV' | rewFunc=='DEV')

df3 <- df3 %>% group_by(dataset,id,run) %>% mutate(prob_lag_sc = scale(lag(probability))) %>% ungroup()

decode_formula <- NULL
decode_formula[[1]] <- formula(~ prob_lag_sc + v_entropy_lag_sc + v_max_wi_lag + run_trial + last_outcome  + (1 | id/run))

splits = c('dataset')
df0 = decode_formula[[1]]
setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
ddf <- mixed_by(df3, outcomes = "rt_swing", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                emmeans_spec = list(
                  #   A = list(outcome='vmPFC_decon',model_name='model1',
                  #            specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
                  W = list(outcome='rt_swing',model_name='model1',
                           specs=formula(~v_max_wi_lag),at=list(v_max_wi_lag=c(-2,0,2))),
                  H = list(outcome='rt_swing', model_name='model1',
                           specs=formula(~v_entropy_lag_sc), at = list(v_entropy_lag_sc=c(-2,0,2))),
                  Y = list(outcome='rt_swing',model_name='model1',
                           specs=formula(~prob_lag_sc), at=list(prob_lag_sc=c(-2,0,2)))
                )
                # emtrends_spec = list(
                #   HTrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:v_entropy_wi),at=list(v_entropy_wi = c(-1.5,0,1.5))),
                #   HVrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:v_max_wi),at=list(v_max_wi = c(-1.5,0,1.5))),
                #   HSrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:score_csv_sc),at=list(score_csv_sc = c(-1.5,0,1.5)))
                # )
)
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf,file=paste0(curr_date,'-All-probability-vmax-rtlag-emtrends-',1,'.Rdata'))

ddq <- ddf$emmeans_list$H %>% filter(v_entropy_lag_sc != 0)
ddq <- ddq %>% mutate(entropy = case_when(v_entropy_lag_sc==-2 ~ 'low entropy',
                                          v_entropy_lag_sc==2 ~ 'high entropy'))

ggplot(data=ddq ,aes(x=entropy,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->')


ddq <- ddf$emmeans_list$W %>% filter(v_max_wi_lag != 0)
ddq <- ddq %>% mutate(value_maximum = case_when(v_max_wi_lag==-2 ~ 'low value max',
                                                v_max_wi_lag==2 ~ 'high value max'))

ggplot(data=ddq ,aes(x=value_maximum,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->')


ddq <- ddf$emmeans_list$Y %>% filter(prob_lag_sc != 0)
ddq <- ddq %>% mutate(probability = case_when(prob_lag_sc==-2 ~ 'low probability',
                                        prob_lag_sc==2 ~ 'high probability'))

ggplot(data=ddq ,aes(x=probability,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->')


decode_formula <- NULL
decode_formula[[1]] <- formula(~ (1+rt_lag_sc | id) + (1|run))

splits = c('dataset')
df0 = decode_formula[[1]]
setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
ddf <- mixed_by(df3, outcomes = "rt_csv", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)#,
                # emmeans_spec = list(
                #   #   A = list(outcome='vmPFC_decon',model_name='model1',
                #   #            specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
                #   W = list(outcome='rt_swing',model_name='model1',
                #            specs=formula(~v_max_wi_lag),at=list(v_max_wi_lag=c(-2,0,2))),
                #   H = list(outcome='rt_swing', model_name='model1',
                #            specs=formula(~v_entropy_lag_sc), at = list(v_entropy_lag_sc=c(-2,0,2))),
                #   Y = list(outcome='rt_swing',model_name='model1',
                #            specs=formula(~prob_lag_sc), at=list(prob_lag_sc=c(-2,0,2)))
                # )
                # emtrends_spec = list(
                #   HTrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:v_entropy_wi),at=list(v_entropy_wi = c(-1.5,0,1.5))),
                #   HVrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:v_max_wi),at=list(v_max_wi = c(-1.5,0,1.5))),
                #   HSrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                #               specs = formula(~rt_lag_sc:score_csv_sc),at=list(score_csv_sc = c(-1.5,0,1.5)))
                # )
)
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf,file=paste0(curr_date,'-All-rt-lag-sc-random-slope-',1,'.Rdata'))

load('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/2025-06-27-All-rt-lag-sc-random-slope-1.Rdata')

ddq <- ddf$coef_df_reml %>% filter(effect == 'ran_vals' & term =='rt_lag_sc') %>% rename(id=level)

# stats mmc-fmri
stats1 <- read_csv('/Users/dnplserv/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats1 <- stats1 %>% select(!dataset)

# stats mmc-meg
stats2 <- read_csv('/Users/dnplserv/clock_analysis/meg/data/mmclock_meg_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats2 <- stats2 %>% select(!dataset)
stats2 <- stats2 %>% mutate(id2 = str_extract(id, "^[0-9]+")) %>% select(!id) %>% rename(id = id2)

# stats Explore
stats3 <- read_csv('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/explore_decay_factorize_selective_psequate_fixedparams_fmri_mfx_sceptic_global_statistics.csv')

stats1 <- stats1 %>% mutate(dataset = "Experiment 1 - fMRI")
stats2 <- stats2 %>% mutate(dataset = "Experiment 1 - MEG Replication")
stats3 <- stats3 %>% mutate(dataset = 'Experiment 2')

statsall <- rbind(stats1,stats2)
statsall <- rbind(statsall,stats3)

ddq <- inner_join(ddq,statsall,by=c('dataset','id'))

ddq <- ddq %>% group_by(dataset) %>% mutate(rt_lag_sc = scale(estimate),tau_ffx_sc = scale(1/beta_transformed_ffx)) %>% ungroup()

ggplot(data=ddq,aes(x=tau_ffx_sc,y=rt_lag_sc)) + geom_point() + facet_wrap(~dataset) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  ggpubr::stat_regline_equation(size = 3) + 
  ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x.npc = "centre",size=3)



################################
### rt ~ rt_lag correlations ###
################################

library(tseries)

df3 <- df3 %>% group_by(dataset,id) %>% mutate(entropy_sc = scale(v_entropy)) %>% ungroup()

decode_formula <- NULL
decode_formula[[1]] <- formula(~ rt_lag_sc*sex + entropy_sc*rt_lag_sc + v_max_wi*rt_lag_sc + rt_lag_sc*run_trial + rt_lag_sc*last_outcome + (1 | id/run))

splits = c('dataset')
df0 = decode_formula[[1]]
setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
ddf <- mixed_by(df3, outcomes = "rt_csv", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                #emmeans_spec = list(
                  #   A = list(outcome='vmPFC_decon',model_name='model1',
                  #            specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
                  # W = list(outcome='rt_csv',model_name='model1',
                  #          specs=formula(~v_max_wi),at=list(v_max_wi=c(-2,0,2))),
                  # H = list(outcome='rt_csv', model_name='model1',
                  #          specs=formula(~v_entropy_wi), at = list(v_entropy_wi=c(-2,0,2)))
                  #   Y = list(outcome='vmPFC_decon',model_name='model1',
                  #            specs=formula(~education_yrs), at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
               # )
                 emtrends_spec = list(
                   HTrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                               specs = formula(~rt_lag_sc:entropy_sc),at=list(entropy_sc = c(-1.5,0,1.5))),
                   HVrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                               specs = formula(~rt_lag_sc:v_max_wi),at=list(v_max_wi = c(-1.5,0,1.5)))
                 )
)
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf,file=paste0(curr_date,'-All-datasets-entropy-vmax-rt-rtlag-emtrends-',1,'.Rdata'))


ddq <- ddf$emtrends_list$HTrt %>% filter(entropy_sc != 0)
ddq <- ddq %>% mutate(entropy = case_when(entropy_sc==-1.5 ~ 'low entropy',
                                          entropy_sc==1.5 ~ 'high entropy'))

ggplot(data=ddq ,aes(x=entropy,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  facet_grid(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->') + scale_y_reverse()


ddq <- ddf$emtrends_list$HVrt %>% filter(v_max_wi != 0)
ddq <- ddq %>% mutate(value_maximum = case_when(v_max_wi==-1.5 ~ 'low value max',
                                                v_max_wi==1.5 ~ 'high value max'))

ggplot(data=ddq ,aes(x=value_maximum,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->') + scale_y_reverse()

ggplot(df3 %>% filter(rewFunc == 'IEV' | rewFunc=='DEV'),aes(x=run_trial,y=rt_swing,color=rewFunc,group=rewFunc)) + 
  geom_smooth(method='loess',span = 0.5) + facet_wrap(~dataset)


df4 <- df3 %>% filter(score_csv != 0)
df4 <- df4 %>% group_by(dataset,id,run) %>% mutate(score_csv_sc = scale(score_csv)) %>% ungroup()

decode_formula <- NULL
decode_formula[[1]] <- formula(~ sex*rt_lag_sc + score_csv_sc*rt_lag_sc + entropy_sc*rt_lag_sc + v_max_wi*rt_lag_sc + run_trial*rt_lag_sc + (1 | id/run))

splits = c('dataset')
df0 = decode_formula[[1]]
setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')
ddf <- mixed_by(df4, outcomes = "rt_csv", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                # emmeans_spec = list(
                #   #   A = list(outcome='vmPFC_decon',model_name='model1',
                #   #            specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
                #   W = list(outcome='rt_swing',model_name='model1',
                #            specs=formula(~v_max_wi_lag),at=list(v_max_wi_lag=c(-2,0,2))),
                #   H = list(outcome='rt_swing', model_name='model1',
                #            specs=formula(~v_entropy_lag_sc), at = list(v_entropy_lag_sc=c(-2,0,2))),
                #   Y = list(outcome='rt_swing',model_name='model1',
                #            specs=formula(~score_csv_wi_lag), at=list(score_csv_wi_lag=c(-2,0,2)))
                # )
                emtrends_spec = list(
                  HTrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                              specs = formula(~rt_lag_sc:entropy_sc),at=list(entropy_sc = c(-1.5,0,1.5))),
                  HVrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                              specs = formula(~rt_lag_sc:v_max_wi),at=list(v_max_wi = c(-1.5,0,1.5))),
                  HSrt = list(outcome='rt_csv',model_name='model1', var = 'rt_lag_sc',
                              specs = formula(~rt_lag_sc:score_csv_sc),at=list(score_csv_sc = c(-1.5,0,1.5)))
                )
)
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf,file=paste0(curr_date,'-All-score_csv-vmax-rtlag-emtrends-',1,'.Rdata'))

ddq <- ddf$emtrends_list$HTrt %>% filter(entropy_sc != 0)
ddq <- ddq %>% mutate(entropy = case_when(entropy_sc==-1.5 ~ 'low entropy',
                                          entropy_sc==1.5 ~ 'high entropy'))

ggplot(data=ddq ,aes(x=entropy,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->')


ddq <- ddf$emtrends_list$HVrt %>% filter(v_max_wi != 0)
ddq <- ddq %>% mutate(value_maximum = case_when(v_max_wi==-1.5 ~ 'low value max',
                                                v_max_wi==1.5 ~ 'high value max'))

ggplot(data=ddq ,aes(x=value_maximum,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->') + scale_y_reverse()


ddq <- ddf$emtrends_list$HSrt %>% filter(score_csv_sc != 0)
ddq <- ddq %>% mutate(score = case_when(score_csv_sc==-1.5 ~ 'low score',
                                        score_csv_sc==1.5 ~ 'high score'))

ggplot(data=ddq ,aes(x=score,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  facet_wrap(~dataset) + geom_errorbar(width=0.5) + ylab('<-- less -- Exploration -- more -->') + scale_y_reverse()



dfrt <- df3 %>% select(dataset,id,run,rt_csv)

# df_acf <- dfrt %>% 
#   group_by(dataset,id,run) %>%
#   nest() %>% 
#   mutate(data = map(data, ~acf(., lag.max=10, type="correlation", plot=F))) %>%
#   mutate(data = map(data, ~as.data.frame(rbind(.x$acf[1,,], .x$acf[2,,], .x$acf[3,,], .x$acf[4,,], .x$acf[5,,], .x$acf[6,,], .x$acf[7,,], .x$acf[8,,], .x$acf[9,,], .x$acf[10,,])))) %>%
#   unnest(data)
# 
# df_acf_lag <- dfrt %>% 
#   group_by(dataset,id,run) %>%
#   nest() %>% 
#   mutate(data = map(data, ~acf(., lag.max=10, type="correlation", plot=F))) %>%
#   mutate(data = map(data, ~as.data.frame(rbind(.x$lag[1,,], .x$lag[1,,], .x$lag[2,,], .x$lag[3,,], .x$lag[4,,], .x$lag[5,,], .x$lag[6,,], .x$lag[7,,], .x$lag[8,,], .x$lag[9,,], .x$lag[10,,])))) %>%
#   unnest(data)
# 
# df_acf <- df_acf %>% arrange(dataset,id,run)
# df_acf_lag <- df_acf_lag %>% arrange(dataset,id,run)

################################################
### Plot changing v_max_wi and v_entropy_wi ####
################################################

ggplot(data=df3 %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV'), aes(x=run_trial,y=v_max,color=rewFunc,group=rewFunc)) + geom_smooth(method='loess',span=0.15) + facet_grid(~dataset)
ggplot(data=df3 %>% filter(rewFunc == 'IEV' | rewFunc == 'DEV'), aes(x=run_trial,y=v_entropy,color=rewFunc,group=rewFunc)) + geom_smooth(method='loess',span=0.15) + facet_grid(~dataset)


################################################
#### plot model stats, regen data frame. #######
################################################

stats_mmcfmri <- read_csv('/Users/dnplserv/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats_mmcfmri <- stats_mmcfmri %>% select(!dataset)
stats_mmcfmri <- stats_mmcfmri %>% mutate(dataset = 'Experiment 1 - fMRI')

stats_mmcmeg <- read_csv('/Users/dnplserv/clock_analysis/meg/data/mmclock_meg_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats_mmcmeg <- stats_mmcmeg %>% select(!dataset)
stats_mmcmeg <- stats_mmcmeg %>% mutate(id2 = str_extract(id, "^[0-9]+")) %>% select(!id) %>% rename(id = id2)
stats_mmcmeg <- stats_mmcmeg %>% mutate(dataset = 'Experiment 1 - MEG Replication')

stats_exp <- read_csv('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/explore_decay_factorize_selective_psequate_fixedparams_fmri_mfx_sceptic_global_statistics.csv')
stats_exp <- stats_exp %>% select(!dataset)
stats_exp <- stats_exp %>% mutate(dataset = 'Experiment 2')
demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
demo <- demo %>% select(!id) %>% select(registration_redcapid,registration_group) %>% rename(id=registration_redcapid,group=registration_group)
stats_exp <- inner_join(stats_exp,demo,by=c('id'))
stats_exp <- stats_exp %>% filter(group=='HC') %>% select(!group)


jstats <- rbind(stats_mmcfmri,stats_mmcmeg)
jstats <- rbind(jstats,stats_exp)

ggplot(data=jstats, aes(x=alpha_transformed_ffx)) + geom_histogram(bins=30) + facet_wrap(~dataset) + xlab('alpha - learning rate')
ggplot(data=jstats, aes(x=beta_transformed_ffx)) + geom_histogram(bins=30) + facet_wrap(~dataset) + xlab('beta - inverse temperature')
ggplot(data=jstats, aes(x=gamma_transformed_ffx)) + geom_histogram(bins=30) + facet_wrap(~dataset) + xlab('gamma - selective maintenance')
ggplot(data=jstats, aes(x=R2)) + geom_histogram(bins=30) + facet_wrap(~dataset) + xlab('model R^2 fits')

ggplot(data=jstats, aes(x=1/beta_transformed_ffx)) + geom_histogram(bins=30) + facet_wrap(~dataset) + xlab('1/beta - temperature')

#############################################
### entropy vs Npeaks and avg prominence ####
#############################################

load('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/mmclock_fmri_peaks.Rdata')
load('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/explore_peaks.Rdata')
load('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/2025-07-vPFC-HC-Behavior.Rdata')


explore_peaks <- explore_peaks %>% filter(!is.infinite(peaks.prominence)) %>% 
  group_by(id,trial) %>% 
  summarize(nP = length(peaks.peak), maxProm = max(peaks.prominence,na.rm=TRUE)) %>% 
  ungroup()

mmclock_fmri_peaks <- mmclock_fmri_peaks %>% filter(!is.infinite(peaks.prominence)) %>% 
  group_by(id,trial) %>% 
  summarize(nP = length(peaks.peak), maxProm = max(peaks.prominence,na.rm=TRUE)) %>% 
  ungroup()

mmclock_fmri_peaks <- mmclock_fmri_peaks %>% mutate(dataset = 'Experiment 1 - fMRI')
explore_peaks <- explore_peaks %>% mutate(dataset = 'Experiment 2')
mmclock_fmri_peaks$id <- as.character(mmclock_fmri_peaks$id)
explore_peaks$id <- as.character(explore_peaks$id)

peaks_fmri <- rbind(mmclock_fmri_peaks,explore_peaks)

df3 <- df3 %>% filter(dataset != 'Experiment 1 - MEG Replication')

df4 <- inner_join(df3,peaks_fmri,by=c('dataset','id','trial'))

df4 <- df4 %>% group_by(dataset,id,run) %>% mutate(zmaxProm = scale(maxProm), qzmaxProm = ntile(zmaxProm,5)) %>% ungroup()

setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R')

pdf('entropy_vs_Npeaks.pdf',height=6,width=8)
gg1 <- ggplot(df4, aes(x=nP,y=v_entropy,group=nP)) + 
  geom_boxplot(notch=TRUE) + facet_wrap(~dataset) + 
  xlab('Number of Peaks') +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 18) +
  # Remove the legend title
  guides(color = guide_legend(title = NULL))
print(gg1)
dev.off()


#########################
### Explore demo ########
#########################

exp <- read_csv('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/explore_clock_with_income.csv') %>% rename(id = registration_redcapid)
exp <- exp %>% filter(registration_group == 'HC')

income <- readxl::read_excel('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/ explore_clock_income.xlsx') %>% rename(id = registration_redcapid)

exp <- inner_join(exp,income,by=c('id'))
exp <- exp %>% filter(registration_group.x == 'HC') %>% select(!registration_group.y) %>% rename(group = registration_group.x)


################################
### Explore model comparison ###
################################


dm <- data.frame(EP = c(0,1,4E-5),MF=c(0.0120,0.8401,0.1479), dMF = c(0.0202,0.0681,0.0659),model = c('Full RL','SCEPTIC','Value + Uncertainty'))
dbor <- 3.6954e-07

setwd('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/')
pdf('Explore HC(N=27) Models.pdf',height=6,width=8)
gg1 <-ggplot(dm, aes(x=MF,y=model,xmin=MF-dMF,xmax=MF+dMF)) + geom_errorbar(width=0.2) + geom_point(size=5) +
  xlab('Model Frequency') +
  theme_minimal(base_size = 18) + annotate("text",x=0.6,y=1,label=paste0('BOR = ',as.character(signif(dbor,1))),size=8,color="black") +
  annotate("rect",xmin=0.4,xmax=0.8,ymin=0.8,ymax=1.2,alpha=0.2) +
  annotate("text",x=0.6,y=2,label="EP=1",size=8,color="black")
print(gg1)
dev.off()
