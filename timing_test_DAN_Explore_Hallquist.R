
library(tidyverse)
library(fmri.pipeline)
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
V1 <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','101','102','103','104','105','106','107','108','109','110','111','112','113','114','115')

# load Explore
rt_exp <- read.csv('/Volumes/Users/Andrew/MEDuSA_data_MMClock_Hallquist/feedback_aligned_transformed_DAN_7Network.csv.gz')
rt_exp <- rt_exp %>% filter(atlas_value %in% V1) %>% 
  select(id,run,trial,evt_time,decon_mean, atlas_value)

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory='~/clock_analysis',dataset='explore')
df <- df %>% select(id,run,trial,iti_prev,rt_csv)

df$id <- as.integer(df$id)
rt_exp <- rt_exp %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
rt_exp <- inner_join(df,rt_exp,by=c('id','run','trial'))
rt_exp <- rt_exp %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                                  trial > 40 & trial <= 80 ~ trial-40,
                                                  trial > 80 & trial <=120 ~ trial-80, 
                                                  trial > 120 & trial <=160 ~ trial-120,
                                                  trial > 160 & trial <=200 ~ trial-160,
                                                  trial > 200 & trial <=240 ~ trial-200))
rt_exp <- rt_exp %>% select(!trial)

decode_formula <- NULL
decode_formula[[1]] = formula(~(1|id))
ncores <- 26
splits = c('evt_time','atlas_value')
setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
df0 <- decode_formula[[1]]
print(df0)
ddf_hallquist <- mixed_by(rt_exp, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf_hallquist,file=paste0(curr_date,'-Explore-Hallquist-Pipeline-timing-test-',1,'.Rdata'))


V1 <- c('4','7','5','9','11','204','206','202','209','211')
# load Explore
rt_bsoc <- read.csv('/Volumes/Users/Andrew/DAN/fixed DAN/feedback_aligned_bsocial_dan.csv.gz')
rt_bsoc <- rt_exp %>% filter(atlas_value %in% V1)

split_ksoc_bsoc <- rt_bsoc %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300])
bsoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240])
ksoc <- rbind(ksoc,221193,221611,220691) # 220691 was bsoc, but was run with > 240 trials
bsoc <- rbind(bsoc,221973,219757,220419,221507,221842,440223)
vmPFC_bsoc <- rt_bsoc %>% filter(id %in% bsoc$id) %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                                                                  trial > 40 & trial <= 80 ~ trial-40,
                                                                                  trial > 80 & trial <=120 ~ trial-80, 
                                                                                  trial > 120 & trial <=160 ~ trial-120,
                                                                                  trial > 160 & trial <=200 ~ trial-160,
                                                                                  trial > 200 & trial <=240 ~ trial-200))

vmPFC_ksoc <- rt_bsoc %>%  filter(id %in% ksoc$id) %>% mutate(run_trial0 = case_when(trial <= 50 ~ trial, 
                                                                                   trial > 50 & trial <= 100 ~ trial-50,
                                                                                   trial > 100 & trial <=150 ~ trial-100, 
                                                                                   trial > 150 & trial <=200 ~ trial-150,
                                                                                   trial > 200 & trial <=250 ~ trial-200,
                                                                                   trial > 250 & trial <=300 ~ trial-250))

rt_bsoc <- rbind(vmPFC_bsoc,vmPFC_ksoc) %>% rename(run_trial = run_trial0)

rt_bsoc <- rt_bsoc %>% mutate(run1 = case_when(run=='run1'~1,run=='run2'~2)) %>% select(!run) %>% rename(run=run1)
rt_bsoc$id <- as.character(rt_bsoc$id)
# select vars of interest
rt_bsoc <- rt_bsoc %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time)

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
behav <- df %>% select(id,scanner_run,trial,run_trial,run_trial0,v_chosen_sc,score_sc,iti_sc,iti_lag_sc,v_max_sc,rt_vmax_sc,
                       rt_lag_sc,rt_vmax_lag_sc,v_entropy_sc,rt_swing_sc,trial_neg_inv_sc,last_outcome,
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev,run_trial0_neg_inv_sc,rewFunc)
behav <- behav %>% rename(run = scanner_run) %>% select(!run_trial) %>% rename(run_trial = run_trial0)

rt_bsoc <- inner_join(rt_bsoc,behav,by=c('id','run','run_trial'))


decode_formula <- NULL
decode_formula[[1]] = formula(~(1|id))
ncores <- 26
splits = c('evt_time','atlas_value')
setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
df0 <- decode_formula[[1]]
print(df0)
ddf_fmriprep <- mixed_by(rt_bsoc, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf_fmriprep,file=paste0(curr_date,'-Bsocial-Fmriprep-Pipeline-timing-test-',1,'.Rdata'))
