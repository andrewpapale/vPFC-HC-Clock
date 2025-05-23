

# libraries we'll need
library(tidyverse)
library(fmri.pipeline)
# set root directory
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
repo_directory <- file.path('/Volumes/Users/Andrew/MEDuSA_data_BSOC')
ncores <- 26
# load mixed_by function for analyses


##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

# load the vmPFC data, filter to within 5s of RT, and select vars of interest
vmPFC <- read_csv(file.path(repo_directory,'clock_aligned_bsocial_vmPFC.csv.gz'))

# calculate some variables
split_ksoc_bsoc <- vmPFC %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300])
bsoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240])
bsoc <- rbind(bsoc,221973,221507,221842,440223)
vmPFC_bsoc <- vmPFC %>% filter(id %in% bsoc$id) %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                                                                  trial > 40 & trial <= 80 ~ trial-40,
                                                                                  trial > 80 & trial <=120 ~ trial-80, 
                                                                                  trial > 120 & trial <=160 ~ trial-120,
                                                                                  trial > 160 & trial <=200 ~ trial-160,
                                                                                  trial > 200 & trial <=240 ~ trial-200))

vmPFC_ksoc <- vmPFC %>%  filter(id %in% ksoc$id) %>% mutate(run_trial0 = case_when(trial <= 50 ~ trial, 
                                                                                   trial > 50 & trial <= 100 ~ trial-50,
                                                                                   trial > 100 & trial <=150 ~ trial-100, 
                                                                                   trial > 150 & trial <=200 ~ trial-150,
                                                                                   trial > 200 & trial <=250 ~ trial-200,
                                                                                   trial > 250 & trial <=300 ~ trial-250))

vmPFC <- rbind(vmPFC_bsoc,vmPFC_ksoc) %>% rename(run_trial = run_trial0)
# network and symmetry group as per the mmclock data:
vmPFC <- vmPFC %>% mutate(network = case_when(atlas_value %in% c(55, 56, 159, 160) ~ 'LIM',
                                              atlas_value %in% c(65,66,67, 170, 171) ~ 'CTR',
                                              atlas_value %in% c(84,86,88,89,161,191,192,194) ~ 'DMN'),
                          symmetry_group = case_when(atlas_value %in% c(65, 66, 170) ~ 1,
                                                     atlas_value %in% c(86, 161) ~ 2,
                                                     atlas_value %in% c(56, 159) ~ 3,
                                                     atlas_value %in% c(84, 191) ~ 4,
                                                     atlas_value %in% c(88, 192) ~ 5,
                                                     atlas_value %in% c(67, 171) ~ 6,
                                                     atlas_value %in% c(89, 194) ~ 7,
                                                     atlas_value %in% c(55, 160) ~ 8))
vmPFC <- vmPFC %>% mutate(run1 = case_when(run=='run1'~1,run=='run2'~2)) %>% select(!run) %>% rename(run=run1)
vmPFC$id <- as.character(vmPFC$id)
# select vars of interest
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
# filter to only event times of interest
vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
rm(vmPFC_ksoc,vmPFC_bsoc)
gc()

# now load in HC data
# load in the hippocampus data, filter to within 5s of RT
# load(file.path(rootdir,'BSOC_HC_clock_TRdiv2.Rdata'))
# hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
# 
# # Compress data from 12 bins to 2 by averaging across anterior 6 bins and posterior 6 bins to create
# # AH and PH averages
# hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>%
#   summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% 
#   ungroup() # 12 -> 2
# 
# # Create a new scaled within-subjects variable (HCwithin) and a between-subjects
# # variable averaged per subject and run (HCbetween)
# hc <- hc %>% group_by(id,run) %>%
#   mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>%
#   ungroup()
# 
# # Merge vmPFC and HC data; remove unscaled within-S variable (decon1)
# Q <- inner_join(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
# Q <- Q %>% select(!decon1)

# Get task behav data
df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
split_ksoc_bsoc <- df %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300]
bsoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240]
bsoc <- rbind(bsoc,221973,221507,221842,440223)
df_bsoc <- df %>% filter(id %in% bsoc) %>% mutate(block = case_when(trial <= 40 ~ 1, 
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
df_bsoc <- df_bsoc %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                  run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                  run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))

df_ksoc <- df %>% filter(id %in% ksoc) %>% mutate(block = case_when(trial <= 50 ~ 1, 
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
df_ksoc <- df_ksoc %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/50.5),
                    run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                    run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))

df <- rbind(df_bsoc,df_ksoc)
#df <- df_bsoc
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
         rt_swing_sc = scale(rt_swing),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% ungroup()

# select only vars of interest and merge into MRI data
behav <- df %>% select(id,scanner_run,trial,run_trial0,run_trial,asc_trial,v_chosen_sc,score_sc,iti_sc,iti_lag_sc,v_max_sc,rt_vmax_sc,
                       rt_lag_sc,rt_vmax_lag_sc,v_entropy_sc,rt_swing_sc,trial_neg_inv_sc,last_outcome,
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev,run_trial0_neg_inv_sc,rewFunc,v_entropy_wi_change_lag,rt_vmax_lag_sc,abs_pe_max_lag_sc,rt_vmax_change_sc)
behav <- behav %>% rename(run = scanner_run) %>% select(!run_trial) %>% rename(run_trial = run_trial0)
Q <- inner_join(behav, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# censor out previous and next trials
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
# Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
# Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
# Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;
# Q$HCbetween[Q$evt_time < -(Q$iti_prev)] = NA;

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
Q <- inner_join(Q,demo2,by=c('id'))
Q <- Q %>% filter(age <= 50)
Q$age <- scale(Q$age)
Q$group <- relevel(factor(Q$group),ref='HC')
#Q <- Q %>% filter(rewFunc == 'IEV')
#Q <- Q %>% filter(group == 'HC')

# scan <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_BSOC/bsoc_ksoc_scan_which.csv')
# scan <- scan %>% select(registration_redcapid,scan_which) %>% rename(id = registration_redcapid, scanner = scan_which)
# scan$id <- as.character(scan$id)
# Q <- inner_join(Q,scan,by='id')
Q <- Q %>% mutate(sex = case_when(sex==1 ~ 'M',
                                    sex==2 ~ 'F'))
Q$sex <- relevel(as.factor(Q$sex),ref='M')
#Q <- Q %>% filter(group=='HC')
#Q <- Q %>% filter(female==0)
# now to add in model fits
#fits <- read_csv('fMRIEmoClock_decay_factorize_selective_psequate_fixedparams_fmri_mfx_sceptic_global_statistics.csv')
#fits <- fits %>% rename(old_id=id)
#fits$id<-gsub("_1","",fits$id)
#Q <- inner_join(Q,fits,by='id')
rm(decode_formula)
decode_formula <- NULL
# decode_formula[[1]] = formula(~  group + age + sex + run_trial0_neg_inv_sc + v_max_wi*sex + rt_lag_sc + last_outcome  + (1|id/run))
# decode_formula[[2]] = formula(~ group + age + sex + run_trial0_neg_inv_sc + v_entropy_wi*sex + rt_lag_sc + last_outcome + (1|id/run))
# decode_formula[[3]] = formula(~ group + age + sex + run_trial0_neg_inv_sc + v_max_wi*sex + v_entropy_wi*sex + rt_lag_sc + last_outcome + (1|id/run))
# decode_formula[[4]] = formula(~ group  + age + sex + run_trial0_neg_inv_sc + v_max_wi*age + rt_lag_sc + last_outcome  + (1|id/run))
# decode_formula[[5]] = formula(~ group  + age + sex + run_trial0_neg_inv_sc + v_entropy_wi*age + rt_lag_sc + last_outcome + (1|id/run))
# decode_formula[[6]] = formula(~ group  + age + sex + run_trial0_neg_inv_sc + v_max_wi*age + v_entropy_wi*age + rt_lag_sc + last_outcome + (1|id/run))
# decode_formula[[7]] = formula(~ group  + age + sex + run_trial0_neg_inv_sc + v_max_wi + v_entropy_wi + rt_lag_sc + last_outcome*sex + (1|id/run))
# decode_formula[[8]] = formula(~ group  + age + sex + run_trial0_neg_inv_sc + v_max_wi + v_entropy_wi + rt_lag_sc + last_outcome*age + (1|id/run))
# decode_formula[[9]] <- formula(~group + v_entropy_wi + (1|id/run))
# decode_formula[[10]] <- formula(~group + v_max_wi + (1|id/run))

decode_formula[[1]] = formula(~ group + sex + v_entropy_wi_change_lag*age + rt_vmax_lag_sc*age + abs_pe_max_lag_sc*age + rt_vmax_change_sc*age +  + rt_lag_sc + iti_lag_sc + (1|id/run))
decode_formula[[2]] = formula(~ group + age + v_entropy_wi_change_lag*sex + rt_vmax_lag_sc*sex + abs_pe_max_lag_sc*sex + rt_vmax_change_sc*sex +  + rt_lag_sc + iti_lag_sc + (1|id/run))


qT2 <- c(-2.62,-0.544,0.372, 0.477)
qT1 <- c(-2.668, -0.12, 0.11, 0.258, 0.288, 0.308, 0.323, 0.348)
splits = c('evt_time','network')
#source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)#,
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
                  # #   HxW = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
                  # #              specs = formula(~v_entropy_wi:wtar),at=list(wtar = c(-1,-0.5,0,0.5,1))),
                  #   T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                  #               specs=formula(~trial_neg_inv_sc:HCwithin),at=list(trial_neg_inv_sc=qT1)),
                  #   T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                  #               specs=formula(~run_trial0_neg_inv_sc:HCwithin),at=list(run_trial0_neg_inv_sc=qT2))
                  # #   H = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
                  # #            specs = formula(~v_entropy_wi),at=list(v_entropy_wi=c(-2,-1,0,1,2))),
                  # #   HxY = list(outcome='vmPFC_decon',model_name='model1',var='v_entropy_wi',
                  # #              specs=formula(~v_entropy_wi:education_yrs),at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
                  # )
  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-Bsocial-vPFC-network-clock-All-RTcorrected-',i,'.Rdata'))
}

