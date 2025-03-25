# Ceci Westbrook, 10/2024
# building models for analysis using Michael Hallquist's mixed_by function
# Based on code by Andrew Papale

# Load required libraries
library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(fmri.pipeline)
library(data.table)
library(parallel)
library(doParallel)
library(lme4)
library(emmeans)
library(sciplot)
library(MplusAutomation)
library(lmerTest)

# set root directory - change for your needs
#rootdir <- '/ix/cladouceur/DNPL'
#repo_directory <- file.path(rootdir,'HC_vPFC_repo/vPFC-HC-Clock')
repo_directory <- "~/clock_analysis"
# load some source code that we'll need to use
#setwd(rootdir)
#source('get_trial_data.R')
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
repo_directory1 <- file.path('/Volumes/Users/Andrew/MEDuSA_data_BSOC')
# load mixed_by function for analyses
#source(file.path(repo_directory,"/mixed_by.R"))
#source(file.path(repo_directory,"/plot_mixed_by_vmPFC_HC.R"))
#source(file.path(rootdir,'/plot_mixed_by_vmPFC_HC_simplified_networksymmetry.R'))
source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')

##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

# load the vmPFC data, filter to within 5s of RT, and select vars of interest
vmPFC <- read_csv(file.path(repo_directory1,'clock_aligned_bsocial_vmPFC.csv.gz'))
# calculate some variables
split_ksoc_bsoc <- vmPFC %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300])
bsoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240])
ksoc <- rbind(ksoc,221193,221611,220691) # 220691 was bsoc, but was run with > 240 trials
bsoc <- rbind(bsoc,221973,219757,220419,221507,221842,440223)
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

# load in the hippocampus data, filter to within 5s of RT
#load('/ix/cladouceur/DNPL/HC_clock_Aug2023.Rdata')
load(file.path(rootdir,'BSOC_HC_clock_TRdiv2.Rdata'))
hc <- hc %>% filter(evt_time > -5 & evt_time < 5)

split_ksoc_bsoc <- hc %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300])
bsoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240])
ksoc <- rbind(ksoc,221193,221611,220691) # 220691 was bsoc, but was run with > 240 trials
bsoc <- rbind(bsoc,221973,219757,220419,221507,221842,440223)
hc_bsoc <- hc %>% filter(id %in% bsoc$id) %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                                                            trial > 40 & trial <= 80 ~ trial-40,
                                                                            trial > 80 & trial <=120 ~ trial-80, 
                                                                            trial > 120 & trial <=160 ~ trial-120,
                                                                            trial > 160 & trial <=200 ~ trial-160,
                                                                            trial > 200 & trial <=240 ~ trial-200))

hc_ksoc <- hc %>%  filter(id %in% ksoc$id) %>% mutate(run_trial0 = case_when(trial <= 50 ~ trial, 
                                                                             trial > 50 & trial <= 100 ~ trial-50,
                                                                             trial > 100 & trial <=150 ~ trial-100, 
                                                                             trial > 150 & trial <=200 ~ trial-150,
                                                                             trial > 200 & trial <=250 ~ trial-200,
                                                                             trial > 250 & trial <=300 ~ trial-250))

hc <- rbind(hc_bsoc,hc_ksoc) %>% select(!run_trial) %>% rename(run_trial = run_trial0)

# Compress data from 12 bins to 2 by averaging across anterior 6 bins and posterior 6 bins to create
# AH and PH averages
hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>%
  summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% 
  ungroup() # 12 -> 2

# Create a new scaled within-subjects variable (HCwithin) and a between-subjects
# variable averaged per subject and run (HCbetween)
hc <- hc %>% group_by(id,run) %>%
  mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>%
  ungroup()

rm(hc_bsoc,hc_ksoc)
gc()

# Merge vmPFC and HC data; remove unscaled within-S variable (decon1)
Q <- inner_join(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
rm(hc,vmPFC)
gc()

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
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev,run_trial0_neg_inv_sc,rewFunc,outcome)
behav <- behav %>% rename(run = scanner_run) %>% select(!run_trial) %>% rename(run_trial = run_trial0)
Q <- inner_join(behav, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# censor out previous and next trials
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCbetween[Q$evt_time < -(Q$iti_prev)] = NA;

# summarize over evt_time and atlas_value and pivot to create new columns with data averaged across network
Q1 <- Q %>% select(!decon1) %>% group_by(id,run,run_trial,network,HC_region) %>% 
  summarize(vmPFC_decon = mean(vmPFC_decon,na.rm=TRUE),HCwithin = mean(HCwithin,na.rm=TRUE), HCbetween = mean(HCbetween,na.rm=TRUE)) %>% 
  ungroup()
Q1 <- Q1 %>% group_by(id,run,run_trial) %>% 
  pivot_wider(values_from=c(vmPFC_decon,HCwithin),names_from = 'network') %>% 
  ungroup()

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
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev,run_trial0_neg_inv_sc,rewFunc,outcome)
behav <- behav %>% rename(run = scanner_run) %>% select(!run_trial) %>% rename(run_trial = run_trial0)
Q1 <- inner_join(behav, Q1, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# add in age and sex variables
# add in age and sex variables
demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(repo_directory1,'2025-02-27-Partial-demo-pull-KSOC.csv'))
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
Q1 <- inner_join(Q1,demo2,by=c('id'))
Q1$female <- ifelse(Q1$sex==1,1,0)
Q1 <- Q1 %>% select(!sex)
Q1$age <- scale(Q1$age)
#Q <- Q %>% filter(group=='HC')
Q1$group <- relevel(factor(Q1$group),ref='HC')

Q1_AH <- Q1 %>% filter(HC_region == "AH")
Q1_PH <- Q1 %>% filter(HC_region == "PH")

# test whether MLMs are still significant

# m1 <- lmer(vmPFC_decon_CTR ~ HCwithin_CTR*female + (1|id), data = Q1_AH)
# summary(m1)
# 
# m2 <- lmer(vmPFC_decon_DMN ~ HCwithin_DMN*female + (1|id), data = Q1_AH)
# summary(m2)
# 
# m3 <- lmer(vmPFC_decon_LIM ~ HCwithin_LIM*female + (1|id), data = Q1_AH)
# summary(m3)
# 
# m4 <- lmer(vmPFC_decon_CTR ~ HCwithin_CTR*female + (1|id), data = Q1_PH)
# summary(m4)
# 
# m5 <- lmer(vmPFC_decon_DMN ~ HCwithin_DMN*female + (1|id), data = Q1_PH)
# summary(m5)
# 
# m6 <- lmer(vmPFC_decon_LIM ~ HCwithin_LIM*female + (1|id), data = Q1_PH)
# summary(m6)

save(Q1_AH,file=file.path(rootdir,'bsocial_HC_vmPFC_clock_AHforMplus.Rdata'))
save(Q1_PH,file=file.path(rootdir,'bsocial_HC_vmPFC_clock_PHforMplus.Rdata'))

setwd(file.path(rootdir))

prepareMplusData(df = Q1_AH, filename = "bsocial_HC_vmPFC_clock_AH_forMplus_taa.dat", dummyCode = c("outcome", "female"), overwrite = TRUE)
prepareMplusData(df = Q1_PH, filename = "bsocial_HC_vmPFC_clock_PH_forMplus_taa.dat", dummyCode = c("outcome", "female"), overwrite = TRUE)

# Alright! Was able to run the random slopes models. Now need to extract random slopes and put them
# back into the dataframe. Can use Michael's MplusAutomation for that.
# library(MplusAutomation)
# 
# list.files(getwd()) # let's see what we got in here
# 
# # let's make a function
# getslopes <- function(HC_region, network){
#   # get filename
#   filename <- paste0("mmclock_get_HC_vmPFC_randslopes_",HC_region,"_",network,".out")
#   renamevar <- paste0(HC_region,"_",network,"_rs")
#   
#   #read in model
#   wm_rs <- readModels(filename)
#   
#   ## save the savedata from the mplus output to an object
#   slopes <- wm_rs$savedata %>% 
#     dplyr::select(ID, RS.Mean) %>% 
#     group_by(ID) %>% slice_head() %>% ungroup() %>%
#     rename(id = ID, !! sym(renamevar) := RS.Mean)
#   
#   return(slopes)
# }
# 
# # ok we ran six models: 3 networks (DMN, CTR, LIM) x 2 HC regions (AH, PH).
# for(netname in c("CTR","DMN","LIM")){
#   Q1_AH <- merge.data.frame(Q1_AH,getslopes("AH",netname),by="id")
#   Q1_PH <- merge.data.frame(Q1_PH,getslopes("PH",netname),by="id")
# }
# 
# # now saving these again
# #prepareMplusData(df = Q1_AH, filename = "mmclock_HC_vmPFC_clock_AH_forMplus_withrs_taa.dat", dummyCode = c("outcome"), overwrite = TRUE)
# #prepareMplusData(df = Q1_PH, filename = "mmclock_HC_vmPFC_clock_PH_forMplus_withrs_taa.dat", dummyCode = c("outcome"), overwrite = TRUE)
# 
