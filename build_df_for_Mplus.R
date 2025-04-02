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
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMClock_MPlus'
# load some source code that we'll need to use
#setwd(rootdir)
#source('get_trial_data.R')

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
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMclock_clock_Aug2023.Rdata')
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)

# load in the hippocampus data, filter to within 5s of RT
#load('/ix/cladouceur/DNPL/HC_clock_Aug2023.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
hc <- hc %>% filter(evt_time > -5 & evt_time < 5)

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

# Merge vmPFC and HC data; remove unscaled within-S variable (decon1)
Q <- inner_join(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
rm(vmPFC, hc)
gc()

# Get task behav data
# (had to download it from UNCDEPENdlab github first)
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')

# select and scale variables of interest
df <- df %>% 
  group_by(id, run) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
         iti_lag_sc = scale(iti_prev),
         v_max_sc = scale(v_max),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         rt_swing_sc = scale(rt_swing))

# select only vars of interest and merge into MRI data
behav <- df %>% select(id,run,trial,run_trial,v_chosen_sc,score_sc,iti_sc,iti_lag_sc,v_max_sc,rt_vmax_sc,
                       rt_lag_sc,rt_vmax_lag_sc,v_entropy_sc,rt_swing_sc,trial_neg_inv_sc,v_entropy_wi_change,outcome,
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev)
Q <- inner_join(behav, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# censor out previous and next trials
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCbetween[Q$evt_time < -(Q$iti_prev)] = NA;

# summarize over evt_time and atlas_value and pivot to create new columns with data averaged across network
Q1 <- Q %>% select(!decon1) %>% filter(evt_time==0) %>% group_by(id,run,run_trial,network,HC_region) %>% 
  summarize(vmPFC_decon = mean(vmPFC_decon,na.rm=TRUE),HCwithin = mean(HCwithin,na.rm=TRUE), HCbetween = mean(HCbetween,na.rm=TRUE)) %>% 
  ungroup()
Q1 <- Q1 %>% group_by(id,run,run_trial) %>% 
  pivot_wider(values_from=c(vmPFC_decon,HCwithin),names_from = 'network') %>% 
  ungroup()

# Get task behav data
# (had to download it from UNCDEPENdlab github first)
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')

# select and scale variables of interest
df <- df %>% 
  group_by(id, run) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
         iti_lag_sc = scale(iti_prev),
         v_max_sc = scale(v_max),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         rt_swing_sc = scale(rt_swing))

# select only vars of interest and merge into MRI data
behav <- df %>% select(id,run,trial,run_trial,v_chosen_sc,score_sc,iti_sc,iti_lag_sc,v_max_sc,rt_vmax_sc,
                       rt_lag_sc,rt_vmax_lag_sc,v_entropy_sc,rt_swing_sc,trial_neg_inv_sc,v_entropy_wi_change,outcome,
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev)
Q1 <- inner_join(behav, Q1, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# add in age and sex variables
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q1 <- inner_join(Q1,demo,by=c('id'))
Q1$female <- relevel(as.factor(Q1$female),ref='0')
Q1$age <- scale(Q1$age)

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

save(Q1_AH,file=file.path(rootdir,'mmclock_HC_vmPFC_clock_evt_time0_AHforMplus.Rdata'))
save(Q1_PH,file=file.path(rootdir,'mmclock_HC_vmPFC_clock_evt_time0_PHforMplus.Rdata'))

setwd(file.path(rootdir))

prepareMplusData(df = Q1_AH, filename = "mmclock_HC_vmPFC_clock_evt_time0_AH_forMplus_taa.dat", dummyCode = c("outcome", "female"), overwrite = TRUE)
prepareMplusData(df = Q1_PH, filename = "mmclock_HC_vmPFC_clock_evt_time0_PH_forMplus_taa.dat", dummyCode = c("outcome", "female"), overwrite = TRUE)

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
# prepareMplusData(df = Q1_AH, filename = "mmclock_HC_vmPFC_clock_AH_forMplus_withrs_taa.dat", dummyCode = c("outcome"), overwrite = TRUE)
# prepareMplusData(df = Q1_PH, filename = "mmclock_HC_vmPFC_clock_PH_forMplus_withrs_taa.dat", dummyCode = c("outcome"), overwrite = TRUE)

