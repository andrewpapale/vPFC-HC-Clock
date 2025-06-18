# Ceci Westbrook, 6/2025
# code for vPFC-HC sex differences paper
# building models for analysis using Michael Hallquist's mixed_by function
# Based on code by Andrew Papale

# modified 2025-06-09 AndyP HC only in Bsocial sample --> All subjects in BSocial sample

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

# set root directory - change for your needs
rootdir <- '/ix/cladouceur/DNPL'
repo_directory <- file.path(rootdir,'HC_vPFC_repo/vPFC-HC-Clock')

# load some source code that we'll need to use
setwd(rootdir)
source('get_trial_data.R')

# load mixed_by function for analyses
source(file.path(repo_directory,"/mixed_by.R"))
source(file.path(repo_directory,"/plot_mixed_by_vmPFC_HC.R"))
source(file.path(rootdir,'/plot_mixed_by_vmPFC_HC_simplified_networksymmetry.R'))

######################################################################################################
######################################################################################################
#########################################   STUDY 1   ################################################
######################################################################################################
######################################################################################################

##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

# load the vmPFC data, filter to within 5s of RT, and select vars of interest
load('MMclock_clock_Aug2023.Rdata')
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)

# load in the hippocampus data, filter to within 5s of RT
load('/ix/cladouceur/DNPL/HC_clock_Aug2023.Rdata')
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
Q <- Q %>% select(!decon1)

# Get task behav data
# (had to download it from UNCDEPENdlab github first)
df <- get_trial_data(repo_directory=rootdir,dataset='mmclock_fmri')

# select and scale variables of interest
df <- df %>% 
  group_by(id, run) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
         iti_lag_sc = scale(iti_prev))

# select only vars of interest and merge into MRI data
behav <- df %>% select(id,run,trial,run_trial,rt_lag_sc,trial_neg_inv_sc,rt_csv_sc,rt_csv,iti_ideal,iti_prev)
Q <- inner_join(behav, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# censor out previous and next trials
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCbetween[Q$evt_time < -(Q$iti_prev)] = NA;

# add in age and sex variables
demo <- read.table(file=file.path(rootdir, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)

save(Q,file=file.path(rootdir,'mmclock_HC_vmPFC_clock.Rdata'))

########################
##### Set up models ####
########################

setwd(paste0(rootdir,'/fmri/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# sex and hippocampal activity only
decode_formula[[1]] = formula(~female*HCwithin + HCbetween + (1 | id/run)) 

# sex, age and hippocampal activity 
decode_formula[[2]] = formula(~female*HCwithin + age*HCwithin + HCbetween + (1 | id/run)) 

# sex, hippocampal activity and control vars (run, trial within run)
decode_formula[[3]] = formula(~female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

# adding age
decode_formula[[4]] = formula(~age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

# with random slope for HCwithin
decode_formula[[5]] = formula(~female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 + HCwithin | id/run)) 

# with age and random slope for HCwithin
decode_formula[[6]] = formula(~age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 + HCwithin | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~female:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  # write to file
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-randslopeHCwithin',i,'.Rdata'))
  
}
#####################
#### plot models ####
#####################

# load in data already run
setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
load(paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-randslopeHCwithin',i,'.Rdata'))

setwd(file.path(rootdir,'fmri/validate_mixed_by_clock_HC_interaction/'))

# Save all plots as pdf - Using Andrew's plot_mixed_by function
plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'MMClock',totest='censored-withrandslopeHCwithin',toalign='clock',
                       toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = 6)

#############################
#############################
#############################
#########  behavior  ########
#############################
#############################
#############################

# using mixed-by: load in fMRI and MEG data and split by dataset
df <- get_trial_data(repo_directory=rootdir,dataset='mmclock_fmri')
df <- df %>% select(outcome,rt_lag, rt_lag_sc, rt_csv_sc,rt_csv,id, run, run_trial, last_outcome,rt_vmax_lag_sc, trial_neg_inv_sc,total_earnings)

df$id <- as.character(df$id)
Q2 <- df;
demo <- read.table(file=file.path(rootdir, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)
Qfmri <- inner_join(Q2,demo,by=c('id'))
Qfmri$age <- scale(Qfmri$age)

df <- get_trial_data(repo_directory=rootdir,dataset='mmclock_meg')
df <- df %>% select(outcome,rt_lag, rt_lag_sc, rt_csv_sc,rt_csv,id, run, run_trial, last_outcome,rt_vmax_lag_sc, trial_neg_inv_sc,total_earnings)

df$id <- as.character(df$id)
Q2 <- df;
demo <- read.table(file=file.path(rootdir, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)
Qmeg <- inner_join(Q2,demo,by=c('id'))
#Qmeg$female <- relevel(as.factor(Qmeg$female),ref='0')
Qmeg$age <- scale(Qmeg$age)

Qmeg <- Qmeg %>% mutate(dataset = 'MEG')
Qfmri <- Qfmri %>% mutate(dataset = 'fMRI')
Q3 <- rbind(Qmeg,Qfmri)

Q3 <- Q3 %>% mutate(sex = case_when(female==0 ~ 'M',
                                    female==1 ~ 'F'))
Q3 <- Q3 %>% filter(rt_csv < 4 & rt_csv > 0.2)
Q3 <- Q3 %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

# look at total earnings by sex
earnings_fmri <- Q3[Q3$dataset=="fMRI",c("id","sex","total_earnings")] %>% distinct(id,sex,total_earnings)
earnings_MEG <- Q3[Q3$dataset=="MEG",c("id","sex","total_earnings")] %>% distinct(id,sex,total_earnings)
earnings_fmri$dataset <- "fMRI"
earnings_MEG$dataset <- "MEG"
earnings_all <- rbind(earnings_fmri,earnings_MEG)

# not significant
t.test(earnings_fmri[earnings_fmri$sex=="M","total_earnings"],earnings_fmri[earnings_fmri$sex=="F","total_earnings"])
t.test(earnings_MEG[earnings_MEG$sex=="M","total_earnings"],earnings_MEG[earnings_MEG$sex=="F","total_earnings"])

# plot:
tab <- earnings_all %>% group_by(sex,dataset) %>% summarise(TotalEarnings = mean(total_earnings),se = se(total_earnings))

g <- ggplot(tab,aes(y=TotalEarnings,x=sex),group=sex)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex)) +
  geom_errorbar(aes(ymin = TotalEarnings - se,ymax = TotalEarnings + se,color=sex),position=position_dodge(width=0.5)) +
  facet_wrap(~dataset) + theme(legend.position = "none")

# set up mixed_by models
setwd(paste0(rootdir,'/fmri/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
rm(decode_formula)
decode_formula <- NULL

# look at sex
decode_formula[[1]] <- formula(~rt_lag_sc*sex + (1 | id/run))

for (j in 1:length(decode_formula)){
  
  ddq <- mixed_by(Q3, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = "dataset",return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:sex), at = list(rt_lag_sc=c(-2,-1,0,1,2)))
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2)))
                  )
  )
  setwd(file.path(rootdir,'fmri/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-Age-Sex-clock-fmri_meg-pred-rt-withrandslope',j,'.Rdata'))
}

# sex by outcome
decode_formula[[2]] <- formula(~rt_lag_sc*reward_lag_rec*sex + (1 | id/run))

# sex only but more control variables
decode_formula[[3]] <- formula(~rt_lag_sc*reward_lag_rec*sex + rt_vmax_lag_sc*trial_neg_inv_sc*sex + rt_lag_sc*trial_neg_inv_sc*sex + (1 | id/run))

# control for age
decode_formula[[4]] <- formula(~age + rt_lag_sc*reward_lag_rec*sex + rt_vmax_lag_sc*trial_neg_inv_sc*sex + rt_lag_sc*trial_neg_inv_sc*sex + (1 | id/run))

for (j in 2:length(decode_formula)){
  
  ddq <- mixed_by(Q3, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = "dataset",return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:sex), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                specs=formula(~rt_lag_sc:reward_lag_rec:sex), at=list(rt_lag_sc=c(-2,-1,0,1,2)))
                    
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                specs=formula(~rt_lag_sc:reward_lag_rec:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2)))
                    )
  )
  setwd(file.path(rootdir,'fmri/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-Age-Sex-clock-fmri_meg-pred-rt-withrandslope',j,'.Rdata'))
}

# plot emtrends for sex
load(paste0(curr_date,'-Age-Sex-clock-fmri_meg-pred-rt-withrandslope3.Rdata'))
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="MEG" & p.value < 0.05)
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="fMRI" & p.value < 0.05)

tab <- as.data.frame(ddq$emtrends_list$RT) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex,dataset) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex),group=sex)
g + geom_point(stat="identity",position=position_dodge(width=0.5), aes(color=sex)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=sex),position=position_dodge(width=0.5)) +
  facet_wrap(~dataset) + scale_y_reverse() + theme(legend.position = "none") +
  ylab("RT Swing (RT ~ previous RT)")

# by sex and last outcome: won't work for model 1
tab <- as.data.frame(ddq$emtrends_list$RTxO) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex,dataset,reward_lag_rec) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
tab$last_outcome <- ifelse(tab$reward_lag_rec==-0.5, "Omission","Reward")

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex),group=dataset)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=dataset)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=dataset),position=position_dodge(width=0.5)) +
  facet_wrap(~last_outcome) + scale_y_reverse()

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=last_outcome),group=dataset)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=dataset)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=dataset),position=position_dodge(width=0.5)) +
  facet_wrap(~sex, labeller = as_labeller(c("F" = "Female","M" = "Male"))) + 
  ylab("RT Swing (RT ~ previous RT)") + xlab("Prior Trial Outcome") + theme(legend.position = "none")

######################################################################################################
######################################################################################################
#########################################   STUDY 2   ################################################
######################################################################################################
######################################################################################################

##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

rootdir <- '/ix/cladouceur/DNPL/BSOC'
repo_directory <- file.path('/ix/cladouceur/DNPL/HC_vPFC_repo/vPFC-HC-Clock')
setwd(rootdir)


# load the vmPFC data, filter to within 5s of RT, and select vars of interest
vmPFC <- read_csv(file.path(rootdir,'clock_aligned_bsocial_vmPFC.csv.gz'))

# calculate some variables
vmPFC <- vmPFC %>% group_by(id, run) %>% mutate(run_trial = trial - min(trial) + 1) %>% ungroup()
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

# now load in HC data
# load in the hippocampus data, filter to within 5s of RT
load(file.path(rootdir,'BSOC_HC_clock_TRdiv2.Rdata'))
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
Q <- Q %>% select(!decon1)

# Get task behav data
df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
split_ksoc_bsoc <- df %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300]
bsoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240]
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
df_bsoc <- df_bsoc %>% mutate(protocol = 'bsocial',
                              run_trial0_c = run_trial0-floor(run_trial0/40.5),
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
behav <- df %>% select(id,scanner_run,trial,run_trial,v_chosen_sc,score_sc,iti_sc,iti_lag_sc,v_max_sc,rt_vmax_sc,
                       rt_lag_sc,rt_vmax_lag_sc,v_entropy_sc,rt_swing_sc,trial_neg_inv_sc,last_outcome,
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev,run_trial0_neg_inv_sc,rewFunc)
behav <- behav %>% rename(run = scanner_run) %>% select(!run_trial) %>% mutate(run_trial = case_when(trial <= 150 ~ trial,
                                                                                                     trial > 150 ~ trial - 150))
Q <- inner_join(behav, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# censor out previous and next trials
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCbetween[Q$evt_time < -(Q$iti_prev)] = NA;

# add in age and sex variables
demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(rootdir,'2025-02-27-Partial-demo-pull-KSOC.csv'))
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
Q$female <- ifelse(Q$sex==1,1,0)
Q <- Q %>% select(!sex)
Q$age <- scale(Q$age)
#Q <- Q %>% filter(group=='HC') # 2025-06-09 AndyP
#Q$group <- relevel(factor(Q$group),ref='HC')

########################
##### Set up models ####
########################

setwd(paste0(rootdir,'/fMRI/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL


# sex and hippocampal activity only
decode_formula[[1]] = formula(~female*HCwithin + HCbetween + (1 | id/run)) 

# sex, age and hippocampal activity 
decode_formula[[2]] = formula(~female*HCwithin + age*HCwithin + HCbetween + (1 | id/run)) 

# sex, hippocampal activity and control vars (run, trial within run)
decode_formula[[3]] = formula(~female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

# adding age
decode_formula[[4]] = formula(~age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

# with random slope for HCwithin
decode_formula[[5]] = formula(~female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 + HCwithin | id/run)) 

# with age and random slope for HCwithin
decode_formula[[6]] = formula(~age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 + HCwithin | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~female:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  # write to file
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study2-randslopeHCwithin',i,'.Rdata'))
  
}

#####################
##### plot models ####
#####################

# load in data already run
setwd(paste0(rootdir,'/fMRI/mixed_by_output/'))
#load('2025-02-18censored-vmPFC-HC-network-clock-BSOC_allmodels_HConly7.Rdata')

# Save all plots as pdf - Using Andrew's plot_mixed_by function
plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'BSOC',totest='allmodels_HConly',toalign='clock',
                       toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = 6)

#############################
#############################
#############################
#########  behavior  ########
#############################
#############################
#############################

# Get task behav data
df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
split_ksoc_bsoc <- df %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300]
bsoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240]
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
df_bsoc <- df_bsoc %>% mutate(protocol = 'bsocial',
                              run_trial0_c = run_trial0-floor(run_trial0/40.5),
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
behav <- df %>% select(id,scanner_run,trial,run_trial,rt_vmax_sc,
                       rt_lag_sc,rt_vmax_lag_sc,trial_neg_inv_sc,last_outcome,total_earnings,
                      rt_csv_sc,rt_csv,run_trial0_neg_inv_sc)
behav <- behav %>% rename(run = scanner_run) %>% select(!run_trial) %>% mutate(run_trial = case_when(trial <= 150 ~ trial,
                                                                                                     trial > 150 ~ trial - 150))
Q <- behav

# add in age and sex variables
demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(rootdir,'2025-02-27-Partial-demo-pull-KSOC.csv'))
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
Q$age <- scale(Q$age)
Q <- Q %>% mutate(sex1 = case_when(sex == 1 ~ 'F', sex == 2 ~ 'M'))
Q <- Q %>% select(!sex) %>% rename(sex=sex1)
#Q <- Q %>% filter(group=='HC') 2025-06-09 AndyP
#Q$group <- relevel(factor(Q$group),ref='HC')
Q <- Q %>% filter(rt_csv < 4 & rt_csv > 0.2)
Q <- Q %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

# look at total earnings by sex
earnings_study2 <- Q[,c("id","sex","total_earnings")] %>% distinct(id,sex,total_earnings)

# not significant
t.test(earnings_study2[earnings_study2$sex=="M","total_earnings"],earnings_study2[earnings_study2$sex=="F","total_earnings"])

# plot:
tab <- earnings_study2 %>% group_by(sex) %>% summarise(TotalEarnings = mean(total_earnings),se = se(total_earnings))

g <- ggplot(tab,aes(y=TotalEarnings,x=sex),group=sex)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex)) +
  geom_errorbar(aes(ymin = TotalEarnings - se,ymax = TotalEarnings + se,color=sex),position=position_dodge(width=0.5))

########################
##### Set up models ####
########################

# set some baseline variables
ncores = 20

# split out by network
rm(decode_formula)
decode_formula <- NULL

# look at sex
decode_formula[[1]] <- formula(~rt_lag_sc*sex + (1 | id/run))

for (j in 1:length(decode_formula)){
  
  ddq <- mixed_by(Q, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:sex), at = list(rt_lag_sc=c(-2,-1,0,1,2)))
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2)))
                  )
  )
  setwd(file.path('/ix/cladouceur/DNPL/fmri/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-Age-Sex-clock-study2-pred-rt-withrandslope',j,'.Rdata'))
}

# sex by outcome
decode_formula[[2]] <- formula(~rt_lag_sc*reward_lag_rec*sex + (1 | id/run))

# sex only but more control variables
decode_formula[[3]] <- formula(~rt_lag_sc*reward_lag_rec*sex + rt_vmax_lag_sc*trial_neg_inv_sc*sex + rt_lag_sc*trial_neg_inv_sc*sex + (1 | id/run))

# control for age
decode_formula[[4]] <- formula(~age + rt_lag_sc*reward_lag_rec*sex + rt_vmax_lag_sc*trial_neg_inv_sc*sex + rt_lag_sc*trial_neg_inv_sc*sex + (1 | id/run))

for (j in 2:length(decode_formula)){
  
  ddq <- mixed_by(Q, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:sex), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                specs=formula(~rt_lag_sc:reward_lag_rec:sex), at=list(rt_lag_sc=c(-2,-1,0,1,2)))
                    
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                specs=formula(~rt_lag_sc:reward_lag_rec:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2)))
                  )
  )
  setwd(file.path('/ix/cladouceur/DNPL/fmri/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-Age-Sex-clock-study2-pred-rt-withrandslope',j,'.Rdata'))
}

# plot emtrends for sex
load(paste0(curr_date,'-Age-Sex-clock-study2-pred-rt-withrandslope3.Rdata'))
ddq$coef_df_reml %>% filter(effect=="fixed" & p.value < 0.05)

tab <- as.data.frame(ddq$emtrends_list$RT) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex),group=sex)
g + geom_point(stat="identity",position=position_dodge(width=0.5), aes(color=sex)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=sex),position=position_dodge(width=0.5)) +
  scale_y_reverse() + theme(legend.position = "none") + ylab("RT Swing (RT ~ previous RT)")

# by sex and last outcome: won't work for model 1
tab <- as.data.frame(ddq$emtrends_list$RTxO) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex,reward_lag_rec) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
tab$last_outcome <- ifelse(tab$reward_lag_rec==-0.5, "Omission","Reward")

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex),group=sex)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=sex),position=position_dodge(width=0.5)) +
  facet_wrap(~last_outcome) + scale_y_reverse()

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=last_outcome),group=sex)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=sex),position=position_dodge(width=0.5)) +
  facet_wrap(~sex, labeller = as_labeller(c("F" = "Female","M" = "Male"))) +
  ylab("RT Swing (RT ~ previous RT)") + xlab("Prior Trial Outcome") + theme(legend.position = "none")
