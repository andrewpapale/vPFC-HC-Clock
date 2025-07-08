# Ceci Westbrook, 6/2025
# code for vPFC-HC sex differences paper
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
library(ggplot2)
library(dplyr)
library(purrr)

# set root directory - change for your needs
rootdir <- '/ix/cladouceur/DNPL'
repo_directory <- file.path(rootdir,'HC_vPFC_repo/vPFC-HC-Clock')

# load some source code that we'll need to use
setwd(rootdir)
source('get_trial_data.R')

# load mixed_by function for analyses
source(file.path(repo_directory,"/mixed_by.R"))
source(file.path(rootdir,"/plot_mixed_by_vmPFC_HC.R"))

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
behav <- df %>% select(id,run,trial,run_trial,rt_lag_sc,trial_neg_inv_sc,iti_ideal,iti_prev,iti_lag_sc,rt_csv,v_entropy_wi,v_max_wi)
Q <- inner_join(behav, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# censor out previous and next trials
Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;
Q$HCbetween[Q$evt_time < -(Q$iti_prev)] = NA;

# add in age and sex variables
demo <- read.csv(file=file.path(rootdir, 'fmri/data/MMC_demog.csv'),header=TRUE)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$sex <- ifelse(Q$female==0,1,
                ifelse(Q$female==1,0,NA))
Q$age <- scale(Q$age)
Q$race <- Q %>% select("AmerIndianAlaskan":"White") %>% rowSums()

#save(Q,file=file.path(rootdir,'mmclock_HC_vmPFC_clock.Rdata'))

########################
##### Set up models ####
########################

setwd(paste0(rootdir,'/fmri/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(basemodel_formula)
basemodel_formula <- NULL

# determine the base model:
# start with sex and hippocampal activity first, then add trial number, iti, HCbetween iteratively and check model fit.

# sex and hippocampal activity only
basemodel_formula[[1]] = formula(~sex*HCwithin + (1 | id/run)) 

# sex, hippocampal activity and control vars (trial, iti)
basemodel_formula[[2]] = formula(~sex*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[3]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + (1 | id/run)) 
basemodel_formula[[4]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[5]] = formula(~sex*HCwithin + iti_lag_sc*HCwithin + (1 | id/run)) 
basemodel_formula[[6]] = formula(~sex*HCwithin + iti_lag_sc*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[7]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + iti_lag_sc*HCwithin + (1 | id/run))
basemodel_formula[[8]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + iti_lag_sc*HCwithin + HCbetween + (1 | id/run))

for(i in 1:length(basemodel_formula)){
  assign(paste0("df",i),basemodel_formula[[i]])
  print(get(paste0("df",i)))
  df0 <- get(paste0("df",i))
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0, split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~sex:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  assign(paste0("ddf",i),ddf)
}

plot_aic_comparison <- function(model_list) {

  # Extract AIC values into a tidy dataframe
  all_df <- imap_dfr(model_list, function(model, name) {
    data.frame(
      Index = seq_len(nrow(model$fit_df)),
      AIC = model$fit_df$AIC,
      Model = name
    )
  })
  
  # Plot
  ggplot(all_df, aes(x = Index, y = AIC, color = Model)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(title = "AIC Values by Model",
         x = "Model Index",
         y = "AIC",
         color = "Model")
}

# can compare models 1, 2, 4, 8; 1, 3, 4, 8; 1, 5, 6, 8; 1, 5, 7, 8; 1, 3, 7, 8; 1, 5, 7, 8 (not all models are nested)

models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf4 = ddf4)
models <- list(ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf4 = ddf4, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf5 = ddf5)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf7 = ddf7)
plot_aic_comparison(models)

# No evidence that including trial or avg hippocampus activity changes model fit. Including ITI worsens it.
# therefore, for the sake of having more control variables accounted for, will include trial num and hippocampal activity.

# Now that we have determined a base model, we can do hypothesis testing:
# rt_lag_sc for RTswings, rt_lag_sc controlled for age, vmax, entropy

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# base effect of sex
decode_formula[[1]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run))

# now adding age
decode_formula[[2]] = formula(~age*HCwithin + sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# age by sex interaction
decode_formula[[3]] = formula(~sex*age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[4]] = formula(~sex*HCwithin + race*HCwithin + Hispanic*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# highest level of education
decode_formula[[5]] = formula(~sex*HCwithin + edu*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding RT lag
decode_formula[[6]] = formula(~sex*HCwithin + rt_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~sex*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of sex with vmax
decode_formula[[8]] = formula(~sex*v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~sex*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of sex with entropy
decode_formula[[10]] = formula(~sex*v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1  | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~sex:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  # write to file
  setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-randslopeHCwithin',i,'.Rdata'))
  
}
#####################
#### plot models ####
#####################

# load in data already run

for(j in 1:length(decode_formula)){
  setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
  load(paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-randslopeHCwithin',j,'.Rdata'))
  
  setwd(file.path(rootdir,'fmri/validate_mixed_by_clock_HC_interaction/'))
  
  # Save all plots as pdf - Using Andrew's plot_mixed_by function
  plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'MMClock',totest='censored-withrandslopeHCwithin',toalign='clock',
                       toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = j)
}

#########################################################################################################
#########################################################################################################

# Now, we can do the same for age!
setwd(paste0(rootdir,'/fmri/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(basemodel_formula)
basemodel_formula <- NULL

# determine the base model:
# start with age and hippocampal activity first, then add trial number, iti, HCbetween iteratively and check model fit.

# sex and hippocampal activity only
basemodel_formula[[1]] = formula(~age*HCwithin + (1 | id/run)) 
# determine the base model:
# start with sex and hippocampal activity first, then add trial number, iti, HCbetween iteratively and check model fit.

# sex and hippocampal activity only
basemodel_formula[[1]] = formula(~sex*HCwithin + (1 | id/run)) 

# sex, hippocampal activity and control vars (trial, iti)
basemodel_formula[[2]] = formula(~age*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[3]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + (1 | id/run)) 
basemodel_formula[[4]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[5]] = formula(~age*HCwithin + iti_lag_sc*HCwithin + (1 | id/run)) 
basemodel_formula[[6]] = formula(~age*HCwithin + iti_lag_sc*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[7]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + iti_lag_sc*HCwithin + (1 | id/run))
basemodel_formula[[8]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + iti_lag_sc*HCwithin + HCbetween + (1 | id/run))

for(i in 1:length(basemodel_formula)){
  assign(paste0("df",i),basemodel_formula[[i]])
  print(get(paste0("df",i)))
  df0 <- get(paste0("df",i))
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0, split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  )
  assign(paste0("ddf",i),ddf)
}

# can compare models 1, 2, 4, 8; 1, 3, 4, 8; 1, 5, 6, 8; 1, 5, 7, 8; 1, 3, 7, 8; 1, 5, 7, 8 (not all models are nested)

models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf4 = ddf4)
models <- list(ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf4 = ddf4, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf5 = ddf5)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf7 = ddf7)
plot_aic_comparison(models)

# Once again: no evidence that including trial or avg hippocampus activity changes model fit. Including ITI worsens it.
# therefore, for the sake of having more control variables accounted for, will include trial num and hippocampal activity.

# Now that we have determined a base model, we can do hypothesis testing:
# rt_lag_sc for RTswings, rt_lag_sc controlled for age, vmax, entropy

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# base effect of age
decode_formula[[1]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run))

# now adding sex
decode_formula[[2]] = formula(~age*HCwithin + sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# age by sex interaction
decode_formula[[3]] = formula(~sex*age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[4]] = formula(~age*HCwithin + race*HCwithin + Hispanic*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# highest level of education
decode_formula[[5]] = formula(~age*HCwithin + level_edum*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding RT lag
decode_formula[[6]] = formula(~age*HCwithin + rt_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~age*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of age with vmax
decode_formula[[8]] = formula(~age*vmax_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~age*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of age with entropy
decode_formula[[10]] = formula(~age*v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~age:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  # write to file
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-agemodels-randslopeHCwithin',i,'.Rdata'))
  
}

#####################
#### plot models ####
#####################

# load in data already run

for(j in 1:length(decode_formula)){
  setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
  load(paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study1-agemodels-randslopeHCwithin',j,'.Rdata'))
  
  setwd(file.path(rootdir,'fmri/validate_mixed_by_clock_HC_interaction/'))
  
  # Save all plots as pdf - Using Andrew's plot_mixed_by function
  plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'MMClock',totest='censored-withrandslopeHCwithin',toalign='clock',
                         toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = j)
  
  setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
}

#############################
#############################
#############################
#########  behavior  ########
#############################
#############################
#############################

# using mixed-by: load in fMRI and MEG data and split by dataset
df <- get_trial_data(repo_directory=rootdir,dataset='mmclock_fmri')
df <- df %>% select(outcome,rt_lag, rt_lag_sc, rt_csv_sc,rt_csv,id, run, run_trial, last_outcome,v_max_wi,v_entropy_wi,trial_neg_inv_sc,total_earnings)

df$id <- as.character(df$id)
Q2 <- df;
demo <- read.csv(file=file.path(rootdir, 'fmri/data/MMC_demog.csv'),header=TRUE)
demo <- demo %>% select(!adult & !scandate)
demo$sex <- ifelse(demo$female==0,1,
                ifelse(demo$female==1,0,NA))
demo$id <- as.character(demo$id)
Qfmri <- inner_join(Q2,demo,by=c('id'))
Qfmri$age <- scale(Qfmri$age)

df <- get_trial_data(repo_directory=rootdir,dataset='mmclock_meg')
df <- df %>% select(outcome,rt_lag, rt_lag_sc, rt_csv_sc,rt_csv,id, run, run_trial, last_outcome,v_max_wi,v_entropy_wi, trial_neg_inv_sc,total_earnings)

df$id <- as.character(df$id)
Q2 <- df;
Q$female <- relevel(as.factor(Q$female),ref='0')

Qmeg <- inner_join(Q2,demo,by=c('id'))
#Qmeg$female <- relevel(as.factor(Qmeg$female),ref='0')
Qmeg$age <- scale(Qmeg$age)

Qmeg <- Qmeg %>% mutate(dataset = 'MEG')
Qfmri <- Qfmri %>% mutate(dataset = 'fMRI')
Q3 <- rbind(Qmeg,Qfmri)
Q3$race <- Q3 %>% select("AmerIndianAlaskan":"White") %>% rowSums()

Q3 <- Q3 %>% filter(rt_csv < 4 & rt_csv > 0.2)
Q3 <- Q3 %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

# look at total earnings by sex
earnings_fmri <- Q3[Q3$dataset=="fMRI",c("id","sex","total_earnings")] %>% distinct(id,sex,total_earnings)
earnings_MEG <- Q3[Q3$dataset=="MEG",c("id","sex","total_earnings")] %>% distinct(id,sex,total_earnings)
earnings_fmri$dataset <- "fMRI"
earnings_MEG$dataset <- "MEG"
earnings_all <- rbind(earnings_fmri,earnings_MEG)

# not significant
t.test(earnings_fmri[earnings_fmri$sex==0,"total_earnings"],earnings_fmri[earnings_fmri$sex==1,"total_earnings"])
t.test(earnings_MEG[earnings_MEG$sex==0,"total_earnings"],earnings_MEG[earnings_MEG$sex==1,"total_earnings"])

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

# sex only
decode_formula[[1]] <- formula(~rt_lag_sc*sex + (1 | id/run))

# matching base model from MEDuSA
decode_formula[[2]] <- formula(~rt_lag_sc*sex + rt_lag_sc*trial_neg_inv_sc + (1 | id/run))

# now adding age
decode_formula[[3]] = formula(~rt_lag_sc*sex + rt_lag_sc*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# age by sex interaction
decode_formula[[4]] = formula(~rt_lag_sc*sex*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[5]] = formula(~rt_lag_sc*sex + rt_lag_sc*race + rt_lag_sc*Hispanic + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# highest level of education
decode_formula[[6]] = formula(~rt_lag_sc*sex + rt_lag_sc*edu + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~rt_lag_sc*sex + rt_lag_sc*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# interaction of sex with vmax
decode_formula[[8]] = formula(~rt_lag_sc*sex*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~rt_lag_sc*sex + rt_lag_sc*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# interaction of sex with entropy
decode_formula[[10]] = formula(~rt_lag_sc*sex*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

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
  save(ddq,file=paste0(curr_date,'-Age-Sex-clock-Study1-fmri-meg-pred-rt-withrandslope',j,'.Rdata'))
}

# plot emtrends for sex
load(paste0(curr_date,'-Age-Sex-clock-mmclock-fmri-meg-pred-rt-withrandslope4.Rdata'))
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="MEG" & p.value < 0.05)
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="fMRI" & p.value < 0.05)

tab <- as.data.frame(ddq$emtrends_list$RT) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex,dataset) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex),group=sex)
g + geom_point(stat="identity",position=position_dodge(width=0.5), aes(color=sex)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=sex),position=position_dodge(width=0.5)) +
  facet_wrap(~dataset) + scale_y_reverse() + theme(legend.position = "none") +
  ylab("RT Swing (RT ~ previous RT)")

# look at differences by last_outcome
rm(decode_formula)
decode_formula <- NULL

# sex only
decode_formula[[1]] <- formula(~rt_lag_sc*sex*reward_lag_rec + (1 | id/run))

# matching base model from MEDuSA
decode_formula[[2]] <- formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run))

# now adding age
decode_formula[[3]] = formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*age*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[4]] = formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*race*reward_lag_rec + rt_lag_sc*Hispanic*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)) 

# highest level of education
decode_formula[[5]] = formula(~rt_lag_sc*sex*reward_lag_rec + rt_lag_sc*edu*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)) 

for (j in 1:length(decode_formula)){
  
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
  save(ddq,file=paste0(curr_date,'-Age-Sex-clock-Study1-fmri-meg-pred-rt-lastoutcome-withrandslope',j,'.Rdata'))
}

load(paste0(curr_date,'-Age-Sex-clock-Study1-fmri-meg-pred-rt-lastoutcome-withrandslope1.Rdata'))
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="MEG" & p.value < 0.05)
ddq$coef_df_reml %>% filter(effect=="fixed" & dataset=="fMRI" & p.value < 0.05)

# by sex and last outcome: won't work for models without outcome
tab <- as.data.frame(ddq$emtrends_list$RTxO) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex,dataset,reward_lag_rec) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
tab$last_outcome <- ifelse(tab$reward_lag_rec==-0.5, "Omission","Reward")
tab$sex1 <- ifelse(tab$sex==0,"Female","Male")

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex1),group=dataset)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=dataset)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=dataset),position=position_dodge(width=0.5)) +
  facet_wrap(~last_outcome) + scale_y_reverse()

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=last_outcome),group=dataset)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=dataset)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=dataset),position=position_dodge(width=0.5)) +
  facet_wrap(~sex1) + 
  ylab("RT Swing (RT ~ previous RT)") + xlab("Prior Trial Outcome")

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
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev)
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
demo$race <- demo %>% select(registration_race___1:registration_race___999) %>% rowSums()
demo$edu <- demo$registration_edu
demo <- demo %>% rename(sex=registration_birthsex,
                        gender=registration_gender,
                        group=registration_group) %>%
  select(id,group,age,sex,gender,race,edu)
demo1$race <- demo1 %>% select(registration_race___1:registration_race___999) %>% rowSums()
demo1$edu <- demo1$registration_edu
demo1 <- demo1 %>% rename(sex=registration_birthsex,
                          gender=registration_gender,
                          group=registration_group) %>%
  select(id,group,age,sex,gender,race,edu)
demo2 <- rbind(demo,demo1)
Q <- inner_join(Q,demo2,by=c('id'))
Q$female <- ifelse(Q$sex==1,1,0)
Q <- Q %>% select(!sex)
# recoding to make female the reference (0) group
Q$sex <- ifelse(Q$female==1,0,
                ifelse(Q$female==0,1,NA))
Q$age <- scale(Q$age)
# in case we want to look only at healthy controls:
#Q <- Q %>% filter(group=='HC')
#Q$group <- relevel(factor(Q$group),ref='HC')

# Remove some participants:
# individuals > 50 years old
Q <- Q[!Q$id %in% c(203521, 219084, 220562, 220886, 221193, 440111, 440151),]

# scanner timing was off for this individual
Q <- Q[!Q$id==221193,]

# this one run failed preprocessing: 440010 run 1
Q <- Q[!Q$id==440010 & !Q$run==1,]


########################
##### Set up models ####
########################

setwd(paste0(rootdir,'/fMRI/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(basemodel_formula)
basemodel_formula <- NULL

# determine the base model:
# start with sex and hippocampal activity first, then add trial number, iti, HCbetween iteratively and check model fit.

# sex and hippocampal activity only
basemodel_formula[[1]] = formula(~sex*HCwithin + (1 | id/run)) 

# sex, hippocampal activity and control vars (trial, iti)
basemodel_formula[[2]] = formula(~sex*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[3]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + (1 | id/run)) 
basemodel_formula[[4]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[5]] = formula(~sex*HCwithin + iti_prev*HCwithin + (1 | id/run)) 
basemodel_formula[[6]] = formula(~sex*HCwithin + iti_prev*HCwithin + HCbetween + (1 | id/run)) 
basemodel_formula[[7]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + iti_prev*HCwithin + (1 | id/run))
basemodel_formula[[8]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + iti_prev*HCwithin + HCbetween + (1 | id/run))

for(i in 1:length(basemodel_formula)){
  assign(paste0("df",i),basemodel_formula[[i]])
  print(get(paste0("df",i)))
  df0 <- get(paste0("df",i))
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0, split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~sex:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  assign(paste0("ddf",i),ddf)
}

plot_aic_comparison <- function(model_list) {
  
  # Extract AIC values into a tidy dataframe
  all_df <- imap_dfr(model_list, function(model, name) {
    data.frame(
      Index = seq_len(nrow(model$fit_df)),
      AIC = model$fit_df$AIC,
      Model = name
    )
  })
  
  # Plot
  ggplot(all_df, aes(x = Index, y = AIC, color = Model)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(title = "AIC Values by Model",
         x = "Model Index",
         y = "AIC",
         color = "Model")
}

# can compare models 1, 2, 4, 8; 1, 3, 4, 8; 1, 5, 6, 8; 1, 5, 7, 8; 1, 3, 7, 8; 1, 5, 7, 8 (not all models are nested)

models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf4 = ddf4)
models <- list(ddf4 = ddf4, ddf8 = ddf8)
models <- list(ddf4 = ddf4, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf2 = ddf2, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf4 = ddf4)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf6 = ddf6)
models <- list(ddf1 = ddf1, ddf5 = ddf5)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7)
models <- list(ddf1 = ddf1, ddf3 = ddf3, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf5 = ddf5, ddf7 = ddf7, ddf8 = ddf8)
models <- list(ddf1 = ddf1, ddf7 = ddf7)
plot_aic_comparison(models)

# same pattern as Study 1, so will use the same base model

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# base effect of sex
decode_formula[[1]] = formula(~sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run))

# now adding age
decode_formula[[2]] = formula(~age*HCwithin + sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# age by sex interaction
decode_formula[[3]] = formula(~sex*age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[4]] = formula(~sex*HCwithin + race*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# highest level of education
decode_formula[[5]] = formula(~sex*HCwithin + edu*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding RT lag
decode_formula[[6]] = formula(~sex*HCwithin + rt_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~sex*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of sex with vmax
decode_formula[[8]] = formula(~sex*v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~sex*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of sex with entropy
decode_formula[[10]] = formula(~sex*v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1  | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~sex:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  # write to file
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study2-randslopeHCwithin',i,'.Rdata'))
  
}
#####################
#### plot models ####
#####################

# load in data already run

for(j in 1:length(decode_formula)){
  setwd(file.path(rootdir,'fMRI/mixed_by_output/'))
  load(paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study2-randslopeHCwithin',j,'.Rdata'))
  
  setwd(file.path(rootdir,'fMRI/validate_mixed_by_clock_HC_interaction/'))
  
  # Save all plots as pdf - Using Andrew's plot_mixed_by function
  plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'BSOC',totest='censored-withrandslopeHCwithin',toalign='clock',
                         toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = j)
}

##############################################################################################################################
##############################################################################################################################
# Repeating for age

setwd(file.path(rootdir,'fMRI/mixed_by_output/'))

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# base effect of age
decode_formula[[1]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run))

# now adding sex
decode_formula[[2]] = formula(~age*HCwithin + sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# age by age interaction
decode_formula[[3]] = formula(~age*sex*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[4]] = formula(~age*HCwithin + race*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# highest level of education
decode_formula[[5]] = formula(~age*HCwithin + edu*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding RT lag
decode_formula[[6]] = formula(~age*HCwithin + rt_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~age*HCwithin + v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of age with vmax
decode_formula[[8]] = formula(~age*v_max_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~age*HCwithin + v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1 | id/run)) 

# interaction of age with entropy
decode_formula[[10]] = formula(~age*v_entropy_wi*HCwithin + trial_neg_inv_sc*HCwithin + HCbetween + (1  | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Sex = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~age:HCwithin), at = list(HCwithin=c(-1.5,1.5)))),
  )
  # write to file
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study2-agemodels-randslopeHCwithin',i,'.Rdata'))
  
}
#####################
#### plot models ####
#####################

# load in data already run

for(j in 1:length(decode_formula)){
  setwd(file.path(rootdir,'fMRI/mixed_by_output/'))
  load(paste0(curr_date,'_censored-vmPFC-HC-network-clock-Study2-agemodels-randslopeHCwithin',j,'.Rdata'))
  
  setwd(file.path(rootdir,'fMRI/validate_mixed_by_clock_HC_interaction/'))
  
  # Save all plots as pdf - Using Andrew's plot_mixed_by function
  plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'BSOC',totest='censored-withrandslopeHCwithin-agemodels',toalign='clock',
                         toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = j)
}

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
behav <- df %>% select(id,scanner_run,trial,run_trial,v_chosen_sc,score_sc,iti_sc,iti_lag_sc,v_max_sc,rt_vmax_sc,
                       rt_lag_sc,rt_vmax_lag_sc,v_entropy_sc,rt_swing_sc,trial_neg_inv_sc,last_outcome,
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,iti_prev,total_earnings)
behav <- behav %>% rename(run = scanner_run) %>% select(!run_trial) %>% mutate(run_trial = case_when(trial <= 150 ~ trial,
                                                                                                     trial > 150 ~ trial - 150))
Q <- behav

# add in age and sex variables
demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(rootdir,'2025-02-27-Partial-demo-pull-KSOC.csv'))
demo$id <- as.character(demo$id)
demo1$id <- as.character(demo1$registration_redcapid)
demo$race <- demo %>% select(registration_race___1:registration_race___999) %>% rowSums()
demo$edu <- demo$registration_edu
demo <- demo %>% rename(gender=registration_gender,
                        group=registration_group) %>%
  select(id,group,age,registration_birthsex,gender,race,edu)
demo1$race <- demo1 %>% select(registration_race___1:registration_race___999) %>% rowSums()
demo1$edu <- demo1$registration_edu
demo1 <- demo1 %>% rename(gender=registration_gender,
                          group=registration_group) %>%
  select(id,group,age,registration_birthsex,gender,race,edu)
demo2 <- rbind(demo,demo1)
demo2 <- demo2 %>% mutate(sex1 = case_when(registration_birthsex == 1 ~ 'F', registration_birthsex == 2 ~ 'M'))
demo2$sex <- ifelse(demo2$registration_birthsex==1,0,
                    ifelse(demo2$registration_birthsex==2,1,NA))
demo2 <- demo2 %>% select(!registration_birthsex)
Q <- inner_join(Q,demo2,by=c('id'))
Q$age <- scale(Q$age)
Q <- Q %>% filter(rt_csv < 4 & rt_csv > 0.2)
Q <- Q %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))
# in case we want to look only at healthy controls:
#Q <- Q %>% filter(group=='HC')
#Q$group <- relevel(factor(Q$group),ref='HC')

# Remove some participants:
# individuals > 50 years old
Q <- Q[!Q$id %in% c(203521, 219084, 220562, 220886, 221193, 440111, 440151),]

# scanner timing was off for this individual
Q <- Q[!Q$id==221193,]

# this one run failed preprocessing: 440010 run 1
Q <- Q[!Q$id==440010 & !Q$run==1,]


# look at total earnings by sex
earnings_study2 <- Q[,c("id","sex1","total_earnings")] %>% distinct(id,sex1,total_earnings)

# not significant
t.test(earnings_study2[earnings_study2$sex1=="M","total_earnings"],earnings_study2[earnings_study2$sex1=="F","total_earnings"])

# plot:
tab <- earnings_study2 %>% group_by(sex1) %>% summarise(TotalEarnings = mean(total_earnings),se = se(total_earnings))

g <- ggplot(tab,aes(y=TotalEarnings,x=sex1),group=sex1)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex1)) +
  geom_errorbar(aes(ymin = TotalEarnings - se,ymax = TotalEarnings + se,color=sex1),position=position_dodge(width=0.5))

########################
##### Set up models ####
########################

# set some baseline variables
ncores = 20

# split out by network
rm(decode_formula)
decode_formula <- NULL

# sex1 only
decode_formula[[1]] <- formula(~rt_lag_sc*sex1 + (1 | id/run))

# matching base model from MEDuSA
decode_formula[[2]] <- formula(~rt_lag_sc*sex1 + rt_lag_sc*trial_neg_inv_sc + (1 | id/run))

# now adding age
decode_formula[[3]] = formula(~rt_lag_sc*sex1 + rt_lag_sc*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# age by sex1 interaction
decode_formula[[4]] = formula(~rt_lag_sc*sex1*age + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[5]] = formula(~rt_lag_sc*sex1 + rt_lag_sc*race + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# highest level of education
decode_formula[[6]] = formula(~rt_lag_sc*sex1 + rt_lag_sc*edu + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding vmax
decode_formula[[7]] = formula(~rt_lag_sc*sex1 + rt_lag_sc*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# interaction of sex1 with vmax
decode_formula[[8]] = formula(~rt_lag_sc*sex1*v_max_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# now adding entropy
decode_formula[[9]] = formula(~rt_lag_sc*sex1 + rt_lag_sc*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

# interaction of sex1 with entropy
decode_formula[[10]] = formula(~rt_lag_sc*sex1*v_entropy_wi + rt_lag_sc*trial_neg_inv_sc + (1 | id/run)) 

for (j in 1:length(decode_formula)){
  
  ddq <- mixed_by(Q, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:sex1), at = list(rt_lag_sc=c(-2,-1,0,1,2)))
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:sex1), at=list(rt_lag_sc = c(-2,-1,0,1,2)))
                  )
  )
  setwd(file.path(rootdir,'/fMRI/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-Age-Sex-clock-study2-pred-rt-withrandslope',j,'.Rdata'))
}

# plot emtrends for sex
load(paste0(curr_date,'-Age-Sex-clock-study2-pred-rt-withrandslope6.Rdata'))
ddq$coef_df_reml %>% filter(effect=="fixed" & p.value < 0.05)

tab <- as.data.frame(ddq$emtrends_list$RT) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex1) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex1),group=sex1)
g + geom_point(stat="identity",position=position_dodge(width=0.5), aes(color=sex1)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=sex1),position=position_dodge(width=0.5)) +
  ylab("RT Swing (RT ~ previous RT)") + scale_y_reverse()

# look at differences by last_outcome
rm(decode_formula)
decode_formula <- NULL

# sex1 only
decode_formula[[1]] <- formula(~rt_lag_sc*sex1*reward_lag_rec + (1 | id/run))

# matching base model from MEDuSA
decode_formula[[2]] <- formula(~rt_lag_sc*sex1*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run))

# now adding age
decode_formula[[3]] = formula(~rt_lag_sc*sex1*reward_lag_rec + rt_lag_sc*age*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)) 

# now adding race and ethnicity
decode_formula[[4]] = formula(~rt_lag_sc*sex1*reward_lag_rec + rt_lag_sc*race*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)) 

# highest level of education
decode_formula[[5]] = formula(~rt_lag_sc*sex1*reward_lag_rec + rt_lag_sc*edu*reward_lag_rec + rt_lag_sc*trial_neg_inv_sc*reward_lag_rec + (1 | id/run)) 

for (j in 1:length(decode_formula)){
  
  ddq <- mixed_by(Q, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:sex1), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                specs=formula(~rt_lag_sc:reward_lag_rec:sex1), at=list(rt_lag_sc=c(-2,-1,0,1,2)))
                    
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:sex1), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                specs=formula(~rt_lag_sc:reward_lag_rec:sex1), at=list(rt_lag_sc = c(-2,-1,0,1,2)))
                  )
  )
  setwd(file.path(rootdir,'fMRI/mixed_by_output/'))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddq,file=paste0(curr_date,'-Age-Sex-clock-Study2-pred-rt-lastoutcome-withrandslope',j,'.Rdata'))
}

load(paste0(curr_date,'-Age-Sex-clock-Study2-pred-rt-lastoutcome-withrandslope5.Rdata'))
ddq$coef_df_reml %>% filter(effect=="fixed" & p.value < 0.05)

# by sex and last outcome: won't work for models without outcome
tab <- as.data.frame(ddq$emtrends_list$RTxO) %>% select(!c(model_name,rhs,outcome)) %>%
  group_by(sex1,reward_lag_rec) %>% summarize(rt_lag_sc.trend = mean(rt_lag_sc.trend), SE = mean(std.error))
tab$last_outcome <- ifelse(tab$reward_lag_rec==-0.5, "Omission","Reward")

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=sex1),group=dataset)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex1)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=sex1),position=position_dodge(width=0.5)) +
  facet_wrap(~last_outcome) + scale_y_reverse()

g <- ggplot(tab,aes(y=rt_lag_sc.trend,x=last_outcome),group=sex1)
g + geom_point(stat="identity",position=position_dodge(width=0.5),aes(color=sex1)) +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE,ymax = rt_lag_sc.trend + SE,color=sex1),position=position_dodge(width=0.5)) +
  facet_wrap(~sex1) + 
  ylab("RT Swing (RT ~ previous RT)") + xlab("Prior Trial Outcome")
