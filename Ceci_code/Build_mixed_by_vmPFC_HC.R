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

# set directories - change for your needs
rootdir <- '/ix/cladouceur/DNPL'
path_to_trial_data <- rootdir
path_to_mixed_by <- rootdir
path_to_data <- rootdir

# load some source code that we'll need to use
setwd(rootdir)
source(paste0(path_to_trial_data,'/get_trial_data.R'))

# load mixed_by function for analyses
source(paste0(path_to_trial_data,"/mixed_by.R"))
source(paste0(rootdir,"/plot_mixed_by_vmPFC_HC.R"))
source(paste0(rootdir,'/plot_mixed_by_vmPFC_HC_simplified_networksymmetry.R'))
##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

# load the vmPFC data, filter to within 5s of RT, and select vars of interest
setwd(path_to_data)
load('MMclock_clock_Aug2023.Rdata')
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)

# load in the hippocampus data, filter to within 5s of RT
load('HC_clock_Aug2023.Rdata')
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
Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
Q <- Q %>% select(!decon1)

# Get task behav data
# (had to download it from UNCDEPENdlab github first)
df <- get_trial_data(repo_directory='/ix/cladouceur/DNPL',dataset='mmclock_fmri')

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
                       rt_lag_sc,v_entropy_sc,rt_swing_sc,trial_neg_inv_sc)
Q <- inner_join(behav, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")

# add in age and sex variables
demo <- read.table(file=file.path('/ix/cladouceur/DNPL', 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)

########################
########################
########################
##### brain models  ####
########################
########################
########################

########################
##### Set up models ####
########################

########################
# by network
########################

setwd(paste0(rootdir,'/fmri/mixed_by_output'))

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# age and hippocampal activity only
decode_formula[[1]] = formula(~age*HCwithin + HCbetween + (1 | id/run)) 

# age, hippocampal activity and control vars (run, trial within run)
decode_formula[[2]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

# adding sex
decode_formula[[3]] = formula(~age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                emmeans_spec = list(Age = list(outcome='vmPFC_decon', model_name='model1', 
                           specs=formula(~age:HCwithin), at = list(age=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))),
                )
  # write to file
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  #save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-',i,'.Rdata'))

}
#####################
##### plot models ####
#####################

# load in data already run
setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
load('2024-11-17-vmPFC-HC-network-clock-3.Rdata')

# plot individual terms to console - using my adapted code
# get terms
unique(ddf$coef_df_reml$term)
plot_mixed_by_vmPFC_HC_simplified_networksymmetry(ddf,totest = "Ceci-models-",behavmodel = "MMClock",toalign='clock',model_iter=4,
                                          term = "HCwithin:female1",toprocess="network-by-HC")

# Save all plots as pdf - Using Andrew's plot_mixed_by function
plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'MMClock',totest='Ceci-models-',toalign='clock',
                       toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = 4)

########################
# by symmetry group
########################

setwd('/ix/cladouceur/DNPL/fmri/mixed_by_output')

# split out by network
splits = c('evt_time','symmetry_group','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# ag ne and hippocampal activity only
decode_formula[[1]] = formula(~age*HCwithin + HCbetween + (1 | id/run)) 

# age, hippocampal activity and control vars (run, trial within run)
decode_formula[[2]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

# adding sex
decode_formula[[3]] = formula(~age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                  emmeans_spec = list(Age = list(outcome='vmPFC_decon', model_name='model1', 
                                                 specs=formula(~age:HCwithin), at = list(age=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))),
  )
  # write to file
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-clock-',i,'.Rdata'))
  
}
#####################
##### plot models ####
#####################
# load in data already run
setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
load('2024-11-17-vmPFC-HC-symmetry-clock-3.Rdata')

# plot individual terms to console - using my adapted code
# get terms
unique(ddf$coef_df_reml$term)
plot_mixed_by_vmPFC_HC_simplified_networksymmetry(ddf,totest = "Ceci-models-",behavmodel = "MMClock",toalign='clock',model_iter=4,
                                          term = "HCwithin:female1",toprocess="symmetry-by-HC")

# Save all plots as pdf - Using Andrew's plot_mixed_by function
plot_mixed_by_vmPFC_HC(ddf=ddf,behavmodel = 'MMClock',totest='Ceci-models-',toalign='clock',
                       toprocess='symmetry-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = 4)

#################################
##### split models by gender ####
#################################

Qf <- Q[Q$female==1,]
Qm <- Q[Q$female==0,]

########################
##### Set up models ####
########################

setwd('/ix/cladouceur/DNPL/fmri/mixed_by_output')

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# age and hippocampal activity only
decode_formula[[1]] = formula(~age*HCwithin + HCbetween + (1 | id/run)) 

# age, hippocampal activity and control vars (run, trial within run)
decode_formula[[2]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  for(gender in c("girls", "boys")){
    print(gender)
    if(gender == "girls"){
      assign("model_to_run",Qf)
      ddf_girls <- mixed_by(model_to_run, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                      emmeans_spec = list(Age = list(outcome='vmPFC_decon', model_name='model1', 
                                                     specs=formula(~age:HCwithin), at = list(age=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))),
      )
      # write to file
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf_girls,file=paste0(curr_date,'-vmPFC-HC-network-clock-',gender,'-',i,'.Rdata'))
    }
    else if(gender == "boys"){
      assign("model_to_run",Qm)
      ddf_boys <- mixed_by(model_to_run, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                      emmeans_spec = list(Age = list(outcome='vmPFC_decon', model_name='model1', 
                                                     specs=formula(~age:HCwithin), at = list(age=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))),
      )
      # write to file
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf_boys,file=paste0(curr_date,'-vmPFC-HC-network-clock-',gender,'-',i,'.Rdata'))

    }
    
  }
  
}

#####################
##### plot models ####
#####################

# load in data already run
setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
load('2024-11-18-vmPFC-HC-network-clock-girls-2.Rdata')

load('2024-11-18-vmPFC-HC-network-clock-boys-2.Rdata')

# plot individual terms to console - using my adapted code
# get terms
unique(ddf_girls$coef_df_reml$term)
plot_mixed_by_vmPFC_HC_simplified_networksymmetry(ddf_girls,totest = "Ceci-models-",behavmodel = "MMClock",toalign='clock',model_iter=4,
                                                  term = "age:HCwithin",toprocess="network-by-HC")

# get terms
unique(ddf_boys$coef_df_reml$term)
plot_mixed_by_vmPFC_HC_simplified_networksymmetry(ddf_boys,totest = "Ceci-models-",behavmodel = "MMClock",toalign='clock',model_iter=4,
                                                  term = "age:HCwithin",toprocess="network-by-HC")

# Save to file using Andrew's plot_mixed_by function

#girls
plot_mixed_by_vmPFC_HC(ddf=ddf_girls,behavmodel = 'MMClock',totest='Ceci-models-girls-',toalign='clock',
                       toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = 4)

#boys
plot_mixed_by_vmPFC_HC(ddf=ddf_boys,behavmodel = 'MMClock',totest='Ceci-models-boys',toalign='clock',
                       toprocess='network-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = 4)


########################
# by symmetry group
########################

setwd('/ix/cladouceur/DNPL/fmri/mixed_by_output')

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','symmetry_group','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# ag ne and hippocampal activity only
decode_formula[[1]] = formula(~age*HCwithin + HCbetween + (1 | id/run)) 

# age, hippocampal activity and control vars (run, trial within run)
decode_formula[[2]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 | id/run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  for(gender in c("girls", "boys")){
    print(gender)
    if(gender == "girls"){
      assign("model_to_run",Qf)
      ddf_girls <- mixed_by(model_to_run, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                            padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                            tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                            emmeans_spec = list(Age = list(outcome='vmPFC_decon', model_name='model1', 
                                                           specs=formula(~age:HCwithin), at = list(age=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))),
      )
      # write to file
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf_girls,file=paste0(curr_date,'-vmPFC-HC-symmetry-clock-',gender,'-',i,'.Rdata'))
    }
    else if(gender == "boys"){
      assign("model_to_run",Qm)
      ddf_boys <- mixed_by(model_to_run, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                           padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                           tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE),
                           emmeans_spec = list(Age = list(outcome='vmPFC_decon', model_name='model1', 
                                                          specs=formula(~age:HCwithin), at = list(age=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))),
      )
      # write to file
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf_boys,file=paste0(curr_date,'-vmPFC-HC-symmetry-clock-',gender,'-',i,'.Rdata'))
    }
    
  }
  
}
#####################
##### plot models ####
#####################

# load in data already run
setwd(paste0(rootdir,'/fmri/mixed_by_output/'))
load('2024-11-18-vmPFC-HC-symmetry-clock-girls-2.Rdata')

load('2024-11-18-vmPFC-HC-symmetry-clock-boys-2.Rdata')

# plot individual terms to console - using my adapted code
# get terms
unique(ddf_girls$coef_df_reml$term)
plot_mixed_by_vmPFC_HC_simplified_networksymmetry(ddf_girls,totest = "Ceci-models-",behavmodel = "MMClock",toalign='clock',model_iter=4,
                                                  term = "age:HCwithin",toprocess="symmetry-by-HC")

# get terms
unique(ddf_boys$coef_df_reml$term)
plot_mixed_by_vmPFC_HC_simplified_networksymmetry(ddf_boys,totest = "Ceci-models-",behavmodel = "MMClock",toalign='clock',model_iter=4,
                                                  term = "age:HCwithin",toprocess="network-by-HC")

# Save to file using Andrew's plot_mixed_by function

#girls
plot_mixed_by_vmPFC_HC(ddf=ddf_girls,behavmodel = 'MMClock',totest='Ceci-models-girls-',toalign='clock',
                       toprocess='symmetry-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = 4)

#boys
plot_mixed_by_vmPFC_HC(ddf=ddf_boys,behavmodel = 'MMClock',totest='Ceci-models-boys',toalign='clock',
                       toprocess='symmetry-by-HC',CTRflag = "FALSE",hc_LorR = 'LR',flipy = 'FALSE',model_iter = 4)

#############################
#############################
#############################
##### brain-to-behavior  ####
#############################
#############################
#############################

########################
##### Set up models ####
########################

########################
# by network
########################

setwd(paste0(rootdir,'/fmri/mixed_by_output'))
#load('2024-11-25-vmPFC-HC-network-clock-withRandslope3.Rdata')
source(paste0(rootdir,'/get_random_slope.R'))

# set some baseline variables
ncores = 20

# split out by network
splits = c('evt_time','network','HC_region')

# store models in a data-frame so they can be bulk run later.
rm(decode_formula)
decode_formula <- NULL

# age and hippocampal activity only
decode_formula[[1]] = formula(~age*HCwithin + HCbetween + (1 + age*HCwithin|id) + (1|run)) 

# age, hippocampal activity and control vars (run, trial within run)
decode_formula[[2]] = formula(~age*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 + age*HCwithin|id) + (1|run)) 

# adding sex
decode_formula[[3]] = formula(~age*HCwithin + female*HCwithin + trial_neg_inv_sc*HCwithin + rt_lag_sc*HCwithin + HCbetween + (1 + age*HCwithin|id) + (1|run)) 

for(i in 1:length(decode_formula)){
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),return_models=TRUE,conf.int=TRUE))
  
  # save for B2B modeling
  assign(paste0("model",i),ddf)
  # write to file
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-withRandslope',i,'.Rdata'))
  
}

# function to get the random slope of the age*HCwithin interaction and combine with dataframe
# sticking with the model without gender for now
Q2 <- get_random_slope(model2,'age:HCwithin')

# to run the B2B model
# define some stuff
splits = c('network','HC_region')
toalign="feedback"
totest <- paste0('rt_csv_sc-',as.character(i))
toprocess <- 'network-by-HC'
behavmodel <- 'compressed'
hc_LorR <- 'LR'

# set up the models
rm(decode_formula)
decode_formula <- NULL
decode_formula[[1]] <- formula(~trial_neg_inv_sc + last_outcome*subj_level_rand_slope + rt_lag_sc*subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope + (1 | id/run))
decode_formula[[2]] <- formula(~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope * trial_neg_inv_sc + (1 | id/run))

# adding random slope for last RT 
decode_formula[[3]] <- formula(~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope * trial_neg_inv_sc + (1 + rt_lag_sc | id/run))

# adding random slope for last RTvmax
decode_formula[[4]] <- formula(~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:subj_level_rand_slope + rt_vmax_lag_sc * subj_level_rand_slope * trial_neg_inv_sc + (1 + rt_vmax_lag_sc | id/run))


# run the models in mixed_by
for (i in 1:length(decode_formula)){
    ddq <- mixed_by(Q2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[i]], split_on = splits,return_models=TRUE,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                    emmeans_spec = list(
                      RT = list(outcome='rt_csv_sc', model_name='model1', 
                                specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                      Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                                  specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                      RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                  specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),        
                      TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                                     specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                      
                    ),
                    emtrends_spec = list(
                      RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                                specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                      Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                  specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                      RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                  specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                      TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
                                     specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                      TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
                                      specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                      TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'subj_level_rand_slope',
                                      specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                    )
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-',toalign,'-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-',i,'.Rdata'))

}


# plotting
library(RColorBrewer)
library(sciplot)

# function for plotting
plot_rand_slopes <- function(data, x, y, facets, xlab, ylab, plotlab){
  g <- ggplot(data, aes(x=x, y=y))
  g + geom_point(aes(color=network)) + 
    geom_errorbar(aes(color=network, ymin = y - std.error, ymax = y + std.error)) + 
    geom_line(aes(color=network)) +
    scale_fill_discrete(name = "Network") + 
    facet_wrap(facets,scales = "free_y") +
    labs(title = plotlab) +
    xlab(xlab) + ylab(ylab)
}

# RT swings
# without rand slope for last RT--significant
load('2024-11-26-vmPFC-HC-network-ranslopes-feedback-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddq$coef_df_reml %>% filter(effect=="fixed" & term=="subj_level_rand_slope:rt_lag_sc" & network=="DMN")

# plotting
toplot <- ddq$emtrends_list$RTxO
#toplot %>% group_by(network,last_outcome,subj_level_rand_slope) %>% summarise(mean = mean(rt_lag_sc.trend),se = se(rt_lag_sc.trend))

plot_rand_slopes(toplot,toplot$subj_level_rand_slope,toplot$rt_lag_sc.trend,formula(~last_outcome+HC_region),
                 "Random slope for age:HCwithin","RT lag scaled","Effect of random slopes on RT")

# including rand slope for last RT--doesn't survive
load('2024-11-26-vmPFC-HC-network-ranslopes-feedback-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-3.Rdata')
ddq$coef_df_reml %>% filter(effect=="fixed" & term=="subj_level_rand_slope:rt_lag_sc" & network=="DMN")

# plotting
toplot <- ddq$emtrends_list$RTxO
#toplot %>% group_by(network,last_outcome,subj_level_rand_slope) %>% summarise(mean = mean(rt_lag_sc.trend),se = se(rt_lag_sc.trend))

plot_rand_slopes(toplot,toplot$subj_level_rand_slope,toplot$rt_lag_sc.trend,formula(~last_outcome+HC_region),
                 "Random slope for age:HCwithin","RT lag scaled","Effect of random slopes on RT")

#RT controlling for RTvmax
# w/o rand slope for last RTvmax--not significant
load('2024-11-26-vmPFC-HC-network-ranslopes-feedback-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddq$coef_df_reml %>% filter(effect=="fixed" & term=="subj_level_rand_slope:rt_vmax_lag_sc" & network=="DMN")

toplot <- ddq$emtrends_list$Vmax
plot_rand_slopes(toplot,toplot$subj_level_rand_slope,toplot$rt_vmax_lag_sc.trend,formula(~HC_region),
                 "Random slope for age:HCwithin","RT Vmax lag scaled","Effect of random slopes on RT")

# w rand slope for last RT -- significant
load('2024-11-26-vmPFC-HC-network-ranslopes-feedback-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-3.Rdata')
ddq$coef_df_reml %>% filter(effect=="fixed" & term=="subj_level_rand_slope:rt_vmax_lag_sc" & network=="DMN")

toplot <- ddq$emtrends_list$Vmax
plot_rand_slopes(toplot,toplot$subj_level_rand_slope,toplot$rt_vmax_lag_sc.trend,formula(~HC_region),
                 "Random slope for age:HCwithin","RT Vmax lag scaled","Effect of random slopes on RT")

# w rand slope for last RTvmax -- doesn't survive
load('2024-11-26-vmPFC-HC-network-ranslopes-feedback-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-4.Rdata')
ddq$coef_df_reml %>% filter(effect=="fixed" & term=="subj_level_rand_slope:rt_vmax_lag_sc" & network=="DMN")

toplot <- ddq$emtrends_list$Vmax
plot_rand_slopes(toplot,toplot$subj_level_rand_slope,toplot$rt_vmax_lag_sc.trend,formula(~HC_region),
                 "Random slope for age:HCwithin","RT Vmax lag scaled","Effect of random slopes on RT")


