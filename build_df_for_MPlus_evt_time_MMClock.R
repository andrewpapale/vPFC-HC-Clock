

# libraries we'll need
library(tidyverse)
library(fmri.pipeline)
library(MplusAutomation)
# set root directory
rootdir1 <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMClock_BSOC_MPlus'
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
repo_directory <- file.path('/Volumes/Users/Andrew/MEDuSA_data_BSOC')
repo_directory1 <- '~/clock_analysis'
ncores <- 26
# load mixed_by function for analyses

##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

load("/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2025-03-05-vmPFC-HC-network-clock-ranslopes-nofixedeffect-noHCbetween-3.Rdata")
Q <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')
Q <- Q %>% rename(id=level)
Q$id <- as.character(Q$id)

Q1ah <- Q %>% filter(HC_region == "AH")
Q1ph <- Q %>% filter(HC_region == "PH")

Qph <- Q1ph %>% select(id,network,evt_time,estimate) %>% group_by(id) %>% 
  pivot_wider(values_from=c(estimate),names_from = c('network','evt_time')) %>% 
  ungroup()  

Qah <- Q1ah %>% select(id,network,evt_time,estimate) %>% group_by(id) %>% 
  pivot_wider(values_from=c(estimate),names_from = c('network','evt_time')) %>% 
  ungroup()  

# Get task behav data
#df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
# Get task behav data
source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory1,dataset='mmclock_fmri')

# select only vars of interest and merge into MRI data
behav <- df %>% select(id,trial,run_trial,
                       rt_lag_sc,rt_vmax_lag_sc,trial_neg_inv_sc,last_outcome,outcome,
                       v_entropy_wi,v_max_wi,rt_csv_sc,rt_csv,iti_ideal,rewFunc,outcome)
behav$id <- as.character(behav$id)

Qah <- inner_join(behav, Qah, by = c("id")) %>% arrange("id","run","trial")
Qph <- inner_join(behav, Qph, by = c("id")) %>% arrange("id","run","trial")

# add in age and sex variables
demo <- read.table(file=file.path(repo_directory1, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)

Qah <- inner_join(Qah,demo,by=c('id'))
Qph <- inner_join(Qph,demo,by=c('id'))

Qah$age <- scale(Qah$age)
Qph$age <- scale(Qph$age)

setwd(file.path(rootdir1))

save(Qah,file=file.path(rootdir1,'MMClock_HC_vmPFC_clock_All_evt_time_AHforMplus.Rdata'))
save(Qph,file=file.path(rootdir1,'MMClock_HC_vmPFC_clock_All_evt_time_PHforMplus.Rdata'))


prepareMplusData(df = Qah, filename = "MMClock_HC_vmPFC_clock_All_evt_time_AH_forMplus_taa.dat", dummyCode = c("outcome", "female"), overwrite = TRUE)
prepareMplusData(df = Qph, filename = "MMClock_HC_vmPFC_clock_All_evt_time_PH_forMplus_taa.dat", dummyCode = c("outcome", "female"), overwrite = TRUE)






