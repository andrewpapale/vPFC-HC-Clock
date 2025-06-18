

# libraries we'll need
library(tidyverse)
library(fmri.pipeline)
library(MplusAutomation)
# set root directory
rootdir1 <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMClock_BSOC_MPlus'
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
repo_directory <- file.path('/Volumes/Users/Andrew/MEDuSA_data_BSOC')
ncores <- 26
# load mixed_by function for analyses

##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2025-04-09-Bsocial-vPFC-HC-network-clock-RTcorrected-ranslopes-1.Rdata')
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
df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
split_ksoc_bsoc <- df %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300])
bsoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240])
bsoc <- rbind(bsoc,221973,221507,221842,440223)
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

Qah <- inner_join(behav, Qah, by = c("id")) %>% arrange("id","run","trial")
Qph <- inner_join(behav, Qph, by = c("id")) %>% arrange("id","run","trial")

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

Qah <- inner_join(Qah,demo2,by=c('id'))
Qph <- inner_join(Qph,demo2,by=c('id'))

Qah$female <- ifelse(Qah$sex==1,1,0)
Qah <- Qah %>% select(!sex)
Qah$age <- scale(Qah$age)

Qph$female <- ifelse(Qph$sex==1,1,0)
Qph <- Qph %>% select(!sex)
Qph$age <- scale(Qph$age)

setwd(file.path(rootdir1))

save(Qah,file=file.path(rootdir,'bsocial_HC_vmPFC_clock_All_evt_time_AHforMplus.Rdata'))
save(Qph,file=file.path(rootdir,'bsocial_HC_vmPFC_clock_All_evt_time_PHforMplus.Rdata'))


prepareMplusData(df = Qah, filename = "bsocial_HC_vmPFC_clock_All_evt_time_AH_forMplus_taa.dat", dummyCode = c("outcome", "female"), overwrite = TRUE)
prepareMplusData(df = Qph, filename = "bsocial_HC_vmPFC_clock_All_evt_time_PH_forMplus_taa.dat", dummyCode = c("outcome", "female"), overwrite = TRUE)






