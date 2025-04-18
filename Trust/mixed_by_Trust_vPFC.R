# 2025-03-03 AndyP
# Trust MEDuSA analysis


library(tidyverse)
library(fmri.pipeline)
library(pracma)

ncores <- 26


vmPFC <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/outcome_aligned_bsocial_vmPFC.csv.gz')
vmPFC <- vmPFC %>% filter(evt_time > -2 & evt_time < 5.4) # need to censor to the highest TR with actual data in it or will get an error in mixed_by
vmPFC <- vmPFC %>% mutate(network = case_when(
  atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
  atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
  atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM'
)) %>% mutate(symmetry_group=case_when(
  atlas_value==67 | atlas_value==171 ~ 6,
  atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 1,
  atlas_value==89 | atlas_value==194 ~ 7,
  atlas_value==88 | atlas_value==192 ~ 5,
  atlas_value==84 | atlas_value==191 ~ 4,
  atlas_value==86 | atlas_value==161 ~ 2,
  atlas_value==55 | atlas_value==160 ~ 8,
  atlas_value==56 | atlas_value==159 ~ 3
))
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean) %>% select(!run)
vmPFC <- vmPFC %>% arrange(id,trial,evt_time)

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Trust/Trust_HC_outcome_TRdiv2.Rdata')
hc <- hc %>% filter(evt_time > -2 & evt_time < 5.4) # need to censor to the highest TR with actual data in it or will get an error in mixed_by
hc <- hc %>% select(!atlas_value & !run)
hc <- hc %>% group_by(id,trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup()
hc <- hc %>% rename(decon_mean=decon1)
hc <- hc %>% group_by(id) %>% mutate(HCwithin = scale(decon_mean),HCbetween=mean(decon_mean,na.rm=TRUE)) %>% ungroup()

df <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/Tim_trust_trialdf_clean.csv')
# there are three subjects for whom the task blew up on the last trial, and they were assigned weird, negative timings. One was already dropped from the trial_df due to a lost nifit (221548). Dropping the last trial for the other 2 subjects manually:
df <- df %>% filter(!(id == 440230 & trial ==192) & !(id == 440124 & trial == 144)) %>% filter(!id==221548)
df <- df %>% select(id,reward,t_decides,p_decides,alpha_transformed,beta_transformed,decision,outcome_duration,pchoice_duration, dchoice_duration,outcome_offset,t_decides,trial,pchoice_onset,trustee,block,pchoice_rt,iti_onset,noresponse,kappaS_transformed,kappaT_bad_transformed,kappaT_good_transformed,kappaT_computer_transformed,PE_mult_z,V_mult_z,reward)
df <- df %>% group_by(id) %>% mutate(iti_duration = lead(iti_onset) - outcome_offset) %>% ungroup()
df <- df %>% arrange(id,trial) %>% group_by(id,trustee,block) %>% mutate(block_trial = seq(1,length(unique(trial)),length.out=length(unique(trial)))) %>% ungroup()
df$trustee <- relevel(as.factor(df$trustee),ref='neutral')
#wrangle trial_type variable
df <- df %>% mutate(p_dec_char = if_else(p_decides == 0, "keep", 
                                                     if_else(p_decides == 1, "share", NA_character_)),
                                trial_type = if_else(!is.na(p_dec_char), 
                                                     paste("i", p_dec_char, "t", t_decides, sep = "_"),
                                                     p_dec_char))
df$trial_type <- as.factor(df$trial_type)
df$trial_type <- relevel(df$trial_type, ref = "i_share_t_share")
df <- df %>% group_by(id,trustee) %>% mutate(tdec_prev = lag(t_decides), pdec_prev = lag(p_decides)) %>% ungroup()

Q <- full_join(vmPFC,hc,by=c('id','trial','evt_time'))
rm(vmPFC,hc)
gc()

Q <- inner_join(Q,df,by=c('id','trial'))
Q <- Q %>% group_by(id) %>% mutate(iti_duration = lead(iti_onset) - outcome_offset) %>% ungroup()

Q$vmPFC_decon[Q$evt_time >  + Q$outcome_duration + Q$iti_duration] = NA
Q$vmPFC_decon[Q$evt_time < -(Q$dchoice_duration + Q$pchoice_duration)] = NA
Q$HCwithin[Q$evt_time >  + Q$outcome_duration + Q$iti_duration] = NA
Q$HCwithin[Q$evt_time < -(Q$dchoice_duration + Q$pchoice_duration)] = NA
Q$HCbetween[Q$evt_time >  + Q$outcome_duration + Q$iti_duration] = NA
Q$HCbetween[Q$evt_time < -(Q$dchoice_duration + Q$pchoice_duration)] = NA

splits = c('network','evt_time','HC_region')
decode_formula <- NULL
decode_formula[[1]] = formula(~ reward*HCwithin + PE_mult_z*HCwithin + V_mult_z*HCwithin + scale(trial) + scale(block) + trustee*HCwithin + scale(iti_duration)*HCwithin + scale(pchoice_duration) + HCbetween + (1|id))
decode_formula[[2]] = formula(~ reward*HCwithin + PE_mult_z*HCwithin + scale(trial)*HCwithin  + scale(block)*HCwithin + trustee*HCwithin + scale(iti_duration)*HCwithin + scale(pchoice_duration) + HCbetween + (1|id))
decode_formula[[3]] = formula(~ reward*HCwithin + scale(alpha_transformed)*HCwithin + scale(beta_transformed)*HCwithin + scale(kappaS_transformed)*HCwithin + scale(kappaT_bad_transformed)*HCwithin + scale(kappaT_good_transformed)*HCwithin + scale(kappaT_computer_transformed)*HCwithin + scale(trial) + scale(block) + scale(iti_duration)*HCwithin + scale(pchoice_duration) + HCbetween + (1|id))

for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-Trust-vmPFC-HC-network-All-RTcorrected-',i,'.Rdata'))
}



vmPFC <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/outcome_aligned_bsocial_vmPFC.csv.gz')
vmPFC <- vmPFC %>% filter(evt_time > -2 & evt_time < 5.4) # need to censor to the highest TR with actual data in it or will get an error in mixed_by
vmPFC <- vmPFC %>% mutate(network = case_when(
  atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
  atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
  atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM'
)) %>% mutate(symmetry_group=case_when(
  atlas_value==67 | atlas_value==171 ~ 6,
  atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 1,
  atlas_value==89 | atlas_value==194 ~ 7,
  atlas_value==88 | atlas_value==192 ~ 5,
  atlas_value==84 | atlas_value==191 ~ 4,
  atlas_value==86 | atlas_value==161 ~ 2,
  atlas_value==55 | atlas_value==160 ~ 8,
  atlas_value==56 | atlas_value==159 ~ 3
))
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean) %>% select(!run)
vmPFC <- vmPFC %>% arrange(id,trial,evt_time)

df <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/Tim_trust_trialdf_clean.csv')
# there are three subjects for whom the task blew up on the last trial, and they were assigned weird, negative timings. One was already dropped from the trial_df due to a lost nifit (221548). Dropping the last trial for the other 2 subjects manually:
df <- df %>% filter(!(id == 440230 & trial ==192) & !(id == 440124 & trial == 144)) %>% filter(!id==221548)
df <- df %>% select(id,reward,alpha_transformed,beta_transformed,decision,outcome_duration,pchoice_duration, dchoice_duration,outcome_offset,t_decides,trial,pchoice_onset,trustee,block,pchoice_rt,iti_onset,noresponse,kappaS_transformed,kappaT_bad_transformed,kappaT_good_transformed,kappaT_computer_transformed,PE_mult_z,V_mult_z,reward)
df <- df %>% group_by(id) %>% mutate(iti_duration = lead(iti_onset) - outcome_offset) %>% ungroup()
df$trustee <- relevel(as.factor(df$trustee),ref='neutral')

Q <- inner_join(vmPFC,df,by=c('id','trial'))
Q <- Q %>% group_by(id) %>% mutate(iti_duration = lead(iti_onset) - outcome_offset) %>% ungroup()

Q$vmPFC_decon[Q$evt_time >  + Q$outcome_duration + Q$iti_duration] = NA
Q$vmPFC_decon[Q$evt_time < -(Q$dchoice_duration + Q$pchoice_duration)] = NA

splits = c('network','evt_time')
decode_formula <- NULL
decode_formula[[1]] = formula(~ reward + PE_mult_z + V_mult_z + scale(trial) + scale(block) + scale(iti_duration) + trustee + scale(pchoice_duration) + (1|id))
decode_formula[[2]] = formula(~ reward + scale(alpha_transformed) + scale(beta_transformed) + scale(kappaS_transformed) + scale(kappaT_bad_transformed) + scale(kappaT_good_transformed) + scale(kappaT_computer_transformed) + scale(trial) + scale(block) + scale(iti_duration) + scale(pchoice_duration) + (1|id))

for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-Trust-vmPFC-network-All-RTcorrected-',i,'.Rdata'))
}

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Trust/Trust_HC_outcome_TRdiv2.Rdata')
hc <- hc %>% filter(evt_time > -2 & evt_time < 5.4) # need to censor to the highest TR with actual data in it or will get an error in mixed_by
hc <- hc %>% select(!atlas_value & !run)
hc <- hc %>% group_by(id,trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup()
hc <- hc %>% rename(decon_mean=decon1)
hc <- hc %>% group_by(id) %>% mutate(HCwithin = scale(decon_mean),HCbetween=mean(decon_mean,na.rm=TRUE)) %>% ungroup()

df <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/Tim_trust_trialdf_clean.csv')
# there are three subjects for whom the task blew up on the last trial, and they were assigned weird, negative timings. One was already dropped from the trial_df due to a lost nifit (221548). Dropping the last trial for the other 2 subjects manually:
df <- df %>% filter(!(id == 440230 & trial ==192) & !(id == 440124 & trial == 144)) %>% filter(!id==221548)
df <- df %>% select(id,reward,alpha_transformed,beta_transformed,decision,outcome_duration,pchoice_duration, dchoice_duration,outcome_offset,t_decides,trial,pchoice_onset,trustee,block,pchoice_rt,iti_onset,noresponse,kappaS_transformed,kappaT_bad_transformed,kappaT_good_transformed,kappaT_computer_transformed,PE_mult_z,V_mult_z,reward)
df <- df %>% group_by(id) %>% mutate(iti_duration = lead(iti_onset) - outcome_offset) %>% ungroup()
df$trustee <- relevel(as.factor(df$trustee),ref='neutral')
df <- df %>% arrange(id,trial) %>% group_by(id,trustee,block) %>% mutate(block_trial = seq(1,length(unique(trial)),length.out=length(unique(trial)))) %>% ungroup()
Q <- inner_join(hc,df,by=c('id','trial'))
rm(hc)
gc()
Q <- Q %>% group_by(id) %>% mutate(iti_duration = lead(iti_onset) - outcome_offset) %>% ungroup()

Q$HCwithin[Q$evt_time >  + Q$outcome_duration + Q$iti_duration] = NA
Q$HCwithin[Q$evt_time < -(Q$dchoice_duration + Q$pchoice_duration)] = NA
Q$HCbetween[Q$evt_time >  + Q$outcome_duration + Q$iti_duration] = NA
Q$HCbetween[Q$evt_time < -(Q$dchoice_duration + Q$pchoice_duration)] = NA

splits = c('HC_region','evt_time')
decode_formula <- NULL
decode_formula[[1]] = formula(~ reward + PE_mult_z + V_mult_z + scale(trial) + scale(block) + trustee + scale(iti_duration) + scale(pchoice_duration) + HCbetween + (1|id))
decode_formula[[2]] = formula(~ reward + PE_mult_z + scale(trial)  + scale(block) + trustee + scale(iti_duration) + scale(pchoice_duration) + HCbetween + (1|id))
decode_formula[[3]] = formula(~ reward + scale(alpha_transformed) + scale(beta_transformed) + scale(kappaS_transformed) + scale(kappaT_bad_transformed) + scale(kappaT_good_transformed) + scale(kappaT_computer_transformed) + scale(trial) + scale(block) + scale(iti_duration) + scale(pchoice_duration) + HCbetween + (1|id))

for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  ddf <- mixed_by(Q, outcomes = "HCwithin", rhs_model_formulae = df0 , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-Trust-HC-network-All-RTcorrected-',i,'.Rdata'))
}
