# 2025-01-17 AndyP
# M1 and V1 timing analysis

library(tidyverse)
library(lmerTest)
library(pracma)
library(fmri.pipeline)
library(lubridate)

### M1 ###

#m1 <- read_csv('/Volumes/Users/Andrew/explore_v1_m1/choice_aligned_explore_clock_Schaefer2018_200Parcel_7Network_V1_M1.csv.gz')
#m1 <- m1 %>% select(!decon_median & !decon_sd)

#source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
#df <- get_trial_data(repo_directory='~/clock_analysis',dataset='explore')
#df <- df %>% select(id,run,trial,iti_prev,rt_csv)

#df$id <- as.integer(df$id)
#m1 <- m1 %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
#m1 <- inner_join(df,m1,by=c('id','run','trial'))
#m1 <- m1 %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                                  # trial > 40 & trial <= 80 ~ trial-40,
                                                  # trial > 80 & trial <=120 ~ trial-80, 
                                                  # trial > 120 & trial <=160 ~ trial-120,
                                                  # trial > 160 & trial <=200 ~ trial-160,
                                                  # trial > 200 & trial <=240 ~ trial-200))
# m1 <- m1 %>% select(!trial)
# 
# scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
# scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = as.vector(t(scale(difftime(scaninfo$scan_date,min(scaninfo$scan_date)))))) %>% select(!registration_redcapid)
# 
# m1 <- inner_join(m1,scaninfo,by='id')
# m1$scan_which <- relevel(as.factor(m1$scan_which),ref='P1')
# m1 <- m1 %>% mutate(experiment = 'Explore')
# m1$scan_which[is.na(m1$scan_which)] = 'P1'
# 
# m1 <- m1 %>% group_by(id) %>% mutate(id_run = paste0(id,'_',run)) %>% ungroup()
# 
# ncores = 26
# decode_formula = formula(~1|(id_run))
# splits = c('evt_time','atlas_value')
# print(decode_formula)
# ddf <- mixed_by(m1, outcomes = "decon_mean", rhs_model_formulae = decode_formula , split_on = splits,
#                 padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                 tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
# 
# ddq <- ddf$coef_df_reml %>% filter(effect=='ran_coefs')
# ddq0 <- ddq %>% filter(atlas_value == 26) # left motor cortex
# m1_tau <- ddq0 %>% group_by(level) %>% arrange(evt_time) %>% summarize(peak = max(estimate,na.rm=TRUE),ipeak = evt_time[which(estimate==peak)]) %>% ungroup()
# hist(m1_tau$ipeak,unique(ddq0$evt_time))
# 
# rm(m1)
# gc()

### V1 ###


v1 <- read_csv('/Volumes/Users/Andrew/explore_v1_m1/clock_aligned_explore_clock_Schaefer2018_200Parcel_7Network_V1_M1.csv.gz')
v1 <- v1 %>% select(!decon_median & !decon_sd)

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory='~/clock_analysis',dataset='explore')
df <- df %>% select(id,run,trial,iti_prev,rt_csv)

df$id <- as.integer(df$id)
v1 <- v1 %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
v1 <- inner_join(df,v1,by=c('id','run','trial'))
v1 <- v1 %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                          trial > 40 & trial <= 80 ~ trial-40,
                                          trial > 80 & trial <=120 ~ trial-80, 
                                          trial > 120 & trial <=160 ~ trial-120,
                                          trial > 160 & trial <=200 ~ trial-160,
                                          trial > 200 & trial <=240 ~ trial-200))
v1 <- v1 %>% select(!trial)

scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = as.vector(t(scale(difftime(scaninfo$scan_date,min(scaninfo$scan_date)))))) %>% select(!registration_redcapid)

v1 <- inner_join(v1,scaninfo,by='id')
v1$scan_which <- relevel(as.factor(v1$scan_which),ref='P1')
v1 <- v1 %>% mutate(experiment = 'Explore')
v1$scan_which[is.na(v1$scan_which)] = 'P1'

v1 <- v1 %>% group_by(id) %>% mutate(id_run = paste0(id,'_',run)) %>% ungroup()

ncores = 26
decode_formula = formula(~1|(id_run))
splits = c('evt_time','atlas_value')
print(decode_formula)
ddf <- mixed_by(v1, outcomes = "decon_mean", rhs_model_formulae = decode_formula , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))

ddv <- ddf$coef_df_reml %>% filter(effect=='ran_coefs')
ddv0 <- ddv %>% filter(atlas_value == 106 | atlas_value == 108 | atlas_value == 7 | atlas_value == 5) # bilateral primary visual cortex
v1_tau <- ddv0 %>% group_by(level,atlas_value) %>% arrange(evt_time) %>% summarize(peak = max(estimate,na.rm=TRUE),ipeak = evt_time[which(estimate==peak)]) %>% ungroup()
hist(v1_tau$ipeak,unique(ddq0$evt_time))

##################
### add BSOC ####
#################

V1 <- c('5','206')
# load BSOC
clock_bsoc <- read.csv('/Volumes/Users/Andrew/DAN/fixed DAN/clock_aligned_bsocial_dan.csv.gz')
clock_bsoc <- clock_bsoc %>% select(id,run,trial,evt_time,decon_mean, atlas_value)
clock_bsoc <- clock_bsoc %>% filter(atlas_value %in% V1) %>% mutate(region = 'V1')

df <- read_csv('/Volumes/Users/Andrew/DAN/bsocial_clock_trial_df.csv')
df <- df %>% select(id,scanner_run,trial,iti_prev,rt_csv)
df <- df %>% rename(run=scanner_run)

df$id <- as.integer(df$id)
clock_bsoc <- clock_bsoc %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
clock_bsoc <- inner_join(df,clock_bsoc,by=c('id','run','trial'))
clock_bsoc <- clock_bsoc %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                                    trial > 40 & trial <= 80 ~ trial-40,
                                                    trial > 80 & trial <=120 ~ trial-80, 
                                                    trial > 120 & trial <=160 ~ trial-120,
                                                    trial > 160 & trial <=200 ~ trial-160,
                                                    trial > 200 & trial <=240 ~ trial-200))


#waiting on BSOC clock scan info
#scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
#scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = as.vector(t(scale(difftime(scaninfo$scan_date,min(scaninfo$scan_date)))))) %>% select(!registration_redcapid)

#rt_bsoc <- inner_join(rt_bsoc,scaninfo,by='id')
#rt_bsoc$scan_which <- relevel(as.factor(rt_bsoc$scan_which),ref='P1')
clock_bsoc <- clock_bsoc %>% mutate(experiment = 'bsocial')
clock_bsoc <- clock_bsoc %>% group_by(id) %>% mutate(id_run = paste0(id,'_',run)) %>% ungroup()


ncores = 26
decode_formula = formula(~1|(id_run))
splits = c('evt_time','atlas_value')
print(decode_formula)
ddf1 <- mixed_by(clock_bsoc, outcomes = "decon_mean", rhs_model_formulae = decode_formula , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))

ddv1 <- ddf1$coef_df_reml %>% filter(effect=='ran_coefs')
ddv1 <- ddv1 %>% filter(atlas_value == 106 | atlas_value == 108 | atlas_value == 7 | atlas_value == 5) # bilateral primary visual cortex
v1_tau2 <- ddv1 %>% group_by(level,atlas_value) %>% arrange(evt_time) %>% summarize(peak = max(estimate,na.rm=TRUE),ipeak = evt_time[which(estimate==peak)]) %>% ungroup()
v1_tau2 <- v1_tau2 %>% rename(id = level)

id0 <- NULL
run0 <- NULL
nT <- nrow(v1_tau2)
for (iT in 1:nT){
  temp <- v1_tau2$id[iT]
  id_temp <- str_split(temp,'_')[[1]][1]
  run_temp <- str_split(temp,'_')[[1]][2]
  id0[iT] <- id_temp
  run0[iT] <- run_temp
}

v1_tau2 <- v1_tau2 %>% mutate(id0=id0,run = run0)
v1_tau2 <- v1_tau2 %>% select(!id) %>% rename(id = id0)



hist(v1_tau2$ipeak,breaks=15)



p1_timing <- read_csv('/Volumes/Users/Andrew/2025-01-24-P1-Clock-Timing.csv')
p1_timing <- p1_timing[-c(90,92),]
p1_timing <- p1_timing %>% mutate(scanner = 'P1')
p2_timing <- read_csv('/Volumes/Users/Andrew/2025-01-22-P2-Clock-Timing.csv')
p2_timing <- p2_timing %>% mutate(scanner = 'P2')
p2_timing <- p2_timing[-c(91,93),]
merged_timing <- rbind(p1_timing,p2_timing)
merged_timing <- merged_timing %>% select(!reldifftiming & !difftiming)
merged_timing <- merged_timing %>% pivot_wider(names_from = 'RunNb', values_from = c('DicomAcquisitionTime','TaskStartTime'))


merged_timing <- merged_timing %>% mutate(at1 = lubridate::dmy_hms(merged_timing$DicomAcquisitionTime_1),at2 = lubridate::dmy_hms(merged_timing$DicomAcquisitionTime_2),ts1 = lubridate::dmy_hms(merged_timing$TaskStartTime_1), ts2 = lubridate::dmy_hms(merged_timing$TaskStartTime_2))
merged_timing <- merged_timing %>% mutate(diff1 = difftime(at1,ts1,units="secs"), diff2 = difftime(at2,ts2,units="secs"))

merged_timing$TR_diff = (merged_timing$diff1 - merged_timing$diff2)/0.6
merged_timing$dod = (merged_timing$diff1 - merged_timing$diff2)
merged_timing <- merged_timing %>% filter(TR_diff < 30 & TR_diff > -30)
merged_timing <- merged_timing %>% rename(id=ParticipantId)

v1_tau <- v1_tau %>% filter(atlas_value == 108 | atlas_value == 5) %>% rename(id = level)
merged_timing$id <- as.character(merged_timing$id)

id0 <- NULL
run0 <- NULL
nT <- nrow(v1_tau)
for (iT in 1:nT){
  temp <- v1_tau$id[iT]
  id_temp <- str_split(temp,'_')[[1]][1]
  run_temp <- str_split(temp,'_')[[1]][2]
  id0[iT] <- id_temp
  run0[iT] <- run_temp
}

v1_tau <- v1_tau %>% mutate(id0=id0,run = run0)
v1_tau <- v1_tau %>% select(!id) %>% rename(id = id0)

v1_tau_merged <- rbind(v1_tau,v1_tau2)

md <- inner_join(v1_tau_merged,merged_timing,by='id')
md <- md %>% mutate(scan_date = date(lubridate::dmy_hms(DicomAcquisitionTime_1)))
md <- md %>% mutate(ddate = as.vector(t(scale(difftime(md$scan_date,min(md$scan_date))))))
