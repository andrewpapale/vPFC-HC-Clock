# 2025-04-14 AndyP
# Plotting random slopes vs. age and AH/PH for devlopmental project

# libraries we'll need
library(tidyverse)
library(fmri.pipeline)
library(MplusAutomation)
# set root directory

repo_directory <- "~/clock_analysis"
ncores <- 26
# load mixed_by function for analyses

##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
model_str <- paste0('-vmPFC-HC-network-','clock','-ranslopes-nofixedeffect-noHCbetween-',3,'.Rdata')
model_str <- Sys.glob(paste0('*',model_str))
load(model_str)

Q <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')

Q <- Q %>% rename(id=level)
Q$id <- as.character(Q$id)

Qah <- Q %>% filter(HC_region == "AH")
Qph <- Q %>% filter(HC_region == "PH")

# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')

Qah <- inner_join(Qah,demo,by=c('id'))
Qph <- inner_join(Qph,demo,by=c('id'))

Qah <- Qah %>% mutate(sex = case_when(female == 1 ~ 'F',
                                      female == 0 ~ 'M'))

Qph <- Qph %>% mutate(sex = case_when(female == 1 ~ 'F',
                                      female == 0 ~ 'M'))

Qah$sex <- as.factor(Qah$sex)
Qph$sex <- as.factor(Qph$sex)

Q <- rbind(Qah,Qph)
Q <- Q %>% mutate(age_bin = case_when(age < 15 ~ '< 15',
                                      age >= 15 & age < 18 ~ '15-18',
                                      age >= 18 & age < 21 ~ '18-21',
                                      age >= 21 & age < 24 ~ '21-24',
                                      age >= 24 & age < 27 ~ '24-27',
                                      age >= 27 ~ '>= 27'))

Q <- Q %>% select(age_bin,HC_region,estimate,sex,network)

Q1 <- Q %>% group_by(age_bin,sex,HC_region,network) %>% summarize(estimate0 = mean(estimate,na.rm=TRUE),sd0 = sd(estimate,na.rm=TRUE), N=n()) %>% ungroup()


#ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex,group=sex)) + geom_line() + facet_wrap(~HC_region)

Q1$age_bin <- factor(Q1$age_bin, levels = c('< 15','15-18','18-21','21-24','24-27','>= 27'))

ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex,group=sex)) + 
  geom_line() + geom_errorbar(aes(ymin=estimate0-sd0/sqrt(N), ymax = estimate0+sd0/sqrt(N))) +
  facet_wrap(network~HC_region)

ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex, group=sex)) + geom_smooth()

# load(file.path(rootdir,'BSOC_HC_clock_TRdiv2.Rdata'))
# hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
# 
# split_ksoc_bsoc <- hc %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
# ksoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300])
# bsoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240])
# bsoc <- rbind(bsoc,221973,221507,221842,440223)
# hc_bsoc <- hc %>% filter(id %in% bsoc$id) %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
#                                                                             trial > 40 & trial <= 80 ~ trial-40,
#                                                                             trial > 80 & trial <=120 ~ trial-80, 
#                                                                             trial > 120 & trial <=160 ~ trial-120,
#                                                                             trial > 160 & trial <=200 ~ trial-160,
#                                                                             trial > 200 & trial <=240 ~ trial-200))
# 
# hc_ksoc <- hc %>%  filter(id %in% ksoc$id) %>% mutate(run_trial0 = case_when(trial <= 50 ~ trial, 
#                                                                              trial > 50 & trial <= 100 ~ trial-50,
#                                                                              trial > 100 & trial <=150 ~ trial-100, 
#                                                                              trial > 150 & trial <=200 ~ trial-150,
#                                                                              trial > 200 & trial <=250 ~ trial-200,
#                                                                              trial > 250 & trial <=300 ~ trial-250))
# 
# hc <- rbind(hc_bsoc,hc_ksoc) %>% select(!run_trial) %>% rename(run_trial = run_trial0)
# 
# # Compress data from 12 bins to 2 by averaging across anterior 6 bins and posterior 6 bins to create
# # AH and PH averages
# hc <- hc %>% group_by(id,run,trial,evt_time,HC_region) %>%
#   summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% 
#   ungroup() # 12 -> 2
# 
# # Create a new scaled within-subjects variable (HCwithin) and a between-subjects
# # variable averaged per subject and run (HCbetween)
# hc <- hc %>% group_by(id,run) %>%
#   mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>%
#   ungroup()
# 
# rm(hc_bsoc,hc_ksoc)
# gc()
# 
# 
# 
# 
