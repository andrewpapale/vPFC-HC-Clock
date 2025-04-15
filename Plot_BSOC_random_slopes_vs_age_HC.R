# 2025-04-14 AndyP
# Plotting random slopes vs. age and AH/PH for devlopmental project

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

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2025-04-09-Bsocial-vPFC-HC-network-clock-RTcorrected-ranslopes-1.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2025-04-09-Bsocial-vPFC-HC-network-clock-RTcorrected-ranslopes-2.Rdata')
Q <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')

Q <- Q %>% rename(id=level)
Q$id <- as.character(Q$id)

Qah <- Q %>% filter(HC_region == "AH")
Qph <- Q %>% filter(HC_region == "PH")

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

Qph$female <- ifelse(Qph$sex==1,1,0)
Qph <- Qph %>% select(!sex)

Qah <- Qah %>% mutate(sex = case_when(female == 1 ~ 'F',
                                      female == 0 ~ 'M'))

Qph <- Qph %>% mutate(sex = case_when(female == 1 ~ 'F',
                                      female == 0 ~ 'M'))

Qah$sex <- as.factor(Qah$sex)
Qph$sex <- as.factor(Qph$sex)

Q <- rbind(Qah,Qph)
Q <- Q %>% mutate(age_bin = case_when(age < 20 ~ '< 20',
                                      age >= 20 & age < 25 ~ '20-25',
                                      age >= 25 & age < 30 ~ '25-30',
                                      age >= 30 & age < 35 ~ '30-35',
                                      age >= 35 & age < 40 ~ '35-40',
                                      age >= 40 & age < 45 ~ '40-45',
                                      age >= 45 & age < 50 ~ '45-50',
                                      age >= 50 ~ '>=50'))

Q <- Q %>% select(id,age_bin,HC_region,estimate,sex,network,evt_time)
#Q <- Q %>% filter(evt_time >=-2 & evt_time <=2)

Q1 <- Q %>% group_by(age_bin,sex,HC_region,network) %>% summarize(estimate0 = mean(estimate,na.rm=TRUE),sd0 = sd(estimate,na.rm=TRUE), N=n()) %>% ungroup()


#ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex,group=sex)) + geom_line() + facet_wrap(~HC_region)

Q1$age_bin <- factor(Q1$age_bin, levels = c('< 20','20-25','25-30','30-35','35-40','40-45','45-50','>=50'))

ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex,group=sex)) + 
  geom_line() + geom_errorbar(aes(ymin=estimate0-sd0/sqrt(N), ymax = estimate0+sd0/sqrt(N))) +
  facet_wrap(network~HC_region)

ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex, group=sex)) + geom_smooth()

load(file.path(rootdir,'BSOC_HC_clock_TRdiv2.Rdata'))
hc <- hc %>% filter(evt_time > -5 & evt_time < 5)

split_ksoc_bsoc <- hc %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300])
bsoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240])
bsoc <- rbind(bsoc,221973,221507,221842,440223)
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
hc <- hc %>% group_by(id,run,trial,evt_time,HC_region) %>%
  summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>%
  ungroup() # 12 -> 2

# Create a new scaled within-subjects variable (HCwithin) and a between-subjects
# variable averaged per subject and run (HCbetween)
hc <- hc %>% group_by(id,run) %>%
  mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>%
  ungroup()

rm(hc_bsoc,hc_ksoc)
gc()

hc0 <- hc %>% group_by(id,HC_region) %>% filter(evt_time > -2 & evt_time < 2) %>% summarize(HCwithin = mean(HCwithin,na.rm=TRUE)) %>% ungroup()

Q2 <- inner_join(Q,hc0,by=c('id','HC_region'))

Q2 <- Q2 %>% mutate(HC_bin = ntile(HCwithin,6))

Q3 <- Q2 %>% group_by(id,HC_bin,sex,network,HC_region) %>% summarize(mE = mean(estimate,na.rm=TRUE)) %>% ungroup()

ggplot(Q3, aes(x=HC_bin,y=mE,color=sex,group=sex)) + geom_smooth() + facet_wrap(network~HC_region,scales='free_y')
ggplot(Q3 %>% filter(network=='LIM' & HC_region=='AH'), aes(x=HC_bin,y=mE,color=sex,group=sex)) + geom_smooth()
