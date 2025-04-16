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

demo <- demo %>% mutate(sex = case_when(female == 1 ~ 'F',
                                          female == 0 ~ 'M'))
demo$sex <- as.factor(demo$sex)

Qah <- inner_join(Qah,demo,by=c('id'))
Qph <- inner_join(Qph,demo,by=c('id'))

Q <- rbind(Qah,Qph)
Q <- Q %>% mutate(age_bin = case_when(age < 15 ~ '< 15',
                                      age >= 15 & age < 18 ~ '15-18',
                                      age >= 18 & age < 21 ~ '18-21',
                                      age >= 21 & age < 24 ~ '21-24',
                                      age >= 24 & age < 27 ~ '24-27',
                                      age >= 27 ~ '>= 27'))

Q <- Q %>% select(id,age_bin,HC_region,estimate,sex,network)
Q <- Q %>% group_by(network,HC_region) %>% mutate(estimate0 = scale(estimate)) %>% ungroup()

Q4 <- Q %>% group_by(id,HC_region,network) %>% summarize(estimate0 = mean(estimate,na.rm=TRUE)) %>% ungroup()
age0 <- demo %>% select(id,age,sex)
Q4 <- inner_join(Q4,age0,by=c('id'))
Q4 <- Q4 %>% group_by(network,HC_region) %>% mutate(estimate0 = scale(estimate0))

ggplot(Q4, aes(x=age,y=estimate0,color=network)) + geom_smooth() + geom_jitter() + facet_wrap(~sex) + xlab('age') + ylab('random slope scaled') + ggtitle('MMClock')


Q1 <- Q %>% group_by(id,age_bin,sex,HC_region,network) %>% summarize(estimate0 = mean(estimate0,na.rm=TRUE),sd0 = sd(estimate,na.rm=TRUE), N=n()) %>% ungroup()


#ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex,group=sex)) + geom_line() + facet_wrap(~HC_region)

Q1$age_bin <- factor(Q1$age_bin, levels = c('< 15','15-18','18-21','21-24','24-27','>= 27'))

# ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex,group=sex)) + 
#   geom_line() + geom_errorbar(aes(ymin=estimate0-sd0/sqrt(N), ymax = estimate0+sd0/sqrt(N))) +
#   facet_wrap(network~HC_region)


ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex, group=sex)) + geom_smooth() + facet_wrap(network~HC_region, scales='free_y') + xlab('age binned') + ylab('estimate scaled') + ggtitle('MMClock')
ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex, group=sex)) + geom_smooth() + xlab('age binned') + ylab('estimate scaled') + ggtitle('MMClock')


#ggplot(Q1, aes(x=age_bin,y=estimate0,color=sex, group=sex)) + geom_smooth() + ylab('age binned') + xlab('estimate scaled') + ggtitle('MMClock')

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
hc <- hc %>% filter(evt_time > -5 & evt_time < 5)

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

hc0 <- hc %>% group_by(id,HC_region) %>% filter(evt_time > -2 & evt_time < 2) %>% summarize(HCwithin0 = mean(HCwithin,na.rm=TRUE)) %>% ungroup()

hc0 <- hc0 %>% group_by(HC_region) %>% mutate(HCwithin0 = scale(HCwithin0))
hc0$id <- as.character(hc0$id)

Q2 <- inner_join(Q,hc0,by=c('id','HC_region'))

Q5 <- Q2 %>% group_by(id,HC_region,network) %>% summarize(estimate0 = mean(estimate,na.rm=TRUE),meanHC = mean(HCwithin0,na.rm=TRUE)) %>% ungroup()
age0 <- demo %>% select(id,age,sex)
Q6 <- inner_join(Q5,age0,by=c('id'))
Q6 <- Q6 %>% group_by(network,HC_region) %>% mutate(estimate0 = scale(estimate0))

ggplot(Q6, aes(x=meanHC,y=estimate0,color=network)) + geom_smooth() + geom_jitter() + facet_wrap(~sex) + xlab('mean scaled HC activity') + ylab('random slope scaled') + ggtitle('MMClock')


Q2 <- Q2 %>% group_by(network,HC_region) %>% mutate(HC_bin = ntile(HCwithin0,6)) %>% ungroup()

Q3 <- Q2 %>% group_by(id,HC_bin,sex,network,HC_region) %>% summarize(estimate = mean(estimate,na.rm=TRUE)) %>% ungroup()

Q3 <- Q3 %>% group_by(network,HC_region) %>% mutate(estimate_scaled = scale(estimate)) %>% ungroup()

ggplot(Q3, aes(x=HC_bin,y=estimate_scaled,color=sex,group=sex)) + geom_smooth() + facet_wrap(network~HC_region,scales='free_y') + xlab('Hippocampal activity binned') + ylab('Scaled Random Slope') + ggtitle('MMClock')
ggplot(Q3 %>% filter(network=='LIM' & HC_region=='AH'), aes(x=HC_bin,y=estimate_scaled,color=sex,group=sex)) + geom_smooth() + xlab('Hippocampal activity binned') + ylab('Scaled Random Slope') + ggtitle('MMClock')
ggplot(Q3, aes(x=HC_bin,y=estimate_scaled,color=sex,group=sex)) + geom_smooth()+ xlab('Hippocampal activity binned') + ylab('Scaled Random Slope') + ggtitle('MMClock - All Networks')

