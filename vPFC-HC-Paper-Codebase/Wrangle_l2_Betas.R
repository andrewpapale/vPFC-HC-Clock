# 2025-09-08 AndyP
# wrangle L2 Betas

# CRC cluster /ix1/adombrovski/
# library(fmri.pipeline)
#L2hwizBetas <- extract_glm_betas_in_mask(gpa,vPFC_mask,extract_l3="none",extract_l2="all", what=c("cope","zstat"),scheduler="slurm",out_dir = "/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-clock-ITI-feedback/")

library(tidyverse)

clock_beta <- read.csv('/Users/andypapale/Downloads/L1m-gPPI.clock/Schaeffer_2018_vmPFC_mask_cope_l2.csv.gz') %>%
  select(id,session,l1_cope_name,mask_value,value) %>% rename(atlas_value=mask_value,beta=value) %>% filter(l1_cope_name == 'EV_clock-PM.ppi' | l1_cope_name == 'EV_clock-PM') %>%
  mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM')) %>%
  mutate(event = 'clock')
feedback_beta <- read.csv('/Users/andypapale/Downloads/L1m-gPPI.feedback/Schaeffer_2018_vmPFC_mask_cope_l2.csv.gz') %>% 
  select(id,session,l1_cope_name,mask_value,value) %>% rename(atlas_value=mask_value,beta=value) %>% filter(l1_cope_name == 'EV_feedback-PM.ppi' | l1_cope_name == 'EV_feedback-PM') %>%
  mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM')) %>%
  mutate(event = 'feedback')
ITI_beta <- read.csv('/Users/andypapale/Downloads/L1m-gPPI.ITI2/Schaeffer_2018_vmPFC_mask_cope_l2.csv.gz') %>%
  select(id,session,l1_cope_name,mask_value,value) %>% rename(atlas_value=mask_value,beta=value) %>% filter(l1_cope_name == 'EV_ITI-PM.ppi' | l1_cope_name == 'EV_ITI-PM') %>%
  mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM')) %>%
  mutate(event = 'ITI')

clock_betavar <- read.csv('/Users/andypapale/Downloads/L1m-gPPI.clock/Schaeffer_2018_vmPFC_mask_varcope_l2.csv.gz') %>%
  select(id,session,l1_cope_name,mask_value,value) %>% rename(atlas_value=mask_value,betavar=value) %>% filter(l1_cope_name == 'EV_clock-PM.ppi' | l1_cope_name == 'EV_clock-PM') %>%
  mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM')) %>%
  mutate(event = 'clock')
feedback_betavar <- read.csv('/Users/andypapale/Downloads/L1m-gPPI.feedback/Schaeffer_2018_vmPFC_mask_varcope_l2.csv.gz') %>% 
  select(id,session,l1_cope_name,mask_value,value) %>% rename(atlas_value=mask_value,betavar=value) %>% filter(l1_cope_name == 'EV_feedback-PM.ppi' | l1_cope_name == 'EV_feedback-PM') %>%
  mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM')) %>%
  mutate(event = 'feedback')
ITI_betavar <- read.csv('/Users/andypapale/Downloads/L1m-gPPI.ITI2/Schaeffer_2018_vmPFC_mask_varcope_l2.csv.gz') %>%
  select(id,session,l1_cope_name,mask_value,value) %>% rename(atlas_value=mask_value,betavar=value) %>% filter(l1_cope_name == 'EV_ITI-PM.ppi' | l1_cope_name == 'EV_ITI-PM') %>%
  mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM')) %>%
  mutate(event = 'ITI')

clock_zstat <- read.csv('/Users/andypapale/Downloads/L1m-gPPI.clock/Schaeffer_2018_vmPFC_mask_zstat_l2.csv.gz') %>%
  select(id,session,l1_cope_name,mask_value,value) %>% rename(atlas_value=mask_value,zstat=value) %>% filter(l1_cope_name == 'EV_clock-PM.ppi' | l1_cope_name == 'EV_clock-PM') %>%
  mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM')) %>%
  mutate(event = 'clock')
feedback_zstat <- read.csv('/Users/andypapale/Downloads/L1m-gPPI.feedback/Schaeffer_2018_vmPFC_mask_zstat_l2.csv.gz') %>% 
  select(id,session,l1_cope_name,mask_value,value) %>% rename(atlas_value=mask_value,zstat=value) %>% filter(l1_cope_name == 'EV_feedback-PM.ppi' | l1_cope_name == 'EV_feedback-PM') %>%
  mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM')) %>%
  mutate(event = 'feedback')
ITI_zstat <- read.csv('/Users/andypapale/Downloads/L1m-gPPI.ITI2/Schaeffer_2018_vmPFC_mask_zstat_l2.csv.gz') %>%
  select(id,session,l1_cope_name,mask_value,value) %>% rename(atlas_value=mask_value,zstat=value) %>% filter(l1_cope_name == 'EV_ITI-PM.ppi' | l1_cope_name == 'EV_ITI-PM') %>%
  mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM')) %>%
  mutate(event = 'ITI')

clock <- inner_join(clock_beta,clock_betavar,by=c('id','session','l1_cope_name','atlas_value'))
clock <- inner_join(clock,clock_zstat,by=c('id','session','l1_cope_name','atlas_value'))

feedback <- inner_join(feedback_beta,feedback_betavar,by=c('id','session','l1_cope_name','atlas_value'))
feedback <- inner_join(feedback,feedback_zstat,by=c('id','session','l1_cope_name','atlas_value'))

ITI <- inner_join(ITI_beta,ITI_betavar,by=c('id','session','l1_cope_name','atlas_value'))
ITI <- inner_join(ITI,ITI_zstat,by=c('id','session','l1_cope_name','atlas_value'))

gPPI_l2_betas <- rbind(clock,feedback,ITI)

