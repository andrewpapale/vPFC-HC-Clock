# 2024-08-20 AndyP
# run acf on vPFC and HC

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

do_MMClock <- TRUE
do_Explore <- TRUE
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26
source("~/fmri.pipeline/R/mixed_by.R")

if (do_MMClock){
  rm(vmPFC)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMclock_clock_Aug2023.Rdata')
  vmPFC <- clock_comb
  #vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(clock_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  vmPFC <- vmPFC %>% group_by(id,run,run_trial,evt_time,network) %>% summarize(decon1 = mean(vmPFC_decon,na.rm=TRUE)) %>% ungroup()
  vmPFC <- vmPFC %>% group_by(id,run) %>% mutate(vmPFC_within = scale(decon1), vmPFC_between = mean(decon1,na.rm=TRUE)) %>% ungroup()
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
  #hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
  #hc <- hc %>% rename(HC_decon = decon_mean)
  #hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
  #hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
  
  #vmPFC <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
  vmPFC <- vmPFC %>% select(!decon1)
  
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
  df <- df %>% 
    group_by(id,run) %>%
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           iti_lag_sc = scale(iti_prev),
           score_sc = scale(score_csv),
           score_lag_sc = scale(lag(score_csv)),
           iti_sc = scale(iti_ideal),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           rt_vmax_sc = scale(rt_vmax),
           ev_sc = scale(ev),
           ev_lag_sc = scale(lag(ev)),
           v_entropy_sc = scale(v_entropy),
           v_entropy_lag_sc = scale(lag(v_entropy)),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
           abs_rt_vmax_change = abs(rt_vmax_change),
           rt_vmax_change_bin = case_when(
             abs_rt_vmax_change < 4/24 ~ 'No Change',
             abs_rt_vmax_change >= 4/24 ~ 'Change'
           ),
           rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld4 = case_when(
             kld4 ==0 ~ NA_real_,
             kld4 >0 ~ log10(kld4)
           )) %>% mutate(log10kld4_lag = case_when(
             kld4_lag==0 ~NA_real_,
             kld4_lag>0 ~ log10(kld4_lag)
           ))
  df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
    rt_csv - rt_lag > 1 ~ 'Longer',
    rt_csv - rt_lag < -1 ~ '0',
    rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
  )))
  df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
    rt_csv - rt_lag > 1 ~ '0',
    rt_csv - rt_lag < -1 ~ 'Shorter',
    rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
  )))
  df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
    rt_csv_sc <= -1 ~ '-1',
    rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
    rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
    rt_csv_sc > 1 ~ '1'
  )))
  df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
    run_trial <= 15 ~ 'Early',
    run_trial > 15 & run_trial < 30 ~ 'Middle',
    run_trial >=30 ~ 'Late',
  )))
  
  df <- df %>% select(id,run,trial,run_trial,trial_neg_inv_sc,iti_ideal,iti_prev,rt_csv,rt_csv_sc,iti_sc,rt_lag_sc,rt_csv_sc,rewFunc,iti_lag_sc)
  
  vmPFC <- inner_join(vmPFC,df,by=c('id','run','run_trial'))
  vmPFC$vmPFC_decon[vmPFC$evt_time > vmPFC$rt_csv + vmPFC$iti_ideal] = NA;
  vmPFC$vmPFC_decon[vmPFC$evt_time < -(vmPFC$iti_prev)] = NA;
  #vmPFC$decon_mean[vmPFC$evt_time > vmPFC$rt_csv + vmPFC$iti_ideal] = NA;
  #vmPFC$decon_mean[vmPFC$evt_time < -(vmPFC$iti_prev)] = NA;  
  
  # test age & sex
  demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
  demo <- demo %>% rename(id=lunaid)
  demo <- demo %>% select(!adult & !scandate)
  vmPFC <- inner_join(vmPFC,demo,by=c('id'))
  vmPFC$female <- relevel(as.factor(vmPFC$female),ref='0')
  vmPFC$age <- scale(vmPFC$age)
  
  #vmPFC <- vmPFC %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
  vmPFC <- vmPFC %>% group_by(network) %>% mutate(vmPFC_between1 = scale(vmPFC_between)) %>% select(!vmPFC_between) %>% rename(vmPFC_between = vmPFC_between1)
  vmPFC1 <- vmPFC %>% group_by(id,run,trial,network) %>% 
    select(id,run,trial,evt_time,vmPFC_within,network) %>% 
    pivot_wider(values_from=vmPFC_within,names_from=evt_time)
  acfPFC <- vmPFC1 %>% 
         group_by(network) %>%
         nest() %>% 
         mutate(data = map(data, ~acf(., lag.max=5, type="correlation", plot=F,na.action=na.pass))) %>%
         mutate(data = map(data, ~as.data.frame(rbind(.x$acf[1,,],.x$acf[2,,],.x$acf[3,,],.x$acf[4,,],.x$acf[5,,])))) %>%
         unnest(data)
  
}