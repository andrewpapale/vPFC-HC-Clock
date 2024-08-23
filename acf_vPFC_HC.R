# 2024-08-20 AndyP
# run acf on vPFC and HC

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(ggplot2)

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
  vmPFC <- vmPFC %>% filter(evt_time > -7 & evt_time < 10)
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
  vmPFC$vmPFC_within[vmPFC$evt_time > vmPFC$rt_csv + vmPFC$iti_ideal] = NA;
  vmPFC$vmPFC_within[vmPFC$evt_time < -(vmPFC$iti_prev)] = NA;
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
  vmPFC <- vmPFC %>% group_by(id,run,trial,network) %>% 
    mutate(vmPFC_lag1 = lag(vmPFC_within,1), 
           vmPFC_lag2 = lag(vmPFC_within,2), 
           vmPFC_lag3 = lag(vmPFC_within,3),
           vmPFC_lag4 = lag(vmPFC_within,4),
           vmPFC_lag5 = lag(vmPFC_within,5),
           vmPFC_lag6 = lag(vmPFC_within,6),
           vmPFC_lead1 = lead(vmPFC_within,1),
           vmPFC_lead2 = lead(vmPFC_within,2),
           vmPFC_lead3 = lead(vmPFC_within,3),
           vmPFC_lead4 = lead(vmPFC_within,4),
           vmPFC_lead5 = lead(vmPFC_within,5),
           vmPFC_lead6 = lead(vmPFC_within,6),
           vmPFC_lead7 = lead(vmPFC_within,7),
           vmPFC_lead8 = lead(vmPFC_within,8),
           vmPFC_lead9 = lead(vmPFC_within,9)
    ) %>% ungroup() %>% filter(evt_time==0)
  
  vmPFC <- vmPFC %>% group_by(network) %>% summarize(cor_lag0 = cor(vmPFC_within,vmPFC_within,use='complete.obs',method='pearson'),
                                                     cor_lag1 = cor(vmPFC_within,vmPFC_lag1,use='complete.obs',method='pearson'),
                                                     cor_lag2 = cor(vmPFC_within,vmPFC_lag2,use='complete.obs',method='pearson'),
                                                     cor_lag3 = cor(vmPFC_within,vmPFC_lag3,use='complete.obs',method='pearson'),
                                                     cor_lag4 = cor(vmPFC_within,vmPFC_lag4,use='complete.obs',method='pearson'),
                                                     cor_lag5 = cor(vmPFC_within,vmPFC_lag5,use='complete.obs',method='pearson'),
                                                     cor_lag6 = cor(vmPFC_within,vmPFC_lag6,use='complete.obs',method='pearson'),
                                                     cor_lead1 = cor(vmPFC_within,vmPFC_lead1,use='complete.obs',method='pearson'),
                                                     cor_lead2 = cor(vmPFC_within,vmPFC_lead2,use='complete.obs',method='pearson'),
                                                     cor_lead3 = cor(vmPFC_within,vmPFC_lead3,use='complete.obs',method='pearson'),
                                                     cor_lead4 = cor(vmPFC_within,vmPFC_lead4,use='complete.obs',method='pearson'),
                                                     cor_lead5 = cor(vmPFC_within,vmPFC_lead5,use='complete.obs',method='pearson'),
                                                     cor_lead6 = cor(vmPFC_within,vmPFC_lead6,use='complete.obs',method='pearson'),
                                                     cor_lead7 = cor(vmPFC_within,vmPFC_lead7,use='complete.obs',method='pearson'),
                                                     cor_lead8 = cor(vmPFC_within,vmPFC_lead8,use='complete.obs',method='pearson'),
                                                     cor_lead9 = cor(vmPFC_within,vmPFC_lead9,use='complete.obs',method='pearson')
  ) %>% ungroup()
  
  # vmPFC1 <- vmPFC %>% group_by(id,run,trial,network) %>% 
  #   select(id,run,trial,evt_time,vmPFC_within,network) %>% 
  #   pivot_wider(values_from=vmPFC_within,names_from=evt_time)
  # acfPFC <- vmPFC1 %>% 
  #        group_by(network) %>% select(!run & !trial & !id) %>%
  #        nest() %>% 
  #        mutate(data = map(data, ~acf(., lag.max=5, type="correlation", plot=F,na.action=na.pass))) %>%
  #        mutate(data = map(data, ~as.data.frame(rbind(.x$acf[1,,],.x$acf[2,,],.x$acf[3,,],.x$acf[4,,],.x$acf[5,,],.x$acf[6,,])))) %>%
  #        unnest(data) %>% ungroup()
  # 
  # acfPFC <- acfPFC %>% rename('-4'=V1,'-3'=V2,'-2'=V3,'-1'=V4,'0'=V5,'1'=V6,'2'=V7,'3'=V8,'4'=V9) %>% 
  #   mutate(lags1 = rep(c(t(rep(0,9)),t(rep(1,9)),t(rep(2,9)),t(rep(3,9)),t(rep(4,9)),t(rep(5,9))),3),lags2 = rep(c(t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9))),3))
  # acfPFC <- acfPFC %>% pivot_longer(cols=c("-4","-3","-2","-1","0","1","2","3","4"))
  # acfPFC_Mmc <- acfPFC %>% rename(acf=value,evt_time=name)
  # acfPFC_Mmc$evt_time <- as.numeric(acfPFC_Mmc$evt_time)
  # vPFC_Mmc0 <- acfPFC_Mmc %>% filter(evt_time==0 & ((lags1==0 & lags2==0) | 
  #                                                     (lags1==1 & lags2==-1) | (lags1==1 & lags2==1) | 
  #                                                     (lags1==2 & lags2==-2) | (lags1==2 & lags2==2) |
  #                                                     (lags1==3 & lags2==-3) | (lags1==3 & lags2==3) |
  #                                                     (lags1==4 & lags2==-4) | (lags1==4 & lags2==4)))
  vmPFC <- vmPFC %>% mutate(HC_region='NA')
  vPFC_Mmc <- vmPFC %>% mutate(dataset='vPFC')
  
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
  hc <- hc %>% filter(evt_time > -7 & evt_time < 10)
  hc <- hc %>% rename(HC_decon = decon_mean)
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(HC_decon,na.rm=TRUE)) %>% ungroup() # 12 -> 2
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
  
  #vmPFC <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
  hc <- hc %>% select(!decon1)
  
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
  
  hc <- inner_join(hc,df,by=c('id','run','run_trial'))
  hc$HCwithin[hc$evt_time > hc$rt_csv + hc$iti_ideal] = NA;
  #vmPFC$HC_decon[vmPFC$evt_time < -(vmPFC$iti_prev)] = NA;
  #vmPFC$decon_mean[vmPFC$evt_time > vmPFC$rt_csv + vmPFC$iti_ideal] = NA;
  hc$HCwithin[hc$evt_time < -(hc$iti_prev)] = NA;  
  
  # test age & sex
  demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
  demo <- demo %>% rename(id=lunaid)
  demo <- demo %>% select(!adult & !scandate)
  hc <- inner_join(hc,demo,by=c('id'))
  hc$female <- relevel(as.factor(hc$female),ref='0')
  hc$age <- scale(hc$age)
  
  #  hc1 <- hc %>% group_by(id,run,trial,HC_region) %>% 
  #   select(id,run,trial,evt_time,HCwithin,HC_region) %>% 
  #   pivot_wider(values_from=HCwithin,names_from=evt_time)
  # acfHC <- hc1 %>% 
  #   group_by(HC_region) %>% select(!run & !trial & !id) %>%
  #   nest() %>% 
  #   mutate(data = map(data, ~acf(., lag.max=5, type="correlation", plot=F,na.action=na.pass))) %>%
  #   mutate(data = map(data, ~as.data.frame(rbind(.x$acf[1,,],.x$acf[2,,],.x$acf[3,,],.x$acf[4,,],.x$acf[5,,],.x$acf[6,,])))) %>%
  #   unnest(data) %>% ungroup()
  # 
  # acfHC <- acfHC %>% rename('-4'=V1,'-3'=V2,'-2'=V3,'-1'=V4,'0'=V5,'1'=V6,'2'=V7,'3'=V8,'4'=V9) %>% 
  #   mutate(lags1 = rep(c(t(rep(0,9)),t(rep(1,9)),t(rep(2,9)),t(rep(3,9)),t(rep(4,9)),t(rep(5,9))),2),lags2 = rep(c(t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9)),t(seq(-4,4,length.out=9))),2))
  # acfHC <- acfHC %>% pivot_longer(cols=c("-4","-3","-2","-1","0","1","2","3","4"))
  # acfHC_Mmc <- acfHC %>% rename(acf=value,evt_time=name)
  # acfHC_Mmc$evt_time <- as.numeric(acfHC_Mmc$evt_time)
  # HC_Mmc0 <- acfHC_Mmc %>% filter(evt_time==0 & ((lags1==0 & lags2==0) | 
  #                                   (lags1==1 & lags2==-1) | (lags1==1 & lags2==1) | 
  #                                   (lags1==2 & lags2==-2) | (lags1==2 & lags2==2) |
  #                                   (lags1==3 & lags2==-3) | (lags1==3 & lags2==3) |
  #                                   (lags1==4 & lags2==-4) | (lags1==4 & lags2==4)))
  
  hc <- hc %>% group_by(id,run,trial,HC_region) %>% 
    mutate(hc_lag1 = lag(HCwithin,1), 
           hc_lag2 = lag(HCwithin,2), 
           hc_lag3 = lag(HCwithin,3),
           hc_lag4 = lag(HCwithin,4),
           hc_lag5 = lag(HCwithin,5),
           hc_lag6 = lag(HCwithin,6),
           hc_lead1 = lead(HCwithin,1),
           hc_lead2 = lead(HCwithin,2),
           hc_lead3 = lead(HCwithin,3),
           hc_lead4 = lead(HCwithin,4),
           hc_lead5 = lead(HCwithin,5),
           hc_lead6 = lead(HCwithin,6),
           hc_lead7 = lead(HCwithin,7),
           hc_lead8 = lead(HCwithin,8),
           hc_lead9 = lead(HCwithin,9)
    ) %>% ungroup() %>% filter(evt_time==0)
  
  hc <- hc %>% group_by(HC_region) %>% summarize(cor_lag0 = cor(HCwithin,HCwithin,use='complete.obs',method='pearson'),
                                                     cor_lag1 = cor(HCwithin,hc_lag1,use='complete.obs',method='pearson'),
                                                     cor_lag2 = cor(HCwithin,hc_lag2,use='complete.obs',method='pearson'),
                                                     cor_lag3 = cor(HCwithin,hc_lag3,use='complete.obs',method='pearson'),
                                                     cor_lag4 = cor(HCwithin,hc_lag4,use='complete.obs',method='pearson'),
                                                     cor_lag5 = cor(HCwithin,hc_lag5,use='complete.obs',method='pearson'),
                                                     cor_lag6 = cor(HCwithin,hc_lag6,use='complete.obs',method='pearson'),
                                                     cor_lead1 = cor(HCwithin,hc_lead1,use='complete.obs',method='pearson'),
                                                     cor_lead2 = cor(HCwithin,hc_lead2,use='complete.obs',method='pearson'),
                                                     cor_lead3 = cor(HCwithin,hc_lead3,use='complete.obs',method='pearson'),
                                                     cor_lead4 = cor(HCwithin,hc_lead4,use='complete.obs',method='pearson'),
                                                     cor_lead5 = cor(HCwithin,hc_lead5,use='complete.obs',method='pearson'),
                                                     cor_lead6 = cor(HCwithin,hc_lead6,use='complete.obs',method='pearson'),
                                                     cor_lead7 = cor(HCwithin,hc_lead7,use='complete.obs',method='pearson'),
                                                     cor_lead8 = cor(HCwithin,hc_lead8,use='complete.obs',method='pearson'),
                                                     cor_lead9 = cor(HCwithin,hc_lead9,use='complete.obs',method='pearson')
  ) %>% ungroup()
  hc <- hc %>% mutate(network='NA')
  HC_Mmc <- hc %>% mutate(dataset='HC')
  
  dq_Mmc <- rbind(vPFC_Mmc,HC_Mmc)
}


if (do_Explore){
  rm(vmPFC)
  load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')
  vmPFC <- md %>% filter(evt_time > -7 & evt_time < 10)
  rm(md)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
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
  vmPFC <- vmPFC %>% group_by(id,run,trial,evt_time,network) %>% summarize(decon1 = mean(vmPFC_decon,na.rm=TRUE)) %>% ungroup()
  vmPFC <- vmPFC %>% group_by(id,run) %>% mutate(vmPFC_within = scale(decon1), vmPFC_between = mean(decon1,na.rm=TRUE)) %>% ungroup()

  
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory='/Volumes/Users/Andrew/MEDuSA_data_Explore',dataset='explore')
  df <- df %>%
    group_by(id,run_number) %>%
    arrange(id, run_number, trial) %>%
    mutate(condition_trial = run_trial-floor(run_trial/40.5)*40,
           condition_trial_neg_inv = -1000 / condition_trial,
           condition_trial_neg_inv_sc = as.vector(scale(condition_trial_neg_inv)),
    ) %>% ungroup()
  df <- df %>%
    group_by(id) %>% 
    mutate(v_entropy_sc_r = scale(v_entropy)) %>% ungroup() %>%
    group_by(id, run) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           score_sc = scale(score_csv),
           iti_sc = scale(iti_ideal),
           iti_lag_sc = lag(iti_sc),
           v_chosen_sc = scale(v_chosen),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           rt_vmax_sc = scale(rt_vmax),
           ev_sc = scale(ev),
           v_entropy_sc = scale(v_entropy),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
           abs_rt_vmax_change = abs(rt_vmax_change),
           rt_vmax_change_bin = case_when(
             abs_rt_vmax_change < 4/24 ~ 'No Change',
             abs_rt_vmax_change >= 4/24 ~ 'Change'
           ),
           v_entropy_wi_change_bin = case_when(
             v_entropy_wi_change < -0.5 ~ 'Decrease',
             v_entropy_wi_change > 0.5  ~ 'Increase',
             v_entropy_wi_change >= -0.5 & v_entropy_wi_change <= 0.5 ~ 'No Change'
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
    run_trial <= 13 ~ 'Early',
    run_trial > 13 & run_trial < 26 ~ 'Middle',
    run_trial >=26 ~ 'Late',
  )))
  df <- df %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                       total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))
  
  df <- df %>% select(total_earnings_split,condition_trial_neg_inv_sc,iti_ideal,rt_lag_sc,iti_lag_sc,iti_prev,iti_sc,v_entropy_wi_change_lag,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
  df$id <- as.character(df$id)
  vmPFC <- full_join(vmPFC,df,by=c('id','run','trial'))
  vmPFC <- vmPFC %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                      trial > 40 & trial <= 80 ~ 2,
                                      trial > 80 & trial <=120 ~ 3, 
                                      trial > 120 & trial <=160 ~ 4,
                                      trial > 160 & trial <=200 ~ 5,
                                      trial > 200 & trial <=240 ~ 6))
  vmPFC <- vmPFC %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                           trial > 40 & trial <= 80 ~ trial-40,
                                           trial > 80 & trial <=120 ~ trial-80, 
                                           trial > 120 & trial <=160 ~ trial-120,
                                           trial > 160 & trial <=200 ~ trial-160,
                                           trial > 200 & trial <=240 ~ trial-200))
  vmPFC <- vmPFC %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                    run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                    run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))
  #Q <- inner_join(Q,hc,by=c('id','run','trial','evt_time'))
  #rm(hc)
  #Q$HC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  vmPFC$vmPFC_within[vmPFC$evt_time > vmPFC$rt_csv + vmPFC$iti_ideal] = NA
  #Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  #Q$HC_decon[Q$evt_time < -Q$iti_prev] = NA
  vmPFC$vmPFC_within[vmPFC$evt_time < -vmPFC$iti_prev] = NA
  vmPFC <- vmPFC %>% arrange(id,run,trial,evt_time)
  
  
  vmPFC$rt_bin <- relevel(as.factor(vmPFC$rt_bin),ref='-0.5')
  vmPFC$trial_bin <- relevel(as.factor(vmPFC$trial_bin),ref='Middle')
  
  demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
  demo <- demo %>% select(registration_redcapid,age,gender,registration_group,wtar,education_yrs)
  demo <- demo %>% rename(id=registration_redcapid,group=registration_group)
  demo$gender <- relevel(as.factor(demo$gender),ref='M')
  demo$age <- scale(demo$age)
  demo$wtar <- scale(demo$wtar)
  demo$education_yrs <- scale(demo$education_yrs)
  
  vmPFC <- merge(demo,vmPFC,by='id')
  #Q <- Q %>% filter(total_earnings_split=='richer')
  vmPFC$group <- relevel(factor(vmPFC$group),ref='HC')
  vmPFC <- vmPFC %>% filter(group!='ATT')
  vmPFC <- vmPFC %>% filter(group=='HC')
  vmPFC <- vmPFC %>% filter(!is.na(rewFunc))
  #Q <- Q %>% filter(trial > 10)
  
  # #Q <- Q %>% group_by(network,HC_region) %>% mutate(HCbetween1 = scale(HCbetween)) %>% select(!HCbetween) %>% rename(HCbetween=HCbetween1)
  # vmPFC <- vmPFC %>% group_by(network) %>% mutate(vmPFC_between1 = scale(vmPFC_between)) %>% select(!vmPFC_between) %>% rename(vmPFC_between = vmPFC_between1)
  # vmPFC1 <- vmPFC %>% group_by(id,run,trial,network) %>% 
  #   select(id,run,trial,evt_time,vmPFC_within,network) %>% 
  #   pivot_wider(values_from=vmPFC_within,names_from=evt_time)
  # acfPFC <- vmPFC1 %>% 
  #   group_by(network) %>% select(!run & !trial & !id) %>%
  #   nest() %>% 
  #   mutate(data = map(data, ~acf(., lag.max=5, type="correlation", plot=F,na.action=na.pass))) %>%
  #   mutate(data = map(data, ~as.data.frame(rbind(.x$acf[1,,],.x$acf[2,,],.x$acf[3,,],.x$acf[4,,],.x$acf[5,,],.x$acf[6,,])))) %>%
  #   unnest(data) %>% ungroup()
  # 
  # acfPFC <- acfPFC %>% rename('-3.6'=V1,'-3'=V2,'-2.4'=V3,'-1.8'=V4,'-1.2'=V5,'-0.6'=V6,'0'=V7,'0.6'=V8,'1.2'=V9,'1.8'=V10,'2.4'=V11,'3'=V12,'3.6'=V13) %>% 
  #   mutate(lags1 = rep(c(t(rep(0,13)),t(rep(1,13)),t(rep(2,13)),t(rep(3,13)),t(rep(4,13)),t(rep(5,13))),3),
  #          lags2 = rep(c(t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1))),3)
  #   )
  # acfPFC <- acfPFC %>% pivot_longer(cols=c("-3.6","-3","-2.4","-1.8","-1.2","-0.6","0","0.6","1.2","1.8","2.4","3","3.6"))
  # acfPFC_Exp <- acfPFC %>% rename(acf=value,evt_time=name)
  # acfPFC_Exp$evt_time <- as.numeric(acfPFC_Exp$evt_time)
  # vPFC_Exp0 <- acfPFC_Exp %>% filter(evt_time==0 & ((lags1==0 & lags2==0) | 
  #                                                     (lags1==1 & lags2==-0.6) | (lags1==1 & lags2==0.6) | 
  #                                                     (lags1==2 & lags2==-1.2) | (lags1==2 & lags2==1.2) |
  #                                                     (lags1==3 & lags2==-1.8) | (lags1==3 & lags2==1.8) |
  #                                                     (lags1==4 & lags2==-2.4) | (lags1==4 & lags2==2.4) |
  #                                                     (lags1==5 & lags2==-3) | (lags1==5 & lags2==3) |
  #                                                     (lags1==6 & lags2==-3.6) | (lags1==6 & lags2==3.6))
  # )
  
  vmPFC <- vmPFC %>% group_by(id,run,trial,network) %>% 
    mutate(vmPFC_lag1 = lag(vmPFC_within,1), 
           vmPFC_lag2 = lag(vmPFC_within,2), 
           vmPFC_lag3 = lag(vmPFC_within,3),
           vmPFC_lag4 = lag(vmPFC_within,4),
           vmPFC_lag5 = lag(vmPFC_within,5),
           vmPFC_lag6 = lag(vmPFC_within,6),
           vmPFC_lag7 = lag(vmPFC_within,7),
           vmPFC_lag8 = lag(vmPFC_within,8),
           vmPFC_lag9 = lag(vmPFC_within,9),
           vmPFC_lag10 = lag(vmPFC_within,10),
           vmPFC_lead1 = lead(vmPFC_within,1),
           vmPFC_lead2 = lead(vmPFC_within,2),
           vmPFC_lead3 = lead(vmPFC_within,3),
           vmPFC_lead4 = lead(vmPFC_within,4),
           vmPFC_lead5 = lead(vmPFC_within,5),
           vmPFC_lead6 = lead(vmPFC_within,6),
           vmPFC_lead7 = lead(vmPFC_within,7),
           vmPFC_lead8 = lead(vmPFC_within,8),
           vmPFC_lead9 = lead(vmPFC_within,9),
           vmPFC_lead10 = lead(vmPFC_within,10),
           vmPFC_lead11 = lead(vmPFC_within,11),
           vmPFC_lead12 = lead(vmPFC_within,12),
           vmPFC_lead13 = lead(vmPFC_within,13),
           vmPFC_lead14 = lead(vmPFC_within,14),
           vmPFC_lead15 = lead(vmPFC_within,15)
           
    ) %>% ungroup() %>% filter(evt_time==0)
  
  vmPFC <- vmPFC %>% group_by(network) %>% summarize(cor_lag0 = cor(vmPFC_within,vmPFC_within,use='complete.obs',method='pearson'),
                                                     cor_lag1 = cor(vmPFC_within,vmPFC_lag1,use='complete.obs',method='pearson'),
                                                     cor_lag2 = cor(vmPFC_within,vmPFC_lag2,use='complete.obs',method='pearson'),
                                                     cor_lag3 = cor(vmPFC_within,vmPFC_lag3,use='complete.obs',method='pearson'),
                                                     cor_lag4 = cor(vmPFC_within,vmPFC_lag4,use='complete.obs',method='pearson'),
                                                     cor_lag5 = cor(vmPFC_within,vmPFC_lag5,use='complete.obs',method='pearson'),
                                                     cor_lag6 = cor(vmPFC_within,vmPFC_lag6,use='complete.obs',method='pearson'),
                                                     cor_lag7 = cor(vmPFC_within,vmPFC_lag7,use='complete.obs',method='pearson'),
                                                     cor_lag8 = cor(vmPFC_within,vmPFC_lag8,use='complete.obs',method='pearson'),
                                                     cor_lag9 = cor(vmPFC_within,vmPFC_lag9,use='complete.obs',method='pearson'),
                                                     cor_lead1 = cor(vmPFC_within,vmPFC_lead1,use='complete.obs',method='pearson'),
                                                     cor_lead2 = cor(vmPFC_within,vmPFC_lead2,use='complete.obs',method='pearson'),
                                                     cor_lead3 = cor(vmPFC_within,vmPFC_lead3,use='complete.obs',method='pearson'),
                                                     cor_lead4 = cor(vmPFC_within,vmPFC_lead4,use='complete.obs',method='pearson'),
                                                     cor_lead5 = cor(vmPFC_within,vmPFC_lead5,use='complete.obs',method='pearson'),
                                                     cor_lead6 = cor(vmPFC_within,vmPFC_lead6,use='complete.obs',method='pearson'),
                                                     cor_lead7 = cor(vmPFC_within,vmPFC_lead7,use='complete.obs',method='pearson'),
                                                     cor_lead8 = cor(vmPFC_within,vmPFC_lead8,use='complete.obs',method='pearson'),
                                                     cor_lead9 = cor(vmPFC_within,vmPFC_lead9,use='complete.obs',method='pearson'),
                                                     cor_lead10 = cor(vmPFC_within,vmPFC_lead10,use='complete.obs',method='pearson'),
                                                     cor_lead11 = cor(vmPFC_within,vmPFC_lead11,use='complete.obs',method='pearson'),
                                                     cor_lead12 = cor(vmPFC_within,vmPFC_lead12,use='complete.obs',method='pearson'),
                                                     cor_lead13 = cor(vmPFC_within,vmPFC_lead13,use='complete.obs',method='pearson'),
                                                     cor_lead14 = cor(vmPFC_within,vmPFC_lead14,use='complete.obs',method='pearson'),
                                                     cor_lead15 = cor(vmPFC_within,vmPFC_lead15,use='complete.obs',method='pearson')
  ) %>% ungroup()
  vmPFC <- vmPFC %>% mutate(HC_region='NA')
  vPFC_Exp <- vmPFC %>% mutate(dataset='vPFC')
  
  
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock.Rdata')
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
  hc <- hc %>% rename(HC_decon = decon_mean)
  hc <- hc %>% group_by(id,run,trial,evt_time,HC_region) %>% summarize(decon1 = mean(HC_decon,na.rm=TRUE)) %>% ungroup() # 12 -> 2
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory='/Volumes/Users/Andrew/MEDuSA_data_Explore',dataset='explore')
  df <- df %>%
    group_by(id,run_number) %>%
    arrange(id, run_number, trial) %>%
    mutate(condition_trial = run_trial-floor(run_trial/40.5)*40,
           condition_trial_neg_inv = -1000 / condition_trial,
           condition_trial_neg_inv_sc = as.vector(scale(condition_trial_neg_inv)),
    ) %>% ungroup()
  df <- df %>%
    group_by(id) %>% 
    mutate(v_entropy_sc_r = scale(v_entropy)) %>% ungroup() %>%
    group_by(id, run) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           score_sc = scale(score_csv),
           iti_sc = scale(iti_ideal),
           iti_lag_sc = lag(iti_sc),
           v_chosen_sc = scale(v_chosen),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           rt_vmax_sc = scale(rt_vmax),
           ev_sc = scale(ev),
           v_entropy_sc = scale(v_entropy),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
           abs_rt_vmax_change = abs(rt_vmax_change),
           rt_vmax_change_bin = case_when(
             abs_rt_vmax_change < 4/24 ~ 'No Change',
             abs_rt_vmax_change >= 4/24 ~ 'Change'
           ),
           v_entropy_wi_change_bin = case_when(
             v_entropy_wi_change < -0.5 ~ 'Decrease',
             v_entropy_wi_change > 0.5  ~ 'Increase',
             v_entropy_wi_change >= -0.5 & v_entropy_wi_change <= 0.5 ~ 'No Change'
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
    run_trial <= 13 ~ 'Early',
    run_trial > 13 & run_trial < 26 ~ 'Middle',
    run_trial >=26 ~ 'Late',
  )))
  df <- df %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                       total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))
  
  df <- df %>% select(total_earnings_split,condition_trial_neg_inv_sc,iti_ideal,rt_lag_sc,iti_lag_sc,iti_prev,iti_sc,v_entropy_wi_change_lag,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
  df$id <- as.character(df$id)
  hc <- full_join(hc,df,by=c('id','run','trial'))
  hc <- hc %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                              trial > 40 & trial <= 80 ~ 2,
                                              trial > 80 & trial <=120 ~ 3, 
                                              trial > 120 & trial <=160 ~ 4,
                                              trial > 160 & trial <=200 ~ 5,
                                              trial > 200 & trial <=240 ~ 6))
  hc <- hc %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                                   trial > 40 & trial <= 80 ~ trial-40,
                                                   trial > 80 & trial <=120 ~ trial-80, 
                                                   trial > 120 & trial <=160 ~ trial-120,
                                                   trial > 160 & trial <=200 ~ trial-160,
                                                   trial > 200 & trial <=240 ~ trial-200))
  hc <- hc %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                            run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                            run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))
  hc$HCwithin[hc$evt_time > hc$rt_csv + hc$iti_ideal] = NA
  #Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  #Q$HC_decon[Q$evt_time < -Q$iti_prev] = NA
  hc$HCwithin[hc$evt_time < -hc$iti_prev] = NA
  hc <- hc %>% arrange(id,run,trial,evt_time)
  
  
  hc$rt_bin <- relevel(as.factor(hc$rt_bin),ref='-0.5')
  hc$trial_bin <- relevel(as.factor(hc$trial_bin),ref='Middle')
  
  demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
  demo <- demo %>% select(registration_redcapid,age,gender,registration_group,wtar,education_yrs)
  demo <- demo %>% rename(id=registration_redcapid,group=registration_group)
  demo$gender <- relevel(as.factor(demo$gender),ref='M')
  demo$age <- scale(demo$age)
  demo$wtar <- scale(demo$wtar)
  demo$education_yrs <- scale(demo$education_yrs)
  
  hc <- merge(demo,hc,by='id')
  #Q <- Q %>% filter(total_earnings_split=='richer')
  hc$group <- relevel(factor(hc$group),ref='HC')
  hc <- hc %>% filter(group!='ATT')
  hc <- hc %>% filter(group=='HC')
  hc <- hc %>% filter(!is.na(rewFunc))
  #Q <- Q %>% filter(trial > 10)
  
  # hc1 <- hc %>% group_by(id,run,trial,HC_region) %>% 
  #   select(id,run,trial,evt_time,HCwithin,HC_region) %>% 
  #   pivot_wider(values_from=HCwithin,names_from=evt_time)
  # acfHC <- hc1 %>% 
  #   group_by(HC_region) %>% select(!run & !trial & !id) %>%
  #   nest() %>% 
  #   mutate(data = map(data, ~acf(., lag.max=5, type="correlation", plot=F,na.action=na.pass))) %>%
  #   mutate(data = map(data, ~as.data.frame(rbind(.x$acf[1,,],.x$acf[2,,],.x$acf[3,,],.x$acf[4,,],.x$acf[5,,],.x$acf[6,,])))) %>%
  #   unnest(data) %>% ungroup()
  # 
  # acfHC <- acfHC %>%rename('-3.6'=V1,'-3'=V2,'-2.4'=V3,'-1.8'=V4,'-1.2'=V5,'-0.6'=V6,'0'=V7,'0.6'=V8,'1.2'=V9,'1.8'=V10,'2.4'=V11,'3'=V12,'3.6'=V13) %>% 
  #   mutate(lags1 = rep(c(t(rep(0,13)),t(rep(1,13)),t(rep(2,13)),t(rep(3,13)),t(rep(4,13)),t(rep(5,13))),2),
  #          lags2 = rep(c(t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1)),t(round(seq(-3.6,3.6,length.out=13),1))),2)
  #   )
  # acfHC <- acfHC %>% pivot_longer(cols=c("-3.6","-3","-2.4","-1.8","-1.2","-0.6","0","0.6","1.2","1.8","2.4","3","3.6"))
  # acfHC_Exp <- acfHC %>% rename(acf=value,evt_time=name)
  # acfHC_Exp$evt_time <- as.numeric(acfHC_Exp$evt_time)
  # HC_Exp0 <- acfHC_Exp %>% filter(evt_time==0 & ((lags1==0 & lags2==0) | 
  #                                                  (lags1==1 & lags2==-0.6) | (lags1==1 & lags2==0.6) | 
  #                                                  (lags1==2 & lags2==-1.2) | (lags1==2 & lags2==1.2) |
  #                                                  (lags1==3 & lags2==-1.8) | (lags1==3 & lags2==1.8) |
  #                                                  (lags1==4 & lags2==-2.4) | (lags1==4 & lags2==2.4) |
  #                                                  (lags1==5 & lags2==-3) | (lags1==5 & lags2==3) |
  #                                                  (lags1==6 & lags2==-3.6) | (lags1==6 & lags2==3.6))
  # )
  hc <- hc %>% group_by(id,run,trial,HC_region) %>% 
    mutate(hc_lag1 = lag(HCwithin,1), 
           hc_lag2 = lag(HCwithin,2), 
           hc_lag3 = lag(HCwithin,3),
           hc_lag4 = lag(HCwithin,4),
           hc_lag5 = lag(HCwithin,5),
           hc_lag6 = lag(HCwithin,6),
           hc_lag7 = lag(HCwithin,7),
           hc_lag8 = lag(HCwithin,8),
           hc_lag9 = lag(HCwithin,9),
           hc_lag10 = lag(HCwithin,10),
           hc_lead1 = lead(HCwithin,1),
           hc_lead2 = lead(HCwithin,2),
           hc_lead3 = lead(HCwithin,3),
           hc_lead4 = lead(HCwithin,4),
           hc_lead5 = lead(HCwithin,5),
           hc_lead6 = lead(HCwithin,6),
           hc_lead7 = lead(HCwithin,7),
           hc_lead8 = lead(HCwithin,8),
           hc_lead9 = lead(HCwithin,9),
           hc_lead10 = lead(HCwithin,10),
           hc_lead11 = lead(HCwithin,11),
           hc_lead12 = lead(HCwithin,12),
           hc_lead13 = lead(HCwithin,13),
           hc_lead14 = lead(HCwithin,14),
           hc_lead15 = lead(HCwithin,15)
           
    ) %>% ungroup() %>% filter(evt_time==0)
  
  hc <- hc %>% group_by(HC_region) %>% summarize(cor_lag0 = cor(HCwithin,HCwithin,use='complete.obs',method='pearson'),
                                                     cor_lag1 = cor(HCwithin,hc_lag1,use='complete.obs',method='pearson'),
                                                     cor_lag2 = cor(HCwithin,hc_lag2,use='complete.obs',method='pearson'),
                                                     cor_lag3 = cor(HCwithin,hc_lag3,use='complete.obs',method='pearson'),
                                                     cor_lag4 = cor(HCwithin,hc_lag4,use='complete.obs',method='pearson'),
                                                     cor_lag5 = cor(HCwithin,hc_lag5,use='complete.obs',method='pearson'),
                                                     cor_lag6 = cor(HCwithin,hc_lag6,use='complete.obs',method='pearson'),
                                                     cor_lag7 = cor(HCwithin,hc_lag7,use='complete.obs',method='pearson'),
                                                     cor_lag8 = cor(HCwithin,hc_lag8,use='complete.obs',method='pearson'),
                                                     cor_lag9 = cor(HCwithin,hc_lag9,use='complete.obs',method='pearson'),
                                                     cor_lead1 = cor(HCwithin,hc_lead1,use='complete.obs',method='pearson'),
                                                     cor_lead2 = cor(HCwithin,hc_lead2,use='complete.obs',method='pearson'),
                                                     cor_lead3 = cor(HCwithin,hc_lead3,use='complete.obs',method='pearson'),
                                                     cor_lead4 = cor(HCwithin,hc_lead4,use='complete.obs',method='pearson'),
                                                     cor_lead5 = cor(HCwithin,hc_lead5,use='complete.obs',method='pearson'),
                                                     cor_lead6 = cor(HCwithin,hc_lead6,use='complete.obs',method='pearson'),
                                                     cor_lead7 = cor(HCwithin,hc_lead7,use='complete.obs',method='pearson'),
                                                     cor_lead8 = cor(HCwithin,hc_lead8,use='complete.obs',method='pearson'),
                                                     cor_lead9 = cor(HCwithin,hc_lead9,use='complete.obs',method='pearson'),
                                                     cor_lead10 = cor(HCwithin,hc_lead10,use='complete.obs',method='pearson'),
                                                     cor_lead11 = cor(HCwithin,hc_lead11,use='complete.obs',method='pearson'),
                                                     cor_lead12 = cor(HCwithin,hc_lead12,use='complete.obs',method='pearson'),
                                                     cor_lead13 = cor(HCwithin,hc_lead13,use='complete.obs',method='pearson'),
                                                     cor_lead14 = cor(HCwithin,hc_lead14,use='complete.obs',method='pearson'),
                                                     cor_lead15 = cor(HCwithin,hc_lead15,use='complete.obs',method='pearson')
  ) %>% ungroup()
  hc <- hc %>% mutate(network='NA')
  HC_Exp <- hc %>% mutate(dataset='HC')
  
  dq_Exp <- rbind(vPFC_Exp,HC_Exp)
}

dq_Mmc1 <- dq_Mmc %>% pivot_longer(cols=starts_with('cor')) %>%
  mutate(lags = case_when(name=='cor_lag0'~0,
                          name=='cor_lag1'~-1,
                          name=='cor_lag2'~-2,
                          name=='cor_lag3'~-3,
                          name=='cor_lag4'~-4,
                          name=='cor_lag5'~-5,
                          name=='cor_lag6'~-6,
                          name=='cor_lead1'~1,
                          name=='cor_lead2'~2,
                          name=='cor_lead3'~3,
                          name=='cor_lead4'~4,
                          name=='cor_lead5'~5,
                          name=='cor_lead6'~6,
                          name=='cor_lead7'~7,
                          name=='cor_lead8'~8,
                          name=='cor_lead9'~9,
                          ))
dq_Exp1 <- dq_Exp %>% pivot_longer(cols=starts_with('cor')) %>%
  mutate(lags = case_when(name=='cor_lag0'~0,
                          name=='cor_lag1'~-1,
                          name=='cor_lag2'~-2,
                          name=='cor_lag3'~-3,
                          name=='cor_lag4'~-4,
                          name=='cor_lag5'~-5,
                          name=='cor_lag6'~-6,
                          name=='cor_lag6'~-7,
                          name=='cor_lag6'~-8,
                          name=='cor_lag6'~-9,
                          name=='cor_lead1'~1,
                          name=='cor_lead2'~2,
                          name=='cor_lead3'~3,
                          name=='cor_lead4'~4,
                          name=='cor_lead5'~5,
                          name=='cor_lead6'~6,
                          name=='cor_lead7'~7,
                          name=='cor_lead8'~8,
                          name=='cor_lead9'~9,
                          name=='cor_lead9'~10,
                          name=='cor_lead9'~11,
                          name=='cor_lead9'~12,
                          name=='cor_lead9'~13,
                          name=='cor_lead9'~14,
                          name=='cor_lead9'~15
  ))


