# 2021-05-25 AndyP
# HC and vmPFC regressor interactions / MEDuSA analyses
# need to redo for RT prediction, due to overinflation of RT's, average within network.


validate_mixed_by_vmPFC_HC <- function(ncores=6,toalign='feedback',toprocess='network',totest='base',behavmodel='compressed',
                                       repo_directory="~/clock_analysis",vmPFC_cache_dir='~/vmPFC/MEDUSA Schaefer Analysis',
                                       HC_cache_dir='~/vmPFC/MEDUSA Schaefer Analysis',
                                       hc_LorR='R',reload=TRUE,replot=TRUE,load_existing=FALSE,ddf_file='',Q=NULL){
  
  
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(pracma)
  library(broom)
  library(grid)
  
  dolmer = 'false'
  #toalign <- 'clock'
  #toprocess <- 'network-by-HC' #'symmetry-by-bin'
  #totest <- 'base' #'anatomy' #'base-reversed' #'expl_code' #nextrt, interactions, base, pe, expl_code
  #behavmodel <- 'compressed' # uncompressed, # valuecomp
  dovif = FALSE
  lagHC <- FALSE
  lagvmPFC <- FALSE
  #hc_LorR <- 'L'
  
  #repo_directory <- "~/clock_analysis"
  #vmPFC_cache_dir = "/Users/andypapale/Box/SCEPTIC_fMRI/vmPFC/cache"
  #HC_cache_dir = "/Users/andypapale/vmPFC/deconvolved_evt_locked_smooth_in_mask_harvardoxford"
  #HC_cache_dir = '/Users/andypapale/vmPFC/MEDUSA Schaefer Analysis'
  #vmPFC_cache_dir = '/Users/andypapale/vmPFC/MEDUSA Schaefer Analysis'
  #gc()
  
  # load MEDUSA deconvolved data
  if (load_existing==FALSE){
    if (reload==TRUE){
      if (strcmp(toalign,'clock')){
        message("Loading HC medusa data from cache: ", HC_cache_dir)
        load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
        hc <- clock_comb
        hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
        rm(clock_comb)
        message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
        load(file.path(vmPFC_cache_dir, 'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
        vmPFC <- clock_comb
        vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
        rm(clock_comb)
      } else if (strcmp(toalign,'feedback')){
        load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
        hc <- fb_comb
        hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
        rm(fb_comb)
        message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
        load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
        vmPFC <- fb_comb
        vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
        rm(fb_comb)
      } else if (strcmp(toalign,'rt') | strcmp(toalign,'rt_vmax')){
        display('currently not working')
      }
      hc <- hc %>% select(id,run,run_trial,decon_mean,evt_time,bin_num,side)
      if (strcmp(hc_LorR,'L')){
        hc <- hc %>% filter(side=='l')
      } else if (strcmp(hc_LorR,'R')){
        hc <- hc %>% filter(side=='r')
      } else if (strcmp(hc_LorR,'LR')){
        # do nothing, summarize will collapse to appropriate size
      }
      if (!strcmp(totest,'anatomy')){
        hc <- hc %>% mutate(
          HC_region = case_when(
            bin_num <= 8 ~ 'AH',
            bin_num >8 ~ 'PH'
          ),
        )
      }
      # change to mean
      if (!strcmp(totest,'anatomy')){
        hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
        if (strcmp(totest,'base') | strcmp(totest,'anatomy') | strcmp(totest,'base-wPE') | strcmp(totest,'simple-model')){ #HC is a predictor variable, scale
          hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
          hc <- hc %>% select(!decon1)
          if (lagHC==TRUE){ # test lagged HC effect on vmPFC
            hc <- hc %>% group_by(id,run,run_trial,evt_time,bin) %>% mutate(HCw2 = lag(HCwithin), HCb2=lag(HCbetween))
            hc <- hc %>% select(!HCwithin & !HCbetween) %>% rename(HCwithin=HCw2,HCbetween=HCb2)
          }
        } else if (strcmp(totest,'base-reversed')){ # HC is an outcome variable for vmPFC activity, do not scale
          hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = decon1, HCbetween = mean(decon1,na.rm=TRUE)) %>% ungroup()
          hc <- hc %>% select(!decon_interp & decon1)
          if (lagHC==TRUE){ # test lagged HC effect on vmPFC
            hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% mutate(HCw2 = lag(HCwithin), HCb2=lag(HCbetween))
            hc <- hc %>% select(!HCwithin & !HCbetween) %>% rename(HCwithin=HCw2,HCbetween=HCb2)
          }
        } else if (strcmp(totest,"nextrt") | strcmp(totest,"rtvmax") | strcmp(totest,"entropy_change") | strcmp(totest,"vmax") | strcmp(totest,"PE") | strcmp(totest,'H')) { # testing behavioral variables, compress
          hc <- hc %>% group_by(id,run) %>% mutate(HCwithin=scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
          hc <- hc %>% select(!decon1)
          if (lagHC==TRUE){ # test lagged HC effect on vmPFC
            hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% mutate(HCw2 = lag(HCwithin), HCb2=lag(HCbetween))
            hc <- hc %>% select(!HCwithin & !HCbetween) %>% rename(HCwithin=HCw2,HCbetween=HCb2)
          }
        } else {
          warning('unknown totest option')
        }
      } else if (strcmp(totest,'anatomy')){
        hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon_mean),HCbetween=mean(decon_mean,na.rm=TRUE)) %>% ungroup()
        hc <- hc %>% select(!decon_mean)
      }
      
      #vmPFC$AP <-vmPFC$'A/P'
      #vmPFC$ML <-vmPFC$'M/L'
      #vmPFC$IS <-vmPFC$'I/S'
      
      # 2021-08-30 AndyP rewrote the vmPFC call
      # test behavior variable, compress
      # predictor variables scale
      vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
      vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
      # if (strcmp(totest,'base') | strcmp(totest,'base-wPE')){ # vmPFC is an outcome variable, do not scale
      #   if (strcmp(toprocess,'network-by-HC')){
      #     vmPFC <- vmPFC %>% mutate(vmPFC_decon = decon_interp)
      #   } else if (strcmp(toprocess,'symmetry-by-HC')){
      #     vmPFC <- vmPFC %>% group_by(id,run,run_trial,evt_time,symmetry_group) %>% summarize(vmPFC_decon = mean(decon_interp,na.rm=TRUE)) %>% ungroup() %>% arrange(id,run,run_trial,evt_time)
      #   }
      #   vmPFC <- vmPFC %>% select(!atlas_value)
      # } else if (strcmp(totest,'base-reversed')){ # vmPFC is a predictor variable for HC activity, scale
      #   vmPFC <- vmPFC %>% group_by(id,run) %>% mutate(vmPFC_decon= scale(decon_interp), vmPFCbetween1=mean(decon_interp,na.rm=TRUE)) %>% ungroup()
      #   if (strcmp(toprocess,'network-by-HC')){
      #     vmPFC <- vmPFC %>% group_by(id,run,run_trial,evt_time,network) %>% summarize(vmPFC_decon_scaled = mean(vmPFC_decon,na.rm=TRUE), vmPFCbetween=mean(vmPFCbetween1,na.rm=TRUE)) %>% ungroup() %>% arrange(id,run,run_trial,evt_time)
      #   } else if (strcmp(toprocess,'symmetry-by-HC')){
      #     vmPFC <- vmPFC %>% group_by(id,run,run_trial,evt_time,symmetry_group) %>% summarize(vmPFC_decon_scaled = mean(vmPFC_decon,na.rm=TRUE), vmPFCbetween=mean(vmPFCbetween1,na.rm=TRUE)) %>% ungroup() %>% arrange(id,run,run_trial,evt_time)
      #   }
      #   vmPFC <- vmPFC %>% select(!atlas_value & !vmPFC_decon & !vmPFCbetween1)
      # } else if  (strcmp(totest,"nextrt") | strcmp(totest,"rtvmax") | strcmp(totest,"entropy_change") | strcmp(totest,"vmax")) { # we are testing behavior variables, compress
      #   vmPFC <- vmPFC %>% group_by(id,run) %>% mutate(vmPFC_decon_scaled1 = scale(decon_interp), vmPFCbetween1=mean(decon_interp,na.rm=TRUE)) %>% ungroup()
      #   vmPFC <- vmPFC %>% mutate(vmPFC_decon_scaled = vmPFC_decon_scaled1)
      #   vmPFC <- vmPFC %>% select(!vmPFC_decon_scaled1 & !vmPFCbetween1)
      # } else {
      #   warning('unknown totest option')
      # }
      vmPFC <- vmPFC %>% mutate(side_vmpfc = as.factor(substr(region, nchar(region), nchar(region))))
      #hc$side_hc <- as.factor(hc$side)
      #hc <- hc %>% select(!side)
      if (strcmp(totest,"nextrt") | strcmp(totest,'PE') | strcmp(totest,'H')){
        vmPFC <- vmPFC %>% group_by(id,run) %>% mutate(vmPFC_decon_scaled1 = scale(vmPFC_decon), vmPFC_between1 = mean(vmPFC_decon,na.rm=TRUE))
        vmPFC <- vmPFC %>% group_by(id,run,run_trial,network,side_vmpfc,evt_time) %>% summarize(vmPFC_decon_scaled = mean(vmPFC_decon_scaled1,na.rm=TRUE),vmPFC_between=mean(vmPFC_between1,na.rm=TRUE))
      }
      
      #vmPFC <- vmPFC %>% filter(evt_time < 8)
      #vmPFC <- vmPFC %>% filter(network == 'C')
      #hc <- hc %>% filter(evt_time < 8)
      
      merged_data <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
      rm(vmPFC)
      rm(hc)
      if (lagHC==TRUE){
        if (strcmp(toalign,"clock")){
          merged_data <- merged_data %>% filter(!evt_time==-1)
        } else if (strcmp(toalign,"feedback")){
          merged_data <- merged_data %>% filter(!evt_time==-1)
        }
      }
      
      # load behavioral data, get lags, merge with signal
      if (!strcmp(totest,'anatomy')){
        if (strcmp(behavmodel,'compressed')){
          source('~/vmPFC/get_trial_data_vmPFC.R')
          df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
          df <- df %>% select(v_entropy,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                              v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                              rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
            group_by(id, run) %>% 
            mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
            mutate(v_chosen_sc = scale(v_chosen),
                   abs_pe_max_sc = scale(abs(pe_max)),
                   abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
                   rt_vmax_sc = scale(rt_vmax),
                   rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
                     kld3 ==0 ~ NA_real_,
                     kld3 >0 ~ log10(kld3)
                   )) %>% mutate(log10kld3_lag = case_when(
                     kld3_lag==0 ~NA_real_,
                     kld3_lag>0 ~ log10(kld3_lag)
                   ))
        } else if (strcmp(behavmodel,'uncompressed')){
          source('~/vmPFC/get_trial_data_vmPFC.R')
          df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
          df <- df %>% select(rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                              v_max_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                              rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,v_entropy_wi_full,kld3_lag,rt_lag_sc,v_entropy_wi_change_full,
                              rt_vmax_full,pe_max_full)
          df <- df %>% arrange(id,run,run_trial)
          df$v_entropy_wi <- df$v_entropy_wi_full # NOTE CHANGE OF VARIABLE
          df$v_entropy_wi_change <- df$v_entropy_wi_change_full # NOTE CHANGE OF VARIABLE
          df <- df %>% group_by(id,run) %>% arrange(id, run, run_trial) %>% 
            mutate(v_chosen_sc = scale(v_chosen),
                   abs_pe_max_sc = scale(abs(pe_max_full)), # NOTE CHANGE OF VARIABLE
                   abs_pe_max_lag_sc = scale(abs(lag(pe_max_full))), # NOTE CHANGE OF VARIABLE
                   rt_vmax_sc = scale(rt_vmax_full)) %>% # NOTE CHANGE OF VARIABLE
            arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
              kld3 ==0 ~ NA_real_,
              kld3 >0 ~ log10(kld3)
            )) %>% mutate(log10kld3_lag = case_when(
              kld3_lag==0 ~NA_real_,
              kld3_lag>0 ~ log10(kld3_lag)
            ))
        }
      }
      Q <- merge(df, merged_data, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
      rm(merged_data)
      
      # # censor prior ITI on clock
      # if (strcmp(toalign,'clock')){
      #   Q$vmPFC_decon[Q$evt_time < -(Q$iti_lag+Q$rt_lag/1000)] <- NA
      #   Q$HCwithin[Q$evt_time < -(Q$iti_lag+Q$rt_lag/1000)] <- NA
      #   Q$HCbetween[Q$evt_time < -(Q$iti_lag+Q$rt_lag/1000)] <- NA
      # }
      
      if (!strcmp(totest,'anatomy')){
        Q <- Q %>% mutate(richorpoor = case_when(total_earnings < median(total_earnings) ~ 'P',
                                                 total_earnings >= median(total_earnings) ~ 'R'))
      }
      #Q$expl_code <- relevel(as.factor(Q$expl_code),ref="Choose-Stay")
      Q$last_outcome <- relevel(as.factor(Q$last_outcome),ref="Omission")
      Q$outcome <- relevel(as.factor(Q$outcome),ref="Omission")
      
      # depends on what you want to test
      #message('change Q select depending on what you want to test')
      #Q <- Q %>% select(vmPFC_decon,trial_neg_inv_sc,HCwithin,rt_csv_sc,v_entropy_wi,v_entropy_wi_change,kld3,v_max_wi,abs_pe_max_sc,abs_pe_max_lag_sc,HCbetween,id,run,side_vmpfc,region,evt_time,network,HC_region,side_hc)
    }
    
    
    if (strcmp(dolmer,'true')){
      #df$rewFunc <- relevel(as.factor(df$rewFunc),ref = "CEV")
      Q$network <- relevel(as.factor(Q$network),ref="L")
      
      ################################
      # clock aligned data -ONLINE  #
      ###############################
      
      
      # DMN
      ph <- Q %>% filter(HC_region=="PH" & online==TRUE) %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfcfp10L")
      ah <- Q %>% filter(HC_region=="AH" & online==TRUE) %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfcfp10L")
      
      # CTRL
      ph <- Q %>% filter(HC_region=="PH" & online==TRUE) %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfc11L")
      ah <- Q %>% filter(HC_region=="AH" & online==TRUE) %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfc11L")
      
      # LIM
      ph <- Q %>% filter(HC_region=="PH" & online==TRUE) %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfcfp10R")
      ah <- Q %>% filter(HC_region=="AH" & online==TRUE) %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfcfp10R")
      
      ############################
      # feedback aligned data   #
      ###########################
      
      # DMN
      ph <- Q %>% filter(HC_region=="PH") %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfcfp10L")
      ah <- Q %>% filter(HC_region=="AH") %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfcfp10L")
      
      # CTRL
      ph <- Q %>% filter(HC_region=="PH") %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfc11L")
      ah <- Q %>% filter(HC_region=="AH") %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfc11L")
      
      # LIM
      ph <- Q %>% filter(HC_region=="PH") %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfcfp10R")
      ah <- Q %>% filter(HC_region=="AH") %>% group_by(id,run,run_trial,network) %>% mutate(mvmPFCw = scale(decon_interp), mHCw = mean(HCwithin,na.rm=TRUE), mvmPFCb = mean(decon_interp,na.rm=TRUE),mHCb = mean(HCbetween,na.rm=TRUE)) %>% ungroup() %>% filter(evt_time==0 & region=="pfcfp10R")
      
      # lm4 <- lmer(decon_interp  ~ HCinterp*network + 
      #               v_entropy_wi*HCinterp*network + 
      #               v_entropy_wi_change*HCinterp*network + 
      #               kld3_lag*HCinterp*network +
      #               v_max_wi*HCinterp*network +
      #               scale(pe_max)*HCinterp*network +
      #               as.factor(evt_time)*network +
      #               (1|region) +
      #               (1|id), merged_data)
      
      
      lm4 <- lmer(vmPFCinterp ~ HCwithin*scale(abs_pe) +
                    HCbetween*scale(abs_pe) +
                    HCwithin*trial_neg_inv_sc +
                    HCbetween*trial_neg_inv_sc +
                    HCwithin*rt_vmax_lag_sc +
                    HCbetween*rt_vmax_lag_sc +
                    HCwithin*v_entropy_wi +
                    HCbetween*v_entropy_wi +
                    as.factor(evt_time)*network +
                    (1|region) +
                    (1|id),Q)
      lm4 <- tidy(lm4)
      
    }
    
    if (strcmp(toalign,"clock")){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction')
    } else if (strcmp(toalign,"feedback")){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC_interaction')
    } else if (strcmp(toalign,"rt")){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_rt_HC_interaction')
    }
    
    if (strcmp(dolmer,'true')){
      fname = paste(totest,"_",toalign, ".Rdata", sep = "")
      save(lm4,file=fname)
    }
    
    # mixed_by call split on, run separate models
    if (strcmp(toprocess,'network')){
    } else if (strcmp(toprocess,'symmetry-by-bin')){
      splits = c("evt_time","symmetry_group","bin_num")
    } else if (strcmp(toprocess,"symmetry_group")){
      splits = c("evt_time","symmetry_group")
    } else if (strcmp(toprocess,"region")){
      splits = c("evt_time","region")
    } else if (strcmp(toprocess,"side1")){
      splits = c("evt_time","side","network")
    } else if (strcmp(toprocess,"side2")){
      splits = c("evt_time","side")
    } else if (strcmp(toprocess,"network-by-HC") | strcmp(toprocess,'network-by-HClag') | strcmp(toprocess,"network-by-HC-wEC") | strcmp(toprocess,"network-by-HC-wO")){
      splits = c("evt_time","network","HC_region")
    } else if (strcmp(toprocess,"region-by-HC")){
      splits = c("evt_time","region","HC_region")
    } else if (strcmp(toprocess,"IS-by-HC")){
      splits = c("evt_time","IS","HC_region")
    } else if (strcmp(toprocess,"symmetry-by-HC")){
      splits = c("evt_time","symmetry_group","HC_region")
    } else if (strcmp(toprocess,'symmetry-group-by-HC-by-outcome')){
      splits = c("evt_time","symmetry_group","HC_region","outcome")
    } else if (strcmp(toprocess,'network-by-HC-by-outcome')){
      splits = c("evt_time","network","HC_region","outcome")
    } else if (strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
      splits = c("evt_time","symmetry_group","HC_region","rewFunc")
    } else if (strcmp(toprocess,"network-by-HC-by-side")){
      splits = c("evt_time","network","HC_region","side_hc")
    }
    if (strcmp(totest,"anatomy")){
      if (strcmp(toprocess,"network-by-bin")){
        if (strcmp(toalign,"feedback")){
          decode_formula = formula(~ HCwithin*outcome + HCbetween +
                                     (1|id/run) + (1| region) + (1|side_hc))
        } else if (strcmp(toalign,"clock")){
          decode_formula = formula(~ HCwithin*last_outcome + HCbetween +
                                     (1|id/run) + (1| region) + (1|side_hc))      
        }
      } else if (strcmp(toprocess,"symmetry-by-bin")){
        if (strcmp(toalign,"feedback")){
          decode_formula = formula(~ HCwithin*outcome + HCbetween +
                                     (1|id/run) + (1|side_hc))
        } else if (strcmp(toalign,"clock")){
          decode_formula = formula(~ HCwithin*last_outcome + HCbetween +
                                     (1|id/run) + (1|side_hc))      
        }
      }
    }
    if (strcmp(totest,"nextrt")){
      if (strcmp(toprocess,"network-by-HC")){
        if (strcmp(toalign,"feedback")){
          decode_formula = formula(~ vmPFC_decon_scaled*rt_csv_sc*HCwithin +  
                                     vmPFC_decon_scaled*rt_vmax_lag_sc*HCwithin + 
                                     vmPFC_decon_scaled*rt_csv_sc*HCwithin*outcome +
                                     HCbetween +  (1|id/run) + (1|region) + (1|side)) # for RT- or feedback-aligned
        } else{
          decode_formula = formula(~ vmPFC_decon_scaled*HCwithin*rt_lag_sc + 
                                     vmPFC_decon_scaled*HCwithin*rt_vmax_lag_sc + 
                                     vmPFC_decon_scaled*rt_lag_sc*HCwithin*last_outcome +
                                     HCbetween + (1|id/run) + (1|region) + (1|side)) # for clock-aligned
        }
      }
      if (strcmp(toprocess,"symmetry-by-HC")){
        if (strcmp(toalign,"feedback")){
          decode_formula = formula(~ vmPFC_decon_scaled*rt_csv_sc*HCwithin +  
                                     vmPFC_decon_scaled*rt_vmax_lag_sc*HCwithin + 
                                     vmPFC_decon_scaled*rt_csv_sc*HCwithin*outcome +
                                     HCbetween +  (1|id/run) + (1|side)) # for RT- or feedback-aligned
        } else{
          decode_formula = formula(~ vmPFC_decon_scaled*HCwithin*rt_lag_sc + 
                                     vmPFC_decon_scaled*HCwithin*rt_vmax_lag_sc + 
                                     vmPFC_decon_scaled*rt_lag_sc*HCwithin*last_outcome +
                                     HCbetween + (1|id/run) + (1|side)) # for clock-aligned
        }
      }
    } 
    if (strcmp(totest,"base-reversed")){
      if (strcmp(toalign,'clock')){
        if (!strcmp(toprocess,"region") & !strcmp(toprocess,"symmetry-by-HC")){
          decode_formula = formula(~ trial_neg_inv_sc*vmPFC_decon_scaled + rt_csv_sc*vmPFC_decon_scaled + 
                                     rt_lag_sc*vmPFC_decon_scaled + rt_vmax_lag_sc*vmPFC_decon_scaled  + 
                                     scale(rt_vmax_change)*vmPFC_decon_scaled + 
                                     v_entropy_wi*vmPFC_decon_scaled + 
                                     v_entropy_wi_change*vmPFC_decon_scaled + 
                                     kld3_lag*vmPFC_decon_scaled  + v_max_wi*vmPFC_decon_scaled  + 
                                     scale(abs_pe)*vmPFC_decon_scaled + last_outcome*vmPFC_decon_scaled + 
                                     vmPFCbetween + (1|id/run) + (1|region))  
        } else {
          decode_formula = formula(~ trial_neg_inv_sc*vmPFC_decon_scaled + rt_csv_sc*vmPFC_decon_scaled +
                                     rt_lag_sc*vmPFC_decon_scaled + rt_vmax_lag_sc*vmPFC_decon_scaled  +
                                     scale(rt_vmax_change)*vmPFC_decon_scaled +
                                     v_entropy_wi*vmPFC_decon_scaled +
                                     v_entropy_wi_change*vmPFC_decon_scaled +
                                     kld3_lag*vmPFC_decon_scaled  + v_max_wi*vmPFC_decon_scaled  +
                                     scale(abs_pe)*vmPFC_decon_scaled + last_outcome*vmPFC_decon_scaled +
                                     vmPFCbetween + (1|id/run))
        }
      } else if (strcmp(toalign,'feedback')) {
        if (!strcmp(toprocess,"region") & !strcmp(toprocess,"symmetry-by-HC")){
          decode_formula = formula(~ trial_neg_inv_sc*vmPFC_decon_scaled + 
                                     expl_code*HCwithin + rt_vmax_lag_sc*vmPFC_decon_scaled  + 
                                     scale(rt_vmax_change)*vmPFC_decon_scaled + 
                                     v_entropy_wi*vmPFC_decon_scaled + 
                                     v_entropy_wi_change*vmPFC_decon_scaled + 
                                     kld3_lag*vmPFC_decon_scaled  + v_max_wi*vmPFC_decon_scaled  + 
                                     scale(abs_pe)*vmPFC_decon_scaled + outcome*vmPFC_decon_scaled + 
                                     vmPFCbetween + (1|id/run) + (1|region))    
        } else {
          decode_formula = formula(~ trial_neg_inv_sc*vmPFC_decon_scaled + rt_csv_sc*vmPFC_decon_scaled + 
                                     rt_lag_sc*vmPFC_decon_scaled + rt_vmax_lag_sc*vmPFC_decon_scaled  + 
                                     scale(rt_vmax_change)*vmPFC_decon_scaled + 
                                     v_entropy_wi*vmPFC_decon_scaled + 
                                     v_entropy_wi_change*vmPFC_decon_scaled + 
                                     kld3_lag*vmPFC_decon_scaled  + v_max_wi*vmPFC_decon_scaled  + 
                                     scale(abs_pe)*vmPFC_decon_scaled + outcome*vmPFC_decon_scaled + 
                                     vmPFCbetween + (1|id/run))       
        }
      }
    }
    if (strcmp(totest,'base')){
      if (strcmp(toprocess,'symmetry-by-HC')){
        if (strcmp(toalign,'clock')){
          decode_formula = formula(~ trial_neg_inv_sc +
                                     rt_lag_sc + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     v_max_wi*HCwithin  + log10kld3 +
                                     last_outcome*HCwithin +
                                     HCbetween + (1|id/run))  
        } else if (strcmp(toalign,'feedback')) {
          decode_formula = formula(~ trial_neg_inv_sc +
                                     rt_csv_sc + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     v_max_wi*HCwithin  + log10kld3 +
                                     outcome*HCwithin +
                                     HCbetween + (1|id/run))    
        }   
      } else if (strcmp(toprocess,'network-by-HC') | strcmp(toprocess,'network-by-HC-by-outcome')){
        if (strcmp(toalign,'clock')){
          decode_formula = formula(~ trial_neg_inv_sc +
                                     rt_lag_sc + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     v_max_wi*HCwithin  + log10kld3 +
                                     last_outcome*HCwithin +
                                     HCbetween + (1|id/run) + (1|region))  
        } else if (strcmp(toalign,'feedback')) {
          decode_formula = formula(~ trial_neg_inv_sc +
                                     rt_csv_sc + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     v_max_wi*HCwithin  + log10kld3 +
                                     outcome*HCwithin +
                                     HCbetween + (1|id/run) + (1|region))    
        }
      } else if (strcmp(toprocess,'symmetry-group-by-HC-by-outcome')){
        if (strcmp(toalign,'clock')){
          decode_formula = formula(~ rt_csv_sc*HCwithin + 
                                     rt_vmax_lag_sc*HCwithin + 
                                     scale(rt_vmax_change)*HCwithin + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     kld3_lag*HCwithin  + v_max_wi*HCwithin  + 
                                     scale(v_chosen)*HCwithin +
                                     HCbetween + (1|id/run))  
        } else if (strcmp(toalign,'feedback')) {
          decode_formula = formula(~ rt_csv_sc*HCwithin + 
                                     rt_vmax_lag_sc*HCwithin  + 
                                     scale(rt_vmax_change)*HCwithin + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     kld3*HCwithin  + v_max_wi*HCwithin  + 
                                     scale(v_chosen)*HCwithin +
                                     HCbetween + (1|id/run))    
        }  
      } else if (strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
        decode_formula = formula(~ rt_csv_sc*HCwithin + 
                                   rt_vmax_lag_sc*HCwithin  + 
                                   scale(rt_vmax_change)*HCwithin + 
                                   v_entropy_wi*HCwithin + 
                                   v_entropy_wi_change*HCwithin + 
                                   kld3*HCwithin  + v_max_wi*HCwithin  + 
                                   scale(v_chosen)*HCwithin + 
                                   scale(abs(pe_max))*HCwithin +
                                   HCbetween + (1|id))  
      } else if(strcmp(toprocess,'network-by-HC-by-side')){
        if (strcmp(toalign,'clock')){
          decode_formula = formula(~ trial_neg_inv_sc*HCwithin +
                                     rt_lag_sc*HCwithin + rt_vmax_lag_sc*HCwithin + 
                                     scale(rt_vmax_change)*HCwithin + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     kld3_lag*HCwithin  + v_max_wi*HCwithin  +
                                     abs_pe_max_lag_sc*HCwithin +
                                     HCbetween + (1|id/run) + (1|side_vmpfc) + (1|region))  
        } else if (strcmp(toalign,'feedback')) {
          decode_formula = formula(~ trial_neg_inv_sc*HCwithin +
                                     rt_csv_sc*HCwithin + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     kld3*HCwithin  + v_max_wi*HCwithin  + 
                                     abs_pe_max_lag_sc*HCwithin +
                                     abs_pe_max_sc*HCwithin +
                                     HCbetween + (1|id/run) + (1|side_vmpfc) + (1|region))    
        }
      }
    } else if (strcmp(totest,'simple-model')){
      decode_formula = formula(~v_entropy_wi*HCwithin + v_max_wi*HCwithin + HCbetween + (1|id/run) + (1|region) + (1|side))
    }
    
    if (strcmp(totest,'base-wPE')){
      if (strcmp(toprocess,'symmetry-by-HC')){
        if (strcmp(toalign,'clock')){
          decode_formula = formula(~ trial_neg_inv_sc +
                                     rt_lag_sc + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     v_max_wi  + log10kld3 +
                                     abs_pe_max_lag_sc*HCwithin +
                                     HCbetween + (1|id/run))  
        } else if (strcmp(toalign,'feedback')) {
          decode_formula = formula(~ trial_neg_inv_sc +
                                     rt_csv_sc + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     v_max_wi*HCwithin  + log10kld3 +
                                     abs_pe_max_sc*HCwithin +
                                     HCbetween + (1|id/run))    
        }   
      } else if (strcmp(toprocess,'network-by-HC')){
        if (strcmp(toalign,'clock')){
          decode_formula = formula(~ rt_csv_sc + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin + 
                                     log10kld3_lag  + v_max_wi*HCwithin  + 
                                     abs_pe_max_lag_sc*HCwithin + 
                                     HCbetween + (1|id/run) + (1|region))  
        } else if (strcmp(toalign,'feedback')) {
          decode_formula = formula(~ rt_csv_sc + 
                                     v_entropy_wi*HCwithin + 
                                     v_entropy_wi_change*HCwithin +
                                     log10kld3  + v_max_wi*HCwithin  + 
                                     abs_pe_max_sc*HCwithin +
                                     HCbetween + (1|id/run) + (1|region)) 
        }
      } 
    } 
    
    if (strcmp(totest,'PE')){
      decode_formula = formula(~ HCwithin*vmPFC_decon_scaled*trial_neg_inv_sc + HCbetween + vmPFC_between + (1|id/run) + (1|side_vmpfc))
    }
    if (strcmp(totest,"H")){
      decode_formula = formula(~ HCwithin*vmPFC_decon_scaled*trial_neg_inv_sc + HCbetween + vmPFC_between + (1|id/run) + (1|side_vmpfc))
    }
    
    source("~/fmri.pipeline/R/mixed_by.R")
    if (strcmp(totest,"nextrt")){
      if (strcmp(toalign,"clock")){
        ddf <- mixed_by(Q, outcomes = "rt_sec", rhs_model_formulae = decode_formula , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
      } else {
        ddf <- mixed_by(Q, outcomes = "rt_next", rhs_model_formulae = decode_formula , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
      }
    }
    if (strcmp(totest,"PE")){
      ddf <- mixed_by(Q, outcomes = "pe_max", rhs_model_formulae = decode_formula , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
    }
    if (strcmp(totest,"H")){
      ddf <- mixed_by(Q, outcomes = "v_entropy_wi", rhs_model_formulae = decode_formula , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
    }
    if (strcmp(totest,"base") | strcmp(totest,'base-wPE') | strcmp(totest,"anatomy")){
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = decode_formula , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
    }
    if (strcmp(totest,'simple-model')){
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = decode_formula , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
    }
    if (strcmp(totest,"rtvmax")){
      if (strcmp(toalign,"clock")){
        ddf <- mixed_by(Q, outcomes = "rt_vmax", rhs_model_formulae = decode_formula , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
      } else {
        ddf <- mixed_by(Q, outcomes = "rt_vmax_lead", rhs_model_formulae = decode_formula , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
      }
    }
    if (strcmp(totest,"entropy_change")){
      if (strcmp(toalign,"clock")){
        ddf <- mixed_by(Q, outcomes = "v_entropy_wi_change", rhs_model_formulae = decode_formula , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
      } else {
        ddf <- mixed_by(Q, outcomes = "v_entropy_wi_change_lead", rhs_model_formulae = decode_formula , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
      }
    }
    if (strcmp(totest,"vmax")){
      if (strcmp(toalign,"clock")){
        ddf <- mixed_by(Q, outcomes = "v_max_wi", rhs_model_formulae = decode_formula , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
      } else {
        ddf <- mixed_by(Q, outcomes = "v_max_wi_lead", rhs_model_formulae = decode_formula , split_on = splits,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
      }
    }
    if (strcmp(totest,"base-reversed")){
      ddf <- mixed_by(Q, outcomes = "HCwithin", rhs_model_formulae = decode_formula , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
    }
    
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-',behavmodel,'-',toalign,'-',totest,'-',toprocess,'-',hc_LorR,'.Rdata')) 
  }
  
  
  if (replot==TRUE){
    
    if (load_existing==TRUE){
      load(ddf_file)
    }
    
    ## Check plots
    message("\nPlotting streams decoding")
    library(viridis)
    library(wesanderson)
    pal = wes_palette("FantasticFox1", 3, type = "discrete")
    epoch_label = paste("Time relative to",toalign, "[s]")
    
    if (strcmp(toprocess,"network-by-outcome") | strcmp(toprocess,"symmetry-group-by-HC-by-outcome") | strcmp(toprocess,'network-by-HC-by-outcome')){
      colnames(ddf$coef_df_reml)[4] <- "outcome1"
    }
    ddg <- ddf
    ddf <- as_tibble(ddf$coef_df_reml)
    ddf$t <- ddf$evt_time
    #if (strcmp(toprocess,"network")){
    ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term, 
                           p_level_fdr = as.factor(case_when(
                             # p_fdr > .1 ~ '0',
                             # p_fdr < .1 & p_fdr > .05 ~ '1',
                             p_fdr > .05 ~ '1',
                             p_fdr < .05 & p_fdr > .01 ~ '2',
                             p_fdr < .01 & p_fdr > .001 ~ '3',
                             p_fdr <.001 ~ '4')))
    ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
    ddf$`p, FDR-corrected` = ddf$p_level_fdr
    #}
    # if (strcmp(toprocess,"region")){
    #   ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term,
    #                          p_level_fdr= as.factor(case_when(
    #                            # p_fdr > .1 ~ '0',
    #                            # p_fdr < .1 & p_fdr > .05 ~ '1',
    #                            p_fdr > .05 ~ '1',
    #                            p_fdr < .05 & p_fdr > .01 ~ '2',
    #                            p_fdr < .01 & p_fdr > .001 ~ '3',
    #                            p_fdr <.001 ~ '4')))
    #   ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
    #   ddf$`p, FDR-corrected` = ddf$p_level_fdr
    # }
    #if ("network" %in% splits) {ddf$visuomotor_grad <- factor(ddf$visuomotor_grad, labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal"))}
    terms <- unique(ddf$term[ddf$effect=="fixed"])
    
    if (strcmp(toprocess,"region")){ # get network for each region
      ddf <- ddf %>% mutate(network=case_when(
        region=="pfc11L" ~ "C",
        region=="pfc11aR" ~ "L",
        region=="pfc11b47R" ~ "C",
        region=="pfc11m13L" ~ "L",
        region=="pfc14m3225L" ~ "D",
        region=="pfc14m3225R" ~ "D",
        region=="pfc14r14c11mL" ~ "L",
        region=="pfc14r14c11mR" ~ "L",
        region=="pfc14rc2L" ~ "C",
        region=="pfc2432R" ~ "D",
        region=="pfc47L" ~ "C",
        region=="pfcd10R" ~ "D",
        region=="pfcfp10L" ~ "D",
        region=="pfcfp10R" ~ "L",
        region=="pfcp3224L" ~ "D",
        region=="pfcrdb10L" ~ "D",
        region=="pfcrl10R" ~ "C"
      )
      )
    }
    
    if (strcmp(toprocess,"symmetry-by-HC") | strcmp(toprocess,'symmetry-by-bin') | strcmp(toprocess,"symmetry-by-HC-by-outcome") | strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc') | strcmp(toprocess,'symmetry-group-by-HC-by-outcome')){
      ddf <- ddf %>% mutate(symmetry_group1=case_when(
        symmetry_group==1 ~ '11/47',
        symmetry_group==2 ~ 'fp10',
        symmetry_group==3 ~ '14rc11m',
        symmetry_group==4 ~ '14m25/32',
        symmetry_group==5 ~ '24/32',
        symmetry_group==6 ~ 'rl9_10',
        symmetry_group==7 ~ 'd10',
        symmetry_group==8 ~ '11/13')) %>% 
        mutate(network_group=case_when(
          symmetry_group==1 ~ 'CTR',
          symmetry_group==2 ~ 'DMN',
          symmetry_group==3 ~ 'LIM',
          symmetry_group==4 ~ 'DMN',
          symmetry_group==5 ~ 'DMN',
          symmetry_group==6 ~ 'CTR',
          symmetry_group==7 ~ 'DMN',
          symmetry_group==8 ~ 'LIM'
        ))
    }
    fills = c('red','red','blue','blue','blue','cyan','green','magenta')
    fills1 = c('green','magenta','red','blue','cyan')
    
    if (strcmp(toprocess,'network-by-HC') | strcmp(toprocess,'network-by-HC-by-side')){
      ddf <- ddf %>% mutate(network1 = 
                              case_when(network=='C'~'2C',
                                        network=='D'~'1D',
                                        network=='L'~'3L')) %>% 
        mutate(network2 = case_when(network=='C'~'CTR',
                                    network=='D'~'DMN',
                                    network=='L'~'LIM'))
    }
    
    if (!strcmp(totest,'anatomy')){
      for (fe in terms) {
        # fe <- terms[1] # test only
        edf <- ddf %>% filter(term == paste(fe) & t < 8) 
        termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
        if (strcmp(toprocess,'symmetry-by-HC') | (strcmp(toprocess,'symmetry-by-bin')) | strcmp(toprocess,'symmetry-group-by-HC-by-outcome') | strcmp(toprocess,'symmetry-by-HC-wPE') | strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
          edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9_10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
          edf$symmetry_group2 <- factor(edf$symmetry_group1, levels = c("14rc11m","11/47",'rl9_10',"fp10","24/32","14m25/32","d10","11/13"))
          edf$network_group = factor(edf$network_group,levels = c('CTR','DMN','LIM'))
        }
        # plot stream gradients
        
        if (!strcmp(toprocess,"symmetry-group-by-HC-by-outcome") & !strcmp(toprocess,"symmetry-group-by-HC-by-rewFunc")){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
        }
        
        if (strcmp(toprocess,"symmetry_group")){
          print(ggplot(edf, aes(t, symmetry_group)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_wrap(~side) +
                  # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } else if (strcmp(toprocess,"network")){
          print(ggplot(edf, aes(t, network2)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } else if (strcmp(toprocess,"region")){
          print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_wrap(~network) +
                  # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } else if (strcmp(toprocess,"side1")){
          print(ggplot(edf, aes(t, side)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_wrap(~network) +
                  # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } else if (strcmp(toprocess,"region")){
          print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_wrap(~network) +
                  # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } else if (strcmp(toprocess,"network-by-HC") | strcmp(toprocess,'network-by-HClag')){
          print(ggplot(edf,aes(t,network2))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_wrap(~HC_region) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } else if (strcmp(toprocess,"network-by-HC-wEC") | strcmp(toprocess,"network-by-HC-wO")) {
          print(ggplot(edf,aes(t,network))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_wrap(~HC_region) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } else if (strcmp(toprocess,"region-by-HC")){
          print(ggplot(edf,aes(t,HC_region))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_grid(~region) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } else if (strcmp(toprocess,"IS-by-HC")){
          print(ggplot(edf,aes(t,IS))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_grid(~HC_region) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } else if (strcmp(toprocess,"symmetry-by-HC")){
          gg1 <- ggplot(edf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
            facet_grid(network_group~HC_region,scales="free_y",space="free_y") +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          
          gg2 <- ggplot_gtable(ggplot_build(gg1))
          stripr <- which(grepl('strip-t', gg2$layout$name))
          k <- 1
          for (i in stripr) {
            j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
            gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
            k <- k+1
          }
          grid.draw(gg2)
          
        } else if (strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
          qdf <- edf %>% filter(HC_region=='PH')
          fname = paste(totest,"_",toalign, "_", toprocess, "_", termstr, '_PH','-',hc_LorR, ".pdf", sep = "")
          
          pdf(fname, width = 9, height = 3.5)
          gg1 <- ggplot(qdf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
            facet_grid(network_group~rewFunc,scales="free_y",space="free_y") +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          
          gg2 <- ggplot_gtable(ggplot_build(gg1))
          stripr <- which(grepl('strip-t', gg2$layout$name))
          k <- 1
          for (i in stripr) {
            j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
            gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
            k <- k+1
          }
          grid.draw(gg2)
          dev.off()
          
          qdf <- edf %>% filter(HC_region=='AH')
          fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr, '_AH','-',hc_LorR, ".pdf", sep = "")
          
          pdf(fname, width = 9, height = 3.5)
          gg1 <- ggplot(qdf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
            facet_grid(network_group~rewFunc,scales="free_y",space="free_y") +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          
          gg2 <- ggplot_gtable(ggplot_build(gg1))
          stripr <- which(grepl('strip-t', gg2$layout$name))
          k <- 1
          for (i in stripr) {
            j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
            gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
            k <- k+1
          }
          grid.draw(gg2)
          dev.off()
          
          
          
          
        } else if (strcmp(toprocess,'symmetry-group-by-HC-by-outcome')){
          qdf <- edf %>% filter(outcome1=='Omission')
          fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr, '_Omission','-',hc_LorR, ".pdf", sep = "")
          
          pdf(fname, width = 9, height = 3.5)
          gg1 <- ggplot(qdf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
            facet_grid(network_group~HC_region,scales="free_y",space="free_y") +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          
          gg2 <- ggplot_gtable(ggplot_build(gg1))
          stripr <- which(grepl('strip-t', gg2$layout$name))
          k <- 1
          for (i in stripr) {
            j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
            gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
            k <- k+1
          }
          grid.draw(gg2)
          dev.off()
          
          qdf <- edf %>% filter(outcome1=='Reward')
          fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr, '_Reward','-',hc_LorR, ".pdf", sep = "")
          
          pdf(fname, width = 9, height = 3.5)
          gg1 <- ggplot(qdf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
            facet_grid(network_group~HC_region,scales="free_y",space="free_y") +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          
          gg2 <- ggplot_gtable(ggplot_build(gg1))
          stripr <- which(grepl('strip-t', gg2$layout$name))
          k <- 1
          for (i in stripr) {
            j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
            gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
            k <- k+1
          }
          grid.draw(gg2)
          dev.off()
        } else if (strcmp(toprocess,'network-by-HC-by-outcome')){
          gg1 <- ggplot(edf,aes(t,network))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
            facet_grid(outcome1~HC_region) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg1)
        } else if (strcmp(toprocess,"network-by-HC-by-side") | strcmp(toprocess,'network-by-HClag')){
          print(ggplot(edf,aes(t,network2))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_wrap(HC_region~side_hc) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        } 
        
        if (!strcmp(toprocess,"symmetry-group-by-HC-by-outcome") & !strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
          dev.off()
        }
        
        
        if (strcmp(toprocess,"network")){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
            geom_line(size = 1) + geom_point() +
            geom_errorbar() +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
        } else if (strcmp(toprocess,"side2")){
          fname = paste(behavmodel,'-',toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=side, size=`p, FDR-corrected`)) +
            geom_line(size = 1) + geom_point() +
            geom_errorbar() +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_color_brewer(palette="Set1") + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
        } else if (strcmp(toprocess,"network-by-HC") | strcmp(toprocess,'network-by-HClag')){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          # gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
          #   geom_line(size = 1) + geom_point() +
          #   geom_errorbar() +
          #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~HC_region) +
          #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          #   scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
          #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          if (all(edf$`p, FDR-corrected`=='p < .001')){
            gg<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network1)) + 
              geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + 
              # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
              geom_line(size = 1) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
              geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
              scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
              #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
              theme_bw(base_size=13) +
              facet_wrap(~HC_region) +
              theme(legend.title = element_blank(),
                    panel.grid.major = element_line(colour = "grey45"), 
                    panel.grid.minor = element_line(colour = "grey45"), 
                    panel.background = element_rect(fill = 'grey40'),
                    axis.title.y = element_text(margin=margin(r=6)),
                    axis.title.x = element_text(margin=margin(t=6)))
          } else {
            gg<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network1)) + 
              geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
              # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
              geom_line(size = 1) + theme(legend.position = "none") +
              geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
              scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
              #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
              theme_bw(base_size=13) +
              facet_wrap(~HC_region) +
              theme(legend.title = element_blank(),
                    panel.grid.major = element_line(colour = "grey45"), 
                    panel.grid.minor = element_line(colour = "grey45"), 
                    panel.background = element_rect(fill = 'grey40'),
                    axis.title.y = element_text(margin=margin(r=6)),
                    axis.title.x = element_text(margin=margin(t=6)))
          }
          print(gg)
          dev.off()
        } else if (strcmp(toprocess,"network-by-HC-wEC") | strcmp(toprocess,"network-by-HC-wO")){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network, size=`p, FDR-corrected`)) +
            geom_line(size = 1) + geom_point() +
            geom_errorbar() +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~HC_region) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_color_brewer(palette="Set1") + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
        } else if (strcmp(toprocess,"region-by-HC")){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=HC_region, size=`p, FDR-corrected`)) +
            geom_line(size = 1) + geom_point() +
            geom_errorbar() +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_grid(~region) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_color_brewer(palette="Set1") + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
        } else if (strcmp(toprocess,"IS-by-HC")){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=as.factor(IS), size=`p, FDR-corrected`)) +
            geom_line(size = 1) + geom_point() +
            geom_errorbar() +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_grid(~HC_region) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_color_brewer(palette="Set1") + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
        } else if (strcmp(toprocess,"symmetry-by-HC")){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=as.factor(HC_region), size=`p, FDR-corrected`)) +
            geom_line(size = 1) + geom_point() +
            geom_errorbar() +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_grid(~symmetry_group1) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_color_brewer(palette="Set1") + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          
          gg2 <- ggplot_gtable(ggplot_build(gg1))
          stripr <- which(grepl('strip-t', gg2$layout$name))
          k <- 1
          for (i in stripr) {
            j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
            gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
            k <- k+1
          }
          grid.draw(gg2)
          #print(gg2)
          dev.off()
        } else if (strcmp(toprocess,"network-by-HC-by-side")){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
            geom_line(size = 1) + geom_point() +
            geom_errorbar() +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(HC_region~side_hc) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
        } 
      }
    } else if (strcmp(totest,'anatomy')){
      
      
      for (fe in terms) {
        edf <- ddf %>% filter(term == paste(fe) & t < 8) 
        termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
        if (strcmp(toprocess,'symmetry-by-HC') | (strcmp(toprocess,'symmetry-by-bin')) | strcmp(toprocess,'symmetry-group-by-HC-by-outcome') | strcmp(toprocess,'symmetry-by-HC-wPE') | strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
          edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9_10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
          edf$symmetry_group2 <- factor(edf$symmetry_group1, levels = c("14rc11m","11/47",'rl9_10',"fp10","24/32","14m25/32","d10","11/13"))
          edf$network_group = factor(edf$network_group,levels = c('CTR','DMN','LIM'))
        }
        
        if (strcmp(toprocess,'region')){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          print(ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
          dev.off()
          
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=HC_region, size=`p, FDR-corrected`)) +
            geom_line(size = 1) + geom_point() +
            geom_errorbar() +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_color_brewer(palette="Set1") + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
        } else if (strcmp(toprocess,'symmetry-by-bin')){
          fname = paste(behavmodel,'-',totest,"_",toalign, toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg1 <- ggplot(edf,aes(t,bin_num))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
            facet_grid(~symmetry_group1,scales="free_y",space="free_y") +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          
          gg2 <- ggplot_gtable(ggplot_build(gg1))
          stripr <- which(grepl('strip-t', gg2$layout$name))
          k <- 1
          for (i in stripr) {
            j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
            gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
            k <- k+1
          }
          grid.draw(gg2)
          dev.off()
          
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=as.factor(bin_num), size=`p, FDR-corrected`)) +
            geom_line(size = 1) + geom_point() +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_color_brewer(palette="Set1") + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
        }
      }
    }
  }
  
  if (dovif==TRUE){
    cl <- detectCores()
    registerDoParallel(6)
    nevt <- unique(Q$evt_time)
    R <- foreach(i=1:length(nevt), .combine=rbind) %dopar% {
      
      Qahl <- Q %>% filter(evt_time==nevt[i] & HC_region=='AH')
      mahl <- lmer(vmPFC_decon~trial_neg_inv_sc + rt_csv_sc + 
                     v_entropy_wi * HCwithin + v_entropy_wi_change * HCwithin + 
                     log10kld3 + v_max_wi * HCwithin + 
                     abs_pe_max_lag_sc * HCwithin + abs_pe_max_sc * HCwithin + 
                     HCbetween + (1 | id/run) + (1 | side_vmpfc) + (1 | region),Qahl)
      Qphl <- Q %>% filter(evt_time==nevt[i] & HC_region=='PH')
      mphl <- lmer(vmPFC_decon~trial_neg_inv_sc * HCwithin + rt_csv_sc * HCwithin + 
                     v_entropy_wi * HCwithin + v_entropy_wi_change * HCwithin + 
                     kld3 * HCwithin + v_max_wi * HCwithin + 
                     abs_pe_max_lag_sc * HCwithin + abs_pe_max_sc * HCwithin + 
                     HCbetween + (1 | id/run) + (1 | side_vmpfc) + (1 | region),Qphl)
      
      Qahr <- Qr %>% filter(evt_time==nevt[i] & HC_region=='AH')
      mahr <- lmer(vmPFC_decon~trial_neg_inv_sc * HCwithin + rt_csv_sc * HCwithin + 
                     v_entropy_wi * HCwithin + v_entropy_wi_change * HCwithin + 
                     kld3 * HCwithin + v_max_wi * HCwithin + 
                     abs_pe_max_lag_sc * HCwithin + abs_pe_max_sc * HCwithin + 
                     HCbetween + (1 | id/run) + (1 | side_vmpfc) + (1 | region),Qahr)
      Qphr <- Qr %>% filter(evt_time==nevt[i] & HC_region=='PH')
      mphr <- lmer(vmPFC_decon~trial_neg_inv_sc * HCwithin + rt_csv_sc * HCwithin + 
                     v_entropy_wi * HCwithin + v_entropy_wi_change * HCwithin + 
                     kld3 * HCwithin + v_max_wi * HCwithin + 
                     abs_pe_max_lag_sc * HCwithin + abs_pe_max_sc * HCwithin + 
                     HCbetween + (1 | id/run) + (1 | side_vmpfc) + (1 | region),Qphr)
      
      # do vif
      vmah <- car::vif(mah)
      vmph <- car::vif(mph)
      
      vif0 <- data.frame(rbind(vmahl,vmphl,vmahr,vmphr))
      vif0 <- vif0 %>% pivot_longer(cols=everything())
      vif0 <- vif0 %>% mutate(side_hc=c(rep('l',18),rep('l',18),rep('r',18),rep('r',18)))
      vif0 <- vif0 %>% mutate(HC_region=c(rep("AH",18),rep("PH",18),rep("AH",18),rep("PH",18)))
      vif0 <- vif0 %>% mutate(evt_time=rep(nevt[i],18*4))
      R <- vif0
      
    }
    
    #terms <- unique(R$name)
    #for (fe in terms) {
    # fe <- terms[1] # test only
    edf <- R %>% filter(name == paste(fe) & evt_time < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    
    if (strcmp(toprocess,"network-by-HC-by-side")){
      fname = paste(totest,"_",toalign, "_line_", toprocess, "_", termstr, "-vif.pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg <- ggplot(R, aes(x=evt_time, y=value)) +
        geom_line(size = 1) + geom_point() +
        geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(HC_region~side_hc) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
        labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      print(gg)
      dev.off()
    } 
    if (strcmp(toprocess,"network-by-HC-by-side")){
      fname = paste(totest,"_",toalign, "_line_", toprocess, "_", termstr, "-vif.pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg <- ggplot(edf, aes(x=evt_time, y=value)) +
        geom_line(size = 1) + geom_point() +
        geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(HC_region~side_hc) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
        labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      print(gg)
      dev.off()
    } 
  }
  return(Q)    
}
