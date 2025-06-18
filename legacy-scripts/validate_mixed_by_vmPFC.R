

validate_mixed_by_vmPFC <- function(ncores=6,toalign='feedback',toprocess='network',totest='base',
                                    behavmodel='compressed',repo_directory="~/clock_analysis",
                                    cache_dir='~/vmPFC/MEDUSA Schaefer Analysis',
                                    reload=TRUE,replot=TRUE,load_existing=FALSE,ddf_file='',d=NULL){
  
  
  library(tidyverse)
  library(lme4)
  library(pracma)
  library(grid)
  
  #toalign <- 'feedback'
  #toprocess <- 'symmetry_group'#'network'
  #toprocess <- 'network'
  #toprocess <- 'symmetry_group_by_rewFunc'
  #totest <-  'base-wPE' #'rt_change' #'rewfun' #'expl_code' #nextrt, interactions, base, pe, expl_code
  #totest <- 'EV' 
  #behavmodel <- 'compressed'
  
  #repo_directory <- "~/clock_analysis"
  #cache_dir = "/Users/andypapale/Box/SCEPTIC_fMRI/vmPFC/cache"
  #cache_dir = '/Users/andypapale/vmPFC/MEDUSA Schaefer Analysis'
  #gc()
  
  # load MEDUSA deconvolved data

if (load_existing==FALSE){
  if (reload==TRUE){
    if (strcmp(toalign,'clock')){
      load(file.path(cache_dir, 'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
      vmPFC <- clock_comb
      vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
      rm(clock_comb)
    } else if (strcmp(toalign,'feedback')){
      load(file.path(cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
      vmPFC <- fb_comb
      vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
      rm(fb_comb)
    } else if (strcmp(toalign,'rt_vmax')){
      load(file.path(cache_dir,  'rt_vmax_vmPFC_Schaefer_tall_ts_1.Rdata'))
      vmPFC <- rt_vmax_comb
      vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
      rm(rt_vmax_comb)
    }
    
    vmPFC <- vmPFC %>% select(id,run,trial,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
    if (strcmp(totest,'base') | strcmp(totest,'base-wEV') | strcmp(totest,'base-wPE')){ # vmPFC is outcome variable, do not scale
      vmPFC <- vmPFC %>% mutate(vmPFCdecon = decon_mean) %>% select(!decon_mean)
    } else if (strcmp(totest,'nextrt') | strcmp(totest,'rt_change') | strcmp(totest,'rtvmax')){ # vmPFC is predictor variable, scale
      vmPFC <- vmPFC %>% group_by(id,run) %>% mutate(vmPFC_decon_scaled = scale(decon_mean), vmPFCbetween1=mean(decon_mean,na.rm=TRUE)) %>% ungroup() %>% arrange(id,run,trial,evt_time) 
      if (strcmp(toprocess,'network')){
        vmPFC <- vmPFC %>% group_by(id,run,trial,evt_time,network) %>% summarize(vmPFC_scale = mean(vmPFC_decon_scaled,na.rm=TRUE), vmPFCbetween=mean(vmPFCbetween1,na.rm=TRUE)) %>% arrange(id,run,trial,evt_time)
      } else if (strcmp(toprocess,'symmetry_group')){
        vmPFC <- vmPFC %>% group_by(id,run,trial,evt_time,symmetry_group) %>% summarize(vmPFC_scale = mean(vmPFC_decon_scaled,na.rm=TRUE), vmPFCbetween=mean(vmPFCbetween1,na.rm=TRUE)) %>% arrange(id,run,trial,evt_time)
      }
    }
    
    
    # load behavioral data, get lags, merge with signal
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
    
    df <- df %>% mutate(richorpoor = case_when(total_earnings < median(total_earnings) ~ 'P',
                                               total_earnings >= median(total_earnings) ~ 'R'))
    
    
    d <- merge(vmPFC, df, by = c("id", "run","run_trial")) %>% arrange('id','run','run_trial','evt_time')
    #d$expl_code <- relevel(as.factor(d$expl_code),ref="Choose-Stay")
    
    # d$AP <-d$'A/P'
    # d$ML <-d$'M/L'
    # d$IS <-d$'I/S'
    
    # add side
    d <- d %>% mutate(side = substr(region, nchar(region), nchar(region)))
    d$rewFunc <- relevel(as.factor(d$rewFunc),ref="CEV")
    d$last_outcome <- relevel(as.factor(d$last_outcome),ref="Omission")
    d$outcome <- relevel(as.factor(d$outcome),ref="Omission")
    d$side <- as.factor(d$side)
    d$symmetry_group <- as.factor(d$symmetry_group)
    d$network <- as.factor(d$network)
    # censor prior ITI on clock
    # if (strcmp(toalign,'clock')){
    #   d$vmPFCdecon[d$evt_time < -(d$iti_lag+d$rt_lag/1000)] <- NA
    # }
  }
  
  # mixed_by call
  if (strcmp(toprocess,'network')){
    splits = c("evt_time","network")
  } else if (strcmp(toprocess,"symmetry_group")){
    splits = c("evt_time","symmetry_group")
  } else if (strcmp(toprocess,"region")){
    splits = c("evt_time","region","network")
  } else if (strcmp(toprocess,"side1")){
    splits = c("evt_time","side","network")
  } else if (strcmp(toprocess,"side2")){
    splits = c("evt_time","side")
  } else if (strcmp(toprocess,"AP")){
    splits = c("evt_time",'AP')
  } else if (strcmp(toprocess,"IS")){
    splits = c("evt_time",'IS')
  } else if (strcmp(toprocess,"ML")){
    splits = c("evt_time","ML")
  } else if (strcmp(toprocess,"network-by-outcome")){
    splits = c("evt_time","network","outcome")
  } else if (strcmp(toprocess,'symmetry_group_by_rewFunc')){
    splits = c("evt_time","symmetry_group","rewFunc")
  } else if (strcmp(toprocess,'symmetry_group_by_outcome')){
    splits = c("evt_time","symmetry_group","outcome")
  }
  
  # take out side as fixed effect
  if (strcmp(totest,"base")){
    if (strcmp(toprocess,"network")){
      if (strcmp(toalign,'clock')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + 
                                   log10kld3_lag  + v_max_wi  + 
                                   last_outcome + (1|id/run) + (1|region))  
      } else if (strcmp(toalign,'feedback')) {
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change +
                                   log10kld3  + v_max_wi  + 
                                   outcome + (1|id/run) + (1|region)) 
      } else if (strcmp(toalign,'rt_vmax')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc +
                                   v_entropy_wi +
                                   v_entropy_wi_change +
                                   log10kld3 + v_max_wi +
                                   (1|id/run) + (1|region) + (1|side))
      }
    } else if (strcmp(toprocess,"symmetry_group")){
      if (strcmp(toalign,'feedback')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + 
                                   log10kld3  + v_max_wi  + outcome +
                                   (1|id/run)) 
      } else if (strcmp(toalign,'clock')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + 
                                   log10kld3_lag  + v_max_wi  +
                                   last_outcome +
                                   (1|id/run)) 
      } else if (strcmp(toalign,'rt_vmax')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc +
                                   v_entropy_wi +
                                   log10kld3 + v_max_wi +
                                   (1|id/run))
      }
    } else if (strcmp(toprocess,'symmetry_group_by_outcome')){
      if (strcmp(toalign,'feedback')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   rt_vmax_sc  + 
                                   rt_vmax_change_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + abs_pe_max_sc +
                                   kld3  + v_max_wi  + 
                                   (1|id/run)) 
      } else if (strcmp(toalign,'clock')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   rt_vmax_lag_sc + 
                                   rt_vmax_change_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + abs_pe_max_lag_sc +
                                   kld3_lag  + v_max_wi  +
                                   (1|id/run)) 
      } else if (strcmp(toalign,'rt_vmax')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc +
                                   rt_vmax_lag_sc + 
                                   v_entropy_wi +
                                   kld3 + v_max_wi +
                                   (1|id/run))
      }
    }
  }
  if (strcmp(toprocess,'symmetry_group_by_rewFunc')){
    if (strcmp(toalign,'feedback')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_vmax_sc + 
                                 rt_vmax_change_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 kld3  + v_max_wi  + outcome +
                                 (1|id/run)) 
    } else if (strcmp(toalign,'clock')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_vmax_lag_sc + 
                                 rt_vmax_change_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 kld3_lag  + v_max_wi +
                                 last_outcome +
                                 (1|id/run)) 
    } else if (strcmp(toalign,'rt_vmax')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc +
                                 rt_vmax_lag_sc + 
                                 v_entropy_wi +
                                 kld3 + v_max_wi +
                                 (1|id/run))
    }
  }
  if (strcmp(totest,"base-wPE")){
    if (strcmp(toprocess,"network")){
      if (strcmp(toalign,'clock')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + 
                                   log10kld3_lag  + v_max_wi  + 
                                   abs_pe_max_lag_sc + (1|id/run) + (1|region))  
      } else if (strcmp(toalign,'feedback')) {
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change +
                                   log10kld3  + v_max_wi  + 
                                   abs_pe_max_sc + (1|id/run) + (1|region)) 
      } else if (strcmp(toalign,'rt_vmax')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc +
                                   v_entropy_wi +
                                   v_entropy_wi_change +
                                   log10kld3 + v_max_wi +
                                   (1|id/run) + (1|region) + (1|side))
      }
    } else if (strcmp(toprocess,"symmetry_group")){
      if (strcmp(toalign,'feedback')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + 
                                   log10kld3  + v_max_wi  + abs_pe_max_sc +
                                   (1|id/run)) 
      } else if (strcmp(toalign,'clock')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + 
                                   log10kld3_lag  + v_max_wi  +
                                   abs_pe_max_lag_sc +
                                   (1|id/run)) 
      } else if (strcmp(toalign,'rt_vmax')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc +
                                   rt_vmax_lag_sc + 
                                   v_entropy_wi +
                                   log10kld3 + v_max_wi +
                                   (1|id/run))
      }
    } else if (strcmp(toprocess,'symmetry_group_by_outcome')){
      if (strcmp(toalign,'feedback')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   rt_vmax_sc  + 
                                   rt_vmax_change_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + abs_pe_max_sc +
                                   kld3  + v_max_wi  + 
                                   (1|id/run)) 
      } else if (strcmp(toalign,'clock')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                   rt_vmax_lag_sc + 
                                   rt_vmax_change_sc + 
                                   v_entropy_wi + 
                                   v_entropy_wi_change + abs_pe_max_lag_sc +
                                   kld3_lag  + v_max_wi  +
                                   (1|id/run)) 
      } else if (strcmp(toalign,'rt_vmax')){
        decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc +
                                   rt_vmax_lag_sc + 
                                   v_entropy_wi +
                                   kld3 + v_max_wi +
                                   (1|id/run))
      }
    }
  }
  if (strcmp(toprocess,'symmetry_group_by_rewFunc')){
    if (strcmp(toalign,'feedback')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_vmax_sc  + 
                                 rt_vmax_change_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 kld3  + v_max_wi  + outcome +
                                 (1|id/run)) 
    } else if (strcmp(toalign,'clock')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_vmax_lag_sc + 
                                 rt_vmax_change_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 kld3_lag  + v_max_wi +
                                 last_outcome +
                                 (1|id/run)) 
    } else if (strcmp(toalign,'rt_vmax')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc +
                                 rt_vmax_lag_sc + 
                                 v_entropy_wi +
                                 kld3 + v_max_wi +
                                 (1|id/run))
    }
  }
  if (strcmp(toprocess,"network-by-outcome")){
    if (strcmp(toalign,'feedback')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_vmax_sc  + 
                                 rt_vmax_change_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 kld3  + v_max_wi  + 
                                 abs_pe_max_sc +
                                 (1|id/run) + (1|region) + (1|side)) 
    } else if (strcmp(toalign,'clock')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_lag_sc + rt_vmax_lag_sc + 
                                 rt_vmax_change_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 kld3_lag  + v_max_wi+
                                 abs_pe_max_sc +
                                 (1|id/run) + (1|region) + (1|side)) 
    }
  }
  
  if (strcmp(totest,"nextrt") | strcmp(totest,"rt_change") | strcmp(totest,'rtvmax')){
    if (strcmp(toalign,"feedback")){
      decode_formula = formula(~ vmPFC_scale*rt_csv_sc +
                                 vmPFC_scale*rt_vmax_lag_sc + 
                                 (1|id/run)) # for RT- or feedback-aligned
    } else{
      decode_formula = formula(~ vmPFC_scale*rt_lag_sc + 
                                 vmPFC_scale*rt_vmax_lag_sc +
                                 (1|id/run)) # for clock-aligned
    }
  }
  
  if (strcmp(totest,'EV')){
    if (strcmp(toprocess,"network")){
      if (strcmp(toalign,'feedback')){
        decode_formula = formula(~ rt_csv_sc + rt_lag_sc +
                                   ev+ trial_neg_inv_sc +
                                   outcome + v_entropy_wi + 
                                   + (1|id/run) + (1|region) + (1|side))
      } else if (strcmp(toalign,'clock')){
        decode_formula = formula(~ rt_lag_sc+ rt_csv_sc +
                                   ev + trial_neg_inv_sc +
                                   last_outcome + v_entropy_wi +
                                   + (1|id/run) + (1|region) + (1|side))
      }
    } else if (strcmp(toprocess,"symmetry-group")){
      if (strcmp(toalign,'feedback')){
        decode_formula = formula(~ rt_csv_sc + rt_lag_sc +
                                   ev+ trial_neg_inv_sc +
                                   outcome + v_entropy_wi +
                                   + (1|id/run))
      } else if (strcmp(toalign,'clock')){
        decode_formula = formula(~ rt_lag_sc+ rt_csv_sc +
                                   ev + trial_neg_inv_sc +
                                   last_outcome + v_entropy_wi +
                                   + (1|id/run))
      }
    }
  }
  
  source("~/fmri.pipeline/R/mixed_by.R")
  if (strcmp(totest,"nextrt")){
    if (strcmp(toalign,"clock")){
      ddf <- mixed_by(d, outcomes = "rt_sec", rhs_model_formulae = decode_formula , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
    } else {
      ddf <- mixed_by(d, outcomes = "rt_next", rhs_model_formulae = decode_formula , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
    }
  }
  if (strcmp(totest,"rtvmax")){
    if (strcmp(toalign,"clock")){
      ddf <- mixed_by(d, outcomes = "rt_vmax", rhs_model_formulae = decode_formula , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
    } else {
      ddf <- mixed_by(d, outcomes = "rt_vmax_lead", rhs_model_formulae = decode_formula , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
    }
  }
  if (strcmp(totest,"base") | (strcmp(totest,"base-wPE")) | strcmp(totest,'EV')){
    ddf <- mixed_by(d,outcomes = "vmPFCdecon", rhs_model_formulae = decode_formula , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3)
  }
  ## Check plots
  if (strcmp(toalign,"clock")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock')
  } else if (strcmp(toalign,"feedback")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback')
  } else if (strcmp(toalign,'rt_vmax')){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_rt_vmax')
  }
  
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-',behavmodel,'-',toalign,'-',totest,'-',toprocess,'.Rdata')) 
}
  if (replot==TRUE){
    
    if (load_existing==TRUE){
      load(ddf_file)
    }
    
    message("\nPlotting streams decoding")
    library(viridis)
    library(wesanderson)
    pal = wes_palette("FantasticFox1", 3, type = "discrete")
    epoch_label = paste("Time relative to",toalign, "[s]")
    
    if (strcmp(toprocess,"network-by-outcome") | strcmp(toprocess,'symmetry_group_by_outcome')){
      colnames(ddf$coef_df_reml)[4] <- "outcome1"
    }
    qdf <- ddf
    ddf <- ddf$coef_df_reml
    ddf$t <- ddf$evt_time
    #if (strcmp(toprocess,"network")){
    ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
    ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
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
    
    if (strcmp(toprocess,"symmetry_group") | strcmp(toprocess,"symmetry_group_by_rewFunc") | strcmp(toprocess,"symmetry_group_by_outcome")){
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
    if (strcmp(toprocess,'network')){
      ddf <- ddf %>% mutate(network1 = 
                              case_when(network=='C'~'2C',
                                        network=='D'~'1D',
                                        network=='L'~'3L')) %>%
        mutate(network2 = case_when(network=='C'~'CTR',
                                    network=='D'~'DMN',
                                    network=='L'~'LIM'))
    }
    
    fills = c('red','red','blue','blue','blue','blue','green','green')
    fills1 = c('red','blue','green')
    ddf$t <- as.numeric(ddf$t)
    for (fe in terms) {
      # fe <- terms[1] # test only
      edf <- ddf %>% filter(term == paste(fe) & t < 8)
      if (strcmp(toprocess,'symmetry_group') | strcmp(toprocess,"symmetry_group_by_rewFunc") | strcmp(toprocess,'symmetry_group_by_outcome')){
        edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9-10',"11-47","d10","24-32","14m25-32","fp10","11-13","14rc11m"))
        edf$symmetry_group2 <- factor(edf$symmetry_group1, levels = c("14rc11m","11-47",'rl9-10',"fp10","24-32","14m25-32","d10","11-13"))
        edf$network_group = factor(edf$network_group,levels = c('CTR','DMN','LIM'))
      }
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      # plot stream gradients
      fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr, ".pdf", sep = "")
      
      pdf(fname, width = 9, height = 3.5)
      
      if (strcmp(toprocess,"network")){
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
        print(ggplot(edf,aes(t,network))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
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
      }  else if (strcmp(toprocess,"region-by-HC")){
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
      } else if (strcmp(toprocess,"symmetry_group")){
        gg1 <- (ggplot(edf,aes(t,symmetry_group2))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_grid(network_group~.,scales="free_y",space="free_y") +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        print(gg1)
      } else if (strcmp(toprocess,"network-by-outcome")){
        gg1 <- (ggplot(edf,aes(t,network))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_grid(~outcome) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        print(gg1)
      } else if (strcmp(toprocess,"symmetry_group_by_rewFunc")){
        gg1 <- (ggplot(edf,aes(t,symmetry_group2))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_grid(~rewFunc) +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        print(gg1)
      } else if (strcmp(toprocess,'symmetry_group_by_outcome')) {
        gg1 <- (ggplot(edf,aes(t,symmetry_group2))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                  facet_grid(network_group~outcome,scales="free_y",space="free_y") +
                  geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                  scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                  labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        print(gg1)
      }
      
      dev.off()
      
      if (strcmp(toprocess,"network")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
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
            #facet_wrap(~HC_region) +
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
            #facet_wrap(~HC_region) +
            theme(legend.title = element_blank(),
                  panel.grid.major = element_line(colour = "grey45"), 
                  panel.grid.minor = element_line(colour = "grey45"), 
                  panel.background = element_rect(fill = 'grey40'),
                  axis.title.y = element_text(margin=margin(r=6)),
                  axis.title.x = element_text(margin=margin(t=6)))
        }
        print(gg)
        dev.off()
      } else if (strcmp(toprocess,"side2")){
        fname = paste(behavmodel,'-',toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
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
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
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
      } else if (strcmp(toprocess,"network-by-HC-wEC") | strcmp(toprocess,"network-by-HC-wO")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
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
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
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
      } else if (strcmp(toprocess,"symmetry_group")){
        #browser()
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, size=`p, FDR-corrected`)) +
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
      } else if (strcmp(toprocess,"network-by-outcome")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network, size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_errorbar() + facet_wrap(~outcome) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_color_brewer(palette="Set1") + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off() 
      } else if (strcmp(toprocess,"symmetry_group_by_rewFunc")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=rewFunc, size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_errorbar() + facet_wrap(~symmetry_group1) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_color_brewer(palette="Set1") + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off() 
      } else if(strcmp(toprocess,'symmetry_group_by_outcome')){
      }
    }
  }
  return(d)
}