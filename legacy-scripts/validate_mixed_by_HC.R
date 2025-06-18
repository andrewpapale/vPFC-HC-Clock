# 2021-09-16 AndyP
# replicating HC analyses


validate_mixed_by_HC <- function(ncores=6,toalign='',toprocess='',totest='',behavmodel='',
                                 repo_directory="~/clock_analysis",HC_cache_dir='~/vmPFC/MEDUSA Schaefer Analysis',
                                 hc_LorR='L',
                                 reload=TRUE,replot=TRUE,load_existing=FALSE,ddf_file='',Q=NULL){
  
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(pracma)
  library(broom)
  library(grid)
  
  dolmer = 'false'
  #toalign <- 'feedback'
  #toprocess <- 'axis' # region
  #totest <- 'base' #'base-reversed' #'expl_code' #nextrt, interactions, base, pe, expl_code
  #behavmodel <- 'compressed' # uncompressed, # valuecomp
  lagHC <- FALSE
  lagvmPFC <- FALSE
  #hc_LorR <- 'L'
  
  #repo_directory <- "~/clock_analysis"
  #HC_cache_dir = "/Users/andypapale/vmPFC/deconvolved_evt_locked_smooth_in_mask_harvardoxford"
  #HC_cache_dir = '/Users/andypapale/vmPFC/MEDUSA Schaefer Analysis'
  #gc()
  force(load_existing)
  force(ncores)
  force(hc_LorR)
if (load_existing==FALSE){
  if (reload==TRUE){
    # load MEDUSA deconvolved data
    if (strcmp(toalign,'clock')){
      message("Loading HC medusa data from cache: ", HC_cache_dir)
      load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
      hc <- clock_comb %>% filter(evt_time > -6 & evt_time < 6)
      rm(clock_comb)
    } else if (strcmp(toalign,'feedback')){
      load(file.path(HC_cache_dir,'feedback_hipp_tall_ts_1.Rdata'))
      hc <- fb_comb
      hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
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
      hc <- hc %>% group_by(id,run,run_trial,evt_time,bin_num) %>% summarize(decon_mean1 = mean(decon_mean,na.rm=TRUE)) 
      hc <- hc %>% rename(decon_mean=decon_mean1)
    }
    hc <- hc %>% mutate(
      HC_region = case_when(
        bin_num <= 8 ~ 'AH',
        bin_num >8 ~ 'PH'
      ),
    )
    if (strcmp(toprocess,'region')){
      hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon_mean1 = mean(decon_mean,na.rm=TRUE))
      hc <- hc %>% rename(decon_mean=decon_mean1)
    }
    
    # load behavioral data, get lags, merge with signal
    if (strcmp(behavmodel,'compressed')){
      source('~/vmPFC/get_trial_data_vmPFC.R')
      df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
      df <- df %>% select(rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                          v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                          rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,pe_max_lag,v_entropy_wi_change) %>% 
        group_by(id, run) %>% arrange(id, run, run_trial) %>% 
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
    } else if (strcmp(behavmodel,'valuecomp')){
      error('not working yet')
      source('~/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
      df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
      df <- df %>% select(run,trial,run_trial,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                          v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                          rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc)
      df <- df %>% arrange(id,run,trial)
      dq <- read.csv('/Users/andypapale/vmPFC/value_compress_results/mmclock_fmri_exp_compress_psequate_variant3_mfx_trial_statistics.csv.gz')
      # need to get into format of compressed model
      nI = unique(dq$id)
      nR = unique(dq$run)
      for (iD in 1:length(nI)){
        iq = sum(dq$id==nI[iD])
        id = sum(df$id==nI[iD])
        if (iq!=id){
          for (iR in 1:length(nR)){
            iq1 = sum(dq$run==nR[iR])
            id1 = sum(df$run==nR[iR])
            if (iq1!=id1){
              dq <- dq %>% filter(!(dq$run==nR[iR] & dq$id==nI[iD]))
            }
          }
        }
      }
      dq$trial <- as.numeric(dq$trial)
      dq <- dq %>% arrange(id,run,trial)
      dq$rt_csv_sc <- df$rt_csv_sc
      dq <- dq %>% select(rt_csv,id, run, trial, pe_max, rt_vmax, score_csv,
                          rewFunc,v_chosen,v_entropy,rewFunc,iti_ideal,v_max,rt_csv_sc) %>%
        group_by(id) %>% 
        mutate(total_earnings=sum(score_csv)) %>% ungroup() %>%
        group_by(id, run) %>% 
        mutate(run_trial=case_when(
          trial >= 1 & trial <= 50 ~ trial,
          trial >= 51 & trial <= 100 ~ trial - 50, #dplyr/rlang has gotten awfully picky about data types!!
          trial >= 101 & trial <= 150 ~ trial - 100,
          trial >= 151 & trial <= 200 ~ trial - 150,
          trial >= 201 & trial <= 250 ~ trial - 200,
          trial >= 251 & trial <= 300 ~ trial - 250,
          trial >= 301 & trial <= 350 ~ trial - 300,
          trial >= 351 & trial <= 400 ~ trial - 350),
          outcome = case_when(
            score_csv > 0 ~ 'Reward',
            score_csv ==0 ~ 'Omission'
          ),
          trial_neg_inv_sc = scale(1/-run_trial),
        ) %>%
        mutate(v_max_wi = scale(v_max),
               v_entropy_wi = scale(v_entropy),
        ) %>%
        mutate(v_entropy_wi_lead = lead(v_entropy_wi),
               v_entropy_wi_change = v_entropy_wi-lag(v_entropy_wi),
               v_entropy_wi_change_lead = lead(v_entropy_wi_change),
               rt_vmax_lead = lead(rt_vmax),
               rt_vmax_lag_sc = lag(scale(rt_vmax)),
               rt_vmax_change = scale(rt_vmax) - lag(rt_vmax),
               pe_max_lag = lag(pe_max),
               expl_code = case_when(
                 rt_csv-lag(rt_csv) <= -800 ~ 'Choose-Sooner',
                 rt_csv-lag(rt_csv)   > -800 & rt_csv-lag(rt_csv) < 800 ~ 'Choose-Stay',
                 rt_csv-lag(rt_csv)   >= 800 ~ 'Choose-Later'),
               rt_lag_sc = lag(rt_csv_sc),
               rt_next = lead(rt_csv)/1000,
               rt_sec = rt_csv/1000,
               v_max_wi_lead = lead(v_max_wi),
               iti_lag = lag(iti_ideal),
               rt_csv_sc = scale(rt_csv),
               last_outcome = lag(outcome)
        ) %>% ungroup() %>% arrange(id,run,run_trial)
      dq$kld3_lag <- df$kld3_lag
      dq$kld3 <- df$kld3
      dq$rt_csv_sc <- df$rt_csv_sc
      dq$rt_lag_sc = df$rt_lag_sc
      df <- dq
    }
    
    Q <- merge(df, hc, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
    rm(hc)
  }
    # # censor prior ITI on clock
    # if (strcmp(toalign,'clock')){
    #   Q$decon_mean[Q$evt_time < -(Q$iti_lag+Q$rt_lag/1000)] <- NA
    # }
  if (strcmp(toprocess,'axis')){
    splits = c('evt_time','bin_num')
  } else if (strcmp(toprocess,'region')){
    splits = c('evt_time','HC_region')
  }
  
  if (strcmp(totest,"nextrt")){
    if (strcmp(toalign,"feedback")){
      decode_formula = formula(~ rt_csv_sc*HCwithin +  
                                 rt_vmax_lag_sc*HCwithin + 
                                 HCbetween +  (1|id/run)) # for RT- or feedback-aligned
    } else{
      decode_formula = formula(~ HCwithin*rt_lag_sc + 
                                 HCwithin*rt_vmax_lag_sc + 
                                 HCbetween + (1|id/run)) # for clock-aligned
    }
  }
  
  if (strcmp(totest,'base')){
    if (strcmp(toalign,'clock')){
      decode_formula = formula(~ trial_neg_inv_sc +
                                 rt_lag_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 v_max_wi  + log10kld3 +
                                 last_outcome +
                                 (1|id/run))  
    } else if (strcmp(toalign,'feedback')) {
      decode_formula = formula(~ trial_neg_inv_sc +
                                 rt_csv_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 v_max_wi + log10kld3 +
                                 outcome +
                                 (1|id/run))    
    }
  }
  if (strcmp(totest,'nokld3')){
    if (strcmp(toalign,'clock')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_lag_sc + rt_vmax_lag_sc + 
                                 scale(rt_vmax_change) + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 scale(kld3_lag) + v_max_wi  + 
                                 last_outcome + (1|id/run))   
    } else if (strcmp(toalign,'feedback')) {
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_lag_sc + rt_vmax_lag_sc  + 
                                 scale(rt_vmax_change) + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 v_max_wi  + 
                                 outcome +  # 
                                 (1|id/run))    
    } 
  }
  if (strcmp(totest,'original')){
    if (strcmp(toalign,'clock')){
      decode_formula = formula(~ scale(rt_csv) + side  + 
                                 scale(rt_vmax_lag_sc)  + 
                                 scale(rt_vmax_change)  + 
                                 v_entropy_wi  + 
                                 v_entropy_wi_change  +
                                 (1|id) + (1|side))
    } else if (strcmp(toalign,'feedback')) {
      decode_formula = formula(~ scale(rt_csv) + side  + 
                                 scale(rt_vmax_lag_sc)  + 
                                 scale(rt_vmax_change)  + 
                                 v_entropy_wi  + 
                                 v_entropy_wi_change  +
                                 (1|id) + (1|side))   
    } 
  }
  if (strcmp(totest,'norandrun')){
    if (strcmp(toalign,'clock')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_lag_sc + rt_vmax_lag_sc + 
                                 scale(rt_vmax_change): + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 scale(kld3_lag) + v_max_wi  + 
                                 last_outcome + (1|id))   
    } else if (strcmp(toalign,'feedback')) {
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_lag_sc + rt_vmax_lag_sc  + 
                                 scale(rt_vmax_change) + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 v_max_wi  + 
                                 outcome +  # 
                                 (1|id))    
    } 
  }
  if (strcmp(totest,'base&side')){
    if (strcmp(toalign,'clock')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + 
                                 rt_lag_sc + rt_vmax_lag_sc + 
                                 scale(rt_vmax_change) + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 v_max_wi  + 
                                 last_outcome + (1|id/run))  
    } else if (strcmp(toalign,'feedback')) {
      decode_formula = formula(~ rt_csv_sc + 
                                 rt_lag_sc + rt_vmax_lag_sc  + 
                                 scale(rt_vmax_change) + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 kld3 + v_max_wi  + 
                                 outcome + 
                                 v_entropy_wi:trial_neg_inv_sc +
                                 v_entropy_wi:kld3 +# 
                                 (1|id/run))    
    } 
  }
  if (strcmp(totest,'base-wPE')){
    if (strcmp(toalign,'clock')){
      decode_formula = formula(~ trial_neg_inv_sc +
                                 rt_lag_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 v_max_wi  + log10kld3 +
                                 abs_pe_max_lag_sc +
                                 (1|id/run))  
    } else if (strcmp(toalign,'feedback')) {
      decode_formula = formula(~ trial_neg_inv_sc +
                                 rt_csv_sc + 
                                 v_entropy_wi + 
                                 v_entropy_wi_change + 
                                 v_max_wi  + log10kld3 +
                                 abs_pe_max_sc +
                                 (1|id/run))    
    }
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
  } else if (strcmp(totest,"base") | strcmp(totest,"nokld3") | strcmp(totest,'original') | strcmp(totest,'norandrun') | strcmp(totest,'base&side') | strcmp(totest,'base-wPE')){
    ddf <- mixed_by(Q, outcomes = "decon_mean", rhs_model_formulae = decode_formula , split_on = splits,
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
  ## Check plots
  if (strcmp(toalign,"clock")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC')
  } else if (strcmp(toalign,"feedback")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC')
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
    
    epoch_label = paste("Time relative to",toalign, "[s]")
    ddq <- ddf
    ddf <- as_tibble(ddf$coef_df_reml)
    ddf$t <- ddf$evt_time
    #ddf$bin_num <- as.factor(ddf$bin_num)
    ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term, 
                           p_level_fdr = as.factor(case_when(
                             # p_fdr > .1 ~ '0',
                             # p_fdr < .1 & p_fdr > .05 ~ '1',
                             p_fdr > .05 ~ '1',
                             p_fdr < .05 & p_fdr > .01 ~ '2',
                             p_fdr < .01 & p_fdr > .001 ~ '3',
                             p_fdr <.001 ~ '4')))
    ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
    ddf$pfdr = ddf$p_level_fdr
    terms <- unique(ddf$term[ddf$effect=="fixed"])
    
    if (strcmp(toprocess,'axis')){
      ddf <- ddf %>% mutate(bin_num1 = case_when(
        bin_num==1 ~ 12,
        bin_num==2 ~ 11,
        bin_num==3 ~ 10,
        bin_num==4 ~ 9,
        bin_num==5 ~ 8,
        bin_num==6 ~ 7,
        bin_num==7 ~ 6,
        bin_num==8 ~ 5,
        bin_num==9 ~ 4,
        bin_num==10 ~ 3,
        bin_num==11 ~ 2,
        bin_num==12 ~ 1
      ))
    }
    library(wesanderson)
    for (fe in terms) {
      if (strcmp(toprocess,'axis')){
        pal = wes_palette("Zissou1", 12, type = "continuous")
      } else if (strcmp(toprocess,'region')){
        pal = wes_palette("Zissou1", 2, type = "continuous")
      }
      edf <- ddf %>% filter(term == paste(fe) & t < 8)
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      
      if (strcmp(toprocess,'region')){
        if (all(edf$`p, FDR-corrected`=='p < .001')){
          fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
            # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
          
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
          #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
          #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
          #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
          #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          pal1 <- palette(c('#F21A00','#3B9AB2'))
          gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
            geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +  + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
            # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
            geom_line(size = 1) + theme(legend.position = "none") +
            geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
            #scale_color_gradientn(colors = pal, guide = 'none') + 
            scale_color_manual(values = pal1,labels=c('AH','PH')) + 
            #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
            theme_bw(base_size=13) +
            theme(legend.title = element_blank(),
                  panel.grid.major = element_line(colour = "grey45"), 
                  panel.grid.minor = element_line(colour = "grey45"), 
                  panel.background = element_rect(fill = 'grey40'),
                  axis.title.y = element_text(margin=margin(r=6)),
                  axis.title.x = element_text(margin=margin(t=6)))
          print(gg)
          dev.off()
        } else {
          fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
            # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          print(gg)
          dev.off()
          
          fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
          pdf(fname, width = 9, height = 3.5)
          # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
          #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
          #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
          #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
          #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
          pal1 <- palette(c('#F21A00','#3B9AB2'))
          gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
            geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
            # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
            geom_line(size = 1) + theme(legend.position = "none") +
            geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
            #scale_color_gradientn(colors = pal, guide = 'none') + 
            scale_color_manual(values = pal1,labels=c('AH','PH')) + 
            #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
            theme_bw(base_size=13) +
            theme(legend.title = element_blank(),
                  panel.grid.major = element_line(colour = "grey45"), 
                  panel.grid.minor = element_line(colour = "grey45"), 
                  panel.background = element_rect(fill = 'grey40'),
                  axis.title.y = element_text(margin=margin(r=6)),
                  axis.title.x = element_text(margin=margin(t=6)))
          print(gg)
          dev.off()
        }
      } else if (strcmp(toprocess,'axis')){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        print(ggplot(edf, aes(t, as.factor(bin_num1))) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
                # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
        #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        gg<-ggplot(edf, aes(x=t, y=estimate,color = as.numeric(bin_num1), group = as.numeric(bin_num1))) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") +
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          scale_color_gradientn(colors = pal, guide = 'none') + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
        print(gg)
        dev.off()
      }
    }
  }
  return(Q)
}