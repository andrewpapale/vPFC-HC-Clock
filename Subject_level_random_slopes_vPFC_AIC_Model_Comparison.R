# 2023-01-03 AndyP
# AIC Comparison for B2B models

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26
toalign <- 'clock'

do_rt_pred_fmri = TRUE
plot_aic = FALSE

if (do_rt_pred_fmri){
  for (i in 1:2){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str <- paste0('-vmPFC-network-',toalign,'-ranslopes-',i,'.Rdata')
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    
    ### Does random slope predict rt_vmax or rt_swing?
    rm(Q2)
    if (strcmp(toalign,'clock')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_sc' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change_lag' & group=='id')
      }
    } else if (strcmp(toalign,'feedback')){
      if (i==1){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_sc' & group=='id')
      } else if (i==2){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_max_wi' & group=='id')
      } else if (i==3){
        qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='v_entropy_wi_change' & group=='id')
      } 
    }
    qdf <- qdf %>% group_by(network) %>% mutate(estimate1 = scale(estimate)) %>% ungroup() %>% select(!estimate) %>% rename(estimate=estimate1)
    qdf <- qdf %>% rename(id=level)
    qdf$id <- as.character(qdf$id)
    qdf <- qdf %>% select(!outcome)
    source('~/vmPFC/get_trial_data_vmPFC.R')
    df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
    df <- df %>% select(rt_vmax_lag,v_max_wi_lag,ev,score_csv,v_max,last_outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                        v_max_wi, iti_ideal, iti_prev,v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,
                        rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
      group_by(id, run) %>% 
      mutate(rt_sec = rt_csv/1000) %>% ungroup() %>%
      mutate(v_chosen_sc = scale(v_chosen),
             abs_pe_max_sc = scale(abs(pe_max)),
             iti_sc = scale(iti_ideal),
             iti_prev_sc = scale(iti_prev),
             score_sc = scale(score_csv),
             pe_max_sc = scale(pe_max),
             pe_max_lag_sc = scale(lag(pe_max)),
             abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
             ev_sc = scale(ev),
             rt_vmax_sc = scale(rt_vmax),
             v_entropy_sc = scale(v_entropy),
             v_max_lag_sc = scale(lag(v_max)),
             v_entropy_wi_change_lag = lag(v_entropy_wi_change),
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
    df <- df %>% group_by(id,run)  %>% mutate(rt_bin = case_when(
      rt_csv_sc <= -1 ~ '-1',
      rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
      rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
      rt_csv_sc > 1 ~ '1'
    ),
    rt_lag_bin = case_when(
      rt_lag_sc <= -1 ~ '-1',
      rt_lag_sc > -1 & rt_lag_sc <= 0 ~ '-0.5',
      rt_lag_sc > 0 & rt_lag_sc <= 1 ~ '0.5',
      rt_lag_sc > 1 ~ '1'
    ))
    df <- df %>% group_by(id,run) %>% mutate(trial_bin = case_when(
      run_trial <= 15 ~ 'Early',
      run_trial > 15 & run_trial < 30 ~ 'Middle',
      run_trial >=30 ~ 'Late',
    ))
    
    Q2 <- merge(qdf,df,by=c('id'))
    Q2$trial_bin <- factor(Q2$trial_bin)
    Q2$trial_bin <- relevel(Q2$trial_bin, ref='Middle')
    Q2 <- Q2 %>% rename(subj_level_rand_slope=estimate)
    
    if (strcmp(toalign,'clock')){
      decode_formula = formula(~ v_entropy_wi + (1|id))
      decode_formula[[1]] = formula(~ trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome + outcome + iti_sc + iti_prev_sc +
                                      rt_lag_sc*subj_level_rand_slope*last_outcome + rt_vmax_lag_sc*subj_level_rand_slope*trial_neg_inv_sc + (1|id/run))
      decode_formula[[2]] = formula(~ trial_neg_inv_sc + rt_lag_bin + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome + outcome + iti_sc + iti_prev_sc +
                                      rt_lag_bin*subj_level_rand_slope*last_outcome + rt_vmax_lag_sc*subj_level_rand_slope*trial_neg_inv_sc + (1|id/run))
      decode_formula[[3]] = formula(~ trial_bin + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome + outcome + iti_sc + iti_prev_sc +
                                      rt_lag_sc*subj_level_rand_slope*last_outcome + rt_vmax_lag_sc*subj_level_rand_slope*trial_bin + (1|id/run))
      decode_formula[[4]] = formula(~ trial_bin + rt_lag_bin + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome + outcome + iti_sc + iti_prev_sc +
                                      rt_lag_bin*subj_level_rand_slope*last_outcome + rt_vmax_lag_sc*subj_level_rand_slope*trial_bin + (1|id/run))
      decode_formula[[5]] = formula(~ trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome + outcome + iti_sc + iti_prev_sc +
                                      rt_lag_sc*subj_level_rand_slope*last_outcome + rt_vmax_lag_sc*subj_level_rand_slope*trial_bin + (1|id/run))
      decode_formula[[6]] = formula(~ trial_bin + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome + outcome + iti_sc + iti_prev_sc +
                                      rt_lag_bin*subj_level_rand_slope*last_outcome + rt_vmax_lag_sc*subj_level_rand_slope*trial_bin + (1|id/run))
    } else if (strcmp(toalign,'feedback')){
      decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome + outcome + iti_sc + iti_prev_sc +
                                 rt_bin*subj_level_rand_slope*outcome + rt_vmax_sc*subj_level_rand_slope*trial_neg_inv_sc + (1|id/run))     
    }
    qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
    qRTV <- quantile(df$rt_vmax_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    qH <- NULL
    qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
    qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    #qH <- c(1,2,3,4)
    qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    qRS <- quantile(Q2$subj_level_rand_slope, c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
    splits = c('network','evt_time')
    source('~/fmri.pipeline/R/mixed_by.R')
    print(i)
    for (j in 5:length(decode_formula)){
      df0 <- decode_formula[[j]]
      if (strcmp(toalign,'clock')){
        ddq <- mixed_by(Q2, outcomes = "rt_csv", rhs_model_formulae = df0, split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed"),conf.int=TRUE)
        )
      } else if (strcmp(toalign,'feedback')){
        ddq <- mixed_by(Q2, outcomes = "rt_next", rhs_model_formulae = decode_formula, split_on = splits,return_models=TRUE,
                        padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                        tidy_args = list(effects=c("fixed"),conf.int=TRUE)
        )  
      }
      if (i==1){
        slrs <- 'entropy'
      } else if (i==2){
        slrs <- 'value'
      }
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddq,file=paste0(curr_date,'-vmPFC-network-ranslopes-',toalign,'-pred-rt_csv_sc-',slrs,'-',j,'.Rdata'))
    }
  }
}

if (plot_aic){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
  ddg <- NULL
  for (j in 1:6){
    load(paste0('2023-01-03-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-entropy-',j,'.Rdata'))
    ddf <- ddq$fit_df
    rm(ddq)
    ddf <- ddf %>% mutate(model = j)
    ddg <- rbind(ddg,ddf)
  }
  
  ddg$model <- as.factor(ddg$model)
  ddg <- ddg %>% mutate(levels = case_when(model==1 ~ 'trial and RT continuous',
                                           model==2 ~ 'RT binned, trial continuous',
                                           model==3 ~ 'trial binned, RT continuous',
                                           model==4 ~ 'RT and trial binned',                                           
                                           model==5 ~ 'trial both, RT continuous',
                                           model==6 ~ 'RT both, trial binned'
                                          ))
  ggplot(ddg, aes(x=evt_time,y=AIC,color=levels,group=levels)) + geom_point(size=5) + geom_line(size=1) + facet_wrap(~network)
  
  
  ddg <- NULL
  for (j in 1:6){
    load(paste0('2023-01-03-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-value-',j,'.Rdata'))
    ddf <- ddq$fit_df
    rm(ddq)
    ddf <- ddf %>% mutate(model = j)
    ddg <- rbind(ddg,ddf)
  }
  
  ddg$model <- as.factor(ddg$model)
  ddg <- ddg %>% mutate(levels = case_when(model==1 ~ 'trial and RT continuous',
                                           model==2 ~ 'RT binned, trial continuous',
                                           model==3 ~ 'trial binned, RT continuous',
                                           model==4 ~ 'RT and trial binned',
                                           model==5 ~ 'trial both, RT continuous',
                                           model==6 ~ 'RT both, trial binned'
                                           ))
  ggplot(ddg, aes(x=evt_time,y=AIC,color=levels,group=levels)) + geom_point(size=5) + geom_line(size=1) + facet_wrap(~network)  
}
