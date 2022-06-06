############
# 2022-05-23 AndyP
# Subject level random slopes

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores <- 26

#### clock ####

rm(Q)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
load(file.path(HC_cache_dir,'clock_hipp_tall_ts_1.Rdata'))
hc <- clock_comb
hc <- hc %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)
hc <- hc %>% mutate(
  HC_region = case_when(
    bin_num <= 8 ~ 'AH',
    bin_num >8 ~ 'PH'
  ),
)
hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
Q <- Q %>% select(!decon1)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld3,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
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
  run_trial <= 10 ~ 'Early',
  run_trial > 10 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)

Q$HCbetween <- scale(Q$HCbetween)

rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (1|id))
decode_formula[[2]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (v_entropy_wi|id))
decode_formula[[3]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (HCwithin|id))
decode_formula[[4]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (v_entropy_wi*HCwithin|id))
decode_formula[[5]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (v_entropy_wi:HCwithin|id))

qHC <- Q %>% filter(atlas_value==55) %>% group_by(evt_time,HC_region) %>% 
  summarize(HC_p2SD = mean(HCwithin,na.rm=TRUE)+2*sd(HCwithin,na.rm=TRUE),
            HC_m2SD = mean(HCwithin,na.rm=TRUE)-2*sd(HCwithin,na.rm=TRUE),
            HC_p1SD = mean(HCwithin,na.rm=TRUE)+1*sd(HCwithin,na.rm=TRUE),
            HC_m1SD = mean(HCwithin,na.rm=TRUE)-1*sd(HCwithin,na.rm=TRUE),
            HC_p05SD = mean(HCwithin,na.rm=TRUE)+0.5*sd(HCwithin,na.rm=TRUE),
            HC_m05SD = mean(HCwithin,na.rm=TRUE)-0.5*sd(HCwithin,na.rm=TRUE))

qHC <- qHC %>% filter(HC_region=='PH')
qHC1 <- NULL
qHC1[1] <- mean(qHC$HC_m2SD,na.rm=TRUE)
qHC1[2] <- mean(qHC$HC_p2SD,na.rm=TRUE)
qHC1[3] <- mean(qHC$HC_m1SD,na.rm=TRUE)
qHC1[4] <- mean(qHC$HC_p1SD,na.rm=TRUE)
qHC1[5] <- mean(qHC$HC_m05SD,na.rm=TRUE)
qHC1[6] <- mean(qHC$HC_p05SD,na.rm=TRUE)

splits = c('evt_time','network','HC_region')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  #qH <- c(1,2,3,4)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emtrends_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=formula(~v_entropy_wi*HCwithin), at = list(v_entropy_wi=qH)),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=formula(~trial_neg_inv_sc*HCwithin), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=formula(~v_max_wi*HCwithin), at=list(v_max_wi=c(qV)))
                  ), 
                  emmeans_spec = list(
                    H_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ v_entropy_wi * HCwithin), at=list(HCwithin=qHC1,v_entropy_wi=qH)),
                    T_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ trial_neg_inv_sc * HCwithin), at=list(HCwithin=qHC1,trial_neg_inv_sc=qT)),
                    V_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ v_max_wi * HCwithin), at=list(HCwithin=qHC1,v_max_wi=qV))
                  )
  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-ranslopes-',i,'.Rdata'))
}



setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
load('2022-06-03-vmPFC-HC-network-clock-ranslopes-4.Rdata')


### Does random slope predict rt_vmax or rt_swing?
rm(Q2)
qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')
qdf <- qdf %>% rename(id=level)
qdf$estimate <- scale(qdf$estimate)
qdf$id <- as.character(qdf$id)
qdf <- qdf %>% select(!outcome)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(id,run,rewFunc,run_trial,rt_swing,rt_vmax,trial_neg_inv_sc,v_entropy_wi,last_outcome,v_max_wi,kld3,rt_csv_sc,v_max_wi_lag,outcome,rt_vmax_lag,rt_lag_sc)
df$id <- as.character(df$id)
Q2 <- inner_join(qdf,df,by='id')
demo$id <- as.character(demo$id)
Q2 <- inner_join(Q2,demo,by=c('id'))
Q2$female <- relevel(as.factor(Q2$female),ref='0')
Q2$age <- scale(Q2$age)

rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula = formula(~ (trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + estimate + last_outcome)^2 +
                                rt_lag_sc:last_outcome:estimate +
                                rt_vmax_lag*trial_neg_inv_sc*estimate +
                                (1 | id/run))
qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
qH <- NULL
qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
#qH <- c(1,2,3,4)
qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
qRS <- quantile(Q2$estimate, c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
splits = c('evt_time','HC_region','network')
source('~/fmri.pipeline/R/mixed_by.R')
ddq <- mixed_by(Q2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula, split_on = splits,return_models=TRUE,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed"),conf.int=TRUE),
                emtrends_spec = list(
                  V = list(outcome='rt_csv_sc', model_name='model1', var='estimate', 
                              specs='v_max_wi_lag', at = list(v_max_wi_lag=qVL)),
                  H = list(outcome='rt_csv_sc', model_name='model1', var='estimate',
                           specs='v_entropy_wi', at=list(v_entropy_wi=qH)),
                  LO = list(outcome='rt_csv_sc', model_name='model1',var='estimate',
                            specs=formula(~last_outcome)),
                  Tr = list(outcome='rt_csv_sc', model_name='model1',var='estimate',
                            specs=formula(~trial_neg_inv_sc), at=list(trial_neg_inv_sc=qT)),
                  RT = list(outcome='rt_csv_sc', model_name='model1',var='estimate',
                            specs=formula(~rt_lag_sc), at=list(rt_lag_sc=qRT)),
                  RTxLO = list(outcome='rt_csv_sc',model_name='model1',var='estimate',
                               specs=formula(~rt_lag_sc:last_outcome),at=list(rt_lag_sc=qRT)),
                  RTvxTR = list(outcome='rt_csv_sc',model_name='model1',var='estimate',
                                specs=formula(~rt_vmax_lag*trial_neg_inv_sc),at=list(trial_neg_inv_sc=qT,rt_vmax_lag=qVL))
                )
                )
setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-pred-rt_csv_sc-','2','.Rdata'))


source('~/vmPFC/plot_subject_level_random_slopes.R')
source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
for (i in 2){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
  model_str <- paste0('-vmPFC-HC-network-ranslopes-pred-rt_csv_sc-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'rt_csv_sc-2'
  #toprocess <- 'symmetry-by-HC'
  toprocess <- 'network-by-HC'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}

#test rt_swing
splits = c('evt_time','HC_region','network')
ddq <- mixed_by(Q2, outcomes = "rt_swing", rhs_model_formulae = decode_formula[[1]] , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                emtrends_spec = list(
                  V_HC = list(outcome='rt_swing', model_name='model1', var='estimate', 
                              specs='v_max_wi', at = list(v_max_wi=qV)),
                  T_HC = list(outcome='rt_swing', model_name='model1', var='estimate', 
                              specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                  O_HC = list(outcome='rt_swing', model_name='model1', var='estimate',
                              specs=c('last_outcome'))
                )
)
setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-pred-rt_swing-',i,'.Rdata'))


source('~/vmPFC/plot_subject_level_random_slopes.R')
source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
  model_str <- paste0('-vmPFC-HC-network-ranslopes-pred-rt_swing-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'rt_swing'
  #toprocess <- 'symmetry-by-HC'
  toprocess <- 'network-by-HC'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}


# replication

### Does random slope predict rt_vmax or rt_swing?
rm(Q2)
qdf <- ddf$coef_df_reml %>% filter(effect=='ran_vals')
qdf$estimate <- scale(qdf$estimate)
qdf <- qdf %>% rename(id=level)
qdf$id <- as.character(qdf$id)
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_meg')
df <- df %>% select(id,run,rewFunc,run_trial,rt_swing,rt_vmax,trial_neg_inv_sc,last_outcome,v_max_wi,kld3,rt_csv_sc)
df$id <- as.character(df$id)
Q2 <- inner_join(qdf,df,by='id')
demo$id <- as.character(demo$id)
Q2 <- inner_join(Q2,demo,by=c('id'))
Q2$female <- relevel(as.factor(Q2$female),ref='0')
Q2$age <- scale(Q2$age)

rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (1|id))
decode_formula[[2]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (v_entropy_wi|id))
decode_formula[[3]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (HCwithin|id))
decode_formula[[4]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (v_entropy_wi*HCwithin|id))
decode_formula[[5]] = formula(~ age*HCwithin + female*HCwithin +  v_max_wi*HCwithin + 
                                trial_neg_inv_sc*HCwithin + 
                                v_entropy_wi*HCwithin + rt_bin + 
                                expl_shorter*HCwithin + expl_longer*HCwithin +  # binary expl_code incr / decr separate variables
                                HCbetween +
                                (v_entropy_wi:HCwithin|id))

qHC <- Q %>% filter(atlas_value==55) %>% group_by(evt_time,HC_region) %>% 
  summarize(HC_p2SD = mean(HCwithin,na.rm=TRUE)+2*sd(HCwithin,na.rm=TRUE),
            HC_m2SD = mean(HCwithin,na.rm=TRUE)-2*sd(HCwithin,na.rm=TRUE),
            HC_p1SD = mean(HCwithin,na.rm=TRUE)+1*sd(HCwithin,na.rm=TRUE),
            HC_m1SD = mean(HCwithin,na.rm=TRUE)-1*sd(HCwithin,na.rm=TRUE),
            HC_p05SD = mean(HCwithin,na.rm=TRUE)+0.5*sd(HCwithin,na.rm=TRUE),
            HC_m05SD = mean(HCwithin,na.rm=TRUE)-0.5*sd(HCwithin,na.rm=TRUE))

qHC <- qHC %>% filter(HC_region=='PH')
qHC1 <- NULL
qHC1[1] <- mean(qHC$HC_m2SD,na.rm=TRUE)
qHC1[2] <- mean(qHC$HC_p2SD,na.rm=TRUE)
qHC1[3] <- mean(qHC$HC_m1SD,na.rm=TRUE)
qHC1[4] <- mean(qHC$HC_p1SD,na.rm=TRUE)
qHC1[5] <- mean(qHC$HC_m05SD,na.rm=TRUE)
qHC1[6] <- mean(qHC$HC_p05SD,na.rm=TRUE)

splits = c('evt_time','network','HC_region')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  qH <- NULL
  qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
  qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
  #qH <- c(1,2,3,4)
  qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
  qV <- quantile(df$v_max_wi,c(0.1,0.9),na.rm=TRUE)
  qdH <- NULL
  qdH[1] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) - 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  qdH[2] <- mean(df$v_entropy_wi_change_lag,na.rm=TRUE) + 2*sd(df$v_entropy_wi_change_lag,na.rm=TRUE)
  ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                  emtrends_spec = list(
                    H_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=formula(~v_entropy_wi*HCwithin), at = list(v_entropy_wi=qH)),
                    T_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin', 
                                specs=formula(~trial_neg_inv_sc*HCwithin), at = list(trial_neg_inv_sc=c(qT))),
                    V_HC = list(outcome='vmPFC_decon', model_name='model1', var='HCwithin',
                                specs=formula(~v_max_wi*HCwithin), at=list(v_max_wi=c(qV)))
                  ), 
                  emmeans_spec = list(
                    H_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ v_entropy_wi * HCwithin), at=list(HCwithin=qHC1,v_entropy_wi=qH)),
                    T_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ trial_neg_inv_sc * HCwithin), at=list(HCwithin=qHC1,trial_neg_inv_sc=qT)),
                    V_HC = list(outcome='vmPFC_decon',model_name='model1',specs=formula(~ v_max_wi * HCwithin), at=list(HCwithin=qHC1,v_max_wi=qV))
                  )
  )
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-rt_csv_sc-ranslopes-replication-',i,'.Rdata'))
}


source('~/vmPFC/plot_subject_level_random_slopes.R')
source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
  model_str <- paste0('-vmPFC-HC-network-clock-rt_csv_sc-ranslopes-replication-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'rt_vmax-replication'
  #toprocess <- 'symmetry-by-HC'
  toprocess <- 'network-by-HC'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}

#test rt_swing
splits = c('evt_time','HC_region','network')
ddq <- mixed_by(Q2, outcomes = "rt_swing", rhs_model_formulae = decode_formula[[1]] , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                emtrends_spec = list(
                  V_HC = list(outcome='rt_swing', model_name='model1', var='estimate', 
                              specs='v_max_wi', at = list(v_max_wi=qV)),
                  T_HC = list(outcome='rt_swing', model_name='model1', var='estimate', 
                              specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=c(qT))),
                  O_HC = list(outcome='rt_swing', model_name='model1', var='estimate',
                              specs=c('last_outcome'))
                )
)
setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddq,file=paste0(curr_date,'-vmPFC-HC-network-ranslopes-pred-rt_swing-replication-',i,'.Rdata'))


source('~/vmPFC/plot_subject_level_random_slopes.R')
source('~/vmPFC/plot_emtrends_subject_level_random_slopes.R')
for (i in 1){
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
  model_str <- paste0('-vmPFC-HC-network-ranslopes-pred-rt_swing-replication-',i,'.Rdata')
  model_str <- Sys.glob(paste0('*',model_str))
  load(model_str)
  model_iter <- i
  totest <- 'rt_swing-replication'
  #toprocess <- 'symmetry-by-HC'
  toprocess <- 'network-by-HC'
  toalign <- 'clock'
  behavmodel <- 'compressed'
  hc_LorR <- 'LR'
  plot_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
  plot_emtrends_subject_level_random_slopes(ddq,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
}
