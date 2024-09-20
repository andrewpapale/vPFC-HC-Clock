# 2024-09-20 AndyP
# subject_level_random_slopes_age




source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% select(rt_vmax_lag,iti_prev,ev,iti_ideal,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,iti_ideal,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
         iti_prev_sc = scale(iti_prev),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_sc = scale(rt_vmax),
         ev_sc = scale(ev),
         v_max_wi_lag = lag(v_max_wi),
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

df$id <- as.character(df$id)
Q2 <- df;
Q2$trial_bin <- factor(Q2$trial_bin)
Q2$trial_bin <- relevel(Q2$trial_bin, ref='Middle')
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)
Qfmri <- inner_join(Q2,demo,by=c('id'))
Qfmri$female <- relevel(as.factor(Qfmri$female),ref='0')
Qfmri$age <- scale(Qfmri$age)


df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_meg')
df <- df %>% select(rt_vmax_lag,ev,score_csv,v_max,outcome,v_entropy,rt_lag,v_entropy_full,v_entropy_wi_full,rt_vmax_full,rt_vmax_change_full,rt_csv_sc,rt_csv,id, run, run_trial, last_outcome, trial_neg_inv_sc,pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi,kld4_lag,kld4,rt_change,total_earnings, rewFunc,rt_csv, pe_max,v_chosen,rewFunc,
                    rt_vmax_lag_sc,rt_vmax_change,outcome,pe_max,kld3_lag,rt_lag_sc,rt_next,v_entropy_wi_change,pe_max_lag) %>% 
  group_by(id, run) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         score_sc = scale(score_csv),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_sc = scale(rt_vmax),
         ev_sc = scale(ev),
         v_max_wi_lag = lag(v_max_wi),
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

df$id <- as.character(df$id)
Q2 <- df;
Q2$trial_bin <- factor(Q2$trial_bin)
Q2$trial_bin <- relevel(Q2$trial_bin, ref='Middle')
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)
Qmeg <- inner_join(Q2,demo,by=c('id'))
Qmeg$female <- relevel(as.factor(Qmeg$female),ref='0')
Qmeg$age <- scale(Qmeg$age)

Qmeg <- Qmeg %>% mutate(dataset = 'MEG')
Qfmri <- Qfmri %>% mutate(dataset = 'fMRI')
Q3 <- rbind(Qmeg,Qfmri)


decode_formula <- NULL
decode_formula[[1]] <- formula(~(rt_lag_sc + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:age + rt_vmax_lag_sc * age * trial_neg_inv_sc + (1 | id/run))
decode_formula[[2]] <- formula(~(trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + subj_level_rand_slope + last_outcome)^2 + rt_lag_sc:last_outcome:age + rt_vmax_lag_sc * trial_neg_inv_sc * age + (1 + rt_vmax_lag_sc + rt_lag_sc | id/run))
qVL <- quantile(df$v_max_wi_lag,c(0.1,0.9),na.rm=TRUE)
qRTV <- quantile(df$rt_vmax_lag,c(0.1,0.9),na.rm=TRUE)
qH <- NULL
qH[1] <- mean(df$v_entropy_wi,na.rm=TRUE) - 2*sd(df$v_entropy_wi,na.rm=TRUE)
qH[2] <- mean(df$v_entropy_wi,na.rm=TRUE) + 2*sd(df$v_entropy_wi,na.rm=TRUE)
qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
qRT <- quantile(df$rt_lag_sc,c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
#qH <- c(1,2,3,4)
qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)
qRS <- quantile(Q$estimate, c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
splits = c('dataset')
source('~/fmri.pipeline/R/mixed_by.R')
print(i)
for (j in 1:length(decode_formula)){
  
  ddq <- mixed_by(Q, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),
                    Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                                specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2),rt_lag_sc=c(-2,-1,0,1,2))),        
                    TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                    
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                    Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                specs=formula(~rt_vmax_lag_sc:subj_level_rand_slope), at = list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                specs=formula(~rt_lag_sc:last_outcome:subj_level_rand_slope), at=list(subj_level_rand_slope=c(-2,-1,0,1,2))),
                    TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
                                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                    TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
                                    specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(subj_level_rand_slope = c(-2,-1,0,1,2),rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                    TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'subj_level_rand_slope',
                                    specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:subj_level_rand_slope), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                  )
  )
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection')
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  if (j==1){
    save(ddq,file=paste0(curr_date,'-Age-vPFC-HC-network-ranslopes-',toalign,'-replication-pred-int-',i,'.Rdata'))
  } else {
    save(ddq,file=paste0(curr_date,'-Age-vPFC-HC-network-ranslopes-',toalign,'-replication-pred-slo-',i,'.Rdata')) 
  } 
  
}

