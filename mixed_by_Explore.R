
ncores = 26

source('/Volumes/Users/Andrew/MEDuSA_data_Explore/get_trial_data_explore.R')

load('/Volumes/Users/Andrew/MEDuSA_data_Explore/fb-vPFC.Rdata')

df <- get_trial_data_explore(repo_directory='/Volumes/Users/Andrew/MEDuSA_data_Explore',censor_early_trials=FALSE)
df <- df %>% rename(run=run_number)
df <- df %>%
  group_by(id) %>% 
  mutate(v_entropy_sc_r = scale(v_entropy)) %>% ungroup() %>%
  group_by(id, run) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
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
df <- df %>% select(iti_ideal,iti_prev,iti_sc,outcome,last_outcome,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_sc,v_max_wi,trial_neg_inv_sc,trial)
df$id <- as.character(df$id)
Q <- full_join(md,df,by=c('id','run','trial'))

Q$decon_mean[Q$evt_time > Q$iti_ideal] = NA
Q <- Q %>% mutate(network = case_when(
  atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
  atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
  atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM'
))
Q <- Q %>% rename(vmPFC_decon = decon_mean)
Q <- Q %>% arrange(id,run,trial,evt_time)
Q <- Q %>% filter(evt_time > -5 & evt_time < 5)

Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')

demo <- read.csv('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_clock_demo.csv')
demo <- demo %>% select(registration_redcapid,age_at_scan,registration_gender,registration_group)
demo <- demo %>% rename(id=registration_redcapid,age=age_at_scan,gender=registration_gender,group=registration_group)
demo$gender <- relevel(as.factor(demo$gender),ref='M')
demo$age <- scale(demo$age)

Q <- merge(demo,Q,by='id')
Q <- Q %>% filter(group=='HC')

rm(decode_formula)
decode_formula <- formula(~ (1|id))
decode_formula[[1]] = formula(~ age + gender + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1|id/run))
decode_formula[[2]] = formula(~ age + gender + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc |id/run))
decode_formula[[3]] = formula(~ age + gender + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | id) + (1 | run))
decode_formula[[4]] = formula(~ age + gender + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | run) + (1| id))

decode_formula[[5]] = formula(~ age + gender + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1|id/run))
decode_formula[[6]] = formula(~ age + gender + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi |id/run))
decode_formula[[7]] = formula(~ age + gender + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | id) + (1 | run))
decode_formula[[8]] = formula(~ age + gender + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | run) + (1| id))

# decode_formula[[1]] = formula(~ v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1|id/run))
# decode_formula[[2]] = formula(~ v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc |id/run))
# decode_formula[[3]] = formula(~ v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | id) + (1 | run))
# decode_formula[[4]] = formula(~ v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | run) + (1| id))
# 
# decode_formula[[5]] = formula(~ v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1|id/run))
# decode_formula[[6]] = formula(~ v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi |id/run))
# decode_formula[[7]] = formula(~ v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | id) + (1 | run))
# decode_formula[[8]] = formula(~ v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | run) + (1| id))


qT <- c(-0.8,0.46)
splits = c('evt_time','network')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  if (i < 5){
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)#,
                    # emmeans_spec = list(
                    #   H = list(outcome='vmPFC_decon', model_name='model1',
                    #            specs=c("v_entropy_sc"), at = list(v_entropy_sc=c(-1.5,1.5))),
                    #   Tr = list(outcome='vmPFC_decon', model_name='model1',
                    #             specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT))
                    # ),
                    # emtrends_spec = list(
                    #   Tr = list(outcome='vmPFC_decon',model_name='model1', var = 'trial_neg_inv_sc',
                    #             specs = formula(~trial_neg_inv_sc),at=list(trial_neg_inv_sc = qT)),
                    #   H = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_sc',
                    #            specs = formula(~v_entropy_sc),at=list(v_entropy_sc=c(-1.5,1.5)))
                    # )
    )
  } else{
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)#,
                    # emmeans_spec = list(
                    #   Tr = list(outcome='vmPFC_decon', model_name='model1', 
                    #             specs=c("trial_neg_inv_sc"), at = list(trial_neg_inv_sc=qT)),
                    #   V = list(outcome='vmPFC_decon', model_name='model1',
                    #            specs=c('v_max_wi'), at=list(v_max_wi=c(-1.5,1.5)))
                    # ),
                    # emtrends_spec = list(
                    #   Tr = list(outcome='vmPFC_decon',model_name='model1', var = 'trial_neg_inv_sc',
                    #             specs = formula(~trial_neg_inv_sc),at=list(trial_neg_inv_sc = qT)),
                    #   V = list(outcome='vmPFC_decon',model_name='model1', var = 'v_max_wi',
                    #            specs = formula(~v_max_wi),at=list(v_max_wi=c(-1.5,1.5)))
                    # )
    )       
  }
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-Explore-HC-vmPFC-network-feedback-',i,'.Rdata'))
}

