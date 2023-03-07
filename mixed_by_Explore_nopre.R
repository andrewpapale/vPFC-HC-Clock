
ncores = 26

source('/Volumes/Users/Andrew/MEDuSA_data_Explore/get_trial_data_explore.R')

load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')

df <- get_trial_data_explore(repo_directory='/Volumes/Users/Andrew/MEDuSA_data_Explore',censor_early_trials=TRUE,trials_to_censor=10)
df <- df %>% rename(run=run_number)
df <- df %>%
  group_by(id) %>% 
  mutate(v_entropy_sc_r = scale(v_entropy)) %>% ungroup() %>%
  group_by(id, run) %>% 
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
         iti_prev_sc = scale(iti_prev),
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
df <- df %>% select(iti_ideal,iti_prev,iti_sc,v_entropy_wi_change,iti_prev_sc,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
df$id <- as.character(df$id)
Q <- full_join(md,df,by=c('id','run','trial'))

Q$decon_mean[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
Q$decon_mean[Q$evt_time < -Q$iti_prev] = NA
Q <- Q %>% mutate(network = case_when(
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
Q <- Q %>% rename(vmPFC_decon = decon_mean)
Q <- Q %>% arrange(id,run,trial,evt_time)
Q <- Q %>% filter(evt_time > -4 & evt_time < 4)


Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')

demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
demo <- demo %>% select(registration_redcapid,age,gender,registration_group,wtar,education_yrs)
demo <- demo %>% rename(id=registration_redcapid,group=registration_group)
demo$gender <- relevel(as.factor(demo$gender),ref='M')
demo$age <- scale(demo$age)
demo$wtar <- scale(demo$wtar)
demo$education_yrs <- scale(demo$education_yrs)

Q <- merge(demo,Q,by='id')
Q <- Q %>% filter(group!='ATT')
Q <- Q %>% mutate(block = case_when(trial <= 40 ~ 1, 
                            trial > 40 & trial <= 80 ~ 2,
                            trial > 80 & trial <=120 ~ 3, 
                            trial > 120 & trial <=160 ~ 4,
                            trial > 160 & trial <=200 ~ 5,
                            trial > 200 & trial <=240 ~ 6))
#Q <- Q %>% filter(group=='HC')

rm(decode_formula)
decode_formula <- NULL
decode_formula[[1]] = formula(~ v_entropy_wi + trial_neg_inv_sc + rt_csv_sc + gender + wtar + education_yrs + age + iti_sc + iti_prev_sc + last_outcome*iti_prev_sc + outcome*rt_csv_sc + (1 |id/run))
#decode_formula[[2]] = formula(~ age + gender + wtar + education_yrs + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc |id/run))
#decode_formula[[3]] = formula(~ age + gender + wtar + education_yrs + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | id) + (1 | run))
#decode_formula[[4]] = formula(~ age + gender + wtar + education_yrs + v_entropy_sc + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | run) + (1| id))

decode_formula[[2]] = formula(~ v_max_wi + trial_neg_inv_sc + rt_csv_sc + gender + wtar + education_yrs + age + iti_sc + iti_prev_sc + last_outcome*iti_prev_sc + outcome*rt_csv_sc + (1 |id/run))
#decode_formula[[6]] = formula(~ age + gender + wtar + education_yrs + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi |id/run))
#decode_formula[[7]] = formula(~ age + gender + wtar + education_yrs + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | id) + (1 | run))
#decode_formula[[8]] = formula(~ age + gender + wtar + education_yrs + v_max_wi + trial_neg_inv_sc + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | run) + (1| id))

decode_formula[[3]] = formula(~age + wtar + education_yrs + v_max_wi + (1|id/run))

# decode_formula[[1]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1|id/run))
# decode_formula[[2]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc |id/run))
# decode_formula[[3]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | id) + (1 | run))
# decode_formula[[4]] = formula(~ age + gender + v_entropy_sc*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome +  (1 + v_entropy_sc | run) + (1| id))
# 
# decode_formula[[5]] = formula(~ age + gender + v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1|id/run))
# decode_formula[[6]] = formula(~ age + gender + v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi |id/run))
# decode_formula[[7]] = formula(~ age + gender + v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | id) + (1 | run))
# decode_formula[[8]] = formula(~ age + gender + v_max_wi*trial_bin + rt_bin + iti_sc + rt_vmax_change_sc + last_outcome + outcome + (1 + v_max_wi | run) + (1| id))

# decode_formula[[1]] = formula(~ v_entropy_sc + (1|id))
# decode_formula[[2]] = formula(~ v_max_wi + (1|id))

qT <- c(-0.8,0.46)
splits = c('evt_time','network')
source("~/fmri.pipeline/R/mixed_by.R")
for (i in 1:length(decode_formula)){
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  df0 <- decode_formula[[i]]
  print(df0)
  if (i < 2){
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                    emmeans_spec = list(
                      A = list(outcome='vmPFC_decon',model_name='model1',
                               specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
                      W = list(outcome='vmPFC_decon',model_name='model1',
                               specs=formula(~wtar),at=list(wtar=c(-1,-0.5,0,0.5,1))),
                      H = list(outcome='vmPFC_decon', model_name='model1',
                               specs=formula(~v_entropy_wi), at = list(v_entropy_wi=c(-2,-1,0,1,2))),
                      Y = list(outcome='vmPFC_decon',model_name='model1',
                               specs=formula(~education_yrs), at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
                    )#,
                    # emtrends_spec = list(
                    #   HxW = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
                    #              specs = formula(~v_entropy_wi:wtar),at=list(wtar = c(-1,-0.5,0,0.5,1))),
                    #   HxA = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
                    #             specs = formula(~v_entropy_wi:age),at=list(age = c(-1,-0.5,0,0.5,1))),
                    #   H = list(outcome='vmPFC_decon',model_name='model1', var = 'v_entropy_wi',
                    #            specs = formula(~v_entropy_wi),at=list(v_entropy_wi=c(-2,-1,0,1,2))),
                    #   HxY = list(outcome='vmPFC_decon',model_name='model1',var='v_entropy_wi',
                    #              specs=formula(~v_entropy_wi:education_yrs),at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
                    # )
    )
  } else{
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                    emmeans_spec = list(
                      A = list(outcome='vmPFC_decon',model_name='model1',
                               specs=formula(~age),at=list(age=c(-1,-0.5,0,0.5,1))),
                      W = list(outcome='vmPFC_decon',model_name='model1',
                               specs=formula(~wtar),at=list(wtar=c(-1,-0.5,0,0.5,1))),
                      V = list(outcome='vmPFC_decon', model_name='model1',
                               specs=formula(~v_max_wi), at = list(v_max_wi=c(-2,-1,0,1,2))),
                      Y = list(outcome='vmPFC_decon',model_name='model1',
                               specs=formula(~education_yrs), at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
                    )#,
                    # emtrends_spec = list(
                    #   HxW = list(outcome='vmPFC_decon',model_name='model1', var = 'v_max_wi',
                    #              specs = formula(~v_max_wi:wtar),at=list(wtar = c(-1,-0.5,0,0.5,1))),
                    #   HxA = list(outcome='vmPFC_decon',model_name='model1', var = 'v_max_wi',
                    #              specs = formula(~v_max_wi:age),at=list(age = c(-1,-0.5,0,0.5,1))),
                    #   HxY = list(outcome='vmPFC_decon',model_name='model1',var='v_max_wi',
                    #              specs=formula(~v_max_wi:education_yrs),at=list(education_yrs=c(-1,-0.5,0,0.5,1)))
                    # )
    )       
  }
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  save(ddf,file=paste0(curr_date,'-Explore-vmPFC-network-clock-nopre-',i,'.Rdata'))
}

