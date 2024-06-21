source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory='/Volumes/Users/Andrew/MEDuSA_data_Explore',dataset='explore')
df <- df %>%
  group_by(id,run_number) %>%
  arrange(id, run_number, trial) %>%
  mutate(condition_trial = run_trial-floor(run_trial/40.5)*40,
         condition_trial_neg_inv = -1000 / condition_trial,
         condition_trial_neg_inv_sc = scale(as.vector(condition_trial_neg_inv)),
  ) %>% ungroup()
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
df <- df %>% select(iti_ideal,condition_trial_neg_inv_sc,rt_vmax_lag_sc, rt_lag_sc,iti_prev,iti_sc,v_entropy_wi_change,iti_prev_sc,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
df$id <- as.integer(df$id)

demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
demo <- demo %>% select(registration_lethality,registration_redcapid,age,gender,registration_group,wtar,education_yrs)
demo <- demo %>% rename(id=registration_redcapid,group=registration_group)
demo$gender <- relevel(as.factor(demo$gender),ref='M')
demo$age <- scale(demo$age)
demo$wtar <- scale(demo$wtar)
demo$education_yrs <- scale(demo$education_yrs)

df1 <- merge(df,demo,by='id')

df1 <- df1 %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                    trial > 40 & trial <= 80 ~ 2,
                                    trial > 80 & trial <=120 ~ 3, 
                                    trial > 120 & trial <=160 ~ 4,
                                    trial > 160 & trial <=200 ~ 5,
                                    trial > 200 & trial <=240 ~ 6))
df1 <- df1 %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                         trial > 40 & trial <= 80 ~ trial-40,
                                         trial > 80 & trial <=120 ~ trial-80, 
                                         trial > 120 & trial <=160 ~ trial-120,
                                         trial > 160 & trial <=200 ~ trial-160,
                                         trial > 200 & trial <=240 ~ trial-200))
df1 <- df1 %>% mutate(run_trial0_c = run_trial0-floor(run_trial0/40.5),
                  run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                  run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))
df1 <- df1 %>% mutate(ATT_3groups = case_when(group=='ATT' & registration_lethality=='hl' ~ 'ATT-HL', group=='ATT' & registration_lethality=='ll' ~ 'ATT-LL', group!='ATT' ~ 'non-ATT'))
df1 <- df1 %>% mutate(groupO = case_when(last_outcome=="Omission" & ATT_3groups=='ATT-HL' ~ 'O-HL', 
                                last_outcome=="Omission" & ATT_3groups=='ATT-LL' ~ 'O-LL', 
                                last_outcome=="Omission" & ATT_3groups=='non-ATT' ~ 'O-NON', 
                                last_outcome=="Reward" & ATT_3groups=='ATT-HL' ~ 'R-HL',
                                last_outcome=="Reward" & ATT_3groups=='ATT-LL' ~ 'R-LL', 
                                last_outcome=="Reward" & ATT_3groups=='non-ATT' ~ 'R-NON'))

df1 <- df1 %>% group_by(id,run) %>% mutate(rt_lag = lag(rt_csv), rt_swing = rt_csv - rt_lag) %>% ungroup()
