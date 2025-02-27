# 2025-02-27 AndyP
# plot mixed by effects for BSOC behavior

library(tidyverse)
library(fmri.pipeline)


# Get task behav data
df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
split_ksoc_bsoc <- df %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300]
bsoc <- split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240]
df_bsoc <- df %>% filter(id %in% bsoc) %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                                                    trial > 40 & trial <= 80 ~ 2,
                                                                    trial > 80 & trial <=120 ~ 3, 
                                                                    trial > 120 & trial <=160 ~ 4,
                                                                    trial > 160 & trial <=200 ~ 5,
                                                                    trial > 200 & trial <=240 ~ 6))
df_bsoc <- df_bsoc %>% mutate(run_trial0 = case_when(trial <= 40 ~ trial, 
                                                     trial > 40 & trial <= 80 ~ trial-40,
                                                     trial > 80 & trial <=120 ~ trial-80, 
                                                     trial > 120 & trial <=160 ~ trial-120,
                                                     trial > 160 & trial <=200 ~ trial-160,
                                                     trial > 200 & trial <=240 ~ trial-200))
df_bsoc <- df_bsoc %>% mutate(protocol = 'bsocial',
                              run_trial0_c = run_trial0-floor(run_trial0/40.5),
                              run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                              run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))

df_ksoc <- df %>% filter(id %in% ksoc) %>% mutate(block = case_when(trial <= 50 ~ 1, 
                                                                    trial > 50 & trial <= 100 ~ 2,
                                                                    trial > 100 & trial <=150 ~ 3, 
                                                                    trial > 150 & trial <=200 ~ 4,
                                                                    trial > 200 & trial <=250 ~ 5,
                                                                    trial > 250 & trial <=300 ~ 6))
df_ksoc <- df_ksoc %>% mutate(run_trial0 = case_when(trial <= 50 ~ trial, 
                                                     trial > 50 & trial <= 100 ~ trial-50,
                                                     trial > 100 & trial <=150 ~ trial-100, 
                                                     trial > 150 & trial <=200 ~ trial-150,
                                                     trial > 200 & trial <=250 ~ trial-200,
                                                     trial > 250 & trial <=300 ~ trial-250))
df_ksoc <- df_ksoc %>% mutate(protocol = 'ksocial',
                              run_trial0_c = run_trial0-floor(run_trial0/50.5),
                              run_trial0_neg_inv = -(1 / run_trial0_c) * 100,
                              run_trial0_neg_inv_sc = as.vector(scale(run_trial0_neg_inv)))

df <- rbind(df_bsoc,df_ksoc)
# select and scale variables of interest
df <- df %>% 
  group_by(id, scanner_run) %>% 
  mutate(id = as.character(id),
         v_chosen_sc = scale(v_chosen),
         score_sc = scale(score_csv),
         iti_sc = scale(iti_ideal),
         iti_lag_sc = scale(iti_prev),
         v_max_sc = scale(v_max),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         rt_swing_sc = scale(rt_swing),
         v_max_wi_lag = lag(v_max_wi)) %>% ungroup()

df <- df %>% select(!run) %>% rename(run = scanner_run) %>% select(!run_trial) %>% mutate(run_trial = case_when(trial <= 150 ~ trial,
                                                                                                     trial > 150 ~ trial - 150))


# add in age and sex variables
demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(repo_directory,'2025-02-27-Partial-demo-pull-KSOC.csv'))
demo$id <- as.character(demo$id)
demo1$id <- as.character(demo1$registration_redcapid)
demo <- demo %>% rename(sex=registration_birthsex,
                        gender=registration_gender,
                        group=registration_group) %>%
  select(id,group,age,sex,gender)
demo1 <- demo1 %>% rename(sex=registration_birthsex,
                          gender=registration_gender,
                          group=registration_group) %>%
  select(id,group,age,sex,gender)
demo2 <- rbind(demo,demo1)

df <- inner_join(df,demo2,by='id')

df <- df %>% filter(group=='HC')
df <- df %>% mutate(sex1 = case_when(sex == 1 ~ 'F', sex == 2 ~ 'M'))
df <- df %>% select(!sex) %>% rename(sex=sex1)
df$sex1 <- relevel(as.factor(df$sex),ref='F')

fits <- read_csv('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL/fMRIEmoClock_decay_factorize_selective_psequate_fixedparams_fmri_mfx_sceptic_global_statistics.csv')
fits$id<-gsub("_1","",fits$id)

df <- inner_join(df, fits,by='id')

df <- df %>% mutate(dataset1='bsocial')
df0 <- df %>% mutate(dataset1 = 'bsocial_copy')
df <- rbind(df,df0)

decode_formula <- NULL
decode_formula[[1]] <- formula(~(rt_lag_sc + sex + last_outcome)^2 + rt_lag_sc:last_outcome:sex + rt_vmax_lag_sc * sex * trial_neg_inv_sc + (1 | id/run))
decode_formula[[2]] <- formula(~(trial_neg_inv_sc + rt_lag_sc + v_max_wi_lag + v_entropy_wi + age + last_outcome)^2 + rt_lag_sc:last_outcome:age + rt_vmax_lag_sc * trial_neg_inv_sc * age + (1 | id/run))
splits = c('dataset1')
for (j in 1:length(decode_formula)){
  
  ddq <- mixed_by(df, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[j]], split_on = splits,return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emmeans_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', 
                              specs=formula(~rt_lag_sc:sex), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
                    Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                                specs=formula(~rt_vmax_lag_sc:sex), at = list(rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',
                                specs=formula(~rt_lag_sc:last_outcome:sex), at=list(rt_lag_sc=c(-2,-1,0,1,2))),        
                    TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:sex), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                    
                  ),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                              specs=formula(~rt_lag_sc:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                                specs=formula(~rt_vmax_lag_sc:sex), at=list(rt_vmax_lag_sc = c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                                specs=formula(~rt_lag_sc:last_outcome:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
                                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:sex), at= list(trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                    TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
                                    specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:sex), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2)))#,
                    # TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'age',
                    #                 specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                  )
  )
  setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection')
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  if (j==1){
    save(ddq,file=paste0(curr_date,'Bsocial-Sex-clock-fmri-pred-rt-int','.Rdata'))
  } else {
    save(ddq,file=paste0(curr_date,'Bsocial-Age-clock-fmri-pred-rt-int','.Rdata')) 
  } 
  
}
