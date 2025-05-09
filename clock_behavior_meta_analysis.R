# 2025-02-28 AndyP
# clock behavior meta-analysis


library(tidyverse)
library(fmri.pipeline)

rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
repo_directory <- file.path('/Volumes/Users/Andrew/MEDuSA_data_BSOC')
repo_directory_mmclock <- "~/clock_analysis"
ncores <- 26

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

df <- df %>% select(!trial_neg_inv_sc) %>% rename(trial_neg_inv_sc = run_trial0_neg_inv_sc)
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

#df <- df %>% filter(group=='HC')
df <- df %>% mutate(sex1 = case_when(sex == 1 ~ 'F', sex == 2 ~ 'M'))
df <- df %>% select(!sex) %>% rename(sex=sex1)
df$sex <- relevel(as.factor(df$sex),ref='F')
df$age <- scale(df$age)
#df$age <- scale(df$age)

fits <- read_csv('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL/fMRIEmoClock_decay_factorize_selective_psequate_fixedparams_fmri_mfx_sceptic_global_statistics.csv')
fits$id<-gsub("_1","",fits$id)
fits <- fits %>% select(!dataset)
df <- inner_join(df, fits,by='id')
df <- df %>% mutate(dataset = 'Bsocial')

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
dfmmc <- get_trial_data(repo_directory=repo_directory_mmclock,dataset='mmclock_fmri')
demo <- read.table(file=file.path(repo_directory_mmclock, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.double(demo$id)
dfmmc <- inner_join(dfmmc,demo,by=c('id'))
dfmmc <- dfmmc %>% mutate(dataset = 'MMClock fMRI')

stats <- read_csv('/Users/dnplserv/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats <- stats %>% select(!dataset)
dfmmc <- inner_join(dfmmc,stats,by='id')
dfmmc <- dfmmc %>% mutate(sex = case_when(female==1 ~ 'F',
                                          female==0 ~ 'M'))
dfmmc$sex <- relevel(as.factor(dfmmc$sex),ref='F')
dfmmc$age <- scale(dfmmc$age)

df <- df %>% select(dataset,id,alpha_ffx,beta_ffx,gamma_ffx,age,sex,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings)
dfmmc <- dfmmc %>% select(dataset,id,alpha_ffx,beta_ffx,gamma_ffx,age,sex,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings)

dfmmc_meg <- get_trial_data(repo_directory=repo_directory_mmclock,dataset='mmclock_meg')
dfmmc_meg <- dfmmc_meg %>% mutate(dataset = 'MMClock MEG')

stats <- read_csv('/Users/dnplserv/clock_analysis/meg/data/mmclock_meg_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv')
stats <- stats %>% select(!dataset)
stats <- stats %>% mutate(id2 = str_extract(id, "^[0-9]+")) %>% select(!id) %>% rename(id = id2)

dfmmc_meg <- inner_join(dfmmc_meg,stats,by=c('id'))

demo <- read.table(file=file.path(repo_directory_mmclock, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
demo$id <- as.character(demo$id)
dfmmc_meg <- inner_join(dfmmc_meg,demo,by='id')

dfmmc_meg <- dfmmc_meg %>% mutate(sex = case_when(female==1 ~ 'F',
                                          female==0 ~ 'M'))
dfmmc_meg$sex <- relevel(as.factor(dfmmc_meg$sex),ref='F')
dfmmc_meg <- dfmmc_meg %>% select(dataset,id,alpha_ffx,beta_ffx,gamma_ffx,age,sex,rt_lag_sc,v_entropy_wi,rt_csv_sc,rt_vmax_lag_sc,trial_neg_inv_sc,run,last_outcome,v_max_wi,total_earnings)



dfmer <- rbind(df,dfmmc)
dfmer1 <- rbind(dfmer,dfmmc_meg)

# splits <- c('dataset')
# decode_formula <- NULL
# decode_formula[[1]] <- formula(~(rt_lag_sc + age + last_outcome)^2 + rt_lag_sc:last_outcome:age + rt_vmax_lag_sc * age * trial_neg_inv_sc + (1 | id/run))
# ddq_age <- mixed_by(dfmer, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[1]], split_on = splits,return_models=TRUE,
#                 padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                 tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
#                 emmeans_spec = list(
#                   RT = list(outcome='rt_csv_sc', model_name='model1', 
#                             specs=formula(~rt_lag_sc:age), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
#                   Vmax = list(outcome='rt_csv_sc', model_name='model1', 
#                               specs=formula(~rt_vmax_lag_sc:age), at = list(rt_vmax_lag_sc=c(-2,-1,0,1,2))),
#                   RTxO = list(outcome='rt_csv_sc',model_name='model1',
#                               specs=formula(~rt_lag_sc:last_outcome:age), at=list(rt_lag_sc=c(-2,-1,0,1,2))),        
#                   TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
#                                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
#                   
#                 ),
#                 emtrends_spec = list(
#                   RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
#                             specs=formula(~rt_lag_sc:age), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
#                   Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
#                               specs=formula(~rt_vmax_lag_sc:age), at=list(rt_vmax_lag_sc = c(-2,-1,0,1,2))),
#                   RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
#                               specs=formula(~rt_lag_sc:last_outcome:age), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
#                   TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
#                                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
#                   TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
#                                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2)))#,
#                   # TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'age',
#                   #                 specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
#                 )
# )
# 
# decode_formula[[1]] <- formula(~(rt_lag_sc + sex + last_outcome)^2 + rt_lag_sc:last_outcome:sex + rt_vmax_lag_sc * sex * trial_neg_inv_sc + (1 | id/run))
# ddq_sex <- mixed_by(dfmer, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[1]], split_on = splits,return_models=TRUE,
#                 padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
#                 tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
#                 emmeans_spec = list(
#                   RT = list(outcome='rt_csv_sc', model_name='model1', 
#                             specs=formula(~rt_lag_sc:sex), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
#                   Vmax = list(outcome='rt_csv_sc', model_name='model1', 
#                               specs=formula(~rt_vmax_lag_sc:sex), at = list(rt_vmax_lag_sc=c(-2,-1,0,1,2))),
#                   RTxO = list(outcome='rt_csv_sc',model_name='model1',
#                               specs=formula(~rt_lag_sc:last_outcome:sex), at=list(rt_lag_sc=c(-2,-1,0,1,2))),        
#                   TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
#                                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:sex), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
#                   
#                 ),
#                 emtrends_spec = list(
#                   RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
#                             specs=formula(~rt_lag_sc:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
#                   Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
#                               specs=formula(~rt_vmax_lag_sc:sex), at=list(rt_vmax_lag_sc = c(-2,-1,0,1,2))),
#                   RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
#                               specs=formula(~rt_lag_sc:last_outcome:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
#                   TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
#                                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:sex), at= list(trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
#                   TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
#                                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:sex), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2)))#,
#                   # TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'age',
#                   #                 specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
#                 )
# )



dfmer0 <- dfmer1 %>% group_by(id,dataset) %>% summarize(alpha_ffx = alpha_ffx[1],beta_ffx=beta_ffx[1],gamma_ffx=gamma_ffx[1],total_earnings=total_earnings[1],age = age[1],sex=sex[1],dataset=dataset[1]) %>% ungroup()

# df$age <- scale(df$age)
# dfmmc$age <- scale(dfmmc$age)
# 
# m_alpha <- lm(data=df, alpha_ffx ~ age*sex)
# m_beta <- lm(data=df, beta_ffx ~ age*sex)
# m_gamma <- lm(data=df, gamma_ffx ~ age*sex)
# m_tE <- lm(data=df, scale(total_earnings) ~ age*sex)
# 
# m_alpha0 <- lm(data=dfmmc, alpha_ffx ~ age*sex)
# m_beta0 <- lm(data=dfmmc, beta_ffx ~ age*sex)
# m_gamma0 <- lm(data=dfmmc, gamma_ffx ~ age*sex)
# m_tE0 <- lm(data=dfmmc, scale(total_earnings) ~ age*sex)
# 
# m_alpha1 <- lm(data=dfmer0, alpha_ffx ~ age*sex*dataset)
# m_beta1 <- lm(data=dfmer0, beta_ffx ~ age*sex*dataset)
# m_gamma1 <- lm(data=dfmer0, gamma_ffx ~ age*sex*dataset)
# m_tE1 <- lm(data=dfmer0, scale(total_earnings) ~ age*sex*dataset)

m_mmc_fmri_alpha <- lm(data = dfmer0 %>% filter(dataset=='MMClock fMRI'), alpha_ffx ~ sex)
m_mmc_meg_alpha <- lm(data = dfmer0 %>% filter(dataset=='MMClock MEG'), alpha_ffx ~ sex)
m_bsoc_fmri_alpha <- lm(data = dfmer0 %>% filter(dataset=='Bsocial'), alpha_ffx ~ sex)
m_meta_alpha <- lmerTest::lmer(data=dfmer0, alpha_ffx ~ sex + (1|dataset),control = lmerControl(optimizer='bobyqa',check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

m_mmc_fmri_beta <- lm(data = dfmer0 %>% filter(dataset=='MMClock fMRI'), beta_ffx ~ age)
m_mmc_meg_beta <- lm(data = dfmer0 %>% filter(dataset=='MMClock MEG'), beta_ffx ~ age)
m_bsoc_fmri_beta <- lm(data = dfmer0 %>% filter(dataset=='Bsocial'), beta_ffx ~ age)
m_meta_beta <- lmerTest::lmer(data=dfmer0, beta_ffx ~ age*sex + (1|dataset),control = lmerControl(optimizer='bobyqa',check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

m_mmc_fmri_gamma <- lm(data = dfmer0 %>% filter(dataset=='MMClock fMRI'), gamma_ffx ~ age)
m_mmc_meg_gamma <- lm(data = dfmer0 %>% filter(dataset=='MMClock MEG'), gamma_ffx ~ age)
m_bsoc_fmri_gamma <- lm(data = dfmer0 %>% filter(dataset=='Bsocial'), gamma_ffx ~ age)
m_meta_gamma <- lmerTest::lmer(data=dfmer0, gamma_ffx ~ age*sex + (1|dataset),control = lmerControl(optimizer='bobyqa',check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

m_mmc_fmri_tE <- lm(data = dfmer0 %>% filter(dataset=='MMClock fMRI'), total_earnings ~ age*sex)
m_mmc_meg_tE <- lm(data = dfmer0 %>% filter(dataset=='MMClock MEG'), total_earnings ~ age*sex)
m_bsoc_fmri_tE <- lm(data = dfmer0 %>% filter(dataset=='Bsocial'), total_earnings ~ age*sex)
m_meta_tE <- lmerTest::lmer(data=dfmer0, total_earnings ~ age*sex + (1|dataset),control = lmerControl(optimizer='bobyqa'))

dfmer2 <- dfmer1 %>% group_by(dataset) %>% mutate(alpha_ffx = scale(alpha_ffx), beta_ffx = scale(beta_ffx), gamma_ffx = scale(gamma_ffx), total_earnings = scale(total_earnings)) %>% ungroup()

splits <- c('dataset')
decode_formula <- NULL
decode_formula[[1]] <- formula(~(rt_lag_sc + alpha_ffx + last_outcome)^2 + rt_lag_sc:last_outcome:alpha_ffx + rt_vmax_lag_sc * alpha_ffx * trial_neg_inv_sc + (1 | id/run))
ddq_alpha <- mixed_by(dfmer2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[1]], split_on = splits,return_models=TRUE,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE)#,
                    # emmeans_spec = list(
                    #   RT = list(outcome='rt_csv_sc', model_name='model1', 
                    #             specs=formula(~rt_lag_sc:alpha_ffx), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
                    #   Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                    #               specs=formula(~rt_vmax_lag_sc:alpha_ffx), at = list(rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                    #   RTxO = list(outcome='rt_csv_sc',model_name='model1',
                    #               specs=formula(~rt_lag_sc:last_outcome:alpha_ffx), at=list(rt_lag_sc=c(-2,-1,0,1,2))),        
                    #   TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                    #                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:alpha_ffx), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                    #   
                    # ),
                    # emtrends_spec = list(
                    #   RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                    #             specs=formula(~rt_lag_sc:alpha_ffx), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    #   Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                    #               specs=formula(~rt_vmax_lag_sc:alpha_ffx), at=list(rt_vmax_lag_sc = c(-2,-1,0,1,2))),
                    #   RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                    #               specs=formula(~rt_lag_sc:last_outcome:alpha_ffx), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    #   TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
                    #                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:alpha_ffx), at= list(trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                    #   TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
                    #                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:alpha_ffx), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2)))#,
                    #   # TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'age',
                    #   #                 specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                    # )
)
decode_formula <- NULL
decode_formula[[1]] <- formula(~(rt_lag_sc + beta_ffx + last_outcome)^2 + rt_lag_sc:last_outcome:beta_ffx + rt_vmax_lag_sc * beta_ffx * trial_neg_inv_sc + (1 | id/run))

ddq_beta <- mixed_by(dfmer2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[1]], split_on = splits,return_models=TRUE,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE)#,
                    # emmeans_spec = list(
                    #   RT = list(outcome='rt_csv_sc', model_name='model1', 
                    #             specs=formula(~rt_lag_sc:beta_ffx), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
                    #   Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                    #               specs=formula(~rt_vmax_lag_sc:beta_ffx), at = list(rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                    #   RTxO = list(outcome='rt_csv_sc',model_name='model1',
                    #               specs=formula(~rt_lag_sc:last_outcome:beta_ffx), at=list(rt_lag_sc=c(-2,-1,0,1,2))),        
                    #   TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                    #                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:beta_ffx), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                    #   
                    # ),
                    # emtrends_spec = list(
                    #   RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                    #             specs=formula(~rt_lag_sc:beta_ffx), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    #   Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                    #               specs=formula(~rt_vmax_lag_sc:beta_ffx), at=list(rt_vmax_lag_sc = c(-2,-1,0,1,2))),
                    #   RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                    #               specs=formula(~rt_lag_sc:last_outcome:beta_ffx), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    #   TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
                    #                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:beta_ffx), at= list(trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                    #   TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
                    #                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:beta_ffx), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2)))#,
                    #   # TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'age',
                    #   #                 specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                    # )
)

decode_formula <- NULL
decode_formula[[1]] <- formula(~(rt_lag_sc + gamma_ffx + last_outcome)^2 + rt_lag_sc:last_outcome:gamma_ffx + rt_vmax_lag_sc * gamma_ffx * trial_neg_inv_sc + (1 | id/run))

ddq_gamma <- mixed_by(dfmer2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[1]], split_on = splits,return_models=TRUE,
                     padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                     tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE)#,
                     # emmeans_spec = list(
                     #   RT = list(outcome='rt_csv_sc', model_name='model1', 
                     #             specs=formula(~rt_lag_sc:gamma_ffx), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
                     #   Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                     #               specs=formula(~rt_vmax_lag_sc:gamma_ffx), at = list(rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                     #   RTxO = list(outcome='rt_csv_sc',model_name='model1',
                     #               specs=formula(~rt_lag_sc:last_outcome:gamma_ffx), at=list(rt_lag_sc=c(-2,-1,0,1,2))),        
                     #   TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                     #                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:gamma_ffx), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                     #   
                     # ),
                     # emtrends_spec = list(
                     #   RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                     #             specs=formula(~rt_lag_sc:gamma_ffx), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                     #   Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                     #               specs=formula(~rt_vmax_lag_sc:gamma_ffx), at=list(rt_vmax_lag_sc = c(-2,-1,0,1,2))),
                     #   RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                     #               specs=formula(~rt_lag_sc:last_outcome:gamma_ffx), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                     #   TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
                     #                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:gamma_ffx), at= list(trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                     #   TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
                     #                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:gamma_ffx), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2)))#,
                     #   # TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'age',
                     #   #                 specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                     # )
)

decode_formula <- NULL
decode_formula[[1]] <- formula(~(rt_lag_sc + total_earnings + last_outcome)^2 + rt_lag_sc:last_outcome:total_earnings + rt_vmax_lag_sc * gamma_ffx * trial_neg_inv_sc + (1 | id/run))

ddq_total_earnings <- mixed_by(dfmer2, outcomes = "rt_csv_sc", rhs_model_formulae = decode_formula[[1]], split_on = splits,return_models=TRUE,
                     padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                     tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE)#,
                     # emmeans_spec = list(
                     #   RT = list(outcome='rt_csv_sc', model_name='model1', 
                     #             specs=formula(~rt_lag_sc:total_earnings), at = list(rt_lag_sc=c(-2,-1,0,1,2))),
                     #   Vmax = list(outcome='rt_csv_sc', model_name='model1', 
                     #               specs=formula(~rt_vmax_lag_sc:total_earnings), at = list(rt_vmax_lag_sc=c(-2,-1,0,1,2))),
                     #   RTxO = list(outcome='rt_csv_sc',model_name='model1',
                     #               specs=formula(~rt_lag_sc:last_outcome:total_earnings), at=list(rt_lag_sc=c(-2,-1,0,1,2))),        
                     #   TrxVmax = list(outcome='rt_csv_sc',model_name='model1',
                     #                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:total_earnings), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                     #   
                     # ),
                     # emtrends_spec = list(
                     #   RT = list(outcome='rt_csv_sc', model_name='model1', var='rt_lag_sc', 
                     #             specs=formula(~rt_lag_sc:total_earnings), at=list(total_earnings = c(-2,-1,0,1,2))),
                     #   Vmax = list(outcome='rt_csv_sc', model_name='model1', var='rt_vmax_lag_sc', 
                     #               specs=formula(~rt_vmax_lag_sc:total_earnings), at=list(total_earnings = c(-2,-1,0,1,2))),
                     #   RTxO = list(outcome='rt_csv_sc',model_name='model1',var='rt_lag_sc',
                     #               specs=formula(~rt_lag_sc:last_outcome:total_earnings), at=list(total_earnings = c(-2,-1,0,1,2))),
                     #   TrxVmax = list(outcome='rt_csv_sc',model_name='model1', var = 'rt_vmax_lag_sc',
                     #                  specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:total_earnings), at= list(trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4),total_earnings = c(-1,-1,0,1,2))),
                     #   TrxVmax1 = list(outcome='rt_csv_sc',model_name='model1', var = 'trial_neg_inv_sc',
                     #                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:total_earnings), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),total_earnings = c(-1,-1,0,1,2)))#,
                     #   # TrxVmax2 = list(outcome='rt_csv_sc',model_name='model1', var = 'age',
                     #   #                 specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:age), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2),trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4)))
                     # )
)
