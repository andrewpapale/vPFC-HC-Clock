# 2025-06-09 AndyP
# plot mixed_by sex effect




load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection/2025-06-09-Sex-clock-fmri_meg-pred-rt-int-1.Rdata')

ddq_mmc <- ddq

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection/2025-06-09-Bsocial-Sex-clock-fmri-pred-rt-int.Rdata')

ddq_bsoc <- ddq


dqRT_bsoc <- ddq_bsoc$emtrends_list$RT %>% filter(dataset1=='bsocial')
dqRT_bsoc <- dqRT_bsoc %>% rename(dataset = dataset1)
dqRT_mmc <- ddq_mmc$emtrends_list$RT

ddqRT <- rbind(dqRT_bsoc,dqRT_mmc)

ddqRT <- ddqRT %>% mutate(dataset1 = case_when(dataset == 'fMRI' ~ 'Experiment 1 - fMRI',
                                               dataset == 'MEG' ~ 'Experiment 1 - MEG',
                                               dataset == 'bsocial' ~ 'Experiment 2')) %>% select(!dataset) %>% rename(dataset = dataset1)


fills <- palette()
fills[1] <- '#E618E6'
fills[2] <- '#181FE6'

setwd('/Volumes/Users/Andrew/2025-Research-Day-Poster-AgeSex')
pdf('RT_swings_by_sex.pdf',height=8,width=12)
gg1 <-ggplot(data=ddqRT, aes(x=sex,y=rt_lag_sc.trend,color=sex,group=sex,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  geom_errorbar(width=0.65,size=2) + geom_point(size=5) + facet_grid(~dataset,scales = 'free_y') + scale_y_reverse() + 
  scale_colour_manual(values=fills) +  
  ylab('<-- less -- Exploration -- more -->') +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=24),
        axis.title.x = element_text(margin=margin(t=6),size=24),
        legend.text = element_text(size=24),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        strip.text = element_text(size=24))
print(gg1)
dev.off()

repo_directory <- "~/clock_analysis"
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
repo_directory1 <- file.path('/Volumes/Users/Andrew/MEDuSA_data_BSOC')

mmcdemo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
mmcdemo <- mmcdemo %>% rename(id=lunaid)
mmcdemo <- mmcdemo %>% select(!adult & !scandate)
mmcdemo$id <- as.character(mmcdemo$id)


demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(repo_directory1,'2025-02-27-Partial-demo-pull-KSOC.csv'))
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
bsocdemo <- rbind(demo,demo1)
bsocdemo <- bsocdemo %>% select(id, age, sex, group) %>% mutate(dataset = 'Experiment 2')

mmcdemo <- mmcdemo %>% mutate(sex = case_when(female==0 ~ 'M',
                                    female==1 ~ 'F'))

mmcdemo <- mmcdemo %>% select(id,age,sex) %>% mutate(dataset = 'Experiment 1 - fMRI', group = 'HC')



source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% group_by(id) %>% summarize(total_earnings = total_earnings[1]) %>% ungroup()
df$id <- as.character(df$id)
mmc_fmri <- inner_join(df,mmcdemo,by=c('id'))
mmc_fmri$total_earnings <- scale(mmc_fmri$total_earnings)

m1 <- lm(data = mmc_fmri, total_earnings ~ age + sex)

df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_meg')
df <- df %>% group_by(id) %>% summarize(total_earnings = total_earnings[1]) %>% ungroup()
df$id <- as.character(df$id)
mmc_meg <- inner_join(df,mmcdemo,by=c('id')) %>% mutate(dataset = 'Experiment 1 - MEG', group = 'HC')
mmc_meg$total_earnings <- scale(mmc_meg$total_earnings)

m2 <- lm(data = mmc_meg, total_earnings ~ age + sex)

df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
df <- df %>% group_by(id) %>% summarize(total_earnings = total_earnings[1]) %>% ungroup()
df$id <- as.character(df$id)
bsoc <- inner_join(df,bsocdemo,by=c('id'))
bsoc <- bsoc %>% mutate(sex1 = case_when(sex == 1 ~ 'F', sex == 2 ~ 'M'))
bsoc <- bsoc %>% select(!sex) %>% rename(sex=sex1)
bsoc$total_earnings <- scale(bsoc$total_earnings)

m3 <- lm(data = bsoc, total_earnings ~ age + sex + group)

demomer <- rbind(mmc_fmri,mmc_meg,bsoc)


setwd('/Volumes/Users/Andrew/2025-Research-Day-Poster-AgeSex')
pdf('total_earnings_by_sex.pdf',height=8,width=12)
gg1 <-ggplot(data=demomer, aes(x=sex,y=total_earnings,color=sex)) +
  geom_violin() + geom_jitter(size=3) + facet_grid(~dataset) + scale_y_reverse() + 
  scale_colour_manual(values=fills) +  
  ylab('Total Earnings (scaled)') +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=24),
        axis.title.x = element_text(margin=margin(t=6),size=24),
        legend.text = element_text(size=24),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        strip.text = element_text(size=24))
print(gg1)
dev.off()


mmcdemo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
mmcdemo <- mmcdemo %>% rename(id=lunaid)
mmcdemo <- mmcdemo %>% select(!adult & !scandate)
mmcdemo$id <- as.character(mmcdemo$id)


demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(repo_directory1,'2025-02-27-Partial-demo-pull-KSOC.csv'))
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
bsocdemo <- rbind(demo,demo1)
bsocdemo <- bsocdemo %>% select(id, age, sex, group) %>% mutate(dataset = 'Experiment 2')

mmcdemo <- mmcdemo %>% mutate(sex = case_when(female==0 ~ 'M',
                                              female==1 ~ 'F'))

mmcdemo <- mmcdemo %>% select(id,age,sex) %>% mutate(dataset = 'Experiment 1 - fMRI', group = 'HC')



source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
df$id <- as.character(df$id)
mmc_fmri <- inner_join(df,mmcdemo,by=c('id'))

df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_meg')
df$id <- as.character(df$id)
mmc_meg <- inner_join(df,mmcdemo,by=c('id')) %>% mutate(dataset = 'Experiment 1 - MEG', group = 'HC')

df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
df$id <- as.character(df$id)
bsoc <- inner_join(df,bsocdemo,by=c('id'))
bsoc <- bsoc %>% mutate(sex1 = case_when(sex == 1 ~ 'F', sex == 2 ~ 'M'))
bsoc <- bsoc %>% select(!sex) %>% rename(sex=sex1)

mmc_fmri <- mmc_fmri %>% select(id,run,sex,age,rt_lag_sc,last_outcome,rt_vmax_lag_sc,trial_neg_inv_sc,rt_csv)
mmc_meg <- mmc_meg %>% select(id,run,sex,age,rt_lag_sc,last_outcome,rt_vmax_lag_sc,trial_neg_inv_sc,rt_csv)
bsoc <- bsoc %>% select(id,run,sex,age,rt_lag_sc,last_outcome,rt_vmax_lag_sc,trial_neg_inv_sc,rt_csv)

mmc_fmri <- mmc_fmri %>% mutate(dataset = 'Experiment 1 - fMRI')
mmc_meg <- mmc_meg %>% mutate(dataset = 'Experiment 1 - MEG')
bsoc <- bsoc %>% mutate(dataset = 'Experiment 2')

dmer <- rbind(mmc_fmri,mmc_meg,bsoc)
dmer <- dmer %>% filter(rt_csv < 4 & rt_csv > 0.2)
dmer <- dmer %>% dplyr::mutate(reward_lag_rec = case_when(last_outcome=="Reward" ~ 0.5, last_outcome=="Omission" ~ -0.5))

dmer$sex <- as.factor(dmer$sex)

decode_formula <- NULL
decode_formula[[1]] <- formula(~rt_lag_sc*reward_lag_rec*sex + rt_vmax_lag_sc*trial_neg_inv_sc*sex + (1 + rt_lag_sc + rt_vmax_lag_sc | id) + (1 | run))


splits = c('dataset')

ddq <- mixed_by(dmer, outcomes = "rt_csv", rhs_model_formulae = decode_formula[[1]], split_on = splits,return_models=TRUE,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals"),conf.int=TRUE),
                  emtrends_spec = list(
                    RT = list(outcome='rt_csv', model_name='model1', var='rt_lag_sc', 
                            specs=formula(~rt_lag_sc:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    Vmax = list(outcome='rt_csv', model_name='model1', var='rt_vmax_lag_sc', 
                              specs=formula(~rt_vmax_lag_sc:sex), at=list(rt_vmax_lag_sc = c(-2,-1,0,1,2))),
                    RTxO = list(outcome='rt_csv',model_name='model1',var='rt_lag_sc',
                                specs=formula(~rt_lag_sc:reward_lag_rec:sex), at=list(rt_lag_sc = c(-2,-1,0,1,2))),
                    TrxVmax = list(outcome='rt_csv',model_name='model1', var = 'rt_vmax_lag_sc',
                                   specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:sex), at= list(trial_neg_inv_sc=c(-0.9,-0.02,0.2,0.34,0.4))),
                    TrxVmax1 = list(outcome='rt_csv',model_name='model1', var = 'trial_neg_inv_sc',
                                    specs=formula(~rt_vmax_lag_sc:trial_neg_inv_sc:sex), at= list(rt_vmax_lag_sc=c(-2,-1,0,1,2)))
                  )
)
 
ddf <- ddq$coef_df_reml %>% filter(effect == 'ran_vals' & term == 'rt_lag_sc') %>% rename(id = level)

mer2 <- inner_join(ddf,demomer,by=c('id','dataset'))

ggplot(data=mer2,aes(x=age,y=estimate,color=sex,group=sex)) + geom_point(size=1) + facet_wrap(~dataset) + geom_smooth(method='lm',formula=y~poly(x,2,raw=TRUE))
ggplot(data=mer2,aes(x=age,y=estimate,color=sex,group=sex)) + geom_point(size=1) + facet_wrap(~dataset) + geom_smooth(method='lm')


ddf <- ddq$coef_df_reml %>% filter(effect == 'ran_vals' & term == 'rt_vmax_lag_sc') %>% rename(id = level)

mer2 <- inner_join(ddf,demomer,by=c('id','dataset'))

mer2 <- mer2 %>% mutate(group1 = case_when(group.y == 'ATT' ~ 'ATT',
                                           group.y == 'DEP' ~ 'NON',
                                           group.y == 'DNA' ~ 'NON',
                                           group.y == 'HC' ~ 'HC',
                                           group.y == 'NON' ~ 'NON'))

ggplot(data=mer2,aes(x=age,y=estimate,color=sex,group=sex)) + geom_point(size=1) + facet_wrap(~dataset) + geom_smooth(method='lm',formula=y~poly(x,2,raw=TRUE))
ggplot(data=mer2,aes(x=age,y=estimate,color=sex,group=sex)) + geom_point(size=1) + facet_wrap(~dataset) + geom_smooth(method='lm')

ddf <- ddq$coef_df_reml %>% filter(effect == 'ran_vals' & term == 'rt_lag_sc') %>% rename(id = level, rt_lag_sc = estimate)
ddg <- ddq$coef_df_reml %>% filter(effect == 'ran_vals' & term == 'rt_vmax_lag_sc') %>% rename(id = level, rt_vmax_lag_sc = estimate)

mer3 <- inner_join(ddf,ddg,by=c('id','dataset'))
mer3 <- inner_join(mer3,demomer,by=c('id','dataset'))

mer3 <- mer3 %>% mutate(group1 = case_when(group == 'ATT' ~ 'ATT',
                                           group == 'DEP' ~ 'NON',
                                           group == 'DNA' ~ 'NON',
                                           group == 'HC' ~ 'HC',
                                           group == 'NON' ~ 'NON'))

ggplot(data=mer3,aes(x=age,y=rt_vmax_lag_sc,color=interaction(group1,sex),group=interaction(group1,sex),shape=group1)) + geom_point(size=1) + facet_wrap(~dataset) + geom_smooth(method='lm')


##################################
### RS HCwithin vs age / sex #####
##################################

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2025-06-16-Bsocial-vPFC-HC-network-clock-All-agesex-HCwithin-RS-2.Rdata')
ddq_bsoc <- ddf$coef_df_reml %>% filter(effect == 'ran_vals' & term == 'HCwithin') %>% filter(group=='id') %>% rename(id = level, HCwithin = estimate)

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2025-06-16-MMClock-vPFC-HC-network-clock-agesex-HCwithin-RS-2.Rdata')
ddq_mmc <- ddf$coef_df_reml %>% filter(effect == 'ran_vals' & term == 'HCwithin') %>% filter(group=='id') %>% rename(id = level, HCwithin = estimate)

ddq_bsoc <- ddq_bsoc %>% mutate(dataset = 'Experiment 2')
ddq_mmc <- ddq_mmc %>% mutate(dataset = 'Experiment 1')


repo_directory <- "~/clock_analysis"
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
repo_directory1 <- file.path('/Volumes/Users/Andrew/MEDuSA_data_BSOC')

mmcdemo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
mmcdemo <- mmcdemo %>% rename(id=lunaid)
mmcdemo <- mmcdemo %>% select(!adult & !scandate)
mmcdemo$id <- as.character(mmcdemo$id)

demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv'))
demo1 <- read_csv(file.path(repo_directory1,'2025-02-27-Partial-demo-pull-KSOC.csv'))
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
bsocdemo <- rbind(demo,demo1)

bsocdemo <- bsocdemo %>% mutate(sex1 = case_when(sex == 1 ~ 'F', sex == 2 ~ 'M'))
bsocdemo <- bsocdemo %>% select(!sex) %>% rename(sex=sex1)

mmcdemo <- mmcdemo %>% mutate(sex = case_when(female==0 ~ 'M',
                                              female==1 ~ 'F'))

ddq_bsoc <- inner_join(ddq_bsoc,bsocdemo,by=c('id')) %>% select(!group.x & !group.y & !gender)
ddq_mmc <- inner_join(ddq_mmc,mmcdemo,by=c('id')) %>% select(!group & !female)

ddqmer <- rbind(ddq_bsoc,ddq_mmc)

ddqmer0 <- ddqmer %>% group_by(id,dataset) %>% summarize(HCwithin_scaled = scale(HCwithin), age=age,sex=sex) %>% ungroup()
ggplot(data=ddqmer0,aes(x=age,y=HCwithin_scaled,color=sex,group=sex)) + facet_wrap(~dataset) + geom_point(size=0.1) + geom_smooth(method = 'lm') + ylim(-0.03,0.03) + ggtitle('AH-DMN slope vs age by sex')

ddqmer0 <- ddqmer %>% group_by(id,dataset) %>% summarize(HCwithin_mean = mean(HCwithin,na.rm=TRUE), age=age,sex=sex) %>% ungroup()
ggplot(data=ddqmer0,aes(x=age,y=HCwithin_mean,color=sex,group=sex)) + facet_wrap(~dataset) + geom_point(size=0.1) + geom_smooth(method = 'lm') + ylim(-0.03,0.03) + ggtitle('AH-DMN slope vs age by sex')
