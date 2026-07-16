# 2026-07-16 AndyP
# Sex and aging paper vmPFC-HC
# Simple script to get the average RT for MMClock (fMRI) and BSOCIAL clock

library(tidyverse)
library(fmri.pipeline)

rootdir <- '/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Sex_differences_Exploration_Paper'
repo_directory <- '/Users/andypapale/Documents/GitHub/clock_analysis/'
setwd(rootdir)

# Get task behav data
df <- read_csv(file.path(rootdir,'bsocial_clock_trial_df.csv'))
split_ksoc_bsoc <- df %>% group_by(id) %>% summarize(maxT = max(trial)) %>% ungroup()
ksoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==300])
bsoc <- data.frame(id = split_ksoc_bsoc$id[split_ksoc_bsoc$maxT==240])
bsoc <- rbind(bsoc,221973,221507,221842,440223)
df_bsoc <- df %>% filter(id %in% bsoc$id) %>% mutate(block = case_when(trial <= 40 ~ 1, 
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

df_ksoc <- df %>% filter(id %in% ksoc$id) %>% mutate(block = case_when(trial <= 50 ~ 1, 
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
         rt_swing_sc = scale(rt_swing)) %>% ungroup()

df <- df %>% mutate(task = 'Study 2')
df <- df %>% select(id,trial,run_trial,task,rt_csv,rewFunc)

demo <- read_csv(file.path(rootdir,'bsoc_clock_N171_dfx_demoonly.csv')) %>% mutate(sex = case_when(registration_birthsex == 1 ~ 'F',
                                                                                                   registration_birthsex == 2 ~ 'M'))
demo <- demo %>% select(id,sex)
demo$id <- as.character(demo$id)
df <- inner_join(df,demo,by='id')

source('/Users/andypapale/Documents/GitHub/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
dfmmc <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri') %>% mutate(task = 'Study 1')
dfmmc <- dfmmc %>% select(id,trial,run_trial,task,rt_csv,rewFunc)

demo <- read_tsv('/Users/andypapale/Documents/GitHub/clock_analysis/fmri/data/mmy3_demographics.tsv')
demo <- demo %>% select(lunaid,female) %>% rename(id = lunaid) %>%
  mutate(sex = case_when(female == 1 ~ 'F',female == 0 ~ 'M')) %>%
  select(!female)

dfmmc <- inner_join(dfmmc,demo,by='id')

df <- rbind(df,dfmmc)

colors <- c('#009b72','#d65b26','#6f6ab0','#e50886')

pdf('Figure_1C_RT_by_Trial.pdf',height=6,width=10)
gg1 <- ggplot(data=df %>% filter(rewFunc %in% c('IEV','DEV')),aes(x=run_trial,y=rt_csv,color=rewFunc,group=interaction(sex,rewFunc),linetype = as.factor(sex))) + geom_smooth(method = "loess",span = 1.25) +
  facet_grid(~task) + theme_minimal() + 
  scale_y_continuous(breaks = c(1.2,1.6,2.0,2.4)) + 
  ylab('Response Time (s)') +
  xlab('Trial') +
  scale_color_manual(values =c('IEV'='#009b72','DEV'='#d65b26','CEV'='#6f6ab0','CEVR'='#e50886')) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 36))
print(gg1)
dev.off()
