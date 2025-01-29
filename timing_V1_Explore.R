library(stringr)
library(fmri.pipeline)
library(tidyverse)
library(pracma)

PPCcaudal <- c('280','281','282','275','277','73','75','77','78','80','81','82')
#V1 <- c('4','7','5','9','11','204','206','202','209','211')


# load Explore
rt_exp <- read.csv('/Volumes/Users/Andrew/DAN/Explore_rt_aligned_444_dan.csv.gz')
rt_exp <- rt_exp %>% filter(atlas_value %in% PPCcaudal) %>% 
  select(id,run,trial,evt_time,decon_mean, atlas_value)

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory='~/clock_analysis',dataset='explore')
df <- df %>% select(id,run,trial,iti_prev,rt_csv)

df$id <- as.integer(df$id)
rt_exp <- rt_exp %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
rt_exp <- inner_join(df,rt_exp,by=c('id','run','trial'))
rt_exp <- rt_exp %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                         trial > 40 & trial <= 80 ~ trial-40,
                                         trial > 80 & trial <=120 ~ trial-80, 
                                         trial > 120 & trial <=160 ~ trial-120,
                                         trial > 160 & trial <=200 ~ trial-160,
                                         trial > 200 & trial <=240 ~ trial-200))
rt_exp <- rt_exp %>% select(!trial)

scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = as.vector(t(scale(difftime(scaninfo$scan_date,min(scaninfo$scan_date)))))) %>% select(!registration_redcapid)

rt_exp <- inner_join(rt_exp,scaninfo,by='id')
rt_exp$scan_which <- relevel(as.factor(rt_exp$scan_which),ref='P1')
rt_exp <- rt_exp %>% mutate(experiment = 'Explore')
rt_exp$scan_which[is.na(rt_exp$scan_which)] = 'P1'
