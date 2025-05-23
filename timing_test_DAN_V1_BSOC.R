library(stringr)
library(fmri.pipeline)
library(tidyverse)
library(pracma)
library(fmri.pipeline)

PPCcaudal <- c('280','281','282','275','277','73','75','77','78','80','81','82')
V1 <- c('7','207')
# load BSOC
rt_bsoc <- read.csv('/Volumes/Users/Andrew/DAN/fixed DAN/feedback_aligned_bsocial_dan.csv.gz')
rt_bsoc <- rt_bsoc %>% select(id,run,trial,evt_time,decon_mean, atlas_value)

df <- read_csv('/Volumes/Users/Andrew/DAN/bsocial_clock_trial_df.csv')
df <- df %>% select(id,scanner_run,trial,iti_prev,rt_bin,iti_ideal,rt_csv)
df <- df %>% rename(run=scanner_run)
df <- df %>% group_by(id,run)

df$id <- as.integer(df$id)
df <- df %>% mutate(iti_ideal_bin = case_when(iti_ideal < 1 ~ '< 1s',
                                              iti_ideal >= 1 & iti_ideal < 3.1 ~ '1-3.1s',
                                              iti_ideal >=3.1 & iti_ideal < 7.1 ~ '3.1-7.1s',
                                              iti_ideal >=7.1 ~ '> 7.1s'))
rt_bsoc <- rt_bsoc %>% mutate(run0 = as.numeric(gsub("run", "", run))) %>% select(!run) %>% rename(run = run0)
rt_bsoc <- inner_join(df,rt_bsoc,by=c('id','run','trial'))
rt_bsoc <- rt_bsoc %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                                  trial > 40 & trial <= 80 ~ trial-40,
                                                  trial > 80 & trial <=120 ~ trial-80, 
                                                  trial > 120 & trial <=160 ~ trial-120,
                                                  trial > 160 & trial <=200 ~ trial-160,
                                                  trial > 200 & trial <=240 ~ trial-200))



#waiting on BSOC clock scan info
#scaninfo <- read_csv('/Volumes/Users/Andrew/v18-2024-12-04/HC_vPFC_Explore_Clock_Scanner_Dates.csv')
#scaninfo <- scaninfo %>% mutate(id = registration_redcapid, ddate = as.vector(t(scale(difftime(scaninfo$scan_date,min(scaninfo$scan_date)))))) %>% select(!registration_redcapid)

#rt_bsoc <- inner_join(rt_bsoc,scaninfo,by='id')
#rt_bsoc$scan_which <- relevel(as.factor(rt_bsoc$scan_which),ref='P1')
rt_bsoc <- rt_bsoc %>% mutate(experiment = 'bsocial')

rtppc <- rt_bsoc %>% filter(atlas_value %in% PPCcaudal) %>% mutate(region = 'PPCcaudal')
rtppc <- rtppc %>% group_by(id) %>% mutate(id_run = paste0(id,'_',run)) %>% ungroup()
rtppc <- rtppc %>% group_by(id,run) %>% mutate(rt_next_bin = lead(rt_bin), rt_next = lead(rt_csv)) %>% ungroup()
rtppc <- rtppc %>% mutate(iti_rt_next = case_when(iti_ideal+rt_next < 2 ~ '< 2s',
                                                      iti_ideal+rt_next >= 2 & iti_ideal+rt_next < 3.4 ~ '2-3.4s',
                                                      iti_ideal+rt_next >=3.4 & iti_ideal+rt_next < 5.8 ~ '3.4-5.8s',
                                                      iti_ideal+rt_next >=5.8 ~ '> 5.8s'))

#   group_by(evt_time,id,run,early_trials) %>% 
#   summarize(mD = mean(decon_mean,na.rm=TRUE)) %>% ungroup()
# 
# rt <- rt %>% group_by(id,run) %>% arrange(evt_time) %>% mutate(peak = max(mD, na.rm=TRUE), ipeak = evt_time[mD==max(mD,na.rm=TRUE)]) %>% ungroup()


gc()


setwd('/Users/dnplserv/DAN')
#rt_comb$scan_which[is.na(rt_comb$scan_which)] = 'P1'

ncores = 26
decode_formula = formula(~1|(id_run))
splits = c('evt_time','rt_bin')
print(decode_formula)
ddf <- mixed_by(rtppc, outcomes = "decon_mean", rhs_model_formulae = decode_formula , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))




PPCcaudal <- c('280','281','282','275','277','73','75','77','78','80','81','82')
V1 <- c('4','7','5','9','11','204','206','202','209','211')

# load MMClock
load('/Users/dnplserv/DAN/MMClock_rt_dan_tall_ts.RData')
rt_comb <- rt_comb %>% 
  filter(atlas_value %in% PPCcaudal) %>% 
  select(id,run,run_trial,evt_time,iti_prev,iti_ideal,decon_interp,rt_csv,atlas_value,rt_bin)
rt_comb <- rt_comb %>% rename(decon_mean = decon_interp)

rt_comb <- rt_comb %>% mutate(iti_ideal_bin = case_when(iti_ideal < 1 ~ '< 1s',
                                              iti_ideal >= 1 & iti_ideal < 3.1 ~ '1-3.1s',
                                              iti_ideal >=3.1 & iti_ideal < 7.1 ~ '3.1-7.1s',
                                              iti_ideal >=7.1 ~ '> 7.1s'))
rt_comb <- rt_comb %>% group_by(id) %>% mutate(id_run = paste0(id,'_',run)) %>% ungroup()
rt_comb <- rt_comb %>% group_by(id,run) %>% mutate(rt_next_bin = lead(rt_bin), rt_next = lead(rt_csv)) %>% ungroup()
rt_comb <- rt_comb %>% mutate(iti_rt_next = case_when(iti_ideal+rt_next < 2 ~ '< 2s',
                                                      iti_ideal+rt_next >= 2 & iti_ideal+rt_next < 3.4 ~ '2-3.4s',
                                                      iti_ideal+rt_next >=3.4 & iti_ideal+rt_next < 5.8 ~ '3.4-5.8s',
                                                      iti_ideal+rt_next >=5.8 ~ '> 5.8s'))
ncores = 26
decode_formula = formula(~1|(id_run))
splits = c('evt_time','iti_rt_next')
print(decode_formula)
ddf1 <- mixed_by(rt_comb, outcomes = "decon_mean", rhs_model_formulae = decode_formula , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
