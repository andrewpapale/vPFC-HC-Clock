V1 <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','101','102','103','104','105','106','107','108','109','110','111','112','113','114','115')

# load Explore
rt_exp <- read.csv('/Volumes/Users/Andrew/MEDuSA_data_MMClock_Hallquist/feedback_aligned_transformed_DAN_7Network.csv.gz')
rt_exp <- rt_exp %>% filter(atlas_value %in% V1) %>% 
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

decode_formula <- NULL
decode_formula[[1]] = formula(~(1|id))
ncores <- 26
splits = c('evt_time','atlas_value')
setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
df0 <- decode_formula[[1]]
print(df0)
ddf <- mixed_by(rt_exp, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf,file=paste0(curr_date,'-Explore-Hallquist-Pipeline-timing-test-',1,'.Rdata'))


V1 <- c('4','7','5','9','11','204','206','202','209','211')
# load Explore
rt_exp <- read.csv('/Volumes/Users/Andrew/DAN/fixed DAN/feedback_aligned_bsocial_dan.csv.gz')
rt_exp <- rt_exp %>% filter(atlas_value %in% V1) %>% 
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

decode_formula <- NULL
decode_formula[[1]] = formula(~(1|id))
ncores <- 26
splits = c('evt_time','atlas_value')
setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
df0 <- decode_formula[[1]]
print(df0)
ddf <- mixed_by(rt_exp, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
save(ddf,file=paste0(curr_date,'-Explore-Fmriprep-Pipeline-timing-test-',1,'.Rdata'))
