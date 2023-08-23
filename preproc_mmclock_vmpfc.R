# preproc MMClock vmPFC

clock_comb <- read_csv('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMClock_vmPFC_clock_Aug2023.csv.gz')
clock_comb <- clock_comb %>% mutate(run1 = case_when(run=='run1'~1,
                                     run=='run2'~2,
                                     run=='run3'~3,
                                     run=='run4'~4,
                                     run=='run5'~5,
                                     run=='run6'~6,
                                     run=='run7'~7,
                                     run=='run8'~8)) %>% select(!run) %>% rename(run=run1)

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory='~/clock_analysis',dataset='mmclock_fmri')
df <- df %>% select(id,run,trial,run_trial)
clock_comb <- inner_join(clock_comb,df,by=c('id','run','trial'))
clock_comb <- clock_comb %>% mutate(network = case_when(
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
save(clock_comb, file='MMclock_clock_Aug2023.Rdata')

fb_comb <- read_csv('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMClock_vmPFC_fb_Aug2023.csv.gz')
fb_comb <- fb_comb %>% mutate(run1 = case_when(run=='run1'~1,
                                     run=='run2'~2,
                                     run=='run3'~3,
                                     run=='run4'~4,
                                     run=='run5'~5,
                                     run=='run6'~6,
                                     run=='run7'~7,
                                     run=='run8'~8)) %>% select(!run) %>% rename(run=run1)
source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory='~/clock_analysis',dataset='mmclock_fmri')
df <- df %>% select(id,run,trial,run_trial)
fb_comb <- inner_join(fb_comb,df,by=c('id','run','trial'))
fb_comb <- fb_comb %>% mutate(network = case_when(
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
save(fb_comb, file='MMclock_fb_Aug2023.Rdata')
