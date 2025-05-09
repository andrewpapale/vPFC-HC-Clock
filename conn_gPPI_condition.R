

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory='/Volumes/Users/Andrew/MEDuSA_data_Explore',dataset='explore')
df <- df %>% filter(id == '202200' & run==1)

df0 <- df %>% select(clock_onset,rt_csv)
df0 <- df0 %>% rename(onsets = clock_onset, durations = rt_csv) %>% mutate(conditions = 1, subject_number = 1, sess= 1)
df0 <- df0[,c(3,4,5,1,2)]

setwd('/Users/dnplserv/gPPI')
write.csv(df0,'clock_onset.csv',row.names=FALSE)


df0 <- df %>% select(feedback_onset)
df0 <- df0 %>% rename(onsets = feedback_onset) %>% mutate(durations = 0.7,conditions = 1, subject_number = 1, sess= 1)
df0 <- df0[,c(3,4,5,1,2)]

setwd('/Users/dnplserv/gPPI')
write.csv(df0,'feedback_onset.csv',row.names=FALSE)
