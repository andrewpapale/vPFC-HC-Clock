# 2025-05-09 AndyP
# Generate AFNI formatted time files for a parametric regressor of entropy

library(tidyverse)

repo_directory <- "~/clock_analysis"

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='explore')
df0 <- df %>% filter(id == '202200' & run == 1) %>% select(clock_onset, rt_csv)

str <- NULL; for (iD in 1:nrow(df0)){str <- paste0(str,paste0(df0$clock_onset[iD],":",df0$rt_csv[iD],' '))}