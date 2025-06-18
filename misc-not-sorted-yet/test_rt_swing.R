library(tidyverse)
library(broom.mixed)
library(data.table)

source('/Users/andypapale/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')

df <- get_trial_data(repo_directory = '~/clock_analysis',dataset='mmclock_fmri')
df <- df %>% group_by(id,run) %>% mutate(rt_swing_sc = scale(rt_swing)) %>% ungroup()
window = seq(from=10,to=50,by=1);

uS = unique(df$id);
uR = 1:8;

rt_swing = NULL;
rt_swing_sc = NULL;
win_len = NULL;
i = 1;
mlist <- list()
rm(df0)
for (iW in window){
  print(iW)
  uT = seq(from=1,to=50,by=iW);
  #for (iS in uS){ 
    #for (iR in 1:8){ 
      for (iT in 1:length(uT)){
        df0 <- df %>% filter((run_trial>=uT[iT] & run_trial < min(uT[iT+1],50,na.rm=TRUE)))
        if (nrow(df0) > 0 && !all(is.na(df0$rt_vmax_lag_sc))){
          #mdf <- broom.mixed::tidy(lm(rt_csv_sc ~ rt_lag_sc + rt_vmax_lag_sc, data= df0));
          mdf <- broom.mixed::tidy(lmerTest::lmer(rt_csv_sc ~ rt_lag_sc + rt_vmax_lag_sc + (1+rt_lag_sc + rt_vmax_lag_sc | id), data= df0));
          mlist[[i]] <- mdf;
          rt_swing <- rbind(rt_swing,median(df0$rt_swing,na.rm=TRUE));
          rt_swing_sc <- rbind(rt_swing_sc,median(df0$rt_swing_sc,na.rm=TRUE));
          win_len <- rbind(win_len, iW);
          i = i + 1;
          rm(df0)
        }
      }
    #}
    #print(iS)
  #}
}
ddf <- rbindlist(mlist);
str(ddf)
