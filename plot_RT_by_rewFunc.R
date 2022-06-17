# 2022-06-15 AndyP
# Plot RT by rewFunc

library(tidyverse)

repo_directory <- '~/clock_analysis/'
pal = palette()
pal[1] <- '#009B72'
pal[2] <- '#D65B26'
pal[3] <- '#6F6AB0'
pal[4] <- '#E50886'

source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_fmri')

dq <- df %>% select(rt_csv,rewFunc,run_trial,v_entropy_wi,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split)
dq <- dq %>% group_by(run_trial,rewFunc) %>% summarize(mRT = mean(rt_csv,na.rm=TRUE), dRT = sd(rt_csv,na.rm=TRUE),N = length(rt_csv)) %>% ungroup()
dq <- dq %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RTbyRewFunc.pdf',height=9,width=12)
gg1 <- ggplot(dq, aes(x=run_trial,y=mRT,color=Contingency,group=Contingency)) + 
  geom_line(size=0.5,linetype=2) + scale_color_manual(values=pal) + 
  geom_smooth(method='loess',span=1) + ylab('Response Time (s)') + xlab('Trial') + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
print(gg1)
dev.off()


dq1 <- df %>% select(rt_swing,rewFunc,run_trial,v_entropy_wi,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split)
dq1 <- dq1 %>% mutate(entropy_bin = ntile(v_entropy_wi,50))
dq1 <- dq1 %>% group_by(rewFunc,entropy_bin) %>% summarize(mRT = mean(rt_swing,na.rm=TRUE), dRT = sd(rt_swing,na.rm=TRUE),N = length(rt_swing)) %>% ungroup()
dq1 <- dq1 %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
dq1 <- dq1 %>% filter(!is.na(entropy_bin))

hbin_edge = quantile(df$v_entropy_wi,probs=seq(0,1,length.out=50),na.rm=TRUE)
Hbin_edge <- as_tibble(hbin_edge) %>% rename(hbin_edge=value)
Hbin_edge <- Hbin_edge %>% mutate(entropy_bin = seq(1,50,length.out=50))

dq1 <- inner_join(dq1,Hbin_edge,by='entropy_bin')

dq1 <- dq1 %>% filter(entropy_bin > 2 & entropy_bin < 48)

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RTbyRewFunc_H.pdf',height=9,width=12)
gg1 <- ggplot(dq1, aes(x=hbin_edge,y=mRT,color=Contingency,group=Contingency)) + 
  geom_line(size=0.5,linetype=2) + scale_color_manual(values=pal) + 
  geom_smooth(method='loess',span=1) + ylab('Change in RT (s)') + xlab('Entropy') + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
print(gg1)
dev.off()



dq1 <- df %>% select(rt_swing,rewFunc,run_trial,v_entropy_wi,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split)
dq1 <- dq1 %>% mutate(entropy_bin = ntile(v_entropy_wi,50))
dq1 <- dq1 %>% group_by(rewFunc,entropy_bin,run_trial) %>% summarize(mRT = mean(rt_swing,na.rm=TRUE), dRT = sd(rt_swing,na.rm=TRUE),N = length(rt_swing)) %>% ungroup()
dq1 <- dq1 %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
dq1 <- dq1 %>% filter(!is.na(entropy_bin))

hbin_edge = quantile(df$v_entropy_wi,probs=seq(0,1,length.out=50),na.rm=TRUE)
Hbin_edge <- as_tibble(hbin_edge) %>% rename(hbin_edge=value)
Hbin_edge <- Hbin_edge %>% mutate(entropy_bin = seq(1,50,length.out=50))


setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_swing_by_RewFunc.pdf',height=9,width=12)
gg1 <- ggplot(dq2, aes(x=run_trial,y=mRT,color=Contingency,group=Contingency)) + 
  geom_line(size=0.5,linetype=2) + scale_color_manual(values=pal) + 
  geom_smooth(method='loess',span=1) + ylab('Change in RT (s)') + xlab('Trial') + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
print(gg1)
dev.off()

rm(dq3)
dq3 <- df %>% select(rt_swing,rt_csv,rewFunc,run_trial,v_entropy_wi,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split,rt_change)
dq3 <- dq3 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 15 ~ 'Early',
  run_trial > 15 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
dq3 <- dq3 %>% mutate(entropy_bin = ntile(v_entropy_wi,3))
dq3 <- dq3 %>% group_by(rewFunc,entropy_bin,run_trial,last_outcome) %>% summarize(mRT = mean(rt_change,na.rm=TRUE), dRT = sd(rt_change,na.rm=TRUE),N = length(rt_change)) %>% ungroup()
dq3 <- dq3 %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
dq3 <- dq3 %>% filter(!is.na(entropy_bin))

hbin_edge = quantile(df$v_entropy_wi,probs=seq(0,1,length.out=3),na.rm=TRUE)
Hbin_edge <- as_tibble(hbin_edge) %>% rename(hbin_edge=value)
Hbin_edge <- Hbin_edge %>% mutate(entropy_bin = seq(1,3,length.out=3))

dq3 <- inner_join(dq3,Hbin_edge,by='entropy_bin')
#dq3 <- dq3 %>% filter(entropy_bin > 2 & entropy_bin < 48)

dq3 <- dq3 %>% filter(!is.na(entropy_bin) & !is.na(last_outcome))

dq3 <- dq3 %>% mutate(entropy_label = case_when(
  entropy_bin == 1 ~ 'low',
  entropy_bin == 2 ~ 'mid',
  entropy_bin == 3 ~ 'high'
))

dq3$entropy_label <- factor(dq3$entropy_label,levels=c('low','mid','high'))
dq3 <- dq3 %>% filter(rewFunc=='IEV' | rewFunc=='DEV')

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_swing_by_RewFunc_byH.pdf',height=9,width=12)
gg1 <- ggplot(dq3, aes(x=run_trial,y=mRT,color=Contingency,group=Contingency)) + 
  #geom_line(size=0.5,linetype=2) + 
  scale_color_manual(values=pal) + 
  geom_smooth(method='loess',span=1) + ylab('Change in RT (s)') + xlab('Trial') + facet_wrap(last_outcome~entropy_label) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
print(gg1)
dev.off()


source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_meg')
rm(dq3)
dq3 <- df %>% select(rt_swing,rt_csv,rewFunc,run_trial,v_entropy_wi,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split,rt_change)
dq3 <- dq3 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 15 ~ 'Early',
  run_trial > 15 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
dq3 <- dq3 %>% mutate(entropy_bin = ntile(v_entropy_wi,3))
dq3 <- dq3 %>% group_by(rewFunc,entropy_bin,run_trial,last_outcome) %>% summarize(mRT = mean(rt_change,na.rm=TRUE), dRT = sd(rt_change,na.rm=TRUE),N = length(rt_change)) %>% ungroup()
dq3 <- dq3 %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
dq3 <- dq3 %>% filter(!is.na(entropy_bin))

hbin_edge = quantile(df$v_entropy_wi,probs=seq(0,1,length.out=3),na.rm=TRUE)
Hbin_edge <- as_tibble(hbin_edge) %>% rename(hbin_edge=value)
Hbin_edge <- Hbin_edge %>% mutate(entropy_bin = seq(1,3,length.out=3))

dq3 <- inner_join(dq3,Hbin_edge,by='entropy_bin')
#dq3 <- dq3 %>% filter(entropy_bin > 2 & entropy_bin < 48)

dq3 <- dq3 %>% filter(!is.na(entropy_bin) & !is.na(last_outcome))

dq3 <- dq3 %>% mutate(entropy_label = case_when(
  entropy_bin == 1 ~ 'low',
  entropy_bin == 2 ~ 'mid',
  entropy_bin == 3 ~ 'high'
))

dq3$entropy_label <- factor(dq3$entropy_label,levels=c('low','mid','high'))
dq3 <- dq3 %>% filter(rewFunc=='IEV' | rewFunc=='DEV')

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_swing_by_RewFunc_byH_rep.pdf',height=9,width=12)
gg1 <- ggplot(dq3, aes(x=run_trial,y=mRT,color=Contingency,group=Contingency)) + 
  #geom_line(size=0.5,linetype=2) + 
  scale_color_manual(values=pal) + 
  geom_smooth(method='loess',span=1) + ylab('Change in RT (s)') + xlab('Trial') + facet_wrap(last_outcome~entropy_label) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
print(gg1)
dev.off()
