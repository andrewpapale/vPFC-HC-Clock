# 2022-06-15 AndyP
# Plot RT by rewFunc

library(tidyverse)
library(ggbeeswarm)

repo_directory <- '~/clock_analysis/'
pal = palette()
pal[1] <- '#009B72'
pal[2] <- '#D65B26'
# pal[3] <- '#6F6AB0'
# pal[4] <- '#E50886'
pal[3] <- '#AA399B'

source('~/vPFC-HC-Clock/get_trial_data_vmPFC.R')
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


dq1 <- df %>% select(rt_swing,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split)
dq1 <- dq1 %>% mutate(entropy_bin = ntile(v_entropy,50))
dq1 <- dq1 %>% group_by(rewFunc,entropy_bin) %>% summarize(mRT = mean(rt_swing,na.rm=TRUE), dRT = sd(rt_swing,na.rm=TRUE),N = length(rt_swing)) %>% ungroup()
dq1 <- dq1 %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
dq1 <- dq1 %>% filter(!is.na(entropy_bin))

hbin_edge = quantile(df$v_entropy,probs=seq(0,1,length.out=50),na.rm=TRUE)
Hbin_edge <- as_tibble(hbin_edge) %>% rename(hbin_edge=value)
Hbin_edge <- Hbin_edge %>% mutate(entropy_bin = seq(1,50,length.out=50))

dq1 <- inner_join(dq1,Hbin_edge,by='entropy_bin')

#dq1 <- dq1 %>% filter(entropy_bin > 2 & entropy_bin < 48)

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

rm(dq4)
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
df1 <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_meg')
rm(dq3)
dq3 <- df1 %>% select(rt_swing,rt_csv,rewFunc,run_trial,v_entropy_wi,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split,rt_change)
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


source('~/vPFC-HC-Clock/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_fmri')
rm(dq3)
dq3 <- df %>% dplyr::select(trial_neg_inv_sc,v_max,v_max_wi,v_max_above_median,rt_swing,rt_csv,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax,id,run,rt_change,rt_vmax)
dq3 <- dq3 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 15 ~ 'Early',
  run_trial > 15 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
dq3$trial_bin <- factor(dq3$trial_bin,levels=c('Early','Middle','Late'),ordered=TRUE)
dq3 <- dq3 %>% group_by(id) %>% mutate(entropy_split = case_when(
   v_entropy < mean(v_entropy,na.rm=TRUE) ~ 'low',
   v_entropy > mean(v_entropy,na.rm=TRUE) ~ 'high'
 )) %>% ungroup()
#dq3 <- dq3 %>% group_by(rewFunc,entropy_bin,trial_bin,last_outcome) %>% summarize(mRT = mean(rt_change,na.rm=TRUE), dRT = sd(rt_change,na.rm=TRUE),N = length(rt_change)) %>% ungroup()
dq3 <- dq3 %>% mutate(Contingency1 = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
#dq3 <- dq3 %>% filter(!is.na(entropy_bin))

#dq3$entropy_label <- factor(dq3$entropy_label, levels=c('low','mid','high'),ordered=TRUE)

#dq3 <- dq3 %>% filter(entropy_label!='mid')
dq3 <- dq3 %>% filter(trial_bin!='Middle')
dq3 <- dq3 %>% mutate(Contingency = case_when(
  Contingency1=='IEV' ~ 'LEA',
  Contingency1=='DEV' ~ 'LEA',
  Contingency1=='CEV' ~ 'UNL',
  Contingency1=='CEVR' ~ 'UNL'
))

dq3$Contingency <- factor(dq3$Contingency,levels=c("LEA","UNL"),ordered=TRUE)
dq3 <- dq3 %>% filter(!is.na(entropy_split))

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_swing_by_RewFunc_byH.pdf',height=9,width=12)
dodge <- position_dodge(width = 0.7)
gg1 <- ggplot(dq3,aes(x=entropy_split,y=rt_swing)) + 
  scale_color_manual(values=pal) + 
  #geom_violin(position=dodge) + 
  coord_cartesian(ylim=c(0,2.5)) +
  geom_boxplot(position=dodge,notch=TRUE,outlier.shape='.',width=0.33) + 
  ylab('Change in RT (s)') + xlab('Entropy') + facet_grid(~trial_bin) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), 
        legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=30))
print(gg1)
dev.off()
 
source('~/vPFC-HC-Clock/get_trial_data_vmPFC.R')
df1 <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_meg')
rm(dq4)
dq4 <- df1 %>% select(v_max_above_median,rt_swing,rt_csv,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax,id,run,rt_change,rt_vmax)
dq4 <- dq4 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 15 ~ 'Early',
  run_trial > 15 & run_trial < 42 ~ 'Middle',
  run_trial >=42 ~ 'Late',
)))
dq4$trial_bin <- factor(dq4$trial_bin,levels=c('Early','Middle','Late'),ordered=TRUE)
dq4 <- dq4 %>% group_by(id) %>% mutate(entropy_split = case_when(
  v_entropy < mean(v_entropy,na.rm=TRUE) ~ 'low',
  v_entropy > mean(v_entropy,na.rm=TRUE) ~ 'high'
)) %>% ungroup()
# dq4 <- dq4 %>% mutate(entropy_label = case_when(
#   v_entropy <= mean(v_entropy,na.rm=TRUE)  - 0.5 * sd(v_entropy,na.rm=TRUE) ~ 'low',
#   v_entropy >= mean(v_entropy,na.rm=TRUE) + 0.5 * sd(v_entropy,na.rm=TRUE) ~ 'high',
#   v_entropy > mean(v_entropy,na.rm=TRUE) - 0.5 * sd(v_entropy,na.rm=TRUE) & v_entropy < mean(v_entropy,na.rm=TRUE) + 0.5 * sd(v_entropy,na.rm=TRUE) ~ 'mid'
# ))
#dq4 <- dq4 %>% group_by(rewFunc,entropy_bin,trial_bin,last_outcome) %>% summarize(mRT = mean(rt_change,na.rm=TRUE), dRT = sd(rt_change,na.rm=TRUE),N = length(rt_change)) %>% ungroup()
dq4 <- dq4 %>% mutate(Contingency1 = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
#dq4 <- dq4 %>% filter(!is.na(entropy_bin))

#dq4$entropy_label <- factor(dq4$entropy_label, levels=c('low','mid','high'),ordered=TRUE)

#dq4 <- dq4 %>% filter(entropy_label!='mid')
dq4 <- dq4 %>% filter(trial_bin!='Middle')
dq4 <- dq4 %>% mutate(Contingency = case_when(
  Contingency1=='IEV' ~ 'LEA',
  Contingency1=='DEV' ~ 'LEA',
  Contingency1=='CEV' ~ 'UNL',
  Contingency1=='CEVR' ~ 'UNL'
))

dq4$Contingency <- factor(dq4$Contingency,levels=c("LEA","UNL"),ordered=TRUE)
dq4 <- dq4 %>% filter(!is.na(entropy_split))

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_swing_by_RewFunc_byH_rep.pdf',height=9,width=12)
dodge <- position_dodge(width = 0.7)
gg1 <- ggplot(dq4,aes(x=entropy_split,y=rt_swing)) + 
  scale_color_manual(values=pal) + 
  #geom_violin(position=dodge) + 
  coord_cartesian(ylim=c(0,2.5)) +
  geom_boxplot(position=dodge,notch=TRUE,outlier.shape='.',width=0.33) + 
  ylab('Change in RT (s)') + xlab('Entropy') + facet_grid(~trial_bin) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), 
        legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=30))
print(gg1)
dev.off()


dq4 <- dq4 %>% select(id,run,run_trial,entropy_split,rt_csv,rt_vmax,trial_bin,last_outcome)
dq3 <- dq3 %>% select(id,run,run_trial,entropy_split,rt_csv,rt_vmax,trial_bin,last_outcome)
dq4$id <- as.double(dq4$id)
dq3 <- dq3 %>% mutate(Experiment='fMRI')
dq4 <- dq4 %>% mutate(Experiment='MEG')
dq5 <- full_join(dq3,dq4)
outliers <- dq5 %>% group_by(Experiment,entropy_split,trial_bin,last_outcome) %>% filter(rt_swing > quantile(rt_swing,probs=0.75,na.rm=TRUE) + 1.5 * IQR(rt_swing,na.rm=TRUE) |
                                                                              rt_swing < quantile(rt_swing,probs=0.25,na.rm=TRUE) - 1.5 * IQR(rt_swing,na.rm=TRUE))
dq5 <- dq5 %>% filter(!is.na(last_outcome))
  
setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_swing_fMRI_MEG.pdf',height=12,width=12)
gg1 <- ggplot(data=dq5,aes(x=trial_bin,y=rt_swing,color=entropy_split)) +
  coord_cartesian(ylim=c(0,1.5)) +
  geom_boxplot(notch=TRUE,width=0.5,outlier.shape=NA) +
  geom_beeswarm(data=outliers,size=0.1,dodge.width=0.5) + 
  ylab('Change in RT (s)') + xlab('Trial') + facet_grid(last_outcome~Experiment) + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), 
        legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=30))
print(gg1)
dev.off() 



dq4 <- dq4 %>% select(id,run,run_trial,v_max_above_median,rt_csv,trial_bin,last_outcome,rt_vmax,rt_swing)
dq3 <- dq3 %>% select(id,run,run_trial,v_max_above_median,rt_csv,trial_bin,last_outcome,rt_vmax,rt_swing)
dq4$id <- as.double(dq4$id)
dq3 <- dq3 %>% mutate(Experiment='fMRI')
dq4 <- dq4 %>% mutate(Experiment='MEG')
dq5 <- full_join(dq3,dq4)
dq5 <- dq5 %>% mutate(convergence = abs(rt_csv-rt_vmax))
#outliers <- dq5 %>% group_by(Experiment,v_max_above_median,trial_bin,last_outcome) %>% filter(convergence > quantile(convergence,probs=0.75,na.rm=TRUE) + 1.5 * IQR(convergence,na.rm=TRUE) |
 #                                                                                               convergence < quantile(convergence,probs=0.25,na.rm=TRUE) - 1.5 * IQR(convergence,na.rm=TRUE))
outliers <- dq5 %>% group_by(Experiment,v_max_above_median,trial_bin,last_outcome) %>% filter(rt_swing > quantile(rt_swing,probs=0.75,na.rm=TRUE) + 1.5 * IQR(rt_swing,na.rm=TRUE) |
                                                                                           rt_swing < quantile(rt_swing,probs=0.25,na.rm=TRUE) - 1.5 * IQR(rt_swing,na.rm=TRUE))
dq5 <- dq5 %>% filter(!is.na(last_outcome))
dq5 <- dq5 %>% filter(!is.na(v_max_above_median))
dq5$v_max_above_median <- factor(dq5$v_max_above_median,levels=c('TRUE','FALSE'))
outliers <- outliers %>% filter(!is.na(last_outcome))
outliers <- outliers %>% filter(!is.na(v_max_above_median))

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RTswing_vmax_fMRI_MEG.pdf',height=12,width=12)
gg1 <- ggplot(data=dq5,aes(x=trial_bin,y=rt_swing,color=v_max_above_median)) +
  coord_cartesian(ylim=c(0,1.5)) +
  geom_boxplot(notch=TRUE,width=0.5,outlier.shape=NA) +
  geom_beeswarm(data=outliers,size=0.1,dodge.width=0.5) + 
  ylab('Change in RT (s)') + xlab('Trial') + facet_grid(last_outcome~Experiment) + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), 
        legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=30))
print(gg1)
dev.off() 
