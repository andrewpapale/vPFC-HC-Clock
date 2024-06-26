# 2022-06-15 AndyP
# Plot RT by rewFunc

library(tidyverse)
library(ggbeeswarm)
library(lmerTest)
library(emmeans)
library(broom)
library(ggplot2)

repo_directory <- '~/clock_analysis/'
pal = palette()
pal[1] <- '#009B72'
pal[2] <- '#D65B26'
pal[3] <- '#6F6AB0'
pal[4] <- '#E50886'
#pal[3] <- '#AA399B'

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df1 <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri') %>% mutate(dataset = 'Experiment 1')
df2 <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_meg') %>% mutate(dataset = 'Experiment 1 - Replication')
df3 <- get_trial_data(repo_directory=repo_directory,dataset='explore') %>% mutate(dataset = 'Experiment 2')

df3 <- df3 %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                      trial > 40 & trial <= 80 ~ 2,
                                      trial > 80 & trial <=120 ~ 3, 
                                      trial > 120 & trial <=160 ~ 4,
                                      trial > 160 & trial <=200 ~ 5,
                                      trial > 200 & trial <=240 ~ 6))
df3 <- df3 %>% mutate(run_trial = case_when(trial <= 40 ~ trial, 
                                          trial > 40 & trial <= 80 ~ trial-40,
                                          trial > 80 & trial <=120 ~ trial-80, 
                                          trial > 120 & trial <=160 ~ trial-120,
                                          trial > 160 & trial <=200 ~ trial-160,
                                          trial > 200 & trial <=240 ~ trial-200))
df3 <- df3 %>% mutate(run_trial_c = run_trial-floor(run_trial/40.5),
                    run_trial_neg_inv = -(1 / run_trial_c) * 100,
                    run_trial_neg_inv_sc = as.vector(scale(run_trial_neg_inv)))

common_cols = intersect(colnames(df1),colnames(df2))
df <- rbind(subset(df1,select = common_cols),subset(df2,select = common_cols))

common_cols = intersect(colnames(df),colnames(df3))
df <- rbind(subset(df,select = common_cols),subset(df3,select = common_cols))

df$rewFunc <- relevel(as.factor(df$rewFunc),order='DEV','IEV','CEV','CEVR')

pdf('RT_by_rewFunc.pdf',height=8,width=26)
ggplot(df,aes(x=run_trial, y=rt_csv, color=rewFunc, group=rewFunc)) +
  geom_smooth(span=1,method='loess',linewidth=3) +
  stat_summary(aes(y=rt_csv,group=rewFunc),fun.y = mean, geom="line",linetype='dashed',size=3) +
  scale_color_manual(values = pal) +
  theme_bw(base_size=13) +  
  ylab('Response Time (s)') +
  xlab('Trial') +
  facet_wrap(~dataset, scales = 'free_x') +
  guides(colour = guide_legend(override.aes = list(size=10,width=5,fill=NA))) +
  theme(legend.title = element_blank(),
        legend.background = element_rect(color='white'),
        axis.title.y = element_text(margin=margin(r=6),size=26),
        axis.title.x = element_text(margin=margin(t=6),size=26),
        legend.text = element_text(size=26),
        axis.text.x = element_text(size=26),
        axis.text.y = element_text(size=26),
        strip.text = element_text(size = 26)
  )
dev.off()



df <- df %>% group_by(id,run) %>% mutate(v_max_lag = lag(v_max)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(v_chosen_lag = lag(v_chosen)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(vD = v_max_lag - v_chosen_lag) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(ev_lag = lag(ev)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(score_lag = lag(score_csv)) %>% ungroup()
df$v_max_lag <- scale(df$v_max_lag)
df$v_chosen_lag <- scale(df$v_chosen_lag)
df$v_entropy <- scale(df$v_entropy)
df$vD <- scale(df$vD)
df$v_max_full <- scale(df$v_max_full)
df$ev_lag <- scale(df$ev_lag)
df$score_lag <- scale(df$score_lag)
# Plot emmeans of mlm model
m1 <- lmer(rt_swing ~ v_entropy*trial_neg_inv_sc*last_outcome + v_max_lag*trial_neg_inv_sc*last_outcome + (1|id/rewFunc),df)

qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)

emH_fmri <- emmeans(m1, outcome = 'rt_swing',specs = c(formula(~ v_entropy*trial_neg_inv_sc*last_outcome)),
                      at = list(trial_neg_inv_sc=qT,v_entropy=c(-1.5,1.5)))
emV_fmri <- emmeans(m1, outcome = 'rt_swing',specs = c(formula(~ v_max_lag*trial_neg_inv_sc*last_outcome)),
               at = list(trial_neg_inv_sc=qT,v_max_lag=c(-1.5,1.5)))

df <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_meg')
df <- df %>% group_by(id,run) %>% mutate(v_max_lag = lag(v_max)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(v_chosen_lag = lag(v_chosen)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(vD = v_max_lag - v_chosen_lag) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(ev_lag = lag(ev)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(score_lag = lag(score_csv)) %>% ungroup()
df$v_max_lag <- scale(df$v_max_lag)
df$v_chosen_lag <- scale(df$v_chosen_lag)
df$v_entropy <- scale(df$v_entropy)
df$vD <- scale(df$vD)
df$v_max_full <- scale(df$v_max_full)
df$ev_lag <- scale(df$ev_lag)
df$score_lag <- scale(df$score_lag)
# Plot emmeans of mlm model
m2 <- lmer(rt_swing ~ v_entropy*trial_neg_inv_sc*last_outcome + v_max_lag*trial_neg_inv_sc*last_outcome + (1|id/rewFunc),df)

qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)

emH_meg <- emmeans(m2, outcome = 'rt_swing',specs = c(formula(~ v_entropy*trial_neg_inv_sc*last_outcome)),
                    at = list(trial_neg_inv_sc=qT,v_entropy=c(-1.5,1.5)))
emV_meg <- emmeans(m2, outcome = 'rt_swing',specs = c(formula(~ v_max_lag*trial_neg_inv_sc*last_outcome)),
                    at = list(trial_neg_inv_sc=qT,v_max_lag=c(-1.5,1.5)))


emH_fmri <- broom::tidy(emH_fmri$`emmeans of v_entropy, trial_neg_inv_sc, last_outcome`)
emH_fmri <- emH_fmri  %>% mutate(Sample = 'fMRI')
emH_meg <- broom::tidy(emH_meg$`emmeans of v_entropy, trial_neg_inv_sc, last_outcome`)
emH_meg <- emH_meg  %>% mutate(Sample = 'MEG')

emV_fmri <- broom::tidy(emV_fmri$`emmeans of v_max_lag, trial_neg_inv_sc, last_outcome`)
emV_fmri <- emV_fmri %>% mutate(Sample = 'fMRI')
emV_meg <- broom::tidy(emV_meg$`emmeans of v_max_lag, trial_neg_inv_sc, last_outcome`)
emV_meg <- emV_meg %>% mutate(Sample = 'MEG')

emH <- rbind(emH_fmri,emH_meg)
emV <- rbind(emV_fmri,emV_meg)

emH <- emH %>% mutate(Entropy = case_when(v_entropy <= -1 ~ 'Low',v_entropy >= 1 ~ 'High'))
emH <- emH %>% mutate(Trial = case_when(trial_neg_inv_sc < 0 ~ 'Early', trial_neg_inv_sc > 0 ~ 'Late'))
emH <- emH %>% mutate(Omission_by_Trial = case_when(Trial=='Early' & last_outcome=='Omission' ~ 'Early Omission',
                      Trial=='Early' & last_outcome=='Reward' ~ 'Early Reward',
                      Trial=='Late' & last_outcome=='Omission' ~ 'Late Omission',
                      Trial=='Late' & last_outcome=='Reward' ~ 'Late Reward'))
emH$Omission_by_Trial <- factor(emH$Omission_by_Trial,levels=c('Early Reward','Late Reward','Early Omission','Late Omission'))

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_swing_by_emmeans_Entropy.pdf',height=12,width=12)
gg1 <- ggplot(emH) + geom_point(aes(x=Omission_by_Trial,y=estimate,color=Entropy),size=8) +
  geom_errorbar(width=0.5,size=1.5,aes(x=Omission_by_Trial,y=estimate,color=Entropy,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_grid(~Sample) + xlab('Trial / Last Outcome') + ylab('Mean RT Swing') +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=40), 
        legend.text=element_text(size=20),legend.title=element_text(size=30),
        legend.spacing.y = unit(1.0,'cm'),
        axis.title.x = element_text(margin = margin(t=20,r=20,b=0,l=0)),
        axis.title.y = element_text(margin = margin(t=0,r=20,b=0,l=0))) +
  theme(axis.text.x = element_text(angle=90,vjust=0.05,hjust=0.5)) +
  theme(strip.text.x = element_text(size=30)) +
  theme(legend.key.size = unit(1.5,"cm")) +
  guides(fill = guide_legend(byrow=TRUE))
print(gg1)
dev.off()

emV <- emV %>% mutate(Value = case_when(v_max_lag <= -1 ~ 'Low',v_max_lag >= 1 ~ 'High'))
emV <- emV %>% mutate(Trial = case_when(trial_neg_inv_sc < 0 ~ 'Early', trial_neg_inv_sc > 0 ~ 'Late'))
emV <- emV %>% mutate(Omission_by_Trial = case_when(Trial=='Early' & last_outcome=='Omission' ~ 'Early Omission',
                                                    Trial=='Early' & last_outcome=='Reward' ~ 'Early Reward',
                                                    Trial=='Late' & last_outcome=='Omission' ~ 'Late Omission',
                                                    Trial=='Late' & last_outcome=='Reward' ~ 'Late Reward'))
emV$Omission_by_Trial <- factor(emV$Omission_by_Trial,levels=c('Early Reward','Late Reward','Early Omission','Late Omission'))
emV$Value <- factor(emV$Value,levels=c('Low','High'))

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_swing_by_emmeans_Value.pdf',height=12,width=12)
gg1 <- ggplot(emV) + geom_point(aes(x=Omission_by_Trial,y=estimate,color=Value),size=8) +
  geom_errorbar(width=0.5,size=1.5,aes(x=Omission_by_Trial,y=estimate,color=Value,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_grid(~Sample) + xlab('Trial / Last Outcome') + ylab('Mean RT Swing') +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=40), 
        legend.text=element_text(size=20),legend.title=element_text(size=30),
        legend.spacing.y = unit(1.0,'cm'),
        axis.title.x = element_text(margin = margin(t=20,r=20,b=0,l=0)),
        axis.title.y = element_text(margin = margin(t=0,r=20,b=0,l=0))) +
  theme(axis.text.x = element_text(angle=90,vjust=0.05,hjust=0.5)) +
  theme(strip.text.x = element_text(size=30)) +
  theme(legend.key.size = unit(1.5,"cm")) +
  guides(fill = guide_legend(byrow=TRUE)) +
  guides(fill=guide_legend(title='Value'))
print(gg1)
dev.off()


# rt_csv - rt_vmax

source('~/vPFC-HC-Clock/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_fmri')
df <- df %>% group_by(id,run) %>% mutate(v_max_lag = lag(v_max)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(v_chosen_lag = lag(v_chosen)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(vD = v_max_lag - v_chosen_lag) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(ev_lag = lag(ev)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(score_lag = lag(score_csv)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(rt_conv = abs(rt_csv - rt_vmax_lag)) %>% ungroup()
df$v_max_lag <- scale(df$v_max_lag)
df$v_chosen_lag <- scale(df$v_chosen_lag)
df$v_entropy <- scale(df$v_entropy)
df$vD <- scale(df$vD)
df$v_max_full <- scale(df$v_max_full)
df$ev_lag <- scale(df$ev_lag)
df$score_lag <- scale(df$score_lag)
# Plot emmeans of mlm model
m1 <- lmer(rt_conv ~ v_entropy*trial_neg_inv_sc*last_outcome + v_max_lag*trial_neg_inv_sc*last_outcome + (1|id/rewFunc),df)

qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)

emH_fmri <- emmeans(m1, outcome = 'rt_conv',specs = c(formula(~ v_entropy*trial_neg_inv_sc*last_outcome)),
                    at = list(trial_neg_inv_sc=qT,v_entropy=c(-1.5,1.5)))
emV_fmri <- emmeans(m1, outcome = 'rt_conv',specs = c(formula(~ v_max_lag*trial_neg_inv_sc*last_outcome)),
                    at = list(trial_neg_inv_sc=qT,v_max_lag=c(-1.5,1.5)))

df <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_meg')
df <- df %>% group_by(id,run) %>% mutate(v_max_lag = lag(v_max)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(v_chosen_lag = lag(v_chosen)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(vD = v_max_lag - v_chosen_lag) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(ev_lag = lag(ev)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(score_lag = lag(score_csv)) %>% ungroup()
df <- df %>% group_by(id,run) %>% mutate(rt_conv = abs(rt_csv - rt_vmax_lag)) %>% ungroup()
df$v_max_lag <- scale(df$v_max_lag)
df$v_chosen_lag <- scale(df$v_chosen_lag)
df$v_entropy <- scale(df$v_entropy)
df$vD <- scale(df$vD)
df$v_max_full <- scale(df$v_max_full)
df$ev_lag <- scale(df$ev_lag)
df$score_lag <- scale(df$score_lag)
# Plot emmeans of mlm model
m2 <- lmer(rt_conv ~ v_entropy*trial_neg_inv_sc*last_outcome + v_max_lag*trial_neg_inv_sc*last_outcome + (1|id/rewFunc),df)

qT <- quantile(df$trial_neg_inv_sc,c(0.1,0.9),na.rm=TRUE)

emH_meg <- emmeans(m2, outcome = 'rt_conv',specs = c(formula(~ v_entropy*trial_neg_inv_sc*last_outcome)),
                   at = list(trial_neg_inv_sc=qT,v_entropy=c(-1.5,1.5)))
emV_meg <- emmeans(m2, outcome = 'rt_conv',specs = c(formula(~ v_max_lag*trial_neg_inv_sc*last_outcome)),
                   at = list(trial_neg_inv_sc=qT,v_max_lag=c(-1.5,1.5)))


emH_fmri <- broom::tidy(emH_fmri$`emmeans of v_entropy, trial_neg_inv_sc, last_outcome`)
emH_fmri <- emH_fmri  %>% mutate(Sample = 'fMRI')
emH_meg <- broom::tidy(emH_meg$`emmeans of v_entropy, trial_neg_inv_sc, last_outcome`)
emH_meg <- emH_meg  %>% mutate(Sample = 'MEG')

emV_fmri <- broom::tidy(emV_fmri$`emmeans of v_max_lag, trial_neg_inv_sc, last_outcome`)
emV_fmri <- emV_fmri %>% mutate(Sample = 'fMRI')
emV_meg <- broom::tidy(emV_meg$`emmeans of v_max_lag, trial_neg_inv_sc, last_outcome`)
emV_meg <- emV_meg %>% mutate(Sample = 'MEG')

emH <- rbind(emH_fmri,emH_meg)
emV <- rbind(emV_fmri,emV_meg)

emH <- emH %>% mutate(Entropy = case_when(v_entropy <= -1 ~ 'Low',v_entropy >= 1 ~ 'High'))
emH <- emH %>% mutate(Trial = case_when(trial_neg_inv_sc < 0 ~ 'Early', trial_neg_inv_sc > 0 ~ 'Late'))
emH <- emH %>% mutate(Omission_by_Trial = case_when(Trial=='Early' & last_outcome=='Omission' ~ 'Early Omission',
                                                    Trial=='Early' & last_outcome=='Reward' ~ 'Early Reward',
                                                    Trial=='Late' & last_outcome=='Omission' ~ 'Late Omission',
                                                    Trial=='Late' & last_outcome=='Reward' ~ 'Late Reward'))
emH$Omission_by_Trial <- factor(emH$Omission_by_Trial,levels=c('Early Reward','Late Reward','Early Omission','Late Omission'))
emH$Entropy <- factor(emH$Entropy, levels = c('High','Low'))

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_conv_by_emmeans_Entropy.pdf',height=12,width=12)
gg1 <- ggplot(emH) + geom_point(aes(x=Omission_by_Trial,y=estimate,color=Entropy),size=8) +
  geom_errorbar(width=0.5,size=1.5,aes(x=Omission_by_Trial,y=estimate,color=Entropy,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_grid(~Sample) + xlab('Trial / Last Outcome') + ylab('Mean RT Convergence') +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=40), 
        legend.text=element_text(size=20),legend.title=element_text(size=30),
        legend.spacing.y = unit(1.0,'cm'),
        axis.title.x = element_text(margin = margin(t=20,r=20,b=0,l=0)),
        axis.title.y = element_text(margin = margin(t=0,r=20,b=0,l=0))) +
  theme(axis.text.x = element_text(angle=90,vjust=0.05,hjust=0.5)) +
  theme(strip.text.x = element_text(size=30)) +
  theme(legend.key.size = unit(1.5,"cm")) +
  guides(fill = guide_legend(byrow=TRUE))
print(gg1)
dev.off()

emV <- emV %>% mutate(Value = case_when(v_max_lag <= -1 ~ 'Low',v_max_lag >= 1 ~ 'High'))
emV <- emV %>% mutate(Trial = case_when(trial_neg_inv_sc < 0 ~ 'Early', trial_neg_inv_sc > 0 ~ 'Late'))
emV <- emV %>% mutate(Omission_by_Trial = case_when(Trial=='Early' & last_outcome=='Omission' ~ 'Early Omission',
                                                    Trial=='Early' & last_outcome=='Reward' ~ 'Early Reward',
                                                    Trial=='Late' & last_outcome=='Omission' ~ 'Late Omission',
                                                    Trial=='Late' & last_outcome=='Reward' ~ 'Late Reward'))
emV$Omission_by_Trial <- factor(emV$Omission_by_Trial,levels=c('Early Reward','Late Reward','Early Omission','Late Omission'))
emV$Value <- factor(emV$Value,levels=c('Low','High'))

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('RT_conv_by_emmeans_Value.pdf',height=12,width=12)
gg1 <- ggplot(emV) + geom_point(aes(x=Omission_by_Trial,y=estimate,color=Value),size=8) +
  geom_errorbar(width=0.5,size=1.5,aes(x=Omission_by_Trial,y=estimate,color=Value,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  facet_grid(~Sample) + xlab('Trial / Last Outcome') + ylab('Mean RT Convergence') +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=40), 
        legend.text=element_text(size=20),legend.title=element_text(size=30),
        legend.spacing.y = unit(1.0,'cm'),
        axis.title.x = element_text(margin = margin(t=20,r=20,b=0,l=0)),
        axis.title.y = element_text(margin = margin(t=0,r=20,b=0,l=0))) +
  theme(axis.text.x = element_text(angle=90,vjust=0.05,hjust=0.5)) +
  theme(strip.text.x = element_text(size=30)) +
  theme(legend.key.size = unit(1.5,"cm")) +
  guides(fill = guide_legend(byrow=TRUE)) +
  guides(fill=guide_legend(title='Value'))
print(gg1)
dev.off()



# dq <- df %>% select(rt_csv,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split)
# dq <- dq %>% group_by(run_trial,rewFunc) %>% summarize(mRT = mean(rt_csv,na.rm=TRUE), dRT = sd(rt_csv,na.rm=TRUE),N = length(rt_csv)) %>% ungroup()
# dq <- dq %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
# 
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
# pdf('RTbyRewFunc.pdf',height=9,width=12)
# gg1 <- ggplot(dq, aes(x=run_trial,y=mRT,color=Contingency,group=Contingency)) + 
#   geom_line(size=0.5,linetype=2) + scale_color_manual(values=pal) + 
#   geom_smooth(method='loess',span=1) + ylab('Response Time (s)') + xlab('Trial') + 
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
# print(gg1)
# dev.off()
# 
# 
# dq1 <- df %>% select(rt_swing,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split)
# dq1 <- dq1 %>% mutate(entropy_bin = ntile(v_entropy,50))
# dq1 <- dq1 %>% group_by(rewFunc,entropy_bin) %>% summarize(mRT = mean(rt_swing,na.rm=TRUE), dRT = sd(rt_swing,na.rm=TRUE),N = length(rt_swing)) %>% ungroup()
# dq1 <- dq1 %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
# dq1 <- dq1 %>% filter(!is.na(entropy_bin))
# 
# hbin_edge = quantile(df$v_entropy,probs=seq(0,1,length.out=50),na.rm=TRUE)
# Hbin_edge <- as_tibble(hbin_edge) %>% rename(hbin_edge=value)
# Hbin_edge <- Hbin_edge %>% mutate(entropy_bin = seq(1,50,length.out=50))
# 
# dq1 <- inner_join(dq1,Hbin_edge,by='entropy_bin')
# 
# #dq1 <- dq1 %>% filter(entropy_bin > 2 & entropy_bin < 48)
# 
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
# pdf('RTbyRewFunc_H.pdf',height=9,width=12)
# gg1 <- ggplot(dq1, aes(x=hbin_edge,y=mRT,color=Contingency,group=Contingency)) + 
#   geom_line(size=0.5,linetype=2) + scale_color_manual(values=pal) + 
#   geom_smooth(method='loess',span=1) + ylab('Change in RT (s)') + xlab('Entropy') + 
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
# print(gg1)
# dev.off()
# 
# 
# 
# dq1 <- df %>% select(rt_swing,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split)
# dq1 <- dq1 %>% mutate(entropy_bin = ntile(v_entropy,50))
# dq1 <- dq1 %>% group_by(rewFunc,entropy_bin,run_trial) %>% summarize(mRT = mean(rt_swing,na.rm=TRUE), dRT = sd(rt_swing,na.rm=TRUE),N = length(rt_swing)) %>% ungroup()
# dq1 <- dq1 %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
# dq1 <- dq1 %>% filter(!is.na(entropy_bin))
# 
# hbin_edge = quantile(df$v_entropy,probs=seq(0,1,length.out=50),na.rm=TRUE)
# Hbin_edge <- as_tibble(hbin_edge) %>% rename(hbin_edge=value)
# Hbin_edge <- Hbin_edge %>% mutate(entropy_bin = seq(1,50,length.out=50))
# 
# 
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
# pdf('RT_swing_by_RewFunc.pdf',height=9,width=12)
# gg1 <- ggplot(dq2, aes(x=run_trial,y=mRT,color=Contingency,group=Contingency)) + 
#   geom_line(size=0.5,linetype=2) + scale_color_manual(values=pal) + 
#   geom_smooth(method='loess',span=1) + ylab('Change in RT (s)') + xlab('Trial') + 
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
# print(gg1)
# dev.off()
# 
# rm(dq4)
# dq3 <- df %>% select(rt_swing,rt_csv,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split,rt_change)
# dq3 <- dq3 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 15 ~ 'Early',
#   run_trial > 15 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# dq3 <- dq3 %>% mutate(entropy_bin = ntile(v_entropy,3))
# dq3 <- dq3 %>% group_by(rewFunc,entropy_bin,run_trial,last_outcome) %>% summarize(mRT = mean(rt_change,na.rm=TRUE), dRT = sd(rt_change,na.rm=TRUE),N = length(rt_change)) %>% ungroup()
# dq3 <- dq3 %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
# dq3 <- dq3 %>% filter(!is.na(entropy_bin))
# 
# hbin_edge = quantile(df$v_entropy,probs=seq(0,1,length.out=3),na.rm=TRUE)
# Hbin_edge <- as_tibble(hbin_edge) %>% rename(hbin_edge=value)
# Hbin_edge <- Hbin_edge %>% mutate(entropy_bin = seq(1,3,length.out=3))
# 
# dq3 <- inner_join(dq3,Hbin_edge,by='entropy_bin')
# #dq3 <- dq3 %>% filter(entropy_bin > 2 & entropy_bin < 48)
# 
# dq3 <- dq3 %>% filter(!is.na(entropy_bin) & !is.na(last_outcome))
# 
# dq3 <- dq3 %>% mutate(entropy_label = case_when(
#   entropy_bin == 1 ~ 'low',
#   entropy_bin == 2 ~ 'mid',
#   entropy_bin == 3 ~ 'high'
# ))
# 
# dq3$entropy_label <- factor(dq3$entropy_label,levels=c('low','mid','high'))
# dq3 <- dq3 %>% filter(rewFunc=='IEV' | rewFunc=='DEV')
# 
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
# pdf('RT_swing_by_RewFunc_byH.pdf',height=9,width=12)
# gg1 <- ggplot(dq3, aes(x=run_trial,y=mRT,color=Contingency,group=Contingency)) + 
#   #geom_line(size=0.5,linetype=2) + 
#   scale_color_manual(values=pal) + 
#   geom_smooth(method='loess',span=1) + ylab('Change in RT (s)') + xlab('Trial') + facet_wrap(last_outcome~entropy_label) +
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
# print(gg1)
# dev.off()
# 
# 
# source('~/vmPFC/get_trial_data_vmPFC.R')
# df1 <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_meg')
# rm(dq3)
# dq3 <- df1 %>% select(rt_swing,rt_csv,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax_lag,id,run,entropy_split,rt_change)
# dq3 <- dq3 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 15 ~ 'Early',
#   run_trial > 15 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# dq3 <- dq3 %>% mutate(entropy_bin = ntile(v_entropy,3))
# dq3 <- dq3 %>% group_by(rewFunc,entropy_bin,run_trial,last_outcome) %>% summarize(mRT = mean(rt_change,na.rm=TRUE), dRT = sd(rt_change,na.rm=TRUE),N = length(rt_change)) %>% ungroup()
# dq3 <- dq3 %>% mutate(Contingency = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
# dq3 <- dq3 %>% filter(!is.na(entropy_bin))
# 
# hbin_edge = quantile(df$v_entropy,probs=seq(0,1,length.out=3),na.rm=TRUE)
# Hbin_edge <- as_tibble(hbin_edge) %>% rename(hbin_edge=value)
# Hbin_edge <- Hbin_edge %>% mutate(entropy_bin = seq(1,3,length.out=3))
# 
# dq3 <- inner_join(dq3,Hbin_edge,by='entropy_bin')
# #dq3 <- dq3 %>% filter(entropy_bin > 2 & entropy_bin < 48)
# 
# dq3 <- dq3 %>% filter(!is.na(entropy_bin) & !is.na(last_outcome))
# 
# dq3 <- dq3 %>% mutate(entropy_label = case_when(
#   entropy_bin == 1 ~ 'low',
#   entropy_bin == 2 ~ 'mid',
#   entropy_bin == 3 ~ 'high'
# ))
# 
# dq3$entropy_label <- factor(dq3$entropy_label,levels=c('low','mid','high'))
# dq3 <- dq3 %>% filter(rewFunc=='IEV' | rewFunc=='DEV')
# 
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
# pdf('RT_swing_by_RewFunc_byH_rep.pdf',height=9,width=12)
# gg1 <- ggplot(dq3, aes(x=run_trial,y=mRT,color=Contingency,group=Contingency)) + 
#   #geom_line(size=0.5,linetype=2) + 
#   scale_color_manual(values=pal) + 
#   geom_smooth(method='loess',span=1) + ylab('Change in RT (s)') + xlab('Trial') + facet_wrap(last_outcome~entropy_label) +
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=40))
# print(gg1)
# dev.off()
# 
# 
# source('~/vPFC-HC-Clock/get_trial_data_vmPFC.R')
# df <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_fmri')
# rm(dq3)
# dq3 <- df %>% dplyr::select(trial_neg_inv_sc,v_max,v_max_wi,v_max_above_median,rt_swing,rt_csv,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax,id,run,rt_change,rt_vmax)
# dq3 <- dq3 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 15 ~ 'Early',
#   run_trial > 15 & run_trial < 30 ~ 'Middle',
#   run_trial >=30 ~ 'Late',
# )))
# dq3$trial_bin <- factor(dq3$trial_bin,levels=c('Early','Middle','Late'),ordered=TRUE)
# dq3 <- dq3 %>% group_by(id) %>% mutate(entropy_split = case_when(
#    v_entropy < mean(v_entropy,na.rm=TRUE) ~ 'low',
#    v_entropy > mean(v_entropy,na.rm=TRUE) ~ 'high'
#  )) %>% ungroup()
# #dq3 <- dq3 %>% group_by(rewFunc,entropy_bin,trial_bin,last_outcome) %>% summarize(mRT = mean(rt_change,na.rm=TRUE), dRT = sd(rt_change,na.rm=TRUE),N = length(rt_change)) %>% ungroup()
# dq3 <- dq3 %>% mutate(Contingency1 = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
# #dq3 <- dq3 %>% filter(!is.na(entropy_bin))
# 
# #dq3$entropy_label <- factor(dq3$entropy_label, levels=c('low','mid','high'),ordered=TRUE)
# 
# #dq3 <- dq3 %>% filter(entropy_label!='mid')
# dq3 <- dq3 %>% filter(trial_bin!='Middle')
# dq3 <- dq3 %>% mutate(Contingency = case_when(
#   Contingency1=='IEV' ~ 'LEA',
#   Contingency1=='DEV' ~ 'LEA',
#   Contingency1=='CEV' ~ 'UNL',
#   Contingency1=='CEVR' ~ 'UNL'
# ))
# 
# dq3$Contingency <- factor(dq3$Contingency,levels=c("LEA","UNL"),ordered=TRUE)
# dq3 <- dq3 %>% filter(!is.na(entropy_split))
# 
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
# pdf('RT_swing_by_RewFunc_byH.pdf',height=9,width=12)
# dodge <- position_dodge(width = 0.7)
# gg1 <- ggplot(dq3,aes(x=entropy_split,y=rt_swing)) + 
#   scale_color_manual(values=pal) + 
#   #geom_violin(position=dodge) + 
#   coord_cartesian(ylim=c(0,2.5)) +
#   geom_boxplot(position=dodge,notch=TRUE,outlier.shape='.',width=0.33) + 
#   ylab('Change in RT (s)') + xlab('Entropy') + facet_grid(~trial_bin) +
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), 
#         legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=30))
# print(gg1)
# dev.off()
#  
# source('~/vPFC-HC-Clock/get_trial_data_vmPFC.R')
# df1 <- get_trial_data_vmPFC(repo_directory = repo_directory,dataset='mmclock_meg')
# rm(dq4)
# dq4 <- df1 %>% select(v_max_above_median,rt_swing,rt_csv,rewFunc,run_trial,v_entropy,last_outcome,outcome,rt_vmax,id,run,rt_change,rt_vmax)
# dq4 <- dq4 %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
#   run_trial <= 15 ~ 'Early',
#   run_trial > 15 & run_trial < 42 ~ 'Middle',
#   run_trial >=42 ~ 'Late',
# )))
# dq4$trial_bin <- factor(dq4$trial_bin,levels=c('Early','Middle','Late'),ordered=TRUE)
# dq4 <- dq4 %>% group_by(id) %>% mutate(entropy_split = case_when(
#   v_entropy < mean(v_entropy,na.rm=TRUE) ~ 'low',
#   v_entropy > mean(v_entropy,na.rm=TRUE) ~ 'high'
# )) %>% ungroup()
# # dq4 <- dq4 %>% mutate(entropy_label = case_when(
# #   v_entropy <= mean(v_entropy,na.rm=TRUE)  - 0.5 * sd(v_entropy,na.rm=TRUE) ~ 'low',
# #   v_entropy >= mean(v_entropy,na.rm=TRUE) + 0.5 * sd(v_entropy,na.rm=TRUE) ~ 'high',
# #   v_entropy > mean(v_entropy,na.rm=TRUE) - 0.5 * sd(v_entropy,na.rm=TRUE) & v_entropy < mean(v_entropy,na.rm=TRUE) + 0.5 * sd(v_entropy,na.rm=TRUE) ~ 'mid'
# # ))
# #dq4 <- dq4 %>% group_by(rewFunc,entropy_bin,trial_bin,last_outcome) %>% summarize(mRT = mean(rt_change,na.rm=TRUE), dRT = sd(rt_change,na.rm=TRUE),N = length(rt_change)) %>% ungroup()
# dq4 <- dq4 %>% mutate(Contingency1 = factor(rewFunc,levels=c('IEV','DEV','CEV','CEVR'),ordered=TRUE))
# #dq4 <- dq4 %>% filter(!is.na(entropy_bin))
# 
# #dq4$entropy_label <- factor(dq4$entropy_label, levels=c('low','mid','high'),ordered=TRUE)
# 
# #dq4 <- dq4 %>% filter(entropy_label!='mid')
# dq4 <- dq4 %>% filter(trial_bin!='Middle')
# dq4 <- dq4 %>% mutate(Contingency = case_when(
#   Contingency1=='IEV' ~ 'LEA',
#   Contingency1=='DEV' ~ 'LEA',
#   Contingency1=='CEV' ~ 'UNL',
#   Contingency1=='CEVR' ~ 'UNL'
# ))
# 
# dq4$Contingency <- factor(dq4$Contingency,levels=c("LEA","UNL"),ordered=TRUE)
# dq4 <- dq4 %>% filter(!is.na(entropy_split))
# 
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
# pdf('RT_swing_by_RewFunc_byH_rep.pdf',height=9,width=12)
# dodge <- position_dodge(width = 0.7)
# gg1 <- ggplot(dq4,aes(x=entropy_split,y=rt_swing)) + 
#   scale_color_manual(values=pal) + 
#   #geom_violin(position=dodge) + 
#   coord_cartesian(ylim=c(0,2.5)) +
#   geom_boxplot(position=dodge,notch=TRUE,outlier.shape='.',width=0.33) + 
#   ylab('Change in RT (s)') + xlab('Entropy') + facet_grid(~trial_bin) +
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), 
#         legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=30))
# print(gg1)
# dev.off()
# 
# 
# dq4 <- dq4 %>% select(id,run,run_trial,entropy_split,rt_csv,rt_vmax,trial_bin,last_outcome)
# dq3 <- dq3 %>% select(id,run,run_trial,entropy_split,rt_csv,rt_vmax,trial_bin,last_outcome)
# dq4$id <- as.double(dq4$id)
# dq3 <- dq3 %>% mutate(Experiment='fMRI')
# dq4 <- dq4 %>% mutate(Experiment='MEG')
# dq5 <- full_join(dq3,dq4)
# outliers <- dq5 %>% group_by(Experiment,entropy_split,trial_bin,last_outcome) %>% filter(rt_swing > quantile(rt_swing,probs=0.75,na.rm=TRUE) + 1.5 * IQR(rt_swing,na.rm=TRUE) |
#                                                                               rt_swing < quantile(rt_swing,probs=0.25,na.rm=TRUE) - 1.5 * IQR(rt_swing,na.rm=TRUE))
# dq5 <- dq5 %>% filter(!is.na(last_outcome))
#   
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
# pdf('RT_swing_fMRI_MEG.pdf',height=12,width=12)
# gg1 <- ggplot(data=dq5,aes(x=trial_bin,y=rt_swing,color=entropy_split)) +
#   coord_cartesian(ylim=c(0,1.5)) +
#   geom_boxplot(notch=TRUE,width=0.5,outlier.shape=NA) +
#   geom_beeswarm(data=outliers,size=0.1,dodge.width=0.5) + 
#   ylab('Change in RT (s)') + xlab('Trial') + facet_grid(last_outcome~Experiment) + 
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), 
#         legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=30))
# print(gg1)
# dev.off() 
# 
# 
# 
# dq4 <- dq4 %>% select(id,run,run_trial,v_max_above_median,rt_csv,trial_bin,last_outcome,rt_vmax,rt_swing)
# dq3 <- dq3 %>% select(id,run,run_trial,v_max_above_median,rt_csv,trial_bin,last_outcome,rt_vmax,rt_swing)
# dq4$id <- as.double(dq4$id)
# dq3 <- dq3 %>% mutate(Experiment='fMRI')
# dq4 <- dq4 %>% mutate(Experiment='MEG')
# dq5 <- full_join(dq3,dq4)
# dq5 <- dq5 %>% mutate(convergence = abs(rt_csv-rt_vmax))
# #outliers <- dq5 %>% group_by(Experiment,v_max_above_median,trial_bin,last_outcome) %>% filter(convergence > quantile(convergence,probs=0.75,na.rm=TRUE) + 1.5 * IQR(convergence,na.rm=TRUE) |
#  #                                                                                               convergence < quantile(convergence,probs=0.25,na.rm=TRUE) - 1.5 * IQR(convergence,na.rm=TRUE))
# outliers <- dq5 %>% group_by(Experiment,v_max_above_median,trial_bin,last_outcome) %>% filter(rt_swing > quantile(rt_swing,probs=0.75,na.rm=TRUE) + 1.5 * IQR(rt_swing,na.rm=TRUE) |
#                                                                                            rt_swing < quantile(rt_swing,probs=0.25,na.rm=TRUE) - 1.5 * IQR(rt_swing,na.rm=TRUE))
# dq5 <- dq5 %>% filter(!is.na(last_outcome))
# dq5 <- dq5 %>% filter(!is.na(v_max_above_median))
# dq5$v_max_above_median <- factor(dq5$v_max_above_median,levels=c('TRUE','FALSE'))
# outliers <- outliers %>% filter(!is.na(last_outcome))
# outliers <- outliers %>% filter(!is.na(v_max_above_median))
# 
# setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
# pdf('RTswing_vmax_fMRI_MEG.pdf',height=12,width=12)
# gg1 <- ggplot(data=dq5,aes(x=trial_bin,y=rt_swing,color=v_max_above_median)) +
#   coord_cartesian(ylim=c(0,1.5)) +
#   geom_boxplot(notch=TRUE,width=0.5,outlier.shape=NA) +
#   geom_beeswarm(data=outliers,size=0.1,dodge.width=0.5) + 
#   ylab('Change in RT (s)') + xlab('Trial') + facet_grid(last_outcome~Experiment) + 
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=40,margin=margin(t=0,r=0,b=10,l=10)), 
#         legend.text=element_text(size=20, margin=margin(t=0, r=10, b=8, l=0)), legend.title=element_text(size=30))
# print(gg1)
# dev.off() 
