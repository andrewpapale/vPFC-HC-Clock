# 202-07-05 AndyP
# Plot emtrends for Fig 6 of vPFC-HC Paper

library(tidyverse)
library(ggplot2)

# 2024-06-27 AndyP
# Plotting emmeans and emtrends for B2B

std_of_subject_level_rand_slope = 1

# model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-int-notimesplit-1.Rdata')
emt_mmclock_fmri <- ddq$emtrends_list
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-network-ranslopes-clock-replication-pred-rt_csv_sc-int-notimesplit-1.Rdata')
emt_mmclock_meg <- ddq$emtrends_list
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-network-ranslopes-clock-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-notimesplit-1.Rdata')
emt_explore <- ddq$emtrends_list


RTxO_mmclock_fmri <- emt_mmclock_fmri$RTxO %>% mutate(dataset = 'Experiment 1 fMRI')
RTxO_mmclock_meg <- emt_mmclock_meg$RTxO %>% mutate(dataset = 'Experiment 1 MEG')
RTxO_explore <- emt_explore$RTxO %>% mutate(dataset = 'Experiment 2')

PFC_mmclock_fmri <- RTxO_mmclock_fmri
PFC_mmclock_meg <- RTxO_mmclock_meg 
PFC_explore <- RTxO_explore 

PFC_merged <- rbind(PFC_mmclock_fmri, PFC_mmclock_meg)
PFC_merged <- rbind(PFC_merged, PFC_explore)

PFC_merged <- PFC_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PFC_merged <- PFC_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PFC_merged$subj_level_rand_slope <- as.factor(PFC_merged$subj_level_rand_slope)
PFC_merged <- PFC_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                     p_level_fdr = as.factor(case_when(
                                       # p_fdr > .1 ~ '0',
                                       # p_fdr < .1 & p_fdr > .05 ~ '1',
                                       p_fdr > .05 ~ '1',
                                       p_fdr < .05 & p_fdr > .01 ~ '2',
                                       p_fdr < .01 & p_fdr > .001 ~ '3',
                                       p_fdr <.001 ~ '4')))

PFC_merged <- PFC_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                     subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



PFC_summary <- PFC_merged %>% 
  group_by(dataset,network,last_outcome,entropy) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 


library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-Pred-Entropy-PFC-emtrends.pdf',height=8,width=10)
ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(last_outcome~dataset,labeller = label_wrap_gen(width=16), scales='free_y') + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  xlab('Network') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()


RT_mmclock_fmri <- emt_mmclock_fmri$RT %>% mutate(dataset = 'Experiment 1 fMRI')
RT_mmclock_meg <- emt_mmclock_meg$RT %>% mutate(dataset = 'Experiment 1 MEG')
RT_explore <- emt_explore$RT %>% mutate(dataset = 'Experiment 2')

PFC_mmclock_fmri <- RT_mmclock_fmri 
PFC_mmclock_meg <- RT_mmclock_meg 
PFC_explore <- RT_explore

PFC_merged <- rbind(PFC_mmclock_fmri, PFC_mmclock_meg)
PFC_merged <- rbind(PFC_merged, PFC_explore)

PFC_merged <- PFC_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PFC_merged <- PFC_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PFC_merged$subj_level_rand_slope <- as.factor(PFC_merged$subj_level_rand_slope)
PFC_merged <- PFC_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                     p_level_fdr = as.factor(case_when(
                                       # p_fdr > .1 ~ '0',
                                       # p_fdr < .1 & p_fdr > .05 ~ '1',
                                       p_fdr > .05 ~ '1',
                                       p_fdr < .05 & p_fdr > .01 ~ '2',
                                       p_fdr < .01 & p_fdr > .001 ~ '3',
                                       p_fdr <.001 ~ '4')))

PFC_merged <- PFC_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                        subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



PFC_summary <- PFC_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-2way-Pred-Entropy-PFC-emtrends.pdf',height=8,width=10)
ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()


RTvmax_mmclock_fmri <- emt_mmclock_fmri$Vmax %>% mutate(dataset = 'Experiment 1 fMRI')
RTvmax_mmclock_meg <- emt_mmclock_meg$Vmax %>% mutate(dataset = 'Experiment 1 MEG')
RTvmax_explore <- emt_explore$Vmax %>% mutate(dataset = 'Experiment 2')

PFC_mmclock_fmri <- RTvmax_mmclock_fmri 
PFC_mmclock_meg <- RTvmax_mmclock_meg 
PFC_explore <- RTvmax_explore

PFC_merged <- rbind(PFC_mmclock_fmri, PFC_mmclock_meg)
PFC_merged <- rbind(PFC_merged, PFC_explore)

PFC_merged <- PFC_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PFC_merged <- PFC_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PFC_merged$subj_level_rand_slope <- as.factor(PFC_merged$subj_level_rand_slope)
PFC_merged <- PFC_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                     p_level_fdr = as.factor(case_when(
                                       # p_fdr > .1 ~ '0',
                                       # p_fdr < .1 & p_fdr > .05 ~ '1',
                                       p_fdr > .05 ~ '1',
                                       p_fdr < .05 & p_fdr > .01 ~ '2',
                                       p_fdr < .01 & p_fdr > .001 ~ '3',
                                       p_fdr <.001 ~ '4')))

PFC_merged <- PFC_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                        subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



PFC_summary <- PFC_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 




library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RTvmax-2way-Pred-Entropy-PFC-emtrends.pdf',height=8,width=10)
ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploitation - more ->') +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()



###################
######  Vmax  #####
###################

# model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-int-notimesplit-2.Rdata')
emt_mmclock_fmri <- ddq$emtrends_list
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-network-ranslopes-clock-replication-pred-rt_csv_sc-int-notimesplit-2.Rdata')
emt_mmclock_meg <- ddq$emtrends_list
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-network-ranslopes-clock-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-notimesplit-2.Rdata')
emt_explore <- ddq$emtrends_list


RTxO_mmclock_fmri <- emt_mmclock_fmri$RTxO %>% mutate(dataset = 'Experiment 1 fMRI')
RTxO_mmclock_meg <- emt_mmclock_meg$RTxO %>% mutate(dataset = 'Experiment 1 MEG')
RTxO_explore <- emt_explore$RTxO %>% mutate(dataset = 'Experiment 2')

PFC_mmclock_fmri <- RTxO_mmclock_fmri
PFC_mmclock_meg <- RTxO_mmclock_meg 
PFC_explore <- RTxO_explore 

PFC_merged <- rbind(PFC_mmclock_fmri, PFC_mmclock_meg)
PFC_merged <- rbind(PFC_merged, PFC_explore)

PFC_merged <- PFC_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PFC_merged <- PFC_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PFC_merged$subj_level_rand_slope <- as.factor(PFC_merged$subj_level_rand_slope)
PFC_merged <- PFC_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PFC_merged <- PFC_merged %>% mutate(Vmax = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                   subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



PFC_summary <- PFC_merged %>% 
  group_by(dataset,network,last_outcome,Vmax) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 


library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-Pred-Vmax-PFC-emtrends.pdf',height=8,width=10)
ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=Vmax)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = Vmax)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(last_outcome~dataset,labeller = label_wrap_gen(width=16), scales='free_y') + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  xlab('Network') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()


RT_mmclock_fmri <- emt_mmclock_fmri$RT %>% mutate(dataset = 'Experiment 1 fMRI')
RT_mmclock_meg <- emt_mmclock_meg$RT %>% mutate(dataset = 'Experiment 1 MEG')
RT_explore <- emt_explore$RT %>% mutate(dataset = 'Experiment 2')

PFC_mmclock_fmri <- RT_mmclock_fmri 
PFC_mmclock_meg <- RT_mmclock_meg 
PFC_explore <- RT_explore

PFC_merged <- rbind(PFC_mmclock_fmri, PFC_mmclock_meg)
PFC_merged <- rbind(PFC_merged, PFC_explore)

PFC_merged <- PFC_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PFC_merged <- PFC_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PFC_merged$subj_level_rand_slope <- as.factor(PFC_merged$subj_level_rand_slope)
PFC_merged <- PFC_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                     p_level_fdr = as.factor(case_when(
                                       # p_fdr > .1 ~ '0',
                                       # p_fdr < .1 & p_fdr > .05 ~ '1',
                                       p_fdr > .05 ~ '1',
                                       p_fdr < .05 & p_fdr > .01 ~ '2',
                                       p_fdr < .01 & p_fdr > .001 ~ '3',
                                       p_fdr <.001 ~ '4')))

PFC_merged <- PFC_merged %>% mutate(Vmax = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                        subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



PFC_summary <- PFC_merged %>% 
  group_by(dataset,network,Vmax) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-2way-Pred-Vmax-PFC-emtrends.pdf',height=8,width=10)
ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=Vmax)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = Vmax)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()


RTvmax_mmclock_fmri <- emt_mmclock_fmri$Vmax %>% mutate(dataset = 'Experiment 1 fMRI')
RTvmax_mmclock_meg <- emt_mmclock_meg$Vmax %>% mutate(dataset = 'Experiment 1 MEG')
RTvmax_explore <- emt_explore$Vmax %>% mutate(dataset = 'Experiment 2')

PFC_mmclock_fmri <- RTvmax_mmclock_fmri 
PFC_mmclock_meg <- RTvmax_mmclock_meg 
PFC_explore <- RTvmax_explore

PFC_merged <- rbind(PFC_mmclock_fmri, PFC_mmclock_meg)
PFC_merged <- rbind(PFC_merged, PFC_explore)

PFC_merged <- PFC_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PFC_merged <- PFC_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PFC_merged$subj_level_rand_slope <- as.factor(PFC_merged$subj_level_rand_slope)
PFC_merged <- PFC_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PFC_merged <- PFC_merged %>% mutate(Vmax = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



PFC_summary <- PFC_merged %>% 
  group_by(dataset,network,Vmax) %>% 
  summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 




library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RTvmax-2way-Pred-Vmax-PFC-emtrends.pdf',height=8,width=10)
ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=Vmax)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = Vmax)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploitation - more ->') +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()


###################
### vPFC-HC #######
###################




std_of_subject_level_rand_slope = 1

# model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
#load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2023-08-23-vmPFC-HC-network-ranslopes-clock-pred-int-1.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-pred-int-1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-1.Rdata')
emt_mmclock_fmri <- ddq$emtrends_list
#load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-02-20-vmPFC-HC-network-ranslopes-clock-replication-pred-int-1.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-replication-pred-int-1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-1.Rdata')
emt_mmclock_meg <- ddq$emtrends_list
#load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-04-27-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-1.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-06-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-1.Rdata')
emt_explore <- ddq$emtrends_list


RTxO_mmclock_fmri <- emt_mmclock_fmri$RTxO %>% mutate(dataset = 'Experiment 1 fMRI')
RTxO_mmclock_meg <- emt_mmclock_meg$RTxO %>% mutate(dataset = 'Experiment 1 MEG')
RTxO_explore <- emt_explore$RTxO %>% mutate(dataset = 'Experiment 2')

PH_mmclock_fmri <- RTxO_mmclock_fmri %>% filter(HC_region == 'PH')
PH_mmclock_meg <- RTxO_mmclock_meg %>% filter(HC_region == 'PH')
PH_explore <- RTxO_explore %>% filter(HC_region == 'PH')

AH_mmclock_fmri <- RTxO_mmclock_fmri %>% filter(HC_region == 'AH')
AH_mmclock_meg <- RTxO_mmclock_meg %>% filter(HC_region == 'AH')
AH_explore <- RTxO_explore %>% filter(HC_region == 'AH')

PH_merged <- rbind(PH_mmclock_fmri,PH_mmclock_meg)
PH_merged <- rbind(PH_merged, PH_explore)

AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
AH_merged <- rbind(AH_merged, AH_explore)


PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PH_merged <- PH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PH_merged <- PH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))


AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
AH_merged <- AH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

AH_merged <- AH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



AH_summary <- AH_merged %>% 
  group_by(dataset,network,last_outcome,entropy) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 

PH_summary <- PH_merged %>% 
  group_by(dataset,network,last_outcome,entropy) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-Pred-Entropy-AH-emtrends.pdf',height=8,width=10)
ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(last_outcome~dataset,labeller = label_wrap_gen(width=16), scales='free_y') + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  xlab('Network') +
  ggtitle('Anterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

pdf('RT-Pred-Entropy-PH-emtrends.pdf',height=8,width=10)
ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(last_outcome~dataset,labeller = label_wrap_gen(width=16), scales='free_y') + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  xlab('Network') +
  ggtitle('Posterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()


RT_mmclock_fmri <- emt_mmclock_fmri$RT %>% mutate(dataset = 'Experiment 1 fMRI')
RT_mmclock_meg <- emt_mmclock_meg$RT %>% mutate(dataset = 'Experiment 1 MEG')
RT_explore <- emt_explore$RT %>% mutate(dataset = 'Experiment 2')

PH_mmclock_fmri <- RT_mmclock_fmri %>% filter(HC_region == 'PH')
PH_mmclock_meg <- RT_mmclock_meg %>% filter(HC_region == 'PH')
PH_explore <- RT_explore %>% filter(HC_region == 'PH')

AH_mmclock_fmri <- RT_mmclock_fmri %>% filter(HC_region == 'AH')
AH_mmclock_meg <- RT_mmclock_meg %>% filter(HC_region == 'AH')
AH_explore <- RT_explore %>% filter(HC_region == 'AH')

PH_merged <- rbind(PH_mmclock_fmri,PH_mmclock_meg)
PH_merged <- rbind(PH_merged, PH_explore)

AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
AH_merged <- rbind(AH_merged, AH_explore)


PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PH_merged <- PH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PH_merged <- PH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))


AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
AH_merged <- AH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

AH_merged <- AH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



AH_summary <- AH_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 

PH_summary <- PH_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-2way-Pred-Entropy-AH-emtrends.pdf',height=8,width=10)
ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  ggtitle('Anterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

pdf('RT-2way-Pred-Entropy-PH-emtrends.pdf',height=8,width=10)
ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  ggtitle('Posterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

# RTdmm$p_level_fdr <- factor(RTdmm$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
# ggplot(RTdmm, aes(x=evt_time,y=rt_lag_sc.trend,color=entropy,group=entropy,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
#   geom_point(size=3) + geom_line(size=1) +
#   geom_errorbar(width=0.5) + scale_y_reverse() + facet_grid(last_outcome~HC_region) +
#   ylab('<--- less --- Exploration --- more --->') + ggtitle('DMN emtrend')
# #ggtitle('Exploration occurs when DMN-HC Connectivity is \n More Modulated by Entropy')


RTvmax_mmclock_fmri <- emt_mmclock_fmri$Vmax %>% mutate(dataset = 'Experiment 1 fMRI')
RTvmax_mmclock_meg <- emt_mmclock_meg$Vmax %>% mutate(dataset = 'Experiment 1 MEG')
RTvmax_explore <- emt_explore$Vmax %>% mutate(dataset = 'Experiment 2')

PH_mmclock_fmri <- RTvmax_mmclock_fmri %>% filter(HC_region == 'PH')
PH_mmclock_meg <- RTvmax_mmclock_meg %>% filter(HC_region == 'PH')
PH_explore <- RTvmax_explore %>% filter(HC_region == 'PH')

AH_mmclock_fmri <- RTvmax_mmclock_fmri %>% filter(HC_region == 'AH')
AH_mmclock_meg <- RTvmax_mmclock_meg %>% filter(HC_region == 'AH')
AH_explore <- RTvmax_explore %>% filter(HC_region == 'AH')

PH_merged <- rbind(PH_mmclock_fmri,PH_mmclock_meg)
PH_merged <- rbind(PH_merged, PH_explore)

AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
AH_merged <- rbind(AH_merged, AH_explore)


PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PH_merged <- PH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PH_merged <- PH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))


AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
AH_merged <- AH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

AH_merged <- AH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



AH_summary <- AH_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 

PH_summary <- PH_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RTvmax-2way-Pred-Entropy-AH-emtrends.pdf',height=8,width=10)
ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploitation - more ->') +
  theme_bw(base_size=13) +
  ggtitle('Anterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

pdf('RTvmax-2way-Pred-Entropy-PH-emtrends.pdf',height=8,width=10)
ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploitation - more ->') +
  theme_bw(base_size=13) +
  ggtitle('Posterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

###################
######  Vmax  #####
###################

# model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
#load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2023-08-23-vmPFC-HC-network-ranslopes-clock-pred-int-2.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-pred-int-2.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-2.Rdata')
emt_mmclock_fmri <- ddq$emtrends_list
#load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-02-20-vmPFC-HC-network-ranslopes-clock-replication-pred-int-2.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-replication-pred-int-2.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-2.Rdata')
emt_mmclock_meg <- ddq$emtrends_list
#load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-04-27-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-2.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-06-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-2.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-2.Rdata')
emt_explore <- ddq$emtrends_list


RTxO_mmclock_fmri <- emt_mmclock_fmri$RTxO %>% mutate(dataset = 'Experiment 1 fMRI')
RTxO_mmclock_meg <- emt_mmclock_meg$RTxO %>% mutate(dataset = 'Experiment 1 MEG')
RTxO_explore <- emt_explore$RTxO %>% mutate(dataset = 'Experiment 2')

PH_mmclock_fmri <- RTxO_mmclock_fmri %>% filter(HC_region == 'PH')
PH_mmclock_meg <- RTxO_mmclock_meg %>% filter(HC_region == 'PH')
PH_explore <- RTxO_explore %>% filter(HC_region == 'PH')

AH_mmclock_fmri <- RTxO_mmclock_fmri %>% filter(HC_region == 'AH')
AH_mmclock_meg <- RTxO_mmclock_meg %>% filter(HC_region == 'AH')
AH_explore <- RTxO_explore %>% filter(HC_region == 'AH')

PH_merged <- rbind(PH_mmclock_fmri,PH_mmclock_meg)
PH_merged <- rbind(PH_merged, PH_explore)

AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
AH_merged <- rbind(AH_merged, AH_explore)


PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PH_merged <- PH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PH_merged <- PH_merged %>% mutate(Vmax = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                   subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))


AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
AH_merged <- AH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

AH_merged <- AH_merged %>% mutate(Vmax = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                   subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



AH_summary <- AH_merged %>% 
  group_by(dataset,network,last_outcome,Vmax) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 

PH_summary <- PH_merged %>% 
  group_by(dataset,network,last_outcome,Vmax) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-Pred-Vmax-AH-emtrends.pdf',height=8,width=10)
ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=Vmax)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = Vmax)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(last_outcome~dataset,labeller = label_wrap_gen(width=16), scales='free_y') + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  xlab('Network') +
  ggtitle('Anterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

pdf('RT-Pred-Vmax-PH-emtrends.pdf',height=8,width=10)
ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=Vmax)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = Vmax)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(last_outcome~dataset,labeller = label_wrap_gen(width=16), scales='free_y') + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  xlab('Network') +
  ggtitle('Posterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()


RT_mmclock_fmri <- emt_mmclock_fmri$RT %>% mutate(dataset = 'Experiment 1 fMRI')
RT_mmclock_meg <- emt_mmclock_meg$RT %>% mutate(dataset = 'Experiment 1 MEG')
RT_explore <- emt_explore$RT %>% mutate(dataset = 'Experiment 2')

PH_mmclock_fmri <- RT_mmclock_fmri %>% filter(HC_region == 'PH')
PH_mmclock_meg <- RT_mmclock_meg %>% filter(HC_region == 'PH')
PH_explore <- RT_explore %>% filter(HC_region == 'PH')

AH_mmclock_fmri <- RT_mmclock_fmri %>% filter(HC_region == 'AH')
AH_mmclock_meg <- RT_mmclock_meg %>% filter(HC_region == 'AH')
AH_explore <- RT_explore %>% filter(HC_region == 'AH')

PH_merged <- rbind(PH_mmclock_fmri,PH_mmclock_meg)
PH_merged <- rbind(PH_merged, PH_explore)

AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
AH_merged <- rbind(AH_merged, AH_explore)


PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PH_merged <- PH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PH_merged <- PH_merged %>% mutate(Vmax = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                   subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))


AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
AH_merged <- AH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

AH_merged <- AH_merged %>% mutate(Vmax = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                   subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



AH_summary <- AH_merged %>% 
  group_by(dataset,network,Vmax) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 

PH_summary <- PH_merged %>% 
  group_by(dataset,network,Vmax) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-2way-Pred-Vmax-AH-emtrends.pdf',height=8,width=10)
ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=Vmax)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = Vmax)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  ggtitle('Anterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

pdf('RT-2way-Pred-Vmax-PH-emtrends.pdf',height=8,width=10)
ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=Vmax)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = Vmax)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  ggtitle('Posterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

RTvmax_mmclock_fmri <- emt_mmclock_fmri$Vmax %>% mutate(dataset = 'Experiment 1 fMRI')
RTvmax_mmclock_meg <- emt_mmclock_meg$Vmax %>% mutate(dataset = 'Experiment 1 MEG')
RTvmax_explore <- emt_explore$Vmax %>% mutate(dataset = 'Experiment 2')

PH_mmclock_fmri <- RTvmax_mmclock_fmri %>% filter(HC_region == 'PH')
PH_mmclock_meg <- RTvmax_mmclock_meg %>% filter(HC_region == 'PH')
PH_explore <- RTvmax_explore %>% filter(HC_region == 'PH')

AH_mmclock_fmri <- RTvmax_mmclock_fmri %>% filter(HC_region == 'AH')
AH_mmclock_meg <- RTvmax_mmclock_meg %>% filter(HC_region == 'AH')
AH_explore <- RTvmax_explore %>% filter(HC_region == 'AH')

PH_merged <- rbind(PH_mmclock_fmri,PH_mmclock_meg)
PH_merged <- rbind(PH_merged, PH_explore)

AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
AH_merged <- rbind(AH_merged, AH_explore)


PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PH_merged <- PH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PH_merged <- PH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))


AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
AH_merged <- AH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

AH_merged <- AH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



AH_summary <- AH_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 

PH_summary <- PH_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RTvmax-2way-Pred-Vmax-AH-emtrends.pdf',height=8,width=10)
ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploitation - more ->') +
  theme_bw(base_size=13) +
  ggtitle('Anterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

pdf('RTvmax-2way-Pred-Vmax-PH-emtrends.pdf',height=8,width=10)
ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploitation - more ->') +
  theme_bw(base_size=13) +
  ggtitle('Posterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()


###################
### HC within #####
###################


# model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
#load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2023-08-23-vmPFC-HC-network-ranslopes-clock-pred-int-3.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-pred-int-3.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-3.Rdata')
emt_mmclock_fmri <- ddq$emtrends_list
#load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-02-20-vmPFC-HC-network-ranslopes-clock-replication-pred-int-3.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-replication-pred-int-3.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-3.Rdata')
emt_mmclock_meg <- ddq$emtrends_list
#load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-04-27-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-3.Rdata')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-06-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-3.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-3.Rdata')
emt_explore <- ddq$emtrends_list


RTxO_mmclock_fmri <- emt_mmclock_fmri$RTxO %>% mutate(dataset = 'Experiment 1 fMRI')
RTxO_mmclock_meg <- emt_mmclock_meg$RTxO %>% mutate(dataset = 'Experiment 1 MEG')
RTxO_explore <- emt_explore$RTxO %>% mutate(dataset = 'Experiment 2')

PH_mmclock_fmri <- RTxO_mmclock_fmri %>% filter(HC_region == 'PH')
PH_mmclock_meg <- RTxO_mmclock_meg %>% filter(HC_region == 'PH')
PH_explore <- RTxO_explore %>% filter(HC_region == 'PH')

AH_mmclock_fmri <- RTxO_mmclock_fmri %>% filter(HC_region == 'AH')
AH_mmclock_meg <- RTxO_mmclock_meg %>% filter(HC_region == 'AH')
AH_explore <- RTxO_explore %>% filter(HC_region == 'AH')

PH_merged <- rbind(PH_mmclock_fmri,PH_mmclock_meg)
PH_merged <- rbind(PH_merged, PH_explore)

AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
AH_merged <- rbind(AH_merged, AH_explore)


PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PH_merged <- PH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PH_merged <- PH_merged %>% mutate(HC = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                 subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))


AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
AH_merged <- AH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

AH_merged <- AH_merged %>% mutate(HC = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                 subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



AH_summary <- AH_merged %>% 
  group_by(dataset,network,last_outcome,HC) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 

PH_summary <- PH_merged %>% 
  group_by(dataset,network,last_outcome,HC) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-Pred-HC-AH-emtrends.pdf',height=8,width=10)
ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=HC)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = HC)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(last_outcome~dataset,labeller = label_wrap_gen(width=16), scales='free_y') + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  xlab('Network') +
  ggtitle('Anterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

pdf('RT-Pred-HC-PH-emtrends.pdf',height=8,width=10)
ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=HC)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = HC)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(last_outcome~dataset,labeller = label_wrap_gen(width=16), scales='free_y') + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  xlab('Network') +
  ggtitle('Posterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()


RT_mmclock_fmri <- emt_mmclock_fmri$RT %>% mutate(dataset = 'Experiment 1 fMRI')
RT_mmclock_meg <- emt_mmclock_meg$RT %>% mutate(dataset = 'Experiment 1 MEG')
RT_explore <- emt_explore$RT %>% mutate(dataset = 'Experiment 2')

PH_mmclock_fmri <- RT_mmclock_fmri %>% filter(HC_region == 'PH')
PH_mmclock_meg <- RT_mmclock_meg %>% filter(HC_region == 'PH')
PH_explore <- RT_explore %>% filter(HC_region == 'PH')

AH_mmclock_fmri <- RT_mmclock_fmri %>% filter(HC_region == 'AH')
AH_mmclock_meg <- RT_mmclock_meg %>% filter(HC_region == 'AH')
AH_explore <- RT_explore %>% filter(HC_region == 'AH')

PH_merged <- rbind(PH_mmclock_fmri,PH_mmclock_meg)
PH_merged <- rbind(PH_merged, PH_explore)

AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
AH_merged <- rbind(AH_merged, AH_explore)


PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PH_merged <- PH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PH_merged <- PH_merged %>% mutate(HC = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                 subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))


AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
AH_merged <- AH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

AH_merged <- AH_merged %>% mutate(HC = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                 subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



AH_summary <- AH_merged %>% 
  group_by(dataset,network,HC) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 

PH_summary <- PH_merged %>% 
  group_by(dataset,network,HC) %>% 
  summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RT-2way-Pred-HC-AH-emtrends.pdf',height=8,width=10)
ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=HC)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = HC)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  ggtitle('Anterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

pdf('RT-2way-Pred-HC-PH-emtrends.pdf',height=8,width=10)
ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=HC)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = HC)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploration - more ->') + scale_y_reverse() +
  theme_bw(base_size=13) +
  ggtitle('Posterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

RTvmax_mmclock_fmri <- emt_mmclock_fmri$Vmax %>% mutate(dataset = 'Experiment 1 fMRI')
RTvmax_mmclock_meg <- emt_mmclock_meg$Vmax %>% mutate(dataset = 'Experiment 1 MEG')
RTvmax_explore <- emt_explore$Vmax %>% mutate(dataset = 'Experiment 2')

PH_mmclock_fmri <- RTvmax_mmclock_fmri %>% filter(HC_region == 'PH')
PH_mmclock_meg <- RTvmax_mmclock_meg %>% filter(HC_region == 'PH')
PH_explore <- RTvmax_explore %>% filter(HC_region == 'PH')

AH_mmclock_fmri <- RTvmax_mmclock_fmri %>% filter(HC_region == 'AH')
AH_mmclock_meg <- RTvmax_mmclock_meg %>% filter(HC_region == 'AH')
AH_explore <- RTvmax_explore %>% filter(HC_region == 'AH')

PH_merged <- rbind(PH_mmclock_fmri,PH_mmclock_meg)
PH_merged <- rbind(PH_merged, PH_explore)

AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
AH_merged <- rbind(AH_merged, AH_explore)


PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
PH_merged <- PH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

PH_merged <- PH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))


AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
AH_merged <- AH_merged %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term, 
                                   p_level_fdr = as.factor(case_when(
                                     # p_fdr > .1 ~ '0',
                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
                                     p_fdr > .05 ~ '1',
                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                     p_fdr <.001 ~ '4')))

AH_merged <- AH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                      subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))



AH_summary <- AH_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 

PH_summary <- PH_merged %>% 
  group_by(dataset,network,entropy) %>% 
  summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
  ungroup() 



library(wesanderson)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]

pdf('RTvmax-2way-Pred-HC-AH-emtrends.pdf',height=8,width=10)
ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploitation - more ->') +
  theme_bw(base_size=13) +
  ggtitle('Anterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

pdf('RTvmax-2way-Pred-HC-PH-emtrends.pdf',height=8,width=10)
ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
  geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
  geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
  facet_grid(~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1,16)) +
  ylab('<- less - Exploitation - more ->') +
  theme_bw(base_size=13) +
  ggtitle('Posterior Hippocampus') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.spacing = unit(0.25,"lines"),
        strip.text = element_text(size=22),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust=0.1, size=22))
dev.off()

