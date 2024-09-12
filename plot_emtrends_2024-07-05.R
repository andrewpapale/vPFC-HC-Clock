# 202-07-05 AndyP
# Plot emtrends for Fig 6 of vPFC-HC Paper

library(tidyverse)
library(ggplot2)

# 2024-06-27 AndyP
# Plotting emmeans and emtrends for B2B
do_vPFC_entropy = TRUE
do_vPFC_vmax = FALSE
do_vPFC_HC_entropy = FALSE
do_vPFC_HC_vmax = FALSE
do_vPFC_HC = FALSE

std_of_subject_level_rand_slope = 1

if (do_vPFC_entropy){
  # model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-int-notimesplit-nofixedeffect-1.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-14-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
  emt_mmclock_fmri <- ddq$emtrends_list
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-14-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
  emt_mmclock_meg <- ddq$emtrends_list
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-network-ranslopes-clock-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-notimesplit-nofixedeffect-1.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-network-ranslopes-clock-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
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
  gg1<- ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
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
  print(gg1)
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
  gg1 <- ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
    geom_point(size=8, position=position_dodge(width=0.5), aes(shape = entropy)) + 
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
  print(gg1)
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
  gg1 <- ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
    geom_point(size=8, position=position_dodge(width=0.5), aes(shape = entropy)) + 
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
  print(gg1)
  dev.off()
  
  
  
  RTvmax_mmclock_fmri <- emt_mmclock_fmri$TrxVmax %>% mutate(dataset = 'Experiment 1 fMRI') %>% 
    filter(trial_neg_inv_sc == -0.9 | trial_neg_inv_sc == 0.4)
  RTvmax_mmclock_meg <- emt_mmclock_meg$TrxVmax %>% mutate(dataset = 'Experiment 1 MEG') %>% 
    filter(trial_neg_inv_sc == -0.9 | trial_neg_inv_sc == 0.4)
  RTvmax_explore <- emt_explore$TrxVmax %>% mutate(dataset = 'Experiment 2') %>%
    rename(trial_neg_inv_sc=run_trial0_neg_inv_sc) %>%
    filter(trial_neg_inv_sc == -0.9 | trial_neg_inv_sc == 0.4)
  
  
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
    group_by(dataset,network,entropy,trial_neg_inv_sc) %>% 
    summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>% 
    ungroup() 
  
  
  
  
  library(wesanderson)
  pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
  pal = palette()
  pal[1] = pal3[2]
  pal[2] = pal3[1]
  pal[3] = pal3[3]
  
  pdf('RTvmax-3way-Pred-Entropy-PFC-emtrends.pdf',height=8,width=10)
  gg1 <- ggplot(PFC_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
    geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) + 
    geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') + 
    facet_grid(trial_neg_inv_sc~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) + 
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
  print(gg1)
  dev.off()
  
}

###################
######  Vmax  #####
###################


if (do_vPFC_vmax){
  # model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-int-notimesplit-nofixedeffect-2.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-2.Rdata')
  emt_mmclock_fmri <- ddq$emtrends_list
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-network-ranslopes-clock-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-2.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-network-ranslopes-clock-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-2.Rdata')
  emt_mmclock_meg <- ddq$emtrends_list
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-network-ranslopes-clock-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-notimesplit-nofixedeffect-2.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-network-ranslopes-clock-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-notimesplit-nofixedeffect-rtvmax_by_trial-2.Rdata')
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
  
}

###################
### vPFC-HC #######
###################

if (do_vPFC_HC_entropy){
  # model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
  #load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2023-08-23-vmPFC-HC-network-ranslopes-clock-pred-int-1.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-pred-int-1.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-nofixedeffect-1.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
  mmclock_fmri <- ddq$coef_df_reml %>% 
    filter(effect=='fixed' & term=="rt_lag_sc:subj_level_rand_slope:last_outcomeReward") %>% 
    mutate(p_fdr_me = padj_fdr_term,p_level_fdr_me = as.factor(case_when(p_fdr_me > .05 ~ '1',p_fdr_me < .05 ~ '2'))) %>%
    select(!rhs & !std.error & !statistic & !p.value & !outcome & !model_name)
  emt_mmclock_fmri <- ddq$emtrends_list
  #load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-02-20-vmPFC-HC-network-ranslopes-clock-replication-pred-int-1.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-replication-pred-int-1.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-nofixedeffect-1.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
  emt_mmclock_meg <- ddq$emtrends_list
  mmclock_meg <- ddq$coef_df_reml %>% 
    filter(effect=='fixed' & term=="rt_lag_sc:subj_level_rand_slope:last_outcomeReward") %>% 
    mutate(p_fdr_me = padj_fdr_term,p_level_fdr_me = as.factor(case_when(p_fdr_me > .05 ~ '1',p_fdr_me < .05 ~ '2'))) %>%
    select(!rhs & !std.error & !statistic & !p.value & !outcome & !model_name)
  #load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-04-27-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-1.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-06-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-1.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-1.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
  emt_explore <- ddq$emtrends_list
  explore <- ddq$coef_df_reml %>% 
    filter(effect=='fixed' & term=="rt_lag_sc:subj_level_rand_slope:last_outcomeReward") %>% 
    mutate(p_fdr_me = padj_fdr_term,p_level_fdr_me = as.factor(case_when(p_fdr_me > .05 ~ '1',p_fdr_me < .05 ~ '2'))) %>%
    select(!rhs & !std.error & !statistic & !p.value & !outcome & !model_name)
  
  
  RTxO_mmclock_fmri <- emt_mmclock_fmri$RTxO %>% mutate(dataset = 'Experiment 1 fMRI')
  RTxO_mmclock_meg <- emt_mmclock_meg$RTxO %>% mutate(dataset = 'Experiment 1 MEG')
  RTxO_explore <- emt_explore$RTxO %>% mutate(dataset = 'Experiment 2')
  
  RTxO_mmclock_fmri <- inner_join(RTxO_mmclock_fmri,mmclock_fmri,by=c('network','HC_region'))
  RTxO_mmclock_meg <- inner_join(RTxO_mmclock_meg,mmclock_meg,by=c('network','HC_region'))
  RTxO_explore <- inner_join(RTxO_explore,explore,by=c('network','HC_region'))
  
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
  #PH_merged <- PH_merged %>% mutate(padj_fdr_term_emt = p.adjust(p.value,method='fdr'))
  PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
  PH_merged <- PH_merged  %>% mutate(p_emt = p.value, 
                                     p_level_emt = as.factor(case_when(
                                       p_emt > .05 ~ '1',
                                       p_emt < .05 ~ '2')))
  
  PH_merged <- PH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                        subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))
  
  
  AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
  #AH_merged <- AH_merged %>% mutate(padj_fdr_term_emt = p.adjust(p.value,method='fdr'))
  AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
  AH_merged <- AH_merged %>% mutate(p_emt = p.value, 
                                    p_level_emt = as.factor(case_when(
                                      p_emt > .05 ~ '1',
                                      p_emt < .05 ~ '2')))
  
  AH_merged <- AH_merged %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                        subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))
  
  
  
  AH_summary <- AH_merged %>% 
    group_by(dataset,network,last_outcome,entropy) %>% 
    summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), significant = as.factor(case_when(p_level_emt==2 & p_level_fdr_me==2 ~ 2, 
                                                                                                                           p_level_emt==1 & p_level_fdr_me==2 ~ 1, 
                                                                                                                           p_level_emt==2 & p_level_fdr_me==1 ~ 1, 
                                                                                                                           p_level_emt==1 & p_level_fdr_me==1 ~ 1))) %>% 
    ungroup() 
  
  AH_summary$significant <- relevel(AH_summary$significant, ref='1')
  
  PH_summary <- PH_merged %>% 
    group_by(dataset,network,last_outcome,entropy) %>% 
    summarize(mean_trend = mean(rt_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), significant = as.factor(case_when(p_level_emt==2 & p_level_fdr_me==2 ~ 2, 
                                                                                                                           p_level_emt==1 & p_level_fdr_me==2 ~ 1, 
                                                                                                                           p_level_emt==2 & p_level_fdr_me==1 ~ 1, 
                                                                                                                           p_level_emt==1 & p_level_fdr_me==1 ~ 1))) %>% 
    ungroup() 
  
  PH_summary$significant <- relevel(PH_summary$significant, ref='1')
  
  library(wesanderson)
  pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
  pal = palette()
  pal[1] = pal3[2]
  pal[2] = pal3[1]
  pal[3] = pal3[3]
  
  pdf('RT-Pred-Entropy-AH-emtrends.pdf',height=8,width=10)
  gg1 <- ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
    geom_point(size=8, position=position_dodge(width=0.5), aes(shape = entropy)) + 
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
  print(gg1)
  dev.off()
  
  pdf('RT-Pred-Entropy-PH-emtrends.pdf',height=8,width=10)
  gg1<- ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) + 
    geom_point(size=8, position=position_dodge(width=0.5), aes(shape = entropy)) + 
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
  print(gg1)
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
  gg1<- ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) +
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
  print(gg1)
  dev.off()

  pdf('RT-2way-Pred-Entropy-PH-emtrends.pdf',height=8,width=10)
  gg1 <- ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) +
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
  print(gg1)
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
  gg1<- ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) +
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
  print(gg1)
  dev.off()

  pdf('RTvmax-2way-Pred-Entropy-PH-emtrends.pdf',height=8,width=10)
  gg1<- ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) +
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
  print(gg1)
  dev.off()
  
  
  
  RTvmax_mmclock_fmri <- emt_mmclock_fmri$TrxVmax %>% mutate(dataset = 'Experiment 1 fMRI') %>% 
    filter(trial_neg_inv_sc == -0.9 | trial_neg_inv_sc == 0.4)
  RTvmax_mmclock_meg <- emt_mmclock_meg$TrxVmax %>% mutate(dataset = 'Experiment 1 MEG') %>% 
    filter(trial_neg_inv_sc == -0.9 | trial_neg_inv_sc == 0.4)
  RTvmax_explore <- emt_explore$TrxVmax %>% mutate(dataset = 'Experiment 2') %>%
    rename(trial_neg_inv_sc=run_trial0_neg_inv_sc) %>%
    filter(trial_neg_inv_sc == -0.9 | trial_neg_inv_sc == 0.4)
  
  
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
    group_by(dataset,network,entropy,trial_neg_inv_sc) %>%
    summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>%
    ungroup()
  
  PH_summary <- PH_merged %>%
    group_by(dataset,network,entropy,trial_neg_inv_sc) %>%
    summarize(mean_trend = mean(rt_vmax_lag_sc.trend,na.rm=TRUE), mean_sd = mean(std.error,na.rm=TRUE), sum_p_ls_001 = sum(p_level_fdr =='4')) %>%
    ungroup()
  
  library(wesanderson)
  pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
  pal = palette()
  pal[1] = pal3[2]
  pal[2] = pal3[1]
  pal[3] = pal3[3]
  
  pdf('RTvmax-3way-Pred-Entropy-AH-emtrends.pdf',height=8,width=10)
  gg1<- ggplot(AH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) +
    geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) +
    geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') +
    facet_grid(trial_neg_inv_sc~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) +
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
  print(gg1)
  dev.off()
  
  pdf('RTvmax-3way-Pred-Entropy-PH-emtrends.pdf',height=8,width=10)
  gg1<- ggplot(PH_summary, aes(x=network,y=mean_trend,ymin = mean_trend-mean_sd,ymax=mean_trend+mean_sd,color=network, group=entropy)) +
    geom_point(size=5, position=position_dodge(width=0.5), aes(shape = entropy)) +
    geom_errorbar(width=0.1, position=position_dodge(width=0.5),color = 'black') +
    facet_grid(trial_neg_inv_sc~dataset,labeller = label_wrap_gen(width=16)) + scale_color_manual(values = pal) +
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
  print(gg1)
  dev.off()  
  
  
  
  
}


if (do_vPFC_HC_vmax){
  
  ###################
  ######  Vmax  #####
  ###################
  
  # model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
  #load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2023-08-23-vmPFC-HC-network-ranslopes-clock-pred-int-2.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-pred-int-2.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-nofixedeffect-2.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-2.Rdata')
  emt_mmclock_fmri <- ddq$emtrends_list
  #load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-02-20-vmPFC-HC-network-ranslopes-clock-replication-pred-int-2.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-replication-pred-int-2.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-nofixedeffect-2.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-2.Rdata')
  emt_mmclock_meg <- ddq$emtrends_list
  #load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-04-27-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-2.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-06-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-2.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-2.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-rtvmax_by_trial-2.Rdata')
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
  
}


if (do_vPFC_HC){
  
  ###################
  ### HC within #####
  ###################
  
  
  # model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
  #load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2023-08-23-vmPFC-HC-network-ranslopes-clock-pred-int-3.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-pred-int-3.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-nofixedeffect-3.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-3.Rdata')
  emt_mmclock_fmri <- ddq$emtrends_list
  AH_3way_mmclock_fmri <- ddq$coef_df_reml %>% filter(HC_region=='AH') %>% 
    mutate(dataset = 'Experiment 1 fMRI') %>%
    filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward') %>%
    select(padj_fdr_term,network,dataset)
  PH_3way_mmclock_fmri <- ddq$coef_df_reml %>% filter(HC_region=='PH') %>%
    mutate(dataset = 'Experiment 1 fMRI') %>%
    filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward') %>%
    select(padj_fdr_term,network,dataset)
  #load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-02-20-vmPFC-HC-network-ranslopes-clock-replication-pred-int-3.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-07-vmPFC-HC-network-ranslopes-clock-replication-pred-int-3.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-nofixedeffect-3.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-3.Rdata')
  emt_mmclock_meg <- ddq$emtrends_list
  AH_3way_mmclock_meg <- ddq$coef_df_reml %>% filter(HC_region=='AH') %>% 
    mutate(dataset = 'Experiment 1 MEG') %>%
    filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward') %>%
    select(padj_fdr_term,network,dataset)
  PH_3way_mmclock_meg <- ddq$coef_df_reml %>% filter(HC_region=='PH') %>%
    mutate(dataset = 'Experiment 1 MEG') %>%
    filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward') %>%
    select(padj_fdr_term,network,dataset)
  #load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2024-04-27-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-3.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-06-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-3.Rdata')
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-3.Rdata')
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-rtvmax_by_trial-3.Rdata')
  emt_explore <- ddq$emtrends_list
  AH_3way_explore <- ddq$coef_df_reml %>% filter(HC_region=='AH') %>% 
    mutate(dataset = 'Experiment 2') %>%
    filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward') %>%
    select(padj_fdr_term,network,dataset)
  PH_3way_explore <- ddq$coef_df_reml %>% filter(HC_region=='PH') %>%
    mutate(dataset = 'Experiment 2') %>%
    filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward') %>%
    select(padj_fdr_term,network,dataset)
  
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
  PH_3way_merged <- rbind(PH_3way_mmclock_fmri,PH_3way_mmclock_meg)
  PH_3way_merged <- rbind(PH_3way_merged, PH_3way_explore)
  
  AH_merged <- rbind(AH_mmclock_fmri, AH_mmclock_meg)
  AH_merged <- rbind(AH_merged, AH_explore)
  AH_3way_merged <- rbind(AH_3way_mmclock_fmri,AH_3way_mmclock_meg)
  AH_3way_merged <- rbind(AH_3way_merged, AH_3way_explore)
  
  PH_merged <- inner_join(PH_merged,PH_3way,by=c('network','dataset'))
  AH_merged <- inner_join(AH_merged,AH_3way,by=c('network','dataset'))
  
  PH_merged <- PH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
  PH_merged <- PH_merged %>% mutate(padj_fdr_term_emt = p.adjust(p.value,method='bonferroni'))
  PH_merged$subj_level_rand_slope <- as.factor(PH_merged$subj_level_rand_slope)
  PH_merged <- PH_merged  %>% mutate(p_fdr = padj_fdr_term_emt, 
                                     p_level_fdr_emt = as.factor(case_when(
                                       # p_fdr > .1 ~ '0',
                                       # p_fdr < .1 & p_fdr > .05 ~ '1',
                                       p_fdr > .05 ~ '1',
                                       p_fdr < .05 & p_fdr > .01 ~ '2',
                                       p_fdr < .01 & p_fdr > .001 ~ '3',
                                       p_fdr <.001 ~ '4')))
  
  PH_merged <- PH_merged %>% mutate(HC = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ '-1 STD',
                                                   subj_level_rand_slope==std_of_subject_level_rand_slope ~ '+1 STD'))
  
  
  AH_merged <- AH_merged %>% filter((subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
  AH_merged <- AH_merged %>% mutate(padj_fdr_term_emt = p.adjust(p.value,method='bonferroni'))
  AH_merged$subj_level_rand_slope <- as.factor(AH_merged$subj_level_rand_slope)
  AH_merged <- AH_merged  %>% mutate(p_fdr = padj_fdr_term_emt, 
                                     p_level_fdr_emt = as.factor(case_when(
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
  
}


#############################################
### Get P-values and statistics for Fig 6 ###
#############################################

##################################
# Exploration Entropy AH/PH    ###
##################################

#model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-nofixedeffect-1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddf1_AH_3way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='AH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward')
ddf1_AH_2way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='AH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
ddf1_PH_3way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='PH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward')
ddf1_PH_2way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='PH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
ddf1_AH_RTvmax <- ddq$coef_df_reml %>% filter(network=='DMN' & HC_region=='AH') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
emt_mmclock_fmri <- ddq$emtrends_list
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-nofixedeffect-1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddf2_AH_3way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='AH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward')
ddf2_AH_2way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='AH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
ddf2_PH_3way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='PH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward')
ddf2_PH_2way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='PH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
ddf2_AH_RTvmax <- ddq$coef_df_reml %>% filter(network=='DMN' & HC_region=='AH') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
emt_mmclock_meg <- ddq$emtrends_list
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddf3_AH_3way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='AH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward')
ddf3_AH_2way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='AH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
ddf3_PH_3way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='PH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward')
ddf3_PH_2way <- ddq$coef_df_reml %>% filter((network=='DMN' | network=='LIM') & HC_region=='PH') %>% filter(term=='rt_lag_sc:subj_level_rand_slope')
ddf3_AH_RTvmax <- ddq$coef_df_reml %>% filter(network=='DMN' & HC_region=='AH') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
emt_explore <- ddq$emtrends_list

#########################
# Entropy vPFC RTvmax ###
#########################

#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-network-ranslopes-clock-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-notimesplit-nofixedeffect-1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddf1_RTvmax1 <- ddq$coef_df_reml %>% filter(network=='DMN') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
ddf1_RTvmax2 <- ddq$coef_df_reml %>% filter(network=='DMN') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc:trial_neg_inv_sc')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-network-ranslopes-clock-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-02-vmPFC-network-ranslopes-clock-replication-pred-rt_csv_sc-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddf2_RTvmax1 <- ddq$coef_df_reml %>% filter(network=='DMN') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
ddf2_RTvmax2 <- ddq$coef_df_reml %>% filter(network=='DMN') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc:trial_neg_inv_sc')
#load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-07-10-vmPFC-network-ranslopes-clock-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-notimesplit-nofixedeffect-1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-network-ranslopes-clock-Explore-pred-rt_csv_sc-int-HConly-trial_mod-trial1-10included-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddf3_RTvmax1 <- ddq$coef_df_reml %>% filter(network=='DMN') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
ddf3_RTvmax2 <- ddq$coef_df_reml %>% filter(network=='DMN') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc:run_trial0_neg_inv_sc')


#######################
# Entropy AH-RTvmax ###
#######################

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddf1_HC_RTvmax1 <- ddq$coef_df_reml %>% filter(HC_region =='AH') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
ddf1_HC_RTvmax2 <- ddq$coef_df_reml %>% filter(HC_region =='AH') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc:trial_neg_inv_sc')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-ranslopes-clock-replication-pred-int-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddf2_HC_RTvmax1 <- ddq$coef_df_reml %>% filter(HC_region =='AH') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
ddf2_HC_RTvmax2 <- ddq$coef_df_reml %>% filter(HC_region =='AH') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc:trial_neg_inv_sc')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-08-01-vmPFC-HC-network-Explore-ranslopes-clock-pred-int-trial_mod-trial1-10included-notimesplit-nofixedeffect-rtvmax_by_trial-1.Rdata')
ddf3_HC_RTvmax1 <- ddq$coef_df_reml %>% filter(HC_region =='AH') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
ddf3_HC_RTvmax2 <- ddq$coef_df_reml %>% filter(HC_region =='AH') %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc:run_trial0_neg_inv_sc')

