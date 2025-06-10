# 2025-06-09 AndyP
# plot mixed_by sex effect




load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection/2025-06-09-Sex-clock-fmri_meg-pred-rt-int-1.Rdata')

ddq_mmc <- ddq

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection/2025-06-09-Bsocial-Sex-clock-fmri-pred-rt-int.Rdata')

ddq_bsoc <- ddq


dqRT_bsoc <- ddq_bsoc$emtrends_list$RT %>% filter(dataset1=='bsocial')
dqRT_bsoc <- dqRT_bsoc %>% rename(dataset = dataset1)
dqRT_mmc <- ddq_mmc$emtrends_list$RT

ddqRT <- rbind(dqRT_bsoc,dqRT_mmc)

ddqRT <- ddqRT %>% mutate(dataset1 = case_when(dataset == 'fMRI' ~ 'Experiment 1 - fMRI',
                                               dataset == 'MEG' ~ 'Experiment 1 - MEG',
                                               dataset == 'bsocial' ~ 'Experiment 2')) %>% select(!dataset) %>% rename(dataset = dataset1)


fills <- palette()
fills[1] <- '#E618E6'
fills[2] <- '#181FE6'


ggplot(data=ddqRT, aes(x=sex,y=rt_lag_sc.trend,color=sex,group=sex,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  geom_errorbar(width=0.65,size=1.2) + geom_point(size=2) + facet_grid(~dataset,scales = 'free_y') + scale_y_reverse() + 
  scale_colour_manual(values=fills) +  
  ylab('<-- less -- Exploration -- more -->') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=22),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22),
        strip.text = element_text(size=22))
