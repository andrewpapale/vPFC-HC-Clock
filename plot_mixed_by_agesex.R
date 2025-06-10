# 2025-06-09 AndyP
# plot mixed_by sex effect




load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection/2024-11-25-Sex-clock-fmri_meg-pred-rt-int-1.Rdata')

ddq_mmc <- ddq

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection/2025-06-09-Bsocial-Sex-clock-fmri-pred-rt-int.Rdata')

ddq_bsoc <- ddq


dqRT_bsoc <- ddq_bsoc$emtrends_list$RTxO %>% filter(dataset1=='bsocial')
dqRT_bsoc <- dqRT_bsoc %>% rename(dataset = dataset1)
dqRT_mmc <- ddq_mmc$emtrends_list$RTxO

ddqRT <- rbind(dqRT_bsoc,dqRT_mmc)

ddqRT <- ddqRT %>% mutate(dataset1 = case_when(dataset == 'fMRI' ~ 'Experiment 1 - fMRI',
                                               dataset == 'MEG' ~ 'Experiment 1 - MEG',
                                               dataset == 'bsocial' ~ 'Experiment 2')) %>% select(!dataset) %>% rename(dataset = dataset1)

ggplot(data=ddqRT, aes(x=sex,y=rt_lag_sc.trend,color=sex,group=sex,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  geom_errorbar() + facet_grid(last_outcome~dataset,scales = 'free_y') + scale_y_reverse()
