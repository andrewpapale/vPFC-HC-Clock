# 2024-06-27 AndyP
# Plotting emmeans and emtrends for B2B

std_of_subject_level_rand_slope = 2

# model iterations 1=HC_within:v_entropy, 2=HC_within:v_max_wi, 3=HC_within (not currently reported)
load('/Volumes/Users/Andrew/v14-vPFC-HC-Figures-and-Models/MLM-Models-v14/2023-08-23-vmPFC-HC-network-ranslopes-clock-pred-int-1.Rdata')

emm <- ddq$emmeans_list
emt <- ddq$emtrends_list

RTxO <- emt$RTxO

RTdmm <- RTxO %>% filter(network=='DMN' & (subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope))
RTdmm <- RTdmm %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
RTdmm$subj_level_rand_slope <- as.factor(RTdmm$subj_level_rand_slope)
RTdmm <- RTdmm  %>% mutate(p_fdr = padj_fdr_term, 
                       p_level_fdr = as.factor(case_when(
                         # p_fdr > .1 ~ '0',
                         # p_fdr < .1 & p_fdr > .05 ~ '1',
                         p_fdr > .05 ~ '1',
                         p_fdr < .05 & p_fdr > .01 ~ '2',
                         p_fdr < .01 & p_fdr > .001 ~ '3',
                         p_fdr <.001 ~ '4')))

RTdmm <- RTdmm %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ 'lower connectivity',
                                              subj_level_rand_slope==std_of_subject_level_rand_slope ~ 'higher connectivity'))

RTdmm$p_level_fdr <- factor(RTdmm$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
ggplot(RTdmm, aes(x=evt_time,y=rt_lag_sc.trend,color=entropy,group=entropy,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error)) + 
  geom_point(size=3) + geom_line(size=1) +
  geom_errorbar(width=0.5) + scale_y_reverse() + facet_grid(last_outcome~HC_region) +
  ylab('<--- less --- Exploration --- more --->') + ggtitle('DMN emtrend')
  #ggtitle('Exploration occurs when DMN-HC Connectivity is \n More Modulated by Entropy')

rm(RTdmm)

RTxO <- emm$RTxO

RTdmm <- RTxO %>% filter(network=='DMN' & (subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope) & HC_region=='AH')
RTdmm <- RTdmm %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
RTdmm$subj_level_rand_slope <- as.factor(RTdmm$subj_level_rand_slope)
RTdmm <- RTdmm  %>% mutate(p_fdr = padj_fdr_term, 
                           p_level_fdr = as.factor(case_when(
                             # p_fdr > .1 ~ '0',
                             # p_fdr < .1 & p_fdr > .05 ~ '1',
                             p_fdr > .05 ~ '1',
                             p_fdr < .05 & p_fdr > .01 ~ '2',
                             p_fdr < .01 & p_fdr > .001 ~ '3',
                             p_fdr <.001 ~ '4')))

RTdmm <- RTdmm %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ 'lower connectivity',
                                              subj_level_rand_slope==std_of_subject_level_rand_slope ~ 'higher connectivity'))

RTdmm$p_level_fdr <- factor(RTdmm$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
ggplot(RTdmm, aes(x=evt_time,y=estimate,color=entropy,group=entropy,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  geom_point(size=3) + geom_line(size=1) +
  geom_errorbar(width=0.5) + scale_y_reverse() + facet_grid(last_outcome~rt_lag_sc) +
  ylab('<--- less --- Exploration --- more --->') + ggtitle('DMN emmeans PH')
#ggtitle('Exploration occurs when DMN-HC Connectivity is \n More Modulated by Entropy')

rm(RTdmm)

RTxO <- emm$RTxO

RTdmm <- RTxO %>% filter(network=='DMN' & (subj_level_rand_slope==-std_of_subject_level_rand_slope | subj_level_rand_slope==std_of_subject_level_rand_slope) & HC_region=='PH')
RTdmm <- RTdmm %>% mutate(padj_fdr_term = p.adjust(p.value,method='bonferroni'))
RTdmm$subj_level_rand_slope <- as.factor(RTdmm$subj_level_rand_slope)
RTdmm <- RTdmm  %>% mutate(p_fdr = padj_fdr_term, 
                           p_level_fdr = as.factor(case_when(
                             # p_fdr > .1 ~ '0',
                             # p_fdr < .1 & p_fdr > .05 ~ '1',
                             p_fdr > .05 ~ '1',
                             p_fdr < .05 & p_fdr > .01 ~ '2',
                             p_fdr < .01 & p_fdr > .001 ~ '3',
                             p_fdr <.001 ~ '4')))

RTdmm <- RTdmm %>% mutate(entropy = case_when(subj_level_rand_slope==-std_of_subject_level_rand_slope ~ 'lower connectivity',
                                              subj_level_rand_slope==std_of_subject_level_rand_slope ~ 'higher connectivity'))

RTdmm$p_level_fdr <- factor(RTdmm$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
ggplot(RTdmm, aes(x=evt_time,y=estimate,color=entropy,group=entropy,ymin=estimate-std.error,ymax=estimate+std.error)) + 
  geom_point(size=3) + geom_line(size=1) +
  geom_errorbar(width=0.5) + scale_y_reverse() + facet_grid(last_outcome~rt_lag_sc) +
  ylab('<--- less --- Exploration --- more --->') + ggtitle('DMN emmeans AH')
#ggtitle('Exploration occurs when DMN-HC Connectivity is \n More Modulated by Entropy')
