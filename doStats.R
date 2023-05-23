# 2023-05-16 AndyP
# doStats Fig 5

library(tidyverse)

# Fig 5A MMClock HC Value Response
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-16-HC-region-clock-2.Rdata')
ddf <- ddf$coef_df_reml
mmclock_hc <- ddf %>% filter(effect=='fixed' & term=='v_max_wi')
m_hc <- mmclock_hc %>% mutate(p_level = case_when(padj_fdr_term <= 0.05 & padj_fdr_term > 0.01 ~ 1,
                                                        padj_fdr_term <= 0.01 & padj_fdr_term > 0.001 ~ 2,
                                                        padj_fdr_term <= 0.001 ~ 3,
                                                        padj_fdr_term > 0.05 ~ 0)) %>%
  group_by(HC_region) %>% summarize(mstat = mean(statistic,na.rm=TRUE), 
                                            mdf = mean(df,na.rm=TRUE),
                                            s05 = sum(p_level==1, na.rm=TRUE),
                                            s01 = sum(p_level==2,na.rm=TRUE),
                                            s001 = sum(p_level==3,na.rm=TRUE),
                                            snP = length(statistic))

# Fig 5B Explore HC Value Response
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-11-Explore-HC-region-clock-HConly-2.Rdata')
ddf <- ddf$coef_df_reml
explore_hc <- ddf %>% filter(effect=='fixed' & term=='v_max_wi')
e_hc <- mmclock_hc %>% mutate(p_level = case_when(padj_fdr_term <= 0.05 & padj_fdr_term > 0.01 ~ 1,
                                                  padj_fdr_term <= 0.01 & padj_fdr_term > 0.001 ~ 2,
                                                  padj_fdr_term <= 0.001 ~ 3,
                                                  padj_fdr_term > 0.05 ~ 0)) %>%
  group_by(HC_region) %>% summarize(mstat = mean(statistic,na.rm=TRUE), 
                                    mdf = mean(df,na.rm=TRUE),
                                    s05 = sum(p_level==1, na.rm=TRUE),
                                    s01 = sum(p_level==2,na.rm=TRUE),
                                    s001 = sum(p_level==3,na.rm=TRUE),
                                    snP = length(statistic))


# Figure 5C MMClock HCwithin
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-09-vmPFC-HC-network-clock-1.Rdata')
ddf <- ddf$coef_df_reml
mmclock_hcwithin <- ddf %>% filter(effect=='fixed' & term=='HCwithin')

m_hc <- mmclock_hcwithin %>% mutate(p_level = case_when(padj_fdr_term <= 0.05 & padj_fdr_term > 0.01 ~ 1,
                                                        padj_fdr_term <= 0.01 & padj_fdr_term > 0.001 ~ 2,
                                                        padj_fdr_term <= 0.001 ~ 3,
                                                        padj_fdr_term > 0.05 ~ 0)) %>%
  group_by(network,HC_region) %>% summarize(mstat = mean(statistic,na.rm=TRUE), 
                                            mdf = mean(df,na.rm=TRUE),
                                            s05 = sum(p_level==1, na.rm=TRUE),
                                            s01 = sum(p_level==2,na.rm=TRUE),
                                            s001 = sum(p_level==3,na.rm=TRUE),
                                            snP = length(statistic))
m_m1 <- aov(estimate ~ network*HC_region + evt_time, mmclock_hcwithin)
TukeyHSD(m_m1)

# Figure 5D Explore HCwithin
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-10-Explore-vPFC-HC-network-clock-HConly-1.Rdata')
ddf <- ddf$coef_df_reml
explore_hcwithin <- ddf %>% filter(effect=='fixed' & term=='HCwithin')

e_hc <- explore_hcwithin %>% mutate(p_level = case_when(padj_fdr_term <= 0.05 & padj_fdr_term > 0.01 ~ 1,
                                                        padj_fdr_term <= 0.01 & padj_fdr_term > 0.001 ~ 2,
                                                        padj_fdr_term <= 0.001 ~ 3,
                                                        padj_fdr_term > 0.05 ~ 0)) %>%
  group_by(network,HC_region) %>% summarize(mstat = mean(statistic,na.rm=TRUE), 
                                            mdf = mean(df,na.rm=TRUE),
                                            s05 = sum(p_level==1, na.rm=TRUE),
                                            s01 = sum(p_level==2,na.rm=TRUE),
                                            s001 = sum(p_level==3,na.rm=TRUE),
                                            snP = length(statistic))
e_m1 <- aov(estimate ~ network*HC_region + evt_time, explore_hcwithin)
TukeyHSD(e_m1)


#mmclock_hcwithin <- mmclock_hcwithin %>% 
#  mutate(dataset = 'MMClock', network1 = case_when(network=='C'~'CTR',network=='D'~'DMN',network=='L'~'LIM')) %>%
#  select(!network) %>% rename(network = network1)
#explore_hcwithin <- explore_hcwithin %>% mutate(dataset = 'Explore')
#c_hc <- rbind(mmclock_hcwithin,explore_hcwithin)
#c_hc$dataset <- as.factor(c_hc$dataset)
#c_hc$HC_region <- as.factor(c_hc$HC_region)
#c_m1 <- aov(estimate ~ network*dataset, c_hc)


# Figure 5E entropy:HCwithin
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-09-vmPFC-HC-network-clock-1.Rdata')
ddf <- ddf$coef_df_reml
mmclock_hcwithin <- ddf %>% filter(effect=='fixed' & term=='HCwithin:v_entropy_wi')
m_hc <- mmclock_hcwithin %>% mutate(p_level = case_when(padj_fdr_term <= 0.05 & padj_fdr_term > 0.01 ~ 1,
                                                        padj_fdr_term <= 0.01 & padj_fdr_term > 0.001 ~ 2,
                                                        padj_fdr_term <= 0.001 ~ 3,
                                                        padj_fdr_term > 0.05 ~ 0)) %>%
  group_by(network,HC_region) %>% summarize(mstat = mean(statistic,na.rm=TRUE), 
                                            mdf = mean(df,na.rm=TRUE),
                                            s05 = sum(p_level==1, na.rm=TRUE),
                                            s01 = sum(p_level==2,na.rm=TRUE),
                                            s001 = sum(p_level==3,na.rm=TRUE),
                                            snP = length(statistic))

# Figure 5F Explore HCwithin
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-10-Explore-vPFC-HC-network-clock-HConly-2.Rdata')
ddf <- ddf$coef_df_reml
explore_hcwithin <- ddf %>% filter(effect=='fixed' & term=='HCwithin:v_entropy_wi')

e_hc <- explore_hcwithin %>% mutate(p_level = case_when(padj_fdr_term <= 0.05 & padj_fdr_term > 0.01 ~ 1,
                                                        padj_fdr_term <= 0.01 & padj_fdr_term > 0.001 ~ 2,
                                                        padj_fdr_term <= 0.001 ~ 3,
                                                        padj_fdr_term > 0.05 ~ 0)) %>%
  group_by(network,HC_region) %>% summarize(mstat = mean(statistic,na.rm=TRUE), 
                                            mdf = mean(df,na.rm=TRUE),
                                            s05 = sum(p_level==1, na.rm=TRUE),
                                            s01 = sum(p_level==2,na.rm=TRUE),
                                            s001 = sum(p_level==3,na.rm=TRUE),
                                            snP = length(statistic))

# Figure 5G entropy:HCwithin
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-09-vmPFC-HC-network-clock-2.Rdata')
ddf <- ddf$coef_df_reml
mmclock_hcwithin <- ddf %>% filter(effect=='fixed' & term=='HCwithin:v_max_wi')
m_hc <- mmclock_hcwithin %>% mutate(p_level = case_when(padj_fdr_term <= 0.05 & padj_fdr_term > 0.01 ~ 1,
                                                        padj_fdr_term <= 0.01 & padj_fdr_term > 0.001 ~ 2,
                                                        padj_fdr_term <= 0.001 ~ 3,
                                                        padj_fdr_term > 0.05 ~ 0)) %>%
  group_by(network,HC_region) %>% summarize(mstat = mean(statistic,na.rm=TRUE), 
                                            mdf = mean(df,na.rm=TRUE),
                                            s05 = sum(p_level==1, na.rm=TRUE),
                                            s01 = sum(p_level==2,na.rm=TRUE),
                                            s001 = sum(p_level==3,na.rm=TRUE),
                                            snP = length(statistic))

# Figure 5H Explore HCwithin
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-10-Explore-vPFC-HC-network-clock-HConly-1.Rdata')
ddf <- ddf$coef_df_reml
explore_hcwithin <- ddf %>% filter(effect=='fixed' & term=='HCwithin:v_max_wi')

e_hc <- explore_hcwithin %>% mutate(p_level = case_when(padj_fdr_term <= 0.05 & padj_fdr_term > 0.01 ~ 1,
                                                        padj_fdr_term <= 0.01 & padj_fdr_term > 0.001 ~ 2,
                                                        padj_fdr_term <= 0.001 ~ 3,
                                                        padj_fdr_term > 0.05 ~ 0)) %>%
  summarize(mstat = mean(statistic,na.rm=TRUE), 
                                            mdf = mean(df,na.rm=TRUE),
                                            s05 = sum(p_level==1, na.rm=TRUE),
                                            s01 = sum(p_level==2,na.rm=TRUE),
                                            s001 = sum(p_level==3,na.rm=TRUE),
                                            snP = length(statistic))

