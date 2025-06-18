# 2023-06-02 AndyP
# Plot b2b emtrends for vPFC Entropy and Vmax random slopes

library(tidyverse)

# load MMClock entropy random slope
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-19-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-both-1.Rdata')


explore <- ddq$coef_df_reml %>% filter(effect=='fixed' & (term=='rt_lag_sc:subj_level_rand_slope' | term=='rt_lag_sc:subj_level_rand_slope:last_outcomeOmission')) # last outcome rewarded/not
exploit <- ddq$coef_df_reml %>% filter(effect=='fixed' & (term=='subj_level_rand_slope:rt_vmax_lag_sc'))

explore <- explore %>% mutate(strategy = case_when(term=='rt_lag_sc:subj_level_rand_slope' ~ 'win-shift',
                                                   term=='rt_lag_sc:subj_level_rand_slope:last_outcomeOmission' ~ 'lose-shift'),
                              last_outcome = case_when(strategy=='win-shift' ~ 'Reward',
                                                       strategy=='lose-shift'~ 'Omission'),
                              p_level_fdr = case_when(
                                padj_fdr_term > .05 ~ 'NS',
                                padj_fdr_term <= .05 & padj_fdr_term > .01 ~ 'p < 0.05',
                                padj_fdr_term <= .01 & padj_fdr_term > .001 ~ 'p < 0.01',
                                padj_fdr_term <=.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
                                padj_fdr_term <=.0001 ~ 'p < 0.0001'
                              ))
exploit <- exploit %>% mutate(p_level_fdr = case_when(
                                padj_fdr_term > .05 ~ 'NS',
                                padj_fdr_term <= .05 & padj_fdr_term > .01 ~ 'p < 0.05',
                                padj_fdr_term <= .01 & padj_fdr_term > .001 ~ 'p < 0.01',
                                padj_fdr_term <=.001 & padj_fdr_term > .0001 ~ 'p < 0.001',
                                padj_fdr_term <=.0001 ~ 'p < 0.0001'
                              ))

explore <- explore %>% select(!std.error)
exploit <- exploit %>% select(!std.error)
exploit$p_level_fdr <- factor(exploit$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
explore$p_level_fdr <- factor(explore$p_level_fdr, levels=c('NS','p < 0.05','p < 0.01','p < 0.001','p < 0.0001'))
emt_explore <- ddq$emtrends_list$RTxO
explore1 <- inner_join(explore,emt_explore, by=c('evt_time','last_outcome','network'))
explore1$subj_level_rand_slope <- as.factor(explore1$subj_level_rand_slope)


ggplot(explore1, aes(x=evt_time,y=rt_lag_sc.trend,color=subj_level_rand_slope,group=subj_level_rand_slope)) + ggh4x::facet_grid2(network~strategy,scales='free',independent='y') +
  geom_point(aes(group=evt_time,size=p_level_fdr,alpha = p_level_fdr)) +
  geom_line(aes(color=subj_level_rand_slope,group=subj_level_rand_slope)) + geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error),width=0,color='black')


explore1 <- explore1 %>% filter(evt_time==0)
ggplot(explore1, aes(x=subj_level_rand_slope,y=rt_lag_sc.trend,group=subj_level_rand_slope)) + geom_point() + facet_grid(network~strategy)

emm_explore = ddq$emmeans_list$RTxO
explore <- explore %>% select(!estimate)
explore2 <- inner_join(explore,emm_explore,by=c('evt_time','last_outcome','network'))
explore2 <- explore2 %>% filter(strategy=='win-shift' & (subj_level_rand_slope==-2 | subj_level_rand_slope==2) & rt_lag_sc==0 & network!='L')
explore2 <- explore2 %>% select(evt_time,strategy,subj_level_rand_slope,network,std.error,p_level_fdr,rt_lag_sc,estimate)
explore2$subj_level_rand_slope <- relevel(as.factor(explore2$subj_level_rand_slope),ref='2')

explore2 <- explore2 %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR',network=='L'~'LIM'),
                                level_of_random_slope = case_when(subj_level_rand_slope=='2' ~ 'high entropy',
                                                                  subj_level_rand_slope=='-2' ~ 'low entropy'))

library(wesanderson)
library(grid)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]
pal1 = palette()
pal1[2] = "#CCCCCC"
pal1[1] = "#333333"
gg1 <- ggplot(explore2, aes(x=evt_time,y=estimate,group=level_of_random_slope,color=level_of_random_slope)) + facet_grid(~network1) + 
  geom_point(position=position_dodge(width=0.5),aes(group=level_of_random_slope,size=p_level_fdr)) + scale_color_manual(values=pal1) +
  geom_line(aes(group=level_of_random_slope,color=level_of_random_slope)) + geom_errorbar(position=position_dodge(width=0.5),aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.1,color='black') +
  ylab('mean RT swing (estimated)') + xlab('time relative to trial onset [s]')
gg2 <- ggplot_gtable(ggplot_build(gg1))
stripr <- which(grepl('strip-t', gg2$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal[k]
  k <- k+1
}
grid.draw(gg2)


emt_exploit <- ddq$emtrends_list$Vmax
exploit1 <- inner_join(exploit,emt_exploit,by=c('evt_time','network'))
exploit1$subj_level_rand_slope <- as.factor(exploit1$subj_level_rand_slope)
ggplot(exploit1, aes(x=evt_time,y=rt_vmax_lag_sc.trend,color=subj_level_rand_slope,group=subj_level_rand_slope)) + facet_grid(~network) +
  geom_point(aes(group=evt_time,size=p_level_fdr,alpha = p_level_fdr)) +
  geom_line(aes(color=subj_level_rand_slope,group=subj_level_rand_slope)) + geom_errorbar(aes(ymin=rt_vmax_lag_sc.trend-std.error,ymax=rt_vmax_lag_sc.trend+std.error),width=0,color='black')
exploit1 <- exploit1 %>% filter(evt_time==0)
ggplot(exploit1, aes(x=subj_level_rand_slope,y=rt_vmax_lag_sc.trend,group=subj_level_rand_slope)) + geom_point() + facet_grid(~network)

emm_exploit = ddq$emmeans_list$Vmax
exploit <- exploit %>% select(!estimate)
exploit2 <- inner_join(exploit,emm_exploit,by=c('evt_time','network'))
exploit2 <- exploit2 %>% filter((subj_level_rand_slope==-1 | subj_level_rand_slope==1) & rt_vmax_lag_sc==0 & network!='L')
exploit2 <- exploit2 %>% select(evt_time,subj_level_rand_slope,network,std.error,p_level_fdr,rt_vmax_lag_sc,estimate)
exploit2$subj_level_rand_slope <- relevel(as.factor(exploit2$subj_level_rand_slope),ref='1')

exploit2 <- exploit2 %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR',network=='L'~'LIM'),
                                level_of_random_slope = case_when(subj_level_rand_slope=='1' ~ 'high entropy',
                                                                  subj_level_rand_slope=='-1' ~ 'low entropy'))

library(wesanderson)
library(grid)
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]
pal1 = palette()
pal1[2] = "#CCCCCC"
pal1[1] = "#333333"
gg1 <- ggplot(exploit2, aes(x=evt_time,y=estimate,group=level_of_random_slope,color=level_of_random_slope)) + facet_grid(~network1) + 
  geom_point(position=position_dodge(width=0.5),aes(group=level_of_random_slope,size=p_level_fdr)) + scale_color_manual(values=pal1) +
  geom_line(aes(group=level_of_random_slope,color=level_of_random_slope)) + geom_errorbar(position=position_dodge(width=0.5),aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.1,color='black') +
  ylab('mean convergence on best RT (estimated)') + xlab('time relative to trial onset [s]')
gg2 <- ggplot_gtable(ggplot_build(gg1))
stripr <- which(grepl('strip-t', gg2$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal[k]
  k <- k+1
}
grid.draw(gg2)


load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-19-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-both-2.Rdata')
