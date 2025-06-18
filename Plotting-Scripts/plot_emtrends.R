#plot_emtrends <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter,slrs_to_test){

library(purrr)
library(stringr)
library(pracma) 
library(stats)
    
  if (strcmp(toalign,'feedback')){
    toalign_str <- 'feedback'
  } else if (strcmp(toalign,'clock')){
    toalign_str <- 'trial start'
  }

  if (strcmp(totest,'replication')){
    totest_str <- 'MMClock out-of-session'
  } else if (strcmp(totest,'explore')){
    totest_str <- 'Explore out-of-sample'
  } else if (strcmp(totest,'MMClock')){
    totest_str <- 'MMClock fMRI'
  }
  
  
  # plot Exploration
  emtR <- ddq$emtrends_list$RTxO
  coefR <- ddq$coef_df_reml %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward') %>% select(!p.value & !statistic & !model_name & !std.error & !outcome & !df & !rhs) 
  coefO <- ddq$coef_df_reml %>% filter(effect=='fixed' & term=='rt_lag_sc:subj_level_rand_slope') %>% select(!p.value & !statistic & !model_name & !std.error & !outcome & !df & !rhs) 
  emtR <- emtR %>% filter(subj_level_rand_slope==-slrs_to_test | subj_level_rand_slope==slrs_to_test)
  emtR$subj_level_rand_slope <- as.factor(emtR$subj_level_rand_slope)
  if (strcmp(toprocess,'network-by-HC')){
    ddz <- merge(emtR,coefR,by=c('evt_time','network','HC_region'))
  } else if (strcmp(toprocess,'network')){
    ddz <- merge(emtR,coefR,by=c('evt_time','network')) %>% mutate(last_outcome = 'Reward')
  }

 ddz <- ddz  %>% 
   mutate(p_level_fdr = case_when(
     padj_fdr_term > .05 ~ '1',
     padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
     padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
     padj_fdr_term <.001 ~ '4'
   ))
 ddz$p_level_fdr <- factor(ddz$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
 ddz$p_level_fdr <- as.factor(ddz$p_level_fdr)
 
 ddz <- ddz %>% group_by(network,HC_region) %>% mutate(padj_fdr_emtrends = stats::p.adjust(p.value, method = 'fdr')) %>% ungroup()
 ddz <- ddz  %>% 
   mutate(p_level_fdr_emtrends = case_when(
     padj_fdr_emtrends > .05 ~ '1',
     padj_fdr_emtrends < .05 & padj_fdr_emtrends > .01 ~ '2',
     padj_fdr_emtrends < .01 & padj_fdr_emtrends > .001 ~ '3',
     padj_fdr_emtrends <.001 ~ '4'
   ))
 ddz$p_level_fdr_emtrends <- factor(ddz$p_level_fdr_emtrends, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
 ddz$p_level_fdr_emtrends <- as.factor(ddz$p_level_fdr_emtrends)
 ggplot(ddz %>% filter(network=='DMN' & HC_region=='AH'), aes(x=evt_time,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color=subj_level_rand_slope, group = subj_level_rand_slope)) + 
   geom_point(aes(size=p_level_fdr,alpha = p_level_fdr)) + 
   geom_line() + geom_errorbar() + facet_grid(~last_outcome) + ggtitle(paste0('DMN-',totest_str, '-coef-pfdr'))
 ggplot(ddz %>% filter(network=='CTR' & HC_region=='AH'), aes(x=evt_time,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color=subj_level_rand_slope, group = subj_level_rand_slope)) + 
   geom_point(aes(size=p_level_fdr,alpha = p_level_fdr)) + 
   geom_line() + geom_errorbar() + facet_grid(~last_outcome) + ggtitle(paste0('CTR-',totest_str, '-coef-pfdr'))
 ggplot(ddz %>% filter(network=='LIM' & HC_region=='AH'), aes(x=evt_time,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color=subj_level_rand_slope, group = subj_level_rand_slope)) + 
   geom_point(aes(size=p_level_fdr,alpha = p_level_fdr)) + 
   geom_line() + geom_errorbar() + facet_grid(~last_outcome) + ggtitle(paste0('LIM-',totest_str, '-coef-pfdr'))
  
 
 ggplot(ddz %>% filter(network=='DMN' & HC_region=='AH'), aes(x=evt_time,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color=subj_level_rand_slope, group = subj_level_rand_slope)) + 
   geom_point(aes(size=p_level_fdr_emtrends,alpha = p_level_fdr_emtrends)) + 
   geom_line() + geom_errorbar() + facet_grid(~last_outcome) + ggtitle(paste0('DMN-',totest_str, '-emtrends-pfdr'))
 ggplot(ddz %>% filter(network=='CTR' & HC_region=='AH'), aes(x=evt_time,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color=subj_level_rand_slope, group = subj_level_rand_slope)) + 
   geom_point(aes(size=p_level_fdr_emtrends,alpha = p_level_fdr_emtrends)) + 
   geom_line() + geom_errorbar() + facet_grid(~last_outcome) + ggtitle(paste0('CTR-',totest_str, '-emtrends-pfdr'))
 ggplot(ddz %>% filter(network=='LIM' & HC_region=='AH'), aes(x=evt_time,y=rt_lag_sc.trend,ymin=rt_lag_sc.trend-std.error,ymax=rt_lag_sc.trend+std.error,color=subj_level_rand_slope, group = subj_level_rand_slope)) + 
   geom_point(aes(size=p_level_fdr_emtrends,alpha = p_level_fdr_emtrends)) + 
   geom_line() + geom_errorbar() + facet_grid(~last_outcome) + ggtitle(paste0('LIM-',totest_str, '-emtrends-pfdr'))
#}

