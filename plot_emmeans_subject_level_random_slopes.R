# 2022-05-10 AndyP
# Plot emmeans_vmPFC_HC.R
# Plot emmeans from model

plot_emmeans_subject_level_random_slopes <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR){
  
  
  if (strcmp(toalign,'feedback')){
    toalign_str <- 'feedback'
  } else if (strcmp(toalign,'clock')){
    toalign_str <- 'trial start'
  }
  
  pal_hc1 <- palette()
  pal_hc1[1] <- '#fc795d'
  pal_hc1[2] <- '#1f53ec'
  
  # pal1090 = palette()
  # pal1090[1] <- '#3863e5'
  # pal1090[2] <- '#f31836'
  
  
  pal1090 = palette()
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1090[1] <- pal[[2]]
  pal1090[2] <- '#7a7745'
  
  
  df <- ddf$coef_df_reml
  
  
  if (strcmp(toprocess,'network-by-HC')){
    
    emt <- ddq$emmeans_list$TrxVmax
    emt <- emt %>% filter(subj_level_rand_slope==min(unique(subj_level_rand_slope)) | subj_level_rand_slope==max(unique(subj_level_rand_slope)))
    emt <- emt %>% mutate(trial_bin=case_when(trial_neg_inv_sc < 0 ~ 'Early',trial_neg_inv_sc > 0 ~ 'Late'))
    #emt <- emt %>% filter(trial_bin != 'Middle')
    emt$levels <- factor(emt$subj_level_rand_slope, labels = c("10'th %ile HC slope","90'th %ile HC slope"))
    emt$levels <- relevel(emt$levels,ref=c("90'th %ile HC slope"))
    
    df0 <- df %>% filter(term=='trial_neg_inv_sc:subj_level_rand_slope:rt_vmax_lag_sc')
    
    Q <- inner_join(emt,df0,by=c('evt_time','network','HC_region'))
    Q <- Q %>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
    Q <- Q  %>% group_by(network1) %>% mutate(padj_BY_term = p.adjust(p.value.y, method = 'bonferroni')) %>% ungroup() %>% 
      mutate(p_level_fdr = as.factor(case_when(
        padj_BY_term > .05 ~ '1',
        padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
        padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
        padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
        padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
        padj_BY_term <.00001 ~ '6'
      )))
    
    
    
    fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'rt_vmax_lag_by_trial','-',hc_LorR, ".pdf", sep = "")
    pdf(fname, width = 9, height = 9)
    gg1 <- ggplot(Q,aes(x=evt_time,y=estimate)) + 
      facet_grid(network~HC_region) +
      geom_point(aes(color=trial_bin,size=as.factor(p_level_fdr),alpha=as.factor(p_level_fdr))) +
      geom_line(aes(color=trial_bin,linetype=as.factor(levels)), size=1) + 
      geom_errorbar(aes(ymin=estimate-std.error.x, ymax=estimate+std.error.x), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
    print(gg1)
    dev.off()
    
    
    emt <- ddq$emmeans_list$LO
    emt <- emt %>% filter(subj_level_rand_slope==min(unique(subj_level_rand_slope)) | subj_level_rand_slope==max(unique(subj_level_rand_slope)))
    emt$levels <- factor(emt$subj_level_rand_slope, labels = c("10'th %ile HC slope","90'th %ile HC slope"))
    emt$levels <- relevel(emt$levels,ref=c("90'th %ile HC slope"))
    
    df0 <- df %>% filter(term=='rt_lag_sc:subj_level_rand_slope' | term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward')
    df0 <- df0 %>% mutate(last_outcome = case_when(term=='rt_lag_sc:subj_level_rand_slope' ~ 'Omission', 
                                                   term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward' ~ 'Reward'))
    
    Q <- inner_join(emt,df0,by=c('evt_time','network','last_outcome','HC_region'))
    Q <- Q %>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
    Q <- Q  %>% group_by(network1) %>% mutate(padj_BY_term = p.adjust(p.value.y, method = 'bonferroni')) %>% ungroup() %>% 
      mutate(p_level_fdr = as.factor(case_when(
        padj_BY_term > .05 ~ '1',
        padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
        padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
        padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
        padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
        padj_BY_term <.00001 ~ '6'
      )))
    
    fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'last_outcome','-',hc_LorR, ".pdf", sep = "")
    pdf(fname, width = 9, height = 9)
    gg1 <- ggplot(Q,aes(x=evt_time,y=estimate)) + 
      facet_grid(network~HC_region) +
      geom_point(aes(color=last_outcome,size=as.factor(p_level_fdr),alpha=as.factor(p_level_fdr))) +
      geom_line(aes(color=last_outcome,linetype=as.factor(levels)), size=1) + 
      geom_errorbar(aes(ymin=estimate-std.error.x, ymax=estimate+std.error.x), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
    
    gg2 <- ggplot_gtable(ggplot_build(gg1))
    stripr <- which(grepl('strip-t', gg2$layout$name))
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
      gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
      k <- k+1
    }
    grid.draw(gg2)
    dev.off()  
    
    
    emt <- ddq$emmeans_list$Vmax
    emt <- emt %>% filter(subj_level_rand_slope==min(unique(subj_level_rand_slope)) | subj_level_rand_slope==max(unique(subj_level_rand_slope)))
    emt$levels <- factor(emt$subj_level_rand_slope, labels = c("10'th %ile HC slope","90'th %ile HC slope"))
    emt$levels <- relevel(emt$levels,ref=c("90'th %ile HC slope"))
    df0 <- df %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
    
    Q <- inner_join(emt,df0,by=c('evt_time','network','HC_region'))
    Q <- Q %>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
    Q <- Q  %>% group_by(network1) %>% mutate(padj_BY_term = p.adjust(p.value.y, method = 'bonferroni')) %>% ungroup() %>% 
      mutate(p_level_fdr = as.factor(case_when(
        padj_BY_term > .05 ~ '1',
        padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
        padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
        padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
        padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
        padj_BY_term <.00001 ~ '6'
      )))
    
    
    #emt <- emt %>% mutate(trial_bin=case_when(trial_neg_inv_sc < 0 ~ 'Early',trial_neg_inv_sc > 0 ~ 'Late'))
    #emt <- emt %>% filter(trial_bin != 'Middle')
    
    fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'rt_vmax_lag','-',hc_LorR, ".pdf", sep = "")
    pdf(fname, width = 9, height = 9)
    gg1 <- ggplot(Q,aes(x=evt_time,y=estimate)) + 
      facet_grid(network1~HC_region) +
      geom_point(aes(color=as.factor(levels),size=as.factor(p_level_fdr),alpha=as.factor(p_level_fdr))) +
      geom_line(aes(color=as.factor(levels),linetype=as.factor(levels)), size=1) + 
      geom_errorbar(aes(ymin=estimate-std.error.x, ymax=estimate+std.error.x), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
    print(gg1)
    dev.off()
    
  } else if (strcmp(toprocess,'network')){
    
    emt <- ddq$emmeans_list$TrxVmax
    emt <- emt %>% filter(subj_level_rand_slope==min(unique(subj_level_rand_slope)) | subj_level_rand_slope==max(unique(subj_level_rand_slope)))
    emt <- emt %>% mutate(trial_bin=case_when(trial_neg_inv_sc < 0 ~ 'Early',trial_neg_inv_sc > 0 ~ 'Late'))
    #emt <- emt %>% filter(trial_bin != 'Middle')
    emt$levels <- factor(emt$subj_level_rand_slope, labels = c("10'th %ile slope","90'th %ile slope"))
    emt$levels <- relevel(emt$levels,ref=c("90'th %ile slope"))
    
    df0 <- df %>% filter(term=='trial_neg_inv_sc:subj_level_rand_slope:rt_vmax_lag_sc')
    
    Q <- inner_join(emt,df0,by=c('evt_time','network'))
    Q <- Q %>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
    Q <- Q  %>% group_by(network1) %>% mutate(padj_BY_term = p.adjust(p.value.y, method = 'bonferroni')) %>% ungroup() %>% 
      mutate(p_level_fdr = as.factor(case_when(
        padj_BY_term > .05 ~ '1',
        padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
        padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
        padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
        padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
        padj_BY_term <.00001 ~ '6'
      )))
    
    
    
    fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'rt_vmax_lag_sc_by_trial','-',hc_LorR, ".pdf", sep = "")
    pdf(fname, width = 9, height = 9)
    gg1 <- ggplot(Q,aes(x=evt_time,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=trial_bin,size=as.factor(p_level_fdr),alpha=as.factor(p_level_fdr))) +
      geom_line(aes(color=trial_bin,linetype=as.factor(levels)), size=1) + 
      geom_errorbar(aes(ymin=estimate-std.error.x, ymax=estimate+std.error.x), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
    print(gg1)
    dev.off()
    
    
    emt <- ddq$emmeans_list$LO
    emt <- emt %>% filter(subj_level_rand_slope==min(unique(subj_level_rand_slope)) | subj_level_rand_slope==max(unique(subj_level_rand_slope)))
    emt$levels <- factor(emt$subj_level_rand_slope, labels = c("10'th %ile slope","90'th %ile slope"))
    emt$levels <- relevel(emt$levels,ref=c("90'th %ile slope"))
    
    df0 <- df %>% filter(term=='rt_lag_sc:subj_level_rand_slope' | term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward')
    df0 <- df0 %>% mutate(last_outcome = case_when(term=='rt_lag_sc:subj_level_rand_slope' ~ 'Omission', 
                                                   term=='rt_lag_sc:subj_level_rand_slope:last_outcomeReward' ~ 'Reward'))
    
    Q <- inner_join(emt,df0,by=c('evt_time','network','last_outcome'))
    Q <- Q %>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
    Q <- Q  %>% group_by(network1) %>% mutate(padj_BY_term = p.adjust(p.value.y, method = 'bonferroni')) %>% ungroup() %>% 
      mutate(p_level_fdr = as.factor(case_when(
        padj_BY_term > .05 ~ '1',
        padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
        padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
        padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
        padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
        padj_BY_term <.00001 ~ '6'
      )))
    
    fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'last_outcome','-',hc_LorR, ".pdf", sep = "")
    pdf(fname, width = 9, height = 9)
    gg1 <- ggplot(Q,aes(x=evt_time,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=last_outcome,size=as.factor(p_level_fdr),alpha=as.factor(p_level_fdr))) +
      geom_line(aes(color=last_outcome,linetype=as.factor(levels)), size=1) + 
      geom_errorbar(aes(ymin=estimate-std.error.x, ymax=estimate+std.error.x), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
    
    gg2 <- ggplot_gtable(ggplot_build(gg1))
    stripr <- which(grepl('strip-t', gg2$layout$name))
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
      gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
      k <- k+1
    }
    grid.draw(gg2)
    dev.off()  
    
    
    emt <- ddq$emmeans_list$Vmax
    emt <- emt %>% filter(subj_level_rand_slope==min(unique(subj_level_rand_slope)) | subj_level_rand_slope==max(unique(subj_level_rand_slope)))
    emt$levels <- factor(emt$subj_level_rand_slope, labels = c("10'th %ile slope","90'th %ile slope"))
    emt$levels <- relevel(emt$levels,ref=c("90'th %ile slope"))
    df0 <- df %>% filter(term=='subj_level_rand_slope:rt_vmax_lag_sc')
    
    Q <- inner_join(emt,df0,by=c('evt_time','network'))
    Q <- Q %>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
    Q <- Q  %>% group_by(network1) %>% mutate(padj_BY_term = p.adjust(p.value.y, method = 'bonferroni')) %>% ungroup() %>% 
      mutate(p_level_fdr = as.factor(case_when(
        padj_BY_term > .05 ~ '1',
        padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
        padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
        padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
        padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
        padj_BY_term <.00001 ~ '6'
      )))
    
    
    #emt <- emt %>% mutate(trial_bin=case_when(trial_neg_inv_sc < 0 ~ 'Early',trial_neg_inv_sc > 0 ~ 'Late'))
    #emt <- emt %>% filter(trial_bin != 'Middle')
    
    fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'rt_vmax_lag_sc','-',hc_LorR, ".pdf", sep = "")
    pdf(fname, width = 9, height = 9)
    gg1 <- ggplot(Q,aes(x=evt_time,y=estimate)) + 
      facet_grid(~network1) +
      geom_point(aes(color=as.factor(levels),size=as.factor(p_level_fdr),alpha=as.factor(p_level_fdr))) +
      geom_line(aes(color=as.factor(levels),linetype=as.factor(levels)), size=1) + 
      geom_errorbar(aes(ymin=estimate-std.error.x, ymax=estimate+std.error.x), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
    print(gg1)
    dev.off()
  } 
  
}  



#   # do v_max_wi
#   
#   emt <- ddf$emtrends_list$V
#   #emt <- emt %>% filter(network=='D')
#   emt$levels <- factor(emt$v_max_wi_lag, labels = c("10'th %ile Value","90'th %ile Value"))
#   
#   fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Value','-',hc_LorR, ".pdf", sep = "")
#   pdf(fname, width = 9, height = 9)
#   gg1 <- ggplot(emt,aes(x=evt_time,y=estimate.trend)) + 
#     facet_grid(network~HC_region) +
#     geom_point(aes(color=levels),size=5) +
#     scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("10th %ile Value","90'th %ile Value")) + 
#     geom_errorbar(aes(ymin=estimate.trend-std.error, ymax=estimate.trend+std.error), width=0.5) +
#     geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#     ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   
#   gg2 <- ggplot_gtable(ggplot_build(gg1))
#   stripr <- which(grepl('strip-t', gg2$layout$name))
#   k <- 1
#   for (i in stripr) {
#     j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#     gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#     k <- k+1
#   }
#   grid.draw(gg2)
#   dev.off()  
#   
#   
#   
#   # rm(emt)
#   # # do trial_neg_inv_sc
#   # emt <- ddf$emtrends_list$Tr
#   # #emt <- emt %>% filter(network=='D')
#   # emt$levels <- factor(emt$trial_neg_inv_sc, labels = c("10'th %ile Trial","90'th %ile Trial"))
#   # 
#   # fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Trial','-',hc_LorR, ".pdf", sep = "")
#   # pdf(fname, width = 9, height = 9)
#   # gg1 <- ggplot(emt,aes(x=evt_time,y=estimate.trend)) + 
#   #   facet_grid(network~HC_region) +
#   #   geom_point(aes(color=levels),size=5) +
#   #   scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("10th %ile Trial","90'th %ile Trial")) + 
#   #   geom_errorbar(aes(ymin=estimate.trend-std.error, ymax=estimate.trend+std.error), width=0.5) +
#   #   geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#   #   ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   # 
#   # gg2 <- ggplot_gtable(ggplot_build(gg1))
#   # stripr <- which(grepl('strip-t', gg2$layout$name))
#   # k <- 1
#   # for (i in stripr) {
#   #   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#   #   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#   #   k <- k+1
#   # }
#   # grid.draw(gg2)
#   # dev.off()  
# 
#   
#   rm(emt)
#   # do last_outcome
#   emt <- ddf$emtrends_list$LO
#   #emt <- emt %>% filter(network=='L')
#   emt$levels <- factor(emt$last_outcome)
#   
#   fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'last_outcome','-',hc_LorR, ".pdf", sep = "")
#   pdf(fname, width = 9, height = 9)
#   gg1 <- ggplot(emt,aes(x=evt_time,y=estimate.trend)) + 
#     facet_grid(network~HC_region) +
#     geom_point(aes(color=levels),size=5) +
#     geom_errorbar(aes(ymin=estimate.trend-std.error, ymax=estimate.trend+std.error), width=0.5) +
#     geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#     ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   
#   gg2 <- ggplot_gtable(ggplot_build(gg1))
#   stripr <- which(grepl('strip-t', gg2$layout$name))
#   k <- 1
#   for (i in stripr) {
#     j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#     gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#     k <- k+1
#   }
#   grid.draw(gg2)
#   dev.off()  
#   
#   rm(emt)
#   
#   # do v_entropy_wi
#   emt <- ddf$emtrends_list$H
#   #emt <- emt %>% filter(network=='L')
#   emt$levels <- factor(emt$v_entropy_wi, labels = c("-2 std entropy","+ 2 std entropy"))
#   
#   fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'v_entropy_wi','-',hc_LorR, ".pdf", sep = "")
#   pdf(fname, width = 9, height = 9)
#   gg1 <- ggplot(emt,aes(x=evt_time,y=estimate.trend)) + 
#     facet_grid(network~HC_region) +
#     geom_point(aes(color=levels),size=5) +
#     scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-2 std entropy","+ 2 std entropy")) + 
#     geom_errorbar(aes(ymin=estimate.trend-std.error, ymax=estimate.trend+std.error), width=0.5) +
#     geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#     ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   
#   gg2 <- ggplot_gtable(ggplot_build(gg1))
#   stripr <- which(grepl('strip-t', gg2$layout$name))
#   k <- 1
#   for (i in stripr) {
#     j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#     gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#     k <- k+1
#   }
#   grid.draw(gg2)
#   dev.off()   
#   
# 
#   
#   
#   # do rt_lag_sc
#   emt <- ddf$emtrends_list$RT
#   #emt <- emt %>% filter(network=='L')
#   emt$levels <- factor(emt$rt_lag_sc, labels = c("10%","25%","50%","75%","90%"))
#   
#   fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'rt_lag_sc','-',hc_LorR, ".pdf", sep = "")
#   pdf(fname, width = 9, height = 9)
#   gg1 <- ggplot(emt,aes(x=evt_time,y=estimate.trend)) + 
#     facet_grid(network~HC_region) +
#     geom_point(aes(color=levels),size=5) +
#     geom_errorbar(aes(ymin=estimate.trend-std.error, ymax=estimate.trend+std.error), width=0.5) +
#     geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#     ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   
#   gg2 <- ggplot_gtable(ggplot_build(gg1))
#   stripr <- which(grepl('strip-t', gg2$layout$name))
#   k <- 1
#   for (i in stripr) {
#     j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#     gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#     k <- k+1
#   }
#   grid.draw(gg2)
#   dev.off()  
#   
#   
#   
#   # do rt_lag_sc * last_outcome
#   emt <- ddf$emtrends_list$RTxLO
#   #emt <- emt %>% filter(network=='L')
#   emt$levels <- factor(emt$rt_lag_sc, labels = c("10%","25%","50%","75%","90%"))
#   fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'rt_lag_sc_by_last_outcome','-',hc_LorR, ".pdf", sep = "")
#   pdf(fname, width = 9, height = 9)
#   gg1 <- ggplot(emt,aes(x=evt_time,y=estimate.trend)) + 
#     facet_grid(network+last_outcome~HC_region) +
#     geom_point(aes(color=levels),size=2.5) +
#     geom_errorbar(aes(ymin=estimate.trend-std.error, ymax=estimate.trend+std.error), width=0.5) +
#     geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#     ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   
#   gg2 <- ggplot_gtable(ggplot_build(gg1))
#   stripr <- which(grepl('strip-t', gg2$layout$name))
#   k <- 1
#   for (i in stripr) {
#     j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#     gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#     k <- k+1
#   }
#   grid.draw(gg2)
#   dev.off() 
# 
#   
#   # do rt_vmax_lag * trial_neg_inv_sc
#   emt <- ddf$emtrends_list$RTvxTR
#   #emt <- emt %>% filter(network=='L')
#   emt$levels <- factor(emt$rt_vmax_lag, labels = c("10%","90%"))
#   
#   fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'rt_vmax_lag_by_trial_bin','-',hc_LorR, ".pdf", sep = "")
#   pdf(fname, width = 9, height = 9)
#   gg1 <- ggplot(emt,aes(x=evt_time,y=estimate.trend)) + 
#     facet_grid(network+rewFunc~HC_region) +
#     geom_point(aes(color=levels),size=2.5) +
#     geom_errorbar(aes(ymin=estimate.trend-std.error, ymax=estimate.trend+std.error), width=0.5) +
#     geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#     ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   
#   gg2 <- ggplot_gtable(ggplot_build(gg1))
#   stripr <- which(grepl('strip-t', gg2$layout$name))
#   k <- 1
#   for (i in stripr) {
#     j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#     gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#     k <- k+1
#   }
#   grid.draw(gg2)
#   dev.off()   
#   
#   
#   
# # do RT_RS  
#   emt <- ddf$emtrends_list$RT_RS
#   #emt <- emt %>% filter(network=='L')
#   emt <- emt %>% filter(estimate==min(unique(estimate)) | estimate==max(unique(estimate)))
#   emt$levels <- factor(emt$estimate, labels = c("10%","90%"))
#   
#   fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'rt_csv_sc-by-estimate','-',hc_LorR, ".pdf", sep = "")
#   pdf(fname, width = 9, height = 9)
#   gg1 <- ggplot(emt,aes(x=evt_time,y=rt_lag_sc.trend)) + 
#     facet_grid(network~HC_region) +
#     geom_point(aes(color=levels),size=2.5) +
#     geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error, ymax=rt_lag_sc.trend+std.error), width=0.5) +
#     geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#     ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   
#   gg2 <- ggplot_gtable(ggplot_build(gg1))
#   stripr <- which(grepl('strip-t', gg2$layout$name))
#   k <- 1
#   for (i in stripr) {
#     j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#     gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#     k <- k+1
#   }
#   grid.draw(gg2)
#   dev.off()     
# 
#   
#   # do RTxLO_RS 
#   emt <- ddf$emtrends_list$RTxLO_RS
#   emt <- emt %>% filter(estimate==min(unique(estimate)) | estimate==max(unique(estimate)))
#   #emt <- emt %>% filter(network=='L')
#   emt$levels <- factor(emt$estimate, labels = c("10%","90%"))
#   
#   fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'rt_csv_sc-by-outcome-by-estimate','-',hc_LorR, ".pdf", sep = "")
#   pdf(fname, width = 9, height = 9)
#   gg1 <- ggplot(emt,aes(x=evt_time,y=rt_lag_sc.trend)) + 
#     facet_grid(network~HC_region) +
#     geom_point(aes(color=levels),shape=as.factor(emt$last_outcome),size=2.5) +
#     geom_line(aes(color=levels,linetype=as.factor(last_outcome)), size=1) + 
#     geom_errorbar(aes(ymin=rt_lag_sc.trend-std.error, ymax=rt_lag_sc.trend+std.error), width=0.5) +
#     geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#     ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   
#   gg2 <- ggplot_gtable(ggplot_build(gg1))
#   stripr <- which(grepl('strip-t', gg2$layout$name))
#   k <- 1
#   for (i in stripr) {
#     j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#     gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#     k <- k+1
#   }
#   grid.draw(gg2)
#   dev.off()    
# 
#   
#   # do HxTrxRT 
#   emt <- ddf$emtrends_list$HxTrxRT
#   emt <- emt %>% filter(network=='D')
#   emt <- emt %>% filter(rt_lag_sc == min(unique(rt_lag_sc)) | rt_lag_sc == max(unique(rt_lag_sc)))
#   emt$levels <- factor(emt$rt_lag_sc, labels = c("10% rt_lag_sc","90% rt_lag_sc"))
#   #emt$trial_bin <- factor(emt$trial_bin, levels=c('Early','Middle','Late'))
#   #emt <- emt %>% filter(trial_bin != 'Middle')
#   emt <- emt %>% filter(HC_region != 'AH')
#   
#   fname = paste('randomslopes','-',behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'rt_csv_sc-by-entropy-by-trial','-',hc_LorR, ".pdf", sep = "")
#   pdf(fname, width = 9, height = 9)
#   gg1 <- ggplot(emt,aes(x=evt_time,y=estimate.trend)) + 
#     facet_grid(levels~rewFunc) +
#     geom_point(shape=as.factor(emt$v_entropy_wi),size=2.5) +
#     geom_line(aes(linetype=as.factor(v_entropy_wi)), size=1) + 
#     geom_errorbar(aes(ymin=estimate.trend-std.error, ymax=estimate.trend+std.error), width=0.5) +
#     geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
#     ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
#   print(gg1)
#   dev.off()
#   # gg2 <- ggplot_gtable(ggplot_build(gg1))
#   # stripr <- which(grepl('strip-t', gg2$layout$name))
#   # k <- 1
#   # for (i in stripr) {
#   #   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#   #   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
#   #   k <- k+1
#   # }
#   # grid.draw(gg2)
#   # dev.off()      
#         
# }
# 
# # gg<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network1)) + 
# #   geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
# #   # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
# #   geom_line(size = 1) + theme(legend.position = "none") +
# #   geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
# #   scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
# #   #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
# #   theme_bw(base_size=13) +
# #   facet_wrap(~HC_region) +
# #   theme(legend.title = element_blank(),
# #         panel.grid.major = element_line(colour = "grey45"), 
# #         panel.grid.minor = element_line(colour = "grey45"), 
# #         panel.background = element_rect(fill = 'grey40'),
# #         axis.title.y = element_text(margin=margin(r=6)),
# #         axis.title.x = element_text(margin=margin(t=6)))
# # }
# # print(gg)
# # dev.off()
# 
# # fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
# # pdf(fname, width = 9, height = 3.5)
# # gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=as.factor(HC_region), size=`p, FDR-corrected`)) +
# #   geom_line(size = 1) + geom_point() +
# #   geom_errorbar() +
# #   geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) + facet_grid(~symmetry_group1) + 
# #   scale_color_manual(values = pal_hc1,labels=c('AH','PH')) + xlab(epoch_label) +
# #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
# # 
# # gg2 <- ggplot_gtable(ggplot_build(gg1))
# # stripr <- which(grepl('strip-t', gg2$layout$name))
# # k <- 1
# # for (i in stripr) {
# #   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
# #   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
# #   k <- k+1
# # }
# # grid.draw(gg2)
# # #print(gg2)
# # dev.off()
