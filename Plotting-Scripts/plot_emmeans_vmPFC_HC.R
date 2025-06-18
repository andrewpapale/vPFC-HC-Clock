# plot_emmeans_vmPFC_HC

plot_emmeans_vmPFC_HC <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR){
  
  library(grid)
  
  if (strcmp(toalign,'feedback')){
    toalign_str <- 'feedback'
  } else if (strcmp(toalign,'clock')){
    toalign_str <- 'trial start'
  }
  
  pal_hc1 <- palette()
  pal_hc1[1] <- '#fc795d'
  pal_hc1[2] <- '#1f53ec'
  
  if (!strcmp(totest,'final-model-testing-')){
  tryCatch({
  pal1090 = palette()
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1090[1] <- pal[[1]]
  pal1090[2] <- '#574c35'
  
  # do v_entropy
  emt <- ddf$emmeans_list$H_HC
  emt <- emt %>% filter(network=='D' & HC_region=='PH' & (HCwithin==sort(unique(emt$HCwithin))[3] | HCwithin==sort(unique(emt$HCwithin))[4]))
  emt$levels <- factor(emt$v_entropy_wi, labels = c("-2 std Entropy","+2 std Entropy"))
  emt <- emt %>% mutate(HCwithin_label = case_when(
    HCwithin < 0 ~ '- 0.5 std HC',
    HCwithin > 0 ~ '+ 0.5 std HC'
  ))
  
  fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Entropy-by-HC','-',hc_LorR, ".pdf", sep = "")
  pdf(fname, width = 9, height = 3.5)
  gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
    facet_wrap(~HCwithin_label) +
    geom_point(aes(color=levels),size=5) +
    scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-2 std Entropy","+2 std Entropy")) + 
    geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
    geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
    ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
  
  # gg2 <- ggplot_gtable(ggplot_build(gg1))
  # stripr <- which(grepl('strip-t', gg2$layout$name))
  # k <- 1
  # for (i in stripr) {
  #   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  #   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
  #   k <- k+1
  # }
  # grid.draw(gg2)
  print(gg1)
  dev.off()
  })
  }
  if (strcmp(totest,'final-model-testing-')){
  tryCatch({
    pal1090 = palette()
    pal = wes_palette("FantasticFox1", 4, type = "discrete")
    
    # do v_entropy
    emt <- ddf$emmeans_list$H_HC
    emt <- emt %>% filter(network=='C',HC_region=='AH' & (HCwithin==sort(unique(emt$HCwithin))[3] | HCwithin==sort(unique(emt$HCwithin))[4]))
    emt <- emt %>% mutate(HCwithin_label = case_when(
      HCwithin < 0 ~ '- 0.5 std HC',
      HCwithin > 0 ~ '+ 0.5 std HC'
    ))
    
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Peaks2-by-HC','-',hc_LorR, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
      facet_wrap(~HCwithin_label) +
      geom_point(aes(color=peaks2),size=5) +
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
    
    # gg2 <- ggplot_gtable(ggplot_build(gg1))
    # stripr <- which(grepl('strip-t', gg2$layout$name))
    # k <- 1
    # for (i in stripr) {
    #   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
    #   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
    #   k <- k+1
    # }
    # grid.draw(gg2)
    print(gg1)
    dev.off()
  })
  }
  rm(emt)
  
  pal1090 = palette()
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1090[1] <- pal[[2]]
  pal1090[2] <- '#7a7745'
  
  # do v_max_wi
  emt <- ddf$emmeans_list$V_HC
  emt <- emt %>% filter(network=='C' & HC_region=='AH' & (HCwithin==sort(unique(emt$HCwithin))[3] | HCwithin==sort(unique(emt$HCwithin))[4]))
  emt$levels <- factor(emt$v_max_wi, labels = c("10'th %ile Value","90'th %ile Value"))
  emt <- emt %>% mutate(HCwithin_label = case_when(
    HCwithin < 0 ~ '- 0.5 std HC',
    HCwithin > 0 ~ '+ 0.5 std HC'
  ))
  
  fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Value-by-HC','-',hc_LorR, ".pdf", sep = "")
  pdf(fname, width = 9, height = 3.5)
  gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
    facet_wrap(~HCwithin_label) +
    geom_point(aes(color=levels),size=5) +
    scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("10th %ile Value","90'th %ile Value")) + 
    geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
    geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
    ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
  
  # gg2 <- ggplot_gtable(ggplot_build(gg1))
  # stripr <- which(grepl('strip-t', gg2$layout$name))
  # k <- 1
  # for (i in stripr) {
  #   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  #   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
  #   k <- k+1
  # }
  # grid.draw(gg2)
  print(gg1)
  dev.off()  
  
  
  rm(emt)
  
  if (strcmp(toalign,'feedback')){
    pal1090 = palette()
    pal = wes_palette("FantasticFox1", 3, type = "discrete")
    pal1090[1] <- pal[[1]]
    pal1090[2] <- '#574c35'
    
    # do v_entropy_wi_change
    emt <- ddf$emmeans_list$dH_HC
    emt <- emt %>% filter(network=='D')
    
    if (strcmp(toalign,'clock')){
      emt$levels <- factor(emt$v_entropy_wi_change_lag, labels = c("-2 std Entropy Change","+2 std Entropy Change"))
    } else if (strcmp(toalign,'feedback')){
      emt$levels <- factor(emt$v_entropy_wi_change, labels = c("-2 std Entropy Change","+2 std Entropy Change"))
    }
    
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy_Change','-',hc_LorR, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=evt_time,y=HCwithin.trend)) + 
      facet_grid(~HC_region,scales='free_y') +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-2 std Entropy Change","+2 std Entropy Change")) + 
      geom_errorbar(aes(ymin=HCwithin.trend-std.error, ymax=HCwithin.trend+std.error), width=0.5) +
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
    grid.draw(gg1)
    dev.off()
    
  } 
  
}