plot_emtrends_vmPFC_HC_entropy_split <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR){
  
  
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
  
  # # do PE
  # emt <- ddf$emtrends_list$PE_HC
  # emt$levels <- factor(emt$pe_max_sc, labels = c("10'th %ile PE","90'th %ile PE"))
  # emt <- emt %>% group_by(evt_time,HC_region,network) %>%
  #   mutate(ds = abs(HCwithin.trend[levels=="10'th %ile PE"]-HCwithin.trend[levels=="90'th %ile PE"])) %>%
  #   ungroup() %>%
  #   mutate(diff_score = case_when(
  #     ds >= mean(ds) + 2.5*std(ds) ~ 2,
  #     ds > mean(ds) + 0.5*std(ds) & ds < mean(ds) + 2.5*std(ds)  ~ 1,
  #     ds < mean(ds) + 0.5*std(ds) ~ 0
  #   ))
  # emt$diff_score <- as.factor(emt$diff_score)
  # 
  # fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'PE','-',hc_LorR, ".pdf", sep = "")
  # pdf(fname, width = 9, height = 3.5)
  # gg1 <- ggplot(emt,aes(x=evt_time,y=HCwithin.trend)) + 
  #   facet_grid(network~HC_region) +
  #   geom_point(aes(size=2,alpha=diff_score,color=levels)) +
  #   scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=levels) + 
  #   geom_errorbar(aes(ymin=HCwithin.trend-std.error, ymax=HCwithin.trend+std.error), width=0.5) +
  #   geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1)
  #   
  # gg2 <- ggplot_gtable(ggplot_build(gg1))
  # stripr <- which(grepl('strip-t', gg2$layout$name))
  # k <- 1
  # for (i in stripr) {
  #   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  #   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
  #   k <- k+1
  # }
  # grid.draw(gg2)
  # dev.off()
  # 
  # 
  # rm(emt)
  
  
  pal1090 = palette()
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1090[1] <- pal[[1]]
  pal1090[2] <- '#574c35'
  
  # do v_entropy
  emt <- ddf$emtrends_list$H_HC
  emt <- emt %>% filter(network=='D')
  emt$levels <- factor(emt$v_entropy_wi, labels = c("-2 std Entropy","+2 std Entropy"))
  
  fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy','-',hc_LorR, ".pdf", sep = "")
  pdf(fname, width = 15, height = 6)
  gg1 <- ggplot(emt,aes(x=evt_time,y=HCwithin.trend)) + 
    facet_wrap(entropy_split~HC_region) +
    geom_point(aes(color=levels),size=5) +
    scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-2 std Entropy","+2 std Entropy")) + 
    geom_errorbar(aes(ymin=HCwithin.trend-std.error, ymax=HCwithin.trend+std.error), width=0.5) +
    geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
    ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
  
  print(gg1)
  dev.off()
  
  
  rm(emt)
  
  pal1090 = palette()
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1090[1] <- pal[[2]]
  pal1090[2] <- '#7a7745'
  
  # do v_max_wi
  emt <- ddf$emtrends_list$V_HC
  emt <- emt %>% filter(network=='C')
  emt$levels <- factor(emt$v_max_wi, labels = c("10'th %ile Value","90'th %ile Value"))
  
  fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Value','-',hc_LorR, ".pdf", sep = "")
  pdf(fname, width = 15, height = 6)
  gg1 <- ggplot(emt,aes(x=evt_time,y=HCwithin.trend)) + 
    facet_wrap(entropy_split~HC_region) +
    geom_point(aes(color=levels),size=5) +
    scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("10th %ile Value","90'th %ile Value")) + 
    geom_errorbar(aes(ymin=HCwithin.trend-std.error, ymax=HCwithin.trend+std.error), width=0.5) +
    geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
    ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
  print(gg1)
  dev.off()
  
  
  
}

