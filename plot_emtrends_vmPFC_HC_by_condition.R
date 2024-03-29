# 2022-05-10 AndyP
# Plot emtrends_vmPFC_HC.R
# Plot emtrends from model

plot_emtrends_vmPFC_HC_by_condition <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR){
  
  
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
  pdf(fname, width = 15, height = 25)
  gg1 <- ggplot(emt,aes(x=evt_time,y=HCwithin.trend)) + 
    facet_grid(rewFunc~HC_region) +
    geom_point(aes(color=levels),size=5) +
    scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-2 std Entropy","+2 std Entropy")) + 
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
  grid.draw(gg2)
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
  pdf(fname, width = 15, height = 25)
  gg1 <- ggplot(emt,aes(x=evt_time,y=HCwithin.trend)) + 
    facet_grid(rewFunc~HC_region) +
    geom_point(aes(color=levels),size=5) +
    scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("10th %ile Value","90'th %ile Value")) + 
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
  grid.draw(gg2)
  dev.off()  
  
  
  rm(emt)
  
  if (strcmp(toalign,'feedback')){
  
  pal1090 = palette()
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1090[1] <- pal[[1]]
  pal1090[2] <- '#574c35'
  
  # do v_entropy_wi_change
  emt <- ddf$emtrends_list$dH_HC
  emt <- emt %>% filter(network=='D')
  
  if (strcmp(toalign,'clock')){
    emt$levels <- factor(emt$v_entropy_wi_change_lag, labels = c("-2 std Entropy Change","+2 std Entropy Change"))
  } else if (strcmp(toalign,'feedback')){
    emt$levels <- factor(emt$v_entropy_wi_change, labels = c("-2 std Entropy Change","+2 std Entropy Change"))
  }
  
  fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy_Change','-',hc_LorR, ".pdf", sep = "")
  pdf(fname, width = 15, height = 25)
  gg1 <- ggplot(emt,aes(x=evt_time,y=HCwithin.trend)) + 
    facet_grid(rewFunc~HC_region) +
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
  grid.draw(gg2)
  dev.off()
  }
}

# gg<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network1)) + 
#   geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
#   # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
#   geom_line(size = 1) + theme(legend.position = "none") +
#   geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
#   scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
#   #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
#   theme_bw(base_size=13) +
#   facet_wrap(~HC_region) +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_line(colour = "grey45"), 
#         panel.grid.minor = element_line(colour = "grey45"), 
#         panel.background = element_rect(fill = 'grey40'),
#         axis.title.y = element_text(margin=margin(r=6)),
#         axis.title.x = element_text(margin=margin(t=6)))
# }
# print(gg)
# dev.off()

# fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
# pdf(fname, width = 9, height = 3.5)
# gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=as.factor(HC_region), size=`p, FDR-corrected`)) +
#   geom_line(size = 1) + geom_point() +
#   geom_errorbar() +
#   geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) + facet_grid(~symmetry_group1) + 
#   scale_color_manual(values = pal_hc1,labels=c('AH','PH')) + xlab(epoch_label) +
#   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
# 
# gg2 <- ggplot_gtable(ggplot_build(gg1))
# stripr <- which(grepl('strip-t', gg2$layout$name))
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
#   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.draw(gg2)
# #print(gg2)
# dev.off()
