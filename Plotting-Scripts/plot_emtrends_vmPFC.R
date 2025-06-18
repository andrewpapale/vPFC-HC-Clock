# 2022-12-07 AndyP
# Plot emtrends_vmPFC.R
# Plot emtrends from model

plot_emtrends_vmPFC <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter){
  
  
  if (strcmp(toalign,'feedback')){
    toalign_str <- 'feedback'
  } else if (strcmp(toalign,'clock')){
    toalign_str <- 'trial start'
  }
  
  
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
  pal1 <- palette()
  pal1[[2]] <- pal[[1]]
  pal1[[1]] <- pal[[2]]
  pal1[[3]] <- pal[[3]]
  pal1090[1] <- pal[[1]]
  pal1090[2] <- '#574c35'
  pal1090[3] <- '#f5b32a'
  
if (!strcmp(behavmodel,'explore')){
if (!strcmp(totest,'online')){
   
    tryCatch({

      # do v_entropy
      emt <- ddf$emtrends_list$H
      emt$levels <- factor(emt$v_entropy_sc, labels = c("-1.5 std Entropy","+1.5 std Entropy"))
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy','-',hc_LorR, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=v_entropy_sc.trend)) + 
        facet_grid(~network) +
        geom_point(aes(color=levels),size=5) +
        scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-1.5 std Entropy","+1.5 std Entropy")) + 
        geom_errorbar(aes(ymin=v_entropy_sc.trend-std.error, ymax=v_entropy_sc.trend+std.error), width=0.5) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
        ylab('') + xlab(toalign_str)
      
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
    })
    
    rm(emt)
  

  
  pal1090 = palette()
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1090[1] <- pal[[2]]
  pal1090[2] <- '#7a7745'
  
  # do v_max_wi
  emt <- ddf$emtrends_list$V
  emt$levels <- factor(emt$v_max_wi, labels = c("-1.5 std Value","+ 1.5 std Value"))
  
  fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Value', ".pdf", sep = "")
  pdf(fname, width = 9, height = 3.5)
  gg1 <- ggplot(emt,aes(x=evt_time,y=v_max_wi.trend)) + 
    facet_grid(~network) +
    geom_point(aes(color=levels),size=5) +
    scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-1.5 std Value","+ 1.5 std Value")) + 
    geom_errorbar(aes(ymin=v_max_wi.trend-std.error, ymax=v_max_wi.trend+std.error), width=0.5) +
    geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
    ylab('') + xlab(toalign_str)
  
  gg2 <- ggplot_gtable(ggplot_build(gg1))
  stripr <- which(grepl('strip-t', gg2$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[k]
    k <- k+1
  }
  grid.draw(gg2)
  dev.off()  
  
  
  

    
    rm(emt)
    
    pal1090 = palette()
    pal = wes_palette("FantasticFox1", 3, type = "discrete")
    pal1090[1] <- pal[[1]]
    pal1090[2] <- '#574c35'
    
    # do v_entropy_wi_change
    emt <- ddf$emtrends_list$dH
    
    if (strcmp(toalign,'clock')){
      emt$levels <- factor(emt$v_entropy_wi_change_lag, labels = c("-0.5 std Entropy Change","+0.5 std Entropy Change"))
      emt$var.trend <- emt$v_entropy_wi_change_lag
    } else if (strcmp(toalign,'feedback')){
      emt$levels <- factor(emt$v_entropy_wi_change, labels = c("-0.5 std Entropy Change","+0.5 std Entropy Change"))
      emt$var.trend <- emt$v_entropy_wi_change
    }
    
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy_Change', ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=evt_time,y=var.trend)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-0.5 std Entropy Change","+0.5 std Entropy Change")) + 
      geom_errorbar(aes(ymin=var.trend-std.error, ymax=var.trend+std.error), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(toalign_str)
    
    gg2 <- ggplot_gtable(ggplot_build(gg1))
    stripr <- which(grepl('strip-t', gg2$layout$name))
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
      gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[k]
      k <- k+1
    }
    grid.draw(gg2)
    dev.off()
  

} else if (strcmp(totest,'online')){
  tryCatch({
    # do v_entropy
    emt <- ddf$emtrends_list$H
    browser()
    emt$levels <- factor(emt$v_entropy_sc, labels = c("-1.5 std Entropy","+1.5 std Entropy"))
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy', ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=online,y=v_entropy_sc.trend)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-1.5 std Entropy","+1.5 std Entropy")) + 
      geom_errorbar(aes(ymin=v_entropy_sc.trend-std.error, ymax=v_entropy_sc.trend+std.error), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(toalign_str)
    
    gg2 <- ggplot_gtable(ggplot_build(gg1))
    stripr <- which(grepl('strip-t', gg2$layout$name))
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
      gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[k]
      k <- k+1
    }
    grid.draw(gg2)
    dev.off()
  })
  
  rm(emt)
  
  
  
  pal1090 = palette()
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  pal1090[1] <- pal[[2]]
  pal1090[2] <- '#7a7745'
  
  # do v_max_wi
  emt <- ddf$emtrends_list$V
  emt$levels <- factor(emt$v_max_wi, labels = c("-1.5 std Value","+ 1.5 std Value"))
  
  fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Value', ".pdf", sep = "")
  pdf(fname, width = 9, height = 3.5)
  gg1 <- ggplot(emt,aes(x=online,y=v_max_wi.trend)) + 
    facet_grid(~network) +
    geom_point(aes(color=levels),size=5) +
    scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-1.5 std Value","+ 1.5 std Value")) + 
    geom_errorbar(aes(ymin=v_max_wi.trend-std.error, ymax=v_max_wi.trend+std.error), width=0.5) +
    geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
    ylab('') + xlab(toalign_str)
  
  gg2 <- ggplot_gtable(ggplot_build(gg1))
  stripr <- which(grepl('strip-t', gg2$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[k]
    k <- k+1
  }
  grid.draw(gg2)
  dev.off()  

    
    rm(emt)
    
    pal1090 = palette()
    pal = wes_palette("FantasticFox1", 3, type = "discrete")
    pal1090[1] <- pal[[1]]
    pal1090[2] <- '#574c35'
    
    # do v_entropy_wi_change
    emt <- ddf$emtrends_list$dH
    
    if (strcmp(toalign,'clock')){
      emt$levels <- factor(emt$v_entropy_wi_change_lag, labels = c("-0.5 std Entropy Change","+0.5 std Entropy Change"))
      emt$var.trend <- emt$v_entropy_wi_change_lag.trend
    } else if (strcmp(toalign,'feedback')){
      emt$levels <- factor(emt$v_entropy_wi_change, labels = c("-0.5 std Entropy Change","+0.5 std Entropy Change"))
      emt$var.trend <- emt$v_entropy_wi_change.trend
    }
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy_Change', ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=online,y=var.trend)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-0.5 std Entropy Change","+0.5 std Entropy Change")) + 
      geom_errorbar(aes(ymin=var.trend-std.error, ymax=var.trend+std.error), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(toalign_str)
    
    gg2 <- ggplot_gtable(ggplot_build(gg1))
    stripr <- which(grepl('strip-t', gg2$layout$name))
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
      gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[k]
      k <- k+1
    }
    grid.draw(gg2)
    dev.off()
  
}
} else if (strcmp(behavmodel,'explore')){
  
  if (strcmp(totest,'entropy')){
    tryCatch({
      emt <- ddf$emtrends_list$HxY
      emt <- emt %>% filter(education_yrs=='-1' | education_yrs=='1')
      emt$education_yrs <- factor(emt$education_yrs)
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy-Education-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=v_entropy_wi.trend)) + 
        facet_grid(~network) +
        geom_point(aes(color=education_yrs),size=5) +
        geom_errorbar(aes(ymin=v_entropy_wi.trend-std.error, ymax=v_entropy_wi.trend+std.error), width=0.5) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
        ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
      
      gg2 <- ggplot_gtable(ggplot_build(gg1))
      stripr <- which(grepl('strip-t', gg2$layout$name))
      k <- 1
      for (i in stripr) {
        j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
        gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[[k]]
        k <- k+1
      }
      grid.draw(gg2)
      dev.off()
      
      rm(emt)
    })
    
    tryCatch({
      emt <- ddf$emtrends_list$HxA
      emt <- emt %>% filter(age=='-1' | age=='1')
      emt$age <- factor(emt$age)
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy-Age-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=v_entropy_wi.trend)) + 
        facet_grid(~network) +
        geom_point(aes(color=age),size=5) +
        geom_errorbar(aes(ymin=v_entropy_wi.trend-std.error, ymax=v_entropy_wi.trend+std.error), width=0.5) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
        ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
      
      gg2 <- ggplot_gtable(ggplot_build(gg1))
      stripr <- which(grepl('strip-t', gg2$layout$name))
      k <- 1
      for (i in stripr) {
        j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
        gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[[k]]
        k <- k+1
      }
      grid.draw(gg2)
      dev.off()
      
      rm(emt)
    })
    
    tryCatch({
      emt <- ddf$emtrends_list$HxW
      emt <- emt %>% filter(wtar=='-1' | wtar=='1')
      emt$wtar <- factor(emt$wtar)
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Entropy-WTAR-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=v_entropy_wi.trend)) + 
        facet_grid(~network) +
        geom_point(aes(color=wtar),size=5) +
        geom_errorbar(aes(ymin=v_entropy_wi.trend-std.error, ymax=v_entropy_wi.trend+std.error), width=0.5) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
        ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
      
      gg2 <- ggplot_gtable(ggplot_build(gg1))
      stripr <- which(grepl('strip-t', gg2$layout$name))
      k <- 1
      for (i in stripr) {
        j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
        gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[[k]]
        k <- k+1
      }
      grid.draw(gg2)
      dev.off()
      
      rm(emt)
    })
  }
  if (strcmp(totest,'value')){

    tryCatch({
      emt <- ddf$emtrends_list$HxY
      emt <- emt %>% filter(education_yrs=='-1' | education_yrs=='1')
      emt$education_yrs <- factor(emt$education_yrs)
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Value-Education-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=v_max_wi.trend)) + 
        facet_grid(~network) +
        geom_point(aes(color=education_yrs),size=5) +
        geom_errorbar(aes(ymin=v_max_wi.trend-std.error, ymax=v_max_wi.trend+std.error), width=0.5) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
        ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
      
      gg2 <- ggplot_gtable(ggplot_build(gg1))
      stripr <- which(grepl('strip-t', gg2$layout$name))
      k <- 1
      for (i in stripr) {
        j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
        gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[[k]]
        k <- k+1
      }
      grid.draw(gg2)
      dev.off()
      
      rm(emt)
    })
    
    
    tryCatch({
      emt <- ddf$emtrends_list$HxA
      emt <- emt %>% filter(age=='-1' | age=='1')
      emt$age <- factor(emt$age)
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Value-Age-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=v_max_wi.trend)) + 
        facet_grid(~network) +
        geom_point(aes(color=age),size=5) +
        geom_errorbar(aes(ymin=v_max_wi.trend-std.error, ymax=v_max_wi.trend+std.error), width=0.5) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
        ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
      
      gg2 <- ggplot_gtable(ggplot_build(gg1))
      stripr <- which(grepl('strip-t', gg2$layout$name))
      k <- 1
      for (i in stripr) {
        j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
        gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[[k]]
        k <- k+1
      }
      grid.draw(gg2)
      dev.off()
      
      rm(emt)
    })
    
    tryCatch({
      emt <- ddf$emtrends_list$HxW
      emt <- emt %>% filter(wtar=='-1' | wtar=='1')
      emt$age <- factor(emt$wtar)
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emtrends_", toprocess, "_", 'Value-WTAR-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=v_max_wi.trend)) + 
        facet_grid(~network) +
        geom_point(aes(color=wtar),size=5) +
        geom_errorbar(aes(ymin=v_max_wi.trend-std.error, ymax=v_max_wi.trend+std.error), width=0.5) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
        ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
      
      gg2 <- ggplot_gtable(ggplot_build(gg1))
      stripr <- which(grepl('strip-t', gg2$layout$name))
      k <- 1
      for (i in stripr) {
        j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
        gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[[k]]
        k <- k+1
      }
      grid.draw(gg2)
      dev.off()
      
      rm(emt)
    })
  }
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
