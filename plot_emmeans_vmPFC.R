# 2022-05-12 AndyP
# Plot emmeans_vmPFC.R
# Plot emmeans from model

plot_emmeans_vmPFC <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter){
  
  library(grid)
  
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
  pal1090[1] <- '#939799'
  pal1090[2] <- '#5A5A5A'
  
  if (!strcmp(behavmodel,'explore' & !grepl('explore',totest,fixed=TRUE))){
  
  if (!strcmp(totest,'online')){
    
    tryCatch({
    # do v_entropy
    emt <- ddf$emmeans_list$H
    emt$levels <- factor(emt$v_entropy_sc, labels = c("-1.5 std Entropy","+1.5 std Entropy"))
    
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Entropy-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-1.5 std Entropy","+1.5 std Entropy")) + 
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
    # do v_max_wi
    emt <- ddf$emmeans_list$V
    emt$levels <- factor(emt$v_max_wi, labels = c("-1.5 std Value","+1.5 std Value"))
    
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Value-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2]),labels= c("-1.5 std Value","+1.5 std Value")) + 
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
    # do v_entropy_wi_change
    emt <- ddf$emmeans_list$dH
    
    if (strcmp(toalign,'feedback')){
      emt$levels <- factor(emt$v_entropy_wi_change)
    }else if (strcmp(toalign,'clock')){
      emt$levels <- factor(emt$v_entropy_wi_change_lag)
    }
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Entropy_Change-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2])) + 
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
    })
    
    rm(emt)
    
    tryCatch({
    # do v_entropy_wi_change
    emt <- ddf$emmeans_list$Tr
    
    emt$levels <- factor(emt$trial_neg_inv_sc)
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Trial-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2])) + 
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
    
    })
    
  } else if (strcmp(totest,'online')){
    
    # do v_entropy
    emt <- ddf$emmeans_list$H
    emt$levels <- factor(emt$v_entropy_sc, labels = c("-1.5 std Entropy","+1.5 std Entropy"))
    if (strcmp(toalign,'feedback')){
      x_lab_str <- ''
    }else if (strcmp(toalign,'clock')){
      x_lab_str <- ''
    }
    if (strcmp(toalign,'feedback')){
      emt$online <- factor(emt$online,levels=c('online_pre','offline_pre','online','offline_post'))
    } else if (strcmp(toalign,'clock')){
      emt$online <- factor(emt$online,levels=c('online_pre','offline_pre','online','offline_post'))
    }
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Entropy-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=online,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c("-1.5 std Entropy","+1.5 std Entropy")) + 
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(x_lab_str) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    
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
    
    # do v_max_wi
    emt <- ddf$emmeans_list$V
    emt$levels <- factor(emt$v_max_wi, labels = c("-1.5 std Value","+1.5 std Value"))
    if (strcmp(toalign,'feedback')){
      x_lab_str <- ''
    }else if (strcmp(toalign,'clock')){
      x_lab_str <- ''
    }
    if (strcmp(toalign,'feedback')){
      emt$online <- factor(emt$online,levels=c('online_pre','offline_pre','online','offline_post'))
    } else if (strcmp(toalign,'clock')){
      emt$online <- factor(emt$online,levels=c('online_pre','offline_pre','online','offline_post'))
    }
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Value-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=online,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2]),labels= c("-1.5 std Value","+1.5 std Value")) + 
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(x_lab_str) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    
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
    
    # do v_entropy_wi_change
    emt <- ddf$emmeans_list$dH
    
    if (strcmp(toalign,'feedback')){
      emt$levels <- factor(emt$v_entropy_wi_change)
      x_lab_str <- ''
    }else if (strcmp(toalign,'clock')){
      emt$levels <- factor(emt$v_entropy_wi_change_lag)
      x_lab_str <- ''
    }
    if (strcmp(toalign,'feedback')){
      emt$online <- factor(emt$online,levels=c('online_pre','offline_pre','online','offline_post'))
    } else if (strcmp(toalign,'clock')){
      emt$online <- factor(emt$online,levels=c('online_pre','offline_pre','online','offline_post'))
    }
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Entropy_Change-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=online,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2])) + 
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(x_lab_str) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    
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
    # do Tr
    emt <- ddf$emmeans_list$Tr
    if (strcmp(toalign,'feedback')){
      x_lab_str <- ''
    }else if (strcmp(toalign,'clock')){
      x_lab_str <- ''
    }
    if (strcmp(toalign,'feedback')){
      emt$online <- factor(emt$online,levels=c('online_pre','offline_pre','online','offline_post'))
    } else if (strcmp(toalign,'clock')){
      emt$online <- factor(emt$online,levels=c('online_pre','offline_pre','online','offline_post'))
    }
    emt$levels <- factor(emt$trial_neg_inv_sc)
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Trial-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=online,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      scale_color_manual(values=c(pal1090[1],pal1090[2])) + 
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
      ylab('') + xlab(x_lab_str) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    
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
  }
  
  
  # # do Trial
  # emt <- ddf$emmeans_list$Tr
  # 
  # emt$levels <- factor(emt$v_entropy_wi_change_lag, labels = c('Trial 5','Trial 45'))
  # 
  # fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Entropy_Change', ".pdf", sep = "")
  # pdf(fname, width = 9, height = 3.5)
  # gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
  #   facet_grid(~network) +
  #   geom_point(aes(color=levels),size=5) +
  #   scale_color_manual(values=c(pal1090[1],pal1090[2]),labels=c('-0.5 std Entropy Change','No Entropy Change','+0.5 std Entropy Change')) + 
  #   geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
  #   geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
  #   ylab('') + xlab(paste0('Time relative to ', toalign_str,' [s]'))
  # 
  # gg2 <- ggplot_gtable(ggplot_build(gg1))
  # stripr <- which(grepl('strip-t', gg2$layout$name))
  # k <- 1
  # for (i in stripr) {
  #   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  #   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal1[[k]]
  #   k <- k+1
  # }
  # grid.draw(gg2)
  # dev.off()
  
  
  } else if (strcmp(behavmodel,'explore')){
    
    if (strcmp(totest,'entropy')){
      
      emt <- ddf$emmeans_list$Y
      emt <- emt %>% filter(education_yrs=='-1' | education_yrs=='1')
      emt$levels <- factor(emt$education_yrs)
      
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Education-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
        facet_grid(~network) +
        geom_point(aes(color=levels),size=5) +
        geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
      
      
      emt <- ddf$emmeans_list$A
      emt <- emt %>% filter(age=='-1' | age=='1')
      emt$levels <- factor(emt$age)
      
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Age-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
        facet_grid(~network) +
        geom_point(aes(color=levels),size=5) +
        geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
      
      emt <- ddf$emmeans_list$W
      emt <- emt %>% filter(wtar=='-1' | wtar=='1')
      emt$levels <- factor(emt$wtar)
      
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'WTAR-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
        facet_grid(~network) +
        geom_point(aes(color=levels),size=5) +
        geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
      
    tryCatch({
    emt <- ddf$emmeans_list$H
    emt <- emt %>% filter(v_entropy_wi=='-2' | v_entropy_wi=='2')
    emt$levels <- factor(emt$v_entropy_wi)
    
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Entropy-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=levels),size=5) +
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
    
      emt <- ddf$emmeans_list$Y
      emt <- emt %>% filter(education_yrs=='-1' | education_yrs=='1')
      emt$levels <- factor(emt$education_yrs)
      
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Education-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
        facet_grid(~network) +
        geom_point(aes(color=levels),size=5) +
        geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
      
      
      emt <- ddf$emmeans_list$A
      emt <- emt %>% filter(age=='-1' | age=='1')
      emt$levels <- factor(emt$age)
      
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Age-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
        facet_grid(~network) +
        geom_point(aes(color=levels),size=5) +
        geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
      
      emt <- ddf$emmeans_list$W
      emt <- emt %>% filter(wtar=='-1' | wtar=='1')
      emt$levels <- factor(emt$wtar)
      
      fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'WTAR-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
        facet_grid(~network) +
        geom_point(aes(color=levels),size=5) +
        geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
      
    tryCatch({
    emt <- ddf$emmeans_list$H
    emt <- emt %>% filter(v_max_wi=='-2' | v_max_wi=='2')
    emt$v_max_wi <- factor(emt$v_max_wi)
    
    fname = paste(behavmodel,'-',totest,"_",toalign, "_emmeans_", toprocess, "_", 'Value-',model_iter, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg1 <- ggplot(emt,aes(x=evt_time,y=estimate)) + 
      facet_grid(~network) +
      geom_point(aes(color=v_max_wi),size=5) +
      geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.5) +
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
