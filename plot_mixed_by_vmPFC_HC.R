plot_mixed_by_vmPFC_HC <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR){
  library(grid)
  ## Check plots
  if (strcmp(toalign,"clock")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction')
  } else if (strcmp(toalign,"feedback")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC_interaction')
  }

  message("\nPlotting streams decoding")
  library(viridis)
  library(wesanderson)
  pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
  pal = palette()
  pal[1] = pal3[2]
  pal[2] = pal3[1]
  pal[3] = pal3[3]
  if (strcmp(toalign,'feedback')){
    toalign_str <- 'feedback'
  } else if (strcmp(toalign,'clock')){
    toalign_str <- 'trial onset'  
  }
  epoch_label = paste("Time relative to",toalign_str, "[s]")
  pal_hc = wes_palette("Royal2", 5, type = "discrete")
  pal_hc1 <- palette()
  pal_hc1[2] <- '#fc795d'
  pal_hc1[1] <- '#1f53ec'
  fills <- palette()
  fills[1] <- pal[2]
  fills[2] <- pal[2]
  fills[3] <- pal[1]
  fills[4] <- pal[1]
  fills[5] <- pal[1]
  fills[6] <- pal[1]
  fills[7] <- pal[3]
  fills[8] <- pal[3]
  
  
  if (strcmp(toprocess,"network-by-outcome") | strcmp(toprocess,"symmetry-group-by-HC-by-outcome") | strcmp(toprocess,'network-by-HC-by-outcome')){
    colnames(ddf$coef_df_reml)[4] <- "outcome1"
  }
  ddg <- ddf
  ddf <- as_tibble(ddf$coef_df_reml)
  ddf$t <- ddf$evt_time
  #if (strcmp(toprocess,"network")){
  ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term, 
                         p_level_fdr = as.factor(case_when(
                           # p_fdr > .1 ~ '0',
                           # p_fdr < .1 & p_fdr > .05 ~ '1',
                           p_fdr > .05 ~ '1',
                           p_fdr < .05 & p_fdr > .01 ~ '2',
                           p_fdr < .01 & p_fdr > .001 ~ '3',
                           p_fdr <.001 ~ '4')))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  
  #ddf$`p, FDR-corrected` = 4;
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  #}
  # if (strcmp(toprocess,"region")){
  #   ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term,
  #                          p_level_fdr= as.factor(case_when(
  #                            # p_fdr > .1 ~ '0',
  #                            # p_fdr < .1 & p_fdr > .05 ~ '1',
  #                            p_fdr > .05 ~ '1',
  #                            p_fdr < .05 & p_fdr > .01 ~ '2',
  #                            p_fdr < .01 & p_fdr > .001 ~ '3',
  #                            p_fdr <.001 ~ '4')))
  #   ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  #   ddf$`p, FDR-corrected` = ddf$p_level_fdr
  # }
  #if ("network" %in% splits) {ddf$visuomotor_grad <- factor(ddf$visuomotor_grad, labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal"))}
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  
  if (strcmp(toprocess,"region")){ # get network for each region
    ddf <- ddf %>% mutate(network=case_when(
      region=="pfc11L" ~ "C",
      region=="pfc11aR" ~ "L",
      region=="pfc11b47R" ~ "C",
      region=="pfc11m13L" ~ "L",
      region=="pfc14m3225L" ~ "D",
      region=="pfc14m3225R" ~ "D",
      region=="pfc14r14c11mL" ~ "L",
      region=="pfc14r14c11mR" ~ "L",
      region=="pfc14rc2L" ~ "C",
      region=="pfc2432R" ~ "D",
      region=="pfc47L" ~ "C",
      region=="pfcd10R" ~ "D",
      region=="pfcfp10L" ~ "D",
      region=="pfcfp10R" ~ "L",
      region=="pfcp3224L" ~ "D",
      region=="pfcrdb10L" ~ "D",
      region=="pfcrl10R" ~ "C"
    )
    )
  }
  

  if (strcmp(toprocess,"symmetry-by-HC") | strcmp(toprocess,'symmetry-by-bin') | strcmp(toprocess,"symmetry-by-HC-by-outcome") | strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc') | strcmp(toprocess,'symmetry-group-by-HC-by-outcome')){
    ddf <- ddf %>% mutate(symmetry_group1=case_when(
      symmetry_group==1 ~ '11/47',
      symmetry_group==2 ~ 'fp10',
      symmetry_group==3 ~ '14rc11m',
      symmetry_group==4 ~ '14m25/32',
      symmetry_group==5 ~ '24/32',
      symmetry_group==6 ~ 'rl9_10',
      symmetry_group==7 ~ 'd10',
      symmetry_group==8 ~ '11/13')) %>% 
      mutate(network_group=case_when(
        symmetry_group==1 ~ 'CTR',
        symmetry_group==2 ~ 'DMN',
        symmetry_group==3 ~ 'LIM',
        symmetry_group==4 ~ 'DMN',
        symmetry_group==5 ~ 'DMN',
        symmetry_group==6 ~ 'CTR',
        symmetry_group==7 ~ 'DMN',
        symmetry_group==8 ~ 'LIM'
      ))
  }
  #fills = c('red','red','blue','blue','blue','cyan','green','magenta')
  #fills1 = c('green','magenta','red','blue','cyan')
  if (strcmp(toprocess,'network-by-HC') | strcmp(toprocess,'network-by-HC-by-side')){
    if (totest!='Explore-'){
      ddf <- ddf %>% mutate(network1 = 
                              case_when(network=='C'~'2C',
                                        network=='D'~'1D',
                                        network=='L'~'3L')) %>% 
        mutate(network2 = case_when(network=='C'~'CTR',
                                    network=='D'~'DMN',
                                    network=='L'~'LIM'))
    } else if (totest=='Explore-'){
      ddf <- ddf %>% mutate(network1 = 
                              case_when(network=='CTR'~'2C',
                                        network=='DMN'~'1D',
                                        network=='LIM'~'3L')) %>% 
        mutate(network2 = case_when(network=='CTR'~'CTR',
                                    network=='DMN'~'DMN',
                                    network=='LIM'~'LIM'))
    }
  }
  
  pal1 = palette()
  pal1[1] = '#ff484d'
  pal1[2] = '#2a52ff'
  
  if (!strcmp(totest,'anatomy')){
    for (fe in terms) {
      # fe <- terms[1] # test only
      edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed' & t < 8) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      if (strcmp(toprocess,'symmetry-by-HC') | (strcmp(toprocess,'symmetry-by-bin')) | strcmp(toprocess,'symmetry-group-by-HC-by-outcome') | strcmp(toprocess,'symmetry-by-HC-wPE') | strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
        edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9_10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
        edf$symmetry_group2 <- factor(edf$symmetry_group1, levels = c("14rc11m","11/47",'rl9_10',"fp10","24/32","14m25/32","d10","11/13"))
        edf$network_group = factor(edf$network_group,levels = c('CTR','DMN','LIM'))
      }
      # plot stream gradients
      
      if (!strcmp(toprocess,"symmetry-group-by-HC-by-outcome") & !strcmp(toprocess,"symmetry-group-by-HC-by-rewFunc")){
        # fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',hc_LorR, '-',model_iter,".pdf", sep = "")
        # pdf(fname, width = 9, height = 3.5)
      }
      
      if (strcmp(toprocess,"symmetry_group")){
        print(ggplot(edf, aes(t, symmetry_group)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_wrap(~side) +
                # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      } else if (strcmp(toprocess,"network")){
        print(ggplot(edf, aes(t, network2)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      } else if (strcmp(toprocess,"region")){
        print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_wrap(~network) +
                # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      } else if (strcmp(toprocess,"side1")){
        print(ggplot(edf, aes(t, side)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_wrap(~network) +
                # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      } else if (strcmp(toprocess,"region")){
        print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_wrap(~network) +
                # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      } else if (strcmp(toprocess,"network-by-HC") | strcmp(toprocess,'network-by-HClag')){
        # gg1 <- ggplot(edf,aes(t,network2))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
        #         facet_wrap(~HC_region) +
        #         geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        #         scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
        #         labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        # gg2 <- ggplot_gtable(ggplot_build(gg1))
        # stripr <- which(grepl('strip-t', gg2$layout$name))
        # k <- 2
        # for (i in stripr) {
        #   j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
        #   gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
        #   k <- k-1
        # }
        # grid.draw(gg2)
        
      } else if (strcmp(toprocess,"network-by-HC-wEC") | strcmp(toprocess,"network-by-HC-wO")) {
        print(ggplot(edf,aes(t,network))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_wrap(~HC_region) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      } else if (strcmp(toprocess,"region-by-HC")){
        print(ggplot(edf,aes(t,HC_region))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_grid(~region) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      } else if (strcmp(toprocess,"IS-by-HC")){
        print(ggplot(edf,aes(t,IS))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_grid(~HC_region) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      } else if (strcmp(toprocess,"symmetry-by-HC")){
        gg1 <- ggplot(edf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
          facet_grid(network_group~HC_region,scales="free_y",space="free_y") +
          geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        
        gg2 <- ggplot_gtable(ggplot_build(gg1))
        stripr <- which(grepl('strip-t', gg2$layout$name))
        k <- 1
        for (i in stripr) {
          j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
          gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_hc1[k]
          k <- k+1
        }
        grid.draw(gg2)
        
      } else if (strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
        qdf <- edf %>% filter(HC_region=='PH')
        fname = paste(totest,"_",toalign, "_", toprocess, "_", termstr, '_PH','-',hc_LorR, ".pdf", sep = "")
        
        pdf(fname, width = 9, height = 3.5)
        gg1 <- ggplot(qdf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
          facet_grid(network_group~rewFunc,scales="free_y",space="free_y") +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        
        gg2 <- ggplot_gtable(ggplot_build(gg1))
        stripr <- which(grepl('strip-t', gg2$layout$name))
        k <- 1
        for (i in stripr) {
          j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
          gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
          k <- k+1
        }
        grid.draw(gg2)
        dev.off()
        
        qdf <- edf %>% filter(HC_region=='AH')
        fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr, '_AH','-',hc_LorR, ".pdf", sep = "")
        
        pdf(fname, width = 9, height = 3.5)
        gg1 <- ggplot(qdf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
          facet_grid(network_group~rewFunc,scales="free_y",space="free_y") +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        
        gg2 <- ggplot_gtable(ggplot_build(gg1))
        stripr <- which(grepl('strip-t', gg2$layout$name))
        k <- 1
        for (i in stripr) {
          j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
          gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
          k <- k+1
        }
        grid.draw(gg2)
        dev.off()
        
        
        
        
      } else if (strcmp(toprocess,'symmetry-group-by-HC-by-outcome')){
        qdf <- edf %>% filter(outcome1=='Omission')
        fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr, '_Omission','-',hc_LorR, ".pdf", sep = "")
        
        pdf(fname, width = 9, height = 3.5)
        gg1 <- ggplot(qdf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
          facet_grid(network_group~HC_region,scales="free_y",space="free_y") +
          geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) +
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        
        gg2 <- ggplot_gtable(ggplot_build(gg1))
        stripr <- which(grepl('strip-t', gg2$layout$name))
        k <- 1
        for (i in stripr) {
          j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
          gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
          k <- k+1
        }
        grid.draw(gg2)
        dev.off()
        
        qdf <- edf %>% filter(outcome1=='Reward')
        fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr, '_Reward','-',hc_LorR, ".pdf", sep = "")
        
        pdf(fname, width = 9, height = 3.5)
        gg1 <- ggplot(qdf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
          facet_grid(network_group~HC_region,scales="free_y",space="free_y") +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        
        gg2 <- ggplot_gtable(ggplot_build(gg1))
        stripr <- which(grepl('strip-t', gg2$layout$name))
        k <- 1
        for (i in stripr) {
          j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
          gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
          k <- k+1
        }
        grid.draw(gg2)
        dev.off()
      } else if (strcmp(toprocess,'network-by-HC-by-outcome')){
        gg1 <- ggplot(edf,aes(t,network))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
          facet_grid(outcome1~HC_region) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg1)
      } else if (strcmp(toprocess,"network-by-HC-by-side") | strcmp(toprocess,'network-by-HClag')){
        print(ggplot(edf,aes(t,network2))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_wrap(HC_region~side_hc) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      } 
      
      if (!strcmp(toprocess,"symmetry-group-by-HC-by-outcome") & !strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
        #dev.off()
      }
      
      
      if (strcmp(toprocess,"network")){
        # fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, '-', model_iter, ".pdf", sep = "")
        # pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network2, size=`p, FDR-corrected`)) +
        #   geom_line(size = 1) + geom_point() +
        #   geom_errorbar() +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        #   scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        # print(gg)
        # dev.off()
      } else if (strcmp(toprocess,"side2")){
        fname = paste(behavmodel,'-',toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=side, size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_errorbar() +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_color_brewer(palette="Set1") + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off()
      } else if (strcmp(toprocess,"network-by-HC") | strcmp(toprocess,'network-by-HClag')){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, '-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
        #   geom_line(size = 1) + geom_point() +
        #   geom_errorbar() +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~HC_region) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        #   scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        if (all(edf$`p, FDR-corrected`=='p < .001')){
          gg1<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network2)) + 
            geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + 
            # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
            geom_line(size = 1) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
            geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
            scale_color_manual(values = pal) + 
            #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
            theme_bw(base_size=13) + ylab('Network Response') + 
            geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
            facet_wrap(~HC_region) + ylab('Network Response') +
            theme(legend.title = element_blank(),
                  panel.grid.major = element_line(colour = "grey45"), 
                  panel.grid.minor = element_line(colour = "grey45"), 
                  panel.background = element_rect(fill = 'grey40'),
                  axis.title.y = element_text(margin=margin(r=6),size=22),
                  axis.title.x = element_text(margin=margin(t=6),size=22),
                  legend.text = element_text(size=22),
                  axis.text.x = element_text(size=22),
                  axis.text.y = element_text(size=22)
            )
        } else {
          gg1<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network2)) + 
            geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
            # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
            geom_line(size = 1) + theme(legend.position = "none") +
            geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
            scale_color_manual(values = pal) + 
            #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
            theme_bw(base_size=13) +
            facet_wrap(~HC_region) + ylab('Network Response') +
            geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
            theme(legend.title = element_blank(),
                  panel.grid.major = element_line(colour = "grey45"), 
                  panel.grid.minor = element_line(colour = "grey45"), 
                  panel.background = element_rect(fill = 'grey40'),
                  axis.title.y = element_text(margin=margin(r=6),size=22),
                  axis.title.x = element_text(margin=margin(t=6),size=22),
                  legend.text = element_text(size=22),
                  axis.text.x = element_text(size=22),
                  axis.text.y = element_text(size=22)
            )
        }
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
        
      } else if (strcmp(toprocess,"network-by-HC-wEC") | strcmp(toprocess,"network-by-HC-wO")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network, size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_errorbar() +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~HC_region) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_color_brewer(palette="Set1") + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off()
      } else if (strcmp(toprocess,"region-by-HC")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=HC_region, size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_errorbar() +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_grid(~region) + 
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_color_brewer(palette="Set1") + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off()
      } else if (strcmp(toprocess,"IS-by-HC")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=as.factor(IS), size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_errorbar() +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_grid(~HC_region) + 
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_color_brewer(palette="Set1") + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off()
      } else if (strcmp(toprocess,"symmetry-by-HC")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, '-',model_iter,".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=as.factor(HC_region), size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_errorbar() +
          geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) + facet_grid(~symmetry_group1) + 
          scale_color_manual(values = pal_hc1,labels=c('Anterior Hippocampus','Posterior Hippocampus')) + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        
        gg2 <- ggplot_gtable(ggplot_build(gg1))
        stripr <- which(grepl('strip-t', gg2$layout$name))
        k <- 1
        for (i in stripr) {
          j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
          gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
          k <- k+1
        }
        grid.draw(gg2)
        #print(gg2)
        dev.off()
      } else if (strcmp(toprocess,"network-by-HC-by-side")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_errorbar() +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(HC_region~side_hc) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off()
      } 
    }
  } else if (strcmp(totest,'anatomy')){
    
    
    for (fe in terms) {
      edf <- ddf %>% filter(term == paste(fe) & t < 8) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      if (strcmp(toprocess,'symmetry-by-HC') | (strcmp(toprocess,'symmetry-by-bin')) | strcmp(toprocess,'symmetry-group-by-HC-by-outcome') | strcmp(toprocess,'symmetry-by-HC-wPE') | strcmp(toprocess,'symmetry-group-by-HC-by-rewFunc')){
        edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9_10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
        edf$symmetry_group2 <- factor(edf$symmetry_group1, levels = c("14rc11m","11/47",'rl9_10',"fp10","24/32","14m25/32","d10","11/13"))
        edf$network_group = factor(edf$network_group,levels = c('CTR','DMN','LIM'))
      }
      
      if (strcmp(toprocess,'region')){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        print(ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
        dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=HC_region, size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_errorbar() +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_color_brewer(palette="Set1") + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off()
      } else if (strcmp(toprocess,'symmetry-by-bin')){
        fname = paste(behavmodel,'-',totest,"_",toalign, toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg1 <- ggplot(edf,aes(t,bin_num))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
          facet_grid(~symmetry_group1,scales="free_y",space="free_y") +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        
        gg2 <- ggplot_gtable(ggplot_build(gg1))
        stripr <- which(grepl('strip-t', gg2$layout$name))
        k <- 1
        for (i in stripr) {
          j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
          gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
          k <- k+1
        }
        grid.draw(gg2)
        dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',hc_LorR, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=as.factor(bin_num), size=`p, FDR-corrected`)) +
          geom_line(size = 1) + geom_point() +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_color_brewer(palette="Set1") + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off()
      }
    }
  }
}