plot_mixed_by_vmPFC_by_trial_bin <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter){
  
  ## Check plots
  if (strcmp(toalign,"clock")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock')
  } else if (strcmp(toalign,"feedback")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback')
  } else if (strcmp(toalign,'rt_vmax')){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_rt_vmax')
  }
  
  message("\nPlotting streams decoding")
  library(wesanderson)
  library(pracma)
  library(ggplot2)
  library(tidyverse)
  library(viridis)
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  if (strcmp(toalign,'feedback')){
    toalign_str <- 'feedback'
  } else if (strcmp(toalign,'clock')){
    toalign_str <- 'trial onset'  
  }
  epoch_label = paste("Time relative to",toalign_str, "[s]")
  # fills for headers for symmetry_group.  Fill by network.
  fills <- palette()
  fills[1] <- pal[2]
  fills[2] <- pal[2]
  fills[3] <- pal[1]
  fills[4] <- pal[1]
  fills[5] <- pal[1]
  fills[6] <- pal[1]
  fills[7] <- pal[3]
  fills[8] <- pal[3]
  
  sym_fill <- palette()
  sym_fill[1] <- '#57007a'
  sym_fill[2] <- '#2d26a4'
  sym_fill[3] <- '#1f515d'
  sym_fill[4] <- '#3b8e1d'
  sym_fill[5] <- '#e2b410'
  sym_fill[6] <- '#ff6769'
  sym_fill[7] <- '#ff6b6d'
  sym_fill[8] <- '#fff5fa'
  # edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9/10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
  sym_fill_rearrange <- palette()
  sym_fill_rearrange[1] <- sym_fill[6]
  sym_fill_rearrange[2] <- sym_fill[1]
  sym_fill_rearrange[3] <- sym_fill[7]
  sym_fill_rearrange[4] <- sym_fill[5]
  sym_fill_rearrange[5] <- sym_fill[4]
  sym_fill_rearrange[6] <- sym_fill[2]
  sym_fill_rearrange[7] <- sym_fill[8]
  sym_fill_rearrange[8] <- sym_fill[3]
  
  #colnames(ddf$coef_df_reml)[3] <- "outcome1"
  qdf <- ddf
  ddf <- ddf$coef_df_reml
  ddf$t <- ddf$evt_time
  #if (strcmp(toprocess,"network")){
  ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
  ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
    # p_fdr > .1 ~ '0',
    # p_fdr < .1 & p_fdr > .05 ~ '1',
    p_fdr > .05 ~ '1',
    p_fdr < .05 & p_fdr > .01 ~ '2',
    p_fdr < .01 & p_fdr > .001 ~ '3',
    p_fdr <.001 ~ '4')))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
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
  
  if (strcmp(toprocess,"symmetry_group") | strcmp(toprocess,"symmetry_group_by_rewFunc") | strcmp(toprocess,"symmetry_group_by_outcome")){
    ddf <- ddf %>% mutate(symmetry_group1=case_when(
      symmetry_group==1 ~ '11/47',
      symmetry_group==2 ~ 'fp10',
      symmetry_group==3 ~ '14rc11m',
      symmetry_group==4 ~ '14m25/32',
      symmetry_group==5 ~ '24/32',
      symmetry_group==6 ~ 'rl9/10',
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
      )) %>%
      mutate(sym_fill1=case_when( # edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9/10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
        symmetry_group==1 ~ sym_fill[1],
        symmetry_group==2 ~ sym_fill[2],
        symmetry_group==3 ~ sym_fill[3],
        symmetry_group==4 ~ sym_fill[4],
        symmetry_group==5 ~ sym_fill[5],
        symmetry_group==6 ~ sym_fill[6],
        symmetry_group==7 ~ sym_fill[7],
        symmetry_group==8 ~ sym_fill[8]
      ))
  }
  if (strcmp(toprocess,'network')){
    ddf <- ddf %>% mutate(network1 = 
                            case_when(network=='C'~'2C',
                                      network=='D'~'1D',
                                      network=='L'~'3L')) %>%
      mutate(network2 = case_when(network=='C'~'CTR',
                                  network=='D'~'DMN',
                                  network=='L'~'LIM'))
  }
  
  #fills = c('red','red','blue','blue','blue','blue','green','green')
  #fills1 = c('red','blue','green')
  ddf$t <- as.numeric(ddf$t)
  for (fe in terms) {
    # fe <- terms[1] # test only
    edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed' & t < 8)
    if (strcmp(toprocess,'symmetry_group') | strcmp(toprocess,"symmetry_group_by_rewFunc") | strcmp(toprocess,'symmetry_group_by_outcome')){
      edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9/10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
      edf$symmetry_group2 <- factor(edf$symmetry_group1, levels = c("14rc11m","11/47",'rl9/10',"fp10","24/32","14m25/32","d10","11/13"))
      # sym_fill[1] <- '#57007a'
      # sym_fill[2] <- '#2d26a4'
      # sym_fill[3] <- '#1f515d'
      # sym_fill[4] <- '#3b8e1d'
      # sym_fill[5] <- '#e2b410'
      # sym_fill[6] <- '#ff6769'
      # sym_fill[7] <- '#ff6b6d'
      # sym_fill[8] <- '#fff5fa'
      edf$sym_fill2 <- factor(edf$sym_fill1, levels = c('#ff6769','#57007a','#ff6b6d','#e2b410','#3b8e1d','#2d26a4','#fff5fa','#1f515d'))
      edf$network_group = factor(edf$network_group,levels = c('CTR','DMN','LIM'))
    }
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    # plot stream gradients
    if (strcmp(toprocess,"network")){
      fname = paste('outcome','-',behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
      pdf(fname, width = 25, height = 15)
      
      if (all(edf$`p, FDR-corrected`=='p < .001')){
        gg<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network1)) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + 
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          facet_wrap(~trial_bin) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
      } else {
        gg<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network1)) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") +
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          facet_wrap(~trial_bin) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
      }
      print(gg)
      dev.off()
    } else if (strcmp(toprocess,"symmetry_group")){
      fname = paste('reFunc','-',behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, '_',model_iter, ".pdf", sep = "")
      pdf(fname, width = 25, height = 15)
      gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, alpha=`p, FDR-corrected`)) +
        geom_point(aes(size=`p, FDR-corrected`,color=sym_fill2)) +
        geom_errorbar(size = 1) +
        scale_color_manual(values = sym_fill_rearrange,labels=sym_fill_rearrange,guide="none") + 
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) + facet_wrap(trial_bin~symmetry_group1) + 
        xlab(epoch_label) +
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
      dev.off()
    }
  }
}
