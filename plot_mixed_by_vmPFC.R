plot_mixed_by_vmPFC <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter){
  ## Check plots
  library(pracma)
  
  if (strcmp(toalign,"clock") & !(strcmp(behavmodel,'explore') | strcmp(behavmodel,'explore-interaction'))){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock')
  } else if (strcmp(toalign,"feedback") & !(strcmp(behavmodel,'explore') | strcmp(behavmodel,'explore-interaction'))){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback')
  } else if (strcmp(toalign,'rt_vmax') & !(strcmp(behavmodel,'explore') | strcmp(behavmodel,'explore-interaction'))){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_rt_vmax')
  } else if (strcmp(toalign,'feedback') & (strcmp(behavmodel,'explore') | strcmp(behavmodel,'explore-interaction'))){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_explore_feedback')
  } else if (strcmp(toalign,'clock') & (strcmp(behavmodel,'explore') | strcmp(behavmodel,'explore-interaction'))){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_explore_clock')
  }

  message("\nPlotting streams decoding")
  library(wesanderson)
  library(pracma)
  library(ggplot2)
  library(tidyverse)
  library(viridis)
  library(grid)
  library(scales)
  pal = wes_palette("FantasticFox1", 3, type = "discrete")
  if (strcmp(toalign,'feedback') & (!strcmp(totest,'online'))){
    toalign_str <- 'feedback'
    epoch_label = paste("Time relative to",toalign_str, "[s]")
  } else if (strcmp(toalign,'clock') & (!strcmp(totest,'online'))){
    toalign_str <- 'trial onset'
    epoch_label = paste("Time relative to",toalign_str, "[s]")
  } else if (strcmp(toalign,'feedback') & strcmp(totest,'online')){
    epoch_label <- ''
  } else if (strcmp(toalign,'clock') & strcmp(totest,'online')){
    epoch_label <- ''
  }
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
  # sym_fill[1] <- '#57007a'
  # sym_fill[2] <- '#2d26a4'
  # sym_fill[3] <- '#1f515d'
  # sym_fill[4] <- '#3b8e1d'
  # sym_fill[5] <- '#e2b410'
  # sym_fill[6] <- '#ff6769'
  # sym_fill[7] <- '#ff6b6d'
  # sym_fill[8] <- '#fff5fa'
  sym_fill[1] <- '#ffe680'
  sym_fill[2] <- '#ffcc00'
  sym_fill[3] <- '#ffb380'
  sym_fill[4] <- '#f58b45'
  sym_fill[5] <- '#c87137'
  sym_fill[6] <- '#aa4400'
  sym_fill[7] <- '#0000ff'
  sym_fill[8] <- '#000080'
  sym_fill_rearranged = palette()
  sym_fill_rearranged[1] <- sym_fill[8]
  sym_fill_rearranged[2] <- sym_fill[7]
  sym_fill_rearranged[3] <- sym_fill[6]
  sym_fill_rearranged[4] <- sym_fill[5]
  sym_fill_rearranged[5] <- sym_fill[4]
  sym_fill_rearranged[6] <- sym_fill[3]
  sym_fill_rearranged[7] <- sym_fill[2]
  sym_fill_rearranged[8] <- sym_fill[1]
    if (strcmp(toprocess,"network-by-outcome") | strcmp(toprocess,'symmetry_group_by_outcome')){
  colnames(ddf$coef_df_reml)[4] <- "outcome1"
}
qdf <- ddf
ddf <- ddf$coef_df_reml
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
        symmetry_group1=='rl9/10' ~ sym_fill[1],
        symmetry_group1=='11/47' ~ sym_fill[2],
        symmetry_group1=='d10' ~ sym_fill[3],
        symmetry_group1=='24/32' ~ sym_fill[4],
        symmetry_group1=='14m25/32' ~ sym_fill[5],
        symmetry_group1=='fp10' ~ sym_fill[6],
        symmetry_group1=='11/13' ~ sym_fill[7],
        symmetry_group1=='14rc11m' ~ sym_fill[8]
      ))
}
if ((strcmp(toprocess,'network') | strcmp(toprocess,'network-by-rewFunc')) & !(strcmp(behavmodel,'explore') | strcmp(behavmodel,'explore-interaction'))){
  ddf <- ddf %>% mutate(network1 = 
                          case_when(network=='C'~'2C',
                                    network=='D'~'1D',
                                    network=='L'~'3L')) %>%
    mutate(network2 = case_when(network=='C'~'CTR',
                                network=='D'~'DMN',
                                network=='L'~'LIM'))
} else if ((strcmp(toprocess,'network') | strcmp(toprocess,'network-by-block') | strcmp(toprocess,'network-by-rewFunc')) & (strcmp(behavmodel,'explore') | strcmp(behavmodel,'Explore') | strcmp(behavmodel,'explore-interaction'))){
  ddf <- ddf %>% mutate(network1 = network, network2 = network)
}

ddf$network2 <- factor(ddf$network2,levels=c('CTR','DMN','LIM'))
pal1 <- palette()
pal1[1] <- pal[2]
pal1[2] <- pal[1]
pal1[3] <- pal[3]
#fills = c('red','red','blue','blue','blue','blue','green','green')
#fills1 = c('red','blue','green')

if (!strcmp(totest,'online')){
  ddf$t <- ddf$evt_time
  ddf$t <- as.numeric(ddf$t)
  for (fe in terms) {
    # fe <- terms[1] # test only
    edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed' & t < 8)
    if (strcmp(toprocess,'symmetry_group') | strcmp(toprocess,"symmetry_group_by_rewFunc") | strcmp(toprocess,'symmetry_group_by_outcome')){
      edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9/10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
      edf$symmetry_group2 <- factor(edf$symmetry_group1, levels = c("14rc11m","11/47",'rl9/10',"fp10","24/32","14m25/32","d10","11/13"))
      # new
      # sym_fill[1] <- '#ffe680'
      # sym_fill[2] <- '#aa4400'
      # sym_fill[3] <- '#000080'
      # sym_fill[4] <- '#c87137'
      # sym_fill[5] <- '#f58b45'
      # sym_fill[6] <- '#ffb380'
      # sym_fill[7] <- '#ff6b6d'
      # sym_fill[8] <- '#0000ff'
      #sym_fill2 <- factor(edf$sym_fill1, levels = c(sym_fill_rearrange[1],sym_fill_rearrange[2],sym_fill_rearrange[3],sym_fill_rearrange[4],sym_fill_rearrange[5],sym_fill_rearrange[6],sym_fill_rearrange[7],sym_fill_rearrange[8]))
      edf$network_group = factor(edf$network_group,levels = c('CTR','DMN','LIM'))
    }
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    # plot stream gradients
    fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'_',model_iter,'_',model_iter,  ".pdf", sep = "")
    
    if (strcmp(toprocess,"network")){
      pdf(fname, width = 9, height = 3.5)
      print(ggplot(edf, aes(t, network2)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
              # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      dev.off()
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
      print(ggplot(edf,aes(t,network))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
              facet_wrap(~HC_region) +
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
    } else if (strcmp(toprocess,"network-by-HC-wEC") | strcmp(toprocess,"network-by-HC-wO")) {
      print(ggplot(edf,aes(t,network))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
              facet_wrap(~HC_region) +
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
    }  else if (strcmp(toprocess,"region-by-HC")){
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
    } else if (strcmp(toprocess,"symmetry_group")){
      pdf(fname, width = 9, height = 3.5)
      gg1 <- (ggplot(edf,aes(t,symmetry_group1))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_grid(network_group~.,scales="free_y",space="free_y") +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      print(gg1)
      dev.off()
    } else if (strcmp(toprocess,"network-by-outcome")){
      gg1 <- (ggplot(edf,aes(t,network))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_grid(~outcome) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      print(gg1)
    } else if (strcmp(toprocess,"symmetry_group_by_rewFunc")){
      gg1 <- (ggplot(edf,aes(t,symmetry_group2))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_grid(~rewFunc) +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      print(gg1)
    } else if (strcmp(toprocess,'symmetry_group_by_outcome')) {
      gg1 <- (ggplot(edf,aes(t,symmetry_group2))+geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
                facet_grid(network_group~outcome,scales="free_y",space="free_y") +
                geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
                scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
                labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      print(gg1)
    }
    
    if (strcmp(toprocess,"network")){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      # gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
      #   geom_line(size = 1) + geom_point() +
      #   geom_errorbar() +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~HC_region) +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
      #   scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
      #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      if (all(edf$`p, FDR-corrected`=='p < .001')){
        gg<-ggplot(edf, aes(x=t, y=estimate,group=network2,color=network2)) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + 
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          scale_color_manual(values = pal1) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          #facet_wrap(~HC_region) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6),size=48),
                axis.title.x = element_text(margin=margin(t=6),size=48),
                legend.text = element_text(size=24),
                axis.text.x = element_text(size=24),
                axis.text.y = element_text(size=24)
          )
      } else {
        gg<-ggplot(edf, aes(x=t, y=estimate,group=network2,color=network2)) + 
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('Network Response') +
          geom_hline(yintercept = 0, lty = 'dashed',color = 'gray', size=1) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") +
          scale_color_manual(values = pal1) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          #facet_wrap(~HC_region) +
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
      print(gg)
      dev.off()
    } else if (strcmp(toprocess,"side2")){
      fname = paste(behavmodel,'-',toalign, "_line_", toprocess, "_", termstr, '_',model_iter, ".pdf", sep = "")
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
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, '_',model_iter, ".pdf", sep = "")
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
    } else if (strcmp(toprocess,"network-by-HC-wEC") | strcmp(toprocess,"network-by-HC-wO")){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
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
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
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
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
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
    } else if (strcmp(toprocess,"symmetry_group")){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, '_',model_iter, ".pdf", sep = "")
      pdf(fname, width = 12, height = 3.5)
      gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, alpha=`p, FDR-corrected`)) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) + facet_grid(~symmetry_group1) + 
        geom_hline(yintercept = 0, lty = 'dashed',color = 'gray', size=1) + 
        geom_point(aes(size=`p, FDR-corrected`,color=sym_fill1)) +
        geom_errorbar(size = 1) + ylab('Region Response') + 
        scale_x_continuous(breaks=c(-3,0,3),labels=c('-3','0','3')) +
        theme(legend.title = element_blank(),
              axis.title.y = element_text(margin=margin(r=6),size=22),
              axis.title.x = element_text(margin=margin(t=6),size=22),
              legend.text = element_text(size=22),
              axis.text.x = element_text(size=22),
              axis.text.y = element_text(size=22)
        ) + 
        scale_color_manual(values = sym_fill_rearranged,guide="none") + 
        xlab(epoch_label) +
        labs(alpha = expression(italic(p)[FDR]))
      
      gg2 <- ggplot_gtable(ggplot_build(gg1))
      stripr <- which(grepl('strip-t', gg2$layout$name))
      k <- 1
      for (i in stripr) {
        j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
        gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
        k <- k+1
      }
      grid.draw(gg2)
      
      # gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, size=`p, FDR-corrected`)) +
      #   geom_line(size = 1) + geom_point() +
      #   geom_errorbar() +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_grid(~symmetry_group1) + 
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
      #   scale_color_brewer(palette="Set1") + xlab(epoch_label) +
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
      dev.off()
    } else if (strcmp(toprocess,"network-by-outcome")){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network, size=`p, FDR-corrected`)) +
        geom_line(size = 1) + geom_point() +
        geom_errorbar() + facet_wrap(~outcome) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        scale_color_brewer(palette="Set1") + xlab(epoch_label) +
        labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      print(gg)
      dev.off() 
    } else if (strcmp(toprocess,"symmetry_group_by_rewFunc")){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=rewFunc, size=`p, FDR-corrected`)) +
        geom_line(size = 1) + geom_point() +
        geom_errorbar() + facet_wrap(~symmetry_group1) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        scale_color_brewer(palette="Set1") + xlab(epoch_label) +
        labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      print(gg)
      dev.off() 
    } else if(strcmp(toprocess,'symmetry_group_by_outcome')){
    } else if (strcmp(toprocess,'atlas_value')){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, '_',model_iter, ".pdf", sep = "")
      pdf(fname, width = 12, height = 3.5)
      gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, alpha=`p, FDR-corrected`)) +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) + facet_grid(~atlas_value) + 
        geom_hline(yintercept = 0, lty = 'dashed',color = 'gray', size=1) + 
        geom_point(aes(size=`p, FDR-corrected`)) +
        geom_errorbar(size = 1) +
        xlab(epoch_label) +
        labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      print(gg1)
      dev.off()
    } else if (strcmp(toprocess,'network-by-block')){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      # gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
      #   geom_line(size = 1) + geom_point() +
      #   geom_errorbar() +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~HC_region) +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
      #   scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
      #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      if (all(edf$`p, FDR-corrected`=='p < .001')){
        gg<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network1)) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + facet_wrap(~block) + 
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          #facet_wrap(~HC_region) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
      } else {
        gg<-ggplot(edf, aes(x=t, y=estimate,group=network1,color=network1)) + facet_wrap(~block) + 
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          geom_hline(yintercept = 0, lty = 'dashed',color = 'gray', size=1) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") +
          scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          #facet_wrap(~HC_region) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
      }
      print(gg)
      dev.off()
    
    } else if (strcmp(toprocess,'atlas_value-by-block')){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
      pdf(fname, width = 36, height = 36)
      # gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
      #   geom_line(size = 1) + geom_point() +
      #   geom_errorbar() +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~HC_region) +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
      #   scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
      #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      if (all(edf$`p, FDR-corrected`=='p < .001')){
        gg<-ggplot(edf, aes(x=t, y=estimate)) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + facet_grid(atlas_value~block) + 
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          #facet_wrap(~HC_region) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
      } else {
        gg<-ggplot(edf, aes(x=t, y=estimate)) + facet_grid(atlas_value~block) + 
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          geom_hline(yintercept = 0, lty = 'dashed',color = 'gray', size=1) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") +
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          #facet_wrap(~HC_region) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
      }
      print(gg)
      dev.off()
    } else if (strcmp(toprocess,'network-by-rewFunc')){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      # gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
      #   geom_line(size = 1) + geom_point() +
      #   geom_errorbar() +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~HC_region) +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
      #   scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
      #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      if (all(edf$`p, FDR-corrected`=='p < .001')){
        gg<-ggplot(edf, aes(x=t, y=estimate,group=network2,color=network2)) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + 
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          scale_color_manual(values = pal1) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          facet_wrap(~rewFunc) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6),size=48),
                axis.title.x = element_text(margin=margin(t=6),size=48),
                legend.text = element_text(size=24),
                axis.text.x = element_text(size=24),
                axis.text.y = element_text(size=24)
          )
      } else {
        gg<-ggplot(edf, aes(x=t, y=estimate,group=network2,color=network2)) + 
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('Network Response') +
          geom_hline(yintercept = 0, lty = 'dashed',color = 'gray', size=1) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") +
          scale_color_manual(values = pal1) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          facet_wrap(~rewFunc) +
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
      print(gg)
      dev.off()
    }
  } 
  
  

  } else if (strcmp(totest,'online')){
  # no evt_time
  for (fe in terms) {
    # fe <- terms[1] # test only
    edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed')
    if (strcmp(toalign,'feedback')){
      edf$online <- factor(edf$online,levels=c('online_pre','offline_pre','online','offline_post'))
    } else if (strcmp(toalign,'clock')){
      edf$online <- factor(edf$online,levels=c('online_pre','offline_pre','online','offline_post'))
    }
    if (strcmp(toprocess,'symmetry_group') | strcmp(toprocess,"symmetry_group_by_rewFunc") | strcmp(toprocess,'symmetry_group_by_outcome')){
      edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9/10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
      edf$symmetry_group2 <- factor(edf$symmetry_group1, levels = c("14rc11m","11/47",'rl9/10',"fp10","24/32","14m25/32","d10","11/13"))
    }
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      # plot stream gradients
      if (strcmp(toprocess,"network")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'_',model_iter,  ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
          gg<-ggplot(edf,aes(x=online,y=estimate,color=network2,fill=network2)) + 
            facet_wrap(~network2) +
            geom_bar(aes(alpha = p_level_fdr),stat='identity') + 
            geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) +
            xlab(epoch_label) + ylab('estimate') +
            scale_color_manual(values = pal1,labels=c('CTR','DMN','LIM')) + guides(color='none') +
            scale_fill_manual(values = pal1,labels=c('CTR','DMN','LIM')) + guides(color='none') +
            #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
            theme_bw(base_size=13) +
            #facet_wrap(~HC_region) +
            theme(legend.title = element_blank(),
                  panel.grid.major = element_line(colour = "grey45"), 
                  panel.grid.minor = element_line(colour = "grey45"), 
                  panel.background = element_rect(fill = 'grey40'),
                  axis.title.y = element_text(margin=margin(r=6)),
                  axis.title.x = element_text(margin=margin(t=6))) + 
            theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
        print(gg)
        dev.off()  
      } else if (strcmp(toprocess,"symmetry_group")){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr, '_',model_iter, ".pdf", sep = "")
        pdf(fname, width = 12, height = 3.5)
        gg1<-ggplot(edf,aes(x=online,y=estimate,color=sym_fill1,fill=sym_fill1)) + 
          facet_wrap(~symmetry_group1) +
          geom_bar(aes(alpha = p_level_fdr),stat='identity') + 
          geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0.2) +
          scale_color_manual(values = sym_fill_rearranged,guide="none") + 
          scale_fill_manual(values = sym_fill_rearranged,guide="none") + 
          geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) + facet_grid(~symmetry_group1) + 
          xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("estimate") + 
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
        
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
}