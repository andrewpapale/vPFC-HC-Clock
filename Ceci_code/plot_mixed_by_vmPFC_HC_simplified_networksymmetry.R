plot_mixed_by_vmPFC_HC_simplified_networksymmetry <- function(ddf,totest,toprocess,behavmodel,toalign,model_iter,term){
  # load required libraries
  library(grid)
  library(pracma)
  library(viridis)
  library(wesanderson)
  
  ## Check if clock- or feedback-aligned, change wd and set string var
  if (strcmp(toalign,"clock")){
    toalign_str <- 'feedback'
  } else if (strcmp(toalign,"feedback")){
    toalign_str <- 'trial onset' 
  }

  # set palette
  pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
  pal = palette()
  pal[1] = pal3[2]
  pal[2] = pal3[1]
  pal[3] = pal3[3]
  pal_hc = wes_palette("Royal2", 5, type = "discrete")
  pal_hc1 <- palette()
  pal_hc1[2] <- '#C9D9F9'
  pal_hc1[1] <- '#818589'
  fills <- palette()
  fills[1] <- pal[2]
  fills[2] <- pal[2]
  fills[3] <- pal[1]
  fills[4] <- pal[1]
  fills[5] <- pal[1]
  fills[6] <- pal[1]
  fills[7] <- pal[3]
  fills[8] <- pal[3]
  
  # get REML coefficients and set up some variables for plotting
  ddg <- ddf # in case we need it
  ddf <- as_tibble(ddf$coef_df_reml)
  ddf$t <- ddf$evt_time
  if (strcmp(toalign,'feedback')){
    ddf <- ddf %>% filter(t > -4)
  }
  if (!all(is.na(ddf$p_level_fdr))){
    ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term, 
                           p_level_fdr = as.factor(case_when(
                             p_fdr > .05 ~ '0',
                             p_fdr < .05 & p_fdr > .01 ~ '0.25',
                             p_fdr < .01 & p_fdr > .001 ~ '0.5',
                             p_fdr <.001 ~ '0.75')))
  } else {
    ddf <- ddf  %>% mutate(p_fdr = p.value, 
                           p_level_fdr = as.factor(case_when(
                             p_fdr > .05 ~ '0',
                             p_fdr < .05 & p_fdr > .01 ~ '0.25',
                             p_fdr < .01 & p_fdr > .001 ~ '0.5',
                             p_fdr <.001 ~ '0.75')))
    warning('p fdr were all nan so deriving p values from p.level')
  }
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('0', '0.25', '0.5', '0.75'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  epoch_label = paste("Time relative to",toalign_str, "[s]")
    # get term we are plotting
    fe = term
    
    print(fe)
    edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed' & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    
  if(toprocess=="network-by-HC"){
    # set networks
    edf <- edf %>% mutate(network1 = 
                            case_when(network=='CTR'~'2C',
                                      network=='DMN'~'1D',
                                      network=='LIM'~'3L')) %>%
      mutate(network2 = case_when(network=='CTR'~'CTR',
                                  network=='DMN'~'DMN',
                                  network=='LIM'~'LIM'))
    
    # do plotting
    if (!all(is.na(edf$`p, FDR-corrected`))){
      if (all(edf$`p, FDR-corrected`=='p < .001')){
        if (fe=='HCwithin'){
          gg1<-ggplot(edf, aes(x=t, y=estimate,group=network2,color=network2))
            gg1 + geom_vline(xintercept = 0, lty = 'dashed', color = '#A9A9A9') +
            geom_point(position=position_dodge(width=0.33),aes(size=p_level_fdr, alpha = p_level_fdr)) + 
            geom_errorbar(position=position_dodge(width=0.33),aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0, color="black") +
            geom_line(position=position_dodge(width=0.33),linewidth = 1) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
            xlab(epoch_label) + ylab('') +
            scale_color_manual(values = pal) + 
            scale_y_reverse() +
            scale_x_continuous(breaks = c(-4,-2,0,2,4)) + 
            theme_bw(base_size=13) + ylab('Response [AU]') + 
            facet_wrap(~HC_region)
        }else {
          gg1<-ggplot(edf, aes(x=t, y=estimate,group=network2,color=network2))
          gg1 +
            geom_hline(yintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) +
            geom_vline(xintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) +
            geom_point(position=position_dodge(width=0.33),aes(size=p_level_fdr, alpha = p_level_fdr)) + 
            geom_errorbar(position=position_dodge(width=0.33),aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0, color="black") +
            geom_line(position=position_dodge(width=0.33),size = 1) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
            xlab(epoch_label) + ylab('') +
            scale_y_reverse() +
            scale_color_manual(values = pal) + 
            scale_x_continuous(breaks = c(-4,-2,0,2,4)) + 
            theme_bw(base_size=13) +  
            facet_wrap(~HC_region) + ylab('Response [AU]')
        }
      } else {
        gg1<-ggplot(edf, aes(x=t, y=estimate,group=network2,color=network2))
        gg1 +
          geom_hline(yintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) +
          geom_vline(xintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) +
          geom_point(position=position_dodge(width=0.33),aes(size=p_level_fdr, alpha = p_level_fdr)) +
          geom_errorbar(position=position_dodge(width=0.33),aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0, color="black") +
          geom_line(position=position_dodge(width=0.33),size = 1) + theme(legend.position = "none") + 
          xlab(epoch_label) + ylab('') +
          scale_color_manual(values = pal) + 
          scale_y_reverse() +
          scale_x_continuous(breaks = c(-4,-2,0,2,4)) + 
          theme_bw(base_size=13) +
          facet_wrap(~HC_region) + ylab('Response [AU]')
  
      }
    }
    
  }
    
    else if(toprocess=="symmetry-by-HC"){
      # get symmetry groups
      edf <- edf %>% mutate(symmetry_group1=case_when(
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
      edf$symmetry_group1 <- factor(edf$symmetry_group1, levels = c('rl9_10',"11/47","d10","24/32","14m25/32","fp10","11/13","14rc11m"))
      edf$symmetry_group2 <- factor(edf$symmetry_group1, levels = c("14rc11m","11/47",'rl9_10',"fp10","24/32","14m25/32","d10","11/13"))
      edf$network_group = factor(edf$network_group,levels = c('CTR','DMN','LIM'))
      
      # do plotting
      gg1 <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=as.factor(HC_region))) 
      gg1 +
        geom_line(size = 1) + geom_point(aes(size= `p, FDR-corrected`)) +
        geom_errorbar() +
        geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 1) + facet_grid(~symmetry_group1) + 
        scale_color_manual(values = pal_hc1,labels=c('Anterior Hippocampus','Posterior Hippocampus')) + xlab(epoch_label) +
        labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
    }
}