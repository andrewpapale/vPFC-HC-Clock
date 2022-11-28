plot_mixed_by_HC <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter){
  
  ## Check plots
  if (strcmp(toalign,"clock")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC')
  } else if (strcmp(toalign,"feedback")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC')
  }
  
  message("\nPlotting streams decoding")
  library(viridis)
  
  epoch_label = paste("Time relative to",toalign, "[s]")
  ddq <- ddf
  ddf <- as_tibble(ddf$coef_df_reml)
  ddf$t <- ddf$evt_time
  #ddf$bin_num <- as.factor(ddf$bin_num)
  ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term, 
                         p_level_fdr = as.factor(case_when(
                           # p_fdr > .1 ~ '0',
                           # p_fdr < .1 & p_fdr > .05 ~ '1',
                           p_fdr > .05 ~ '1',
                           p_fdr < .05 & p_fdr > .01 ~ '2',
                           p_fdr < .01 & p_fdr > .001 ~ '3',
                           p_fdr <.001 ~ '4')))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$pfdr = ddf$p_level_fdr
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  
  if (strcmp(toprocess,'axis')){
    ddf <- ddf %>% mutate(bin_num1 = case_when(
      bin_num==1 ~ 12,
      bin_num==2 ~ 11,
      bin_num==3 ~ 10,
      bin_num==4 ~ 9,
      bin_num==5 ~ 8,
      bin_num==6 ~ 7,
      bin_num==7 ~ 6,
      bin_num==8 ~ 5,
      bin_num==9 ~ 4,
      bin_num==10 ~ 3,
      bin_num==11 ~ 2,
      bin_num==12 ~ 1
    ))
  }
  library(wesanderson)
  for (fe in terms) {
    if (strcmp(toprocess,'axis')){
      pal = wes_palette("Zissou1", 12, type = "continuous")
    } else if (strcmp(toprocess,'region')){
      pal = wes_palette("Zissou1", 2, type = "continuous")
    }
    edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed' & t < 8)
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    
    if (strcmp(toprocess,'region')){
      if (all(edf$p_level_fdr=='p < .001')){
        fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
          # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
        #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        pal1 <- palette(c('#F21A00','#3B9AB2'))
        gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") +
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          #scale_color_gradientn(colors = pal, guide = 'none') + 
          scale_color_manual(values = pal1,labels=c('AH','PH')) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
        print(gg)
        dev.off()
      } else {
        fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
          # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
          geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
          labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        print(gg)
        dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
        #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        pal1 <- palette(c('#F21A00','#3B9AB2'))
        gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
          # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
          geom_line(size = 1) + theme(legend.position = "none") +
          geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
          #scale_color_gradientn(colors = pal, guide = 'none') + 
          scale_color_manual(values = pal1,labels=c('AH','PH')) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          theme(legend.title = element_blank(),
                panel.grid.major = element_line(colour = "grey45"), 
                panel.grid.minor = element_line(colour = "grey45"), 
                panel.background = element_rect(fill = 'grey40'),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
        print(gg)
        dev.off()
      }
    } else if (strcmp(toprocess,'axis')){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      print(ggplot(edf, aes(t, as.factor(bin_num1))) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
              # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      dev.off()
      
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
      #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
      #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      gg<-ggplot(edf, aes(x=t, y=estimate,color = as.numeric(bin_num1), group = as.numeric(bin_num1))) + 
        geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
        # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
        geom_line(size = 1) + theme(legend.position = "none") +
        geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
        scale_color_gradientn(colors = pal, guide = 'none') + 
        #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
        theme_bw(base_size=13) +
        theme(legend.title = element_blank(),
              panel.grid.major = element_line(colour = "grey45"), 
              panel.grid.minor = element_line(colour = "grey45"), 
              panel.background = element_rect(fill = 'grey40'),
              axis.title.y = element_text(margin=margin(r=6)),
              axis.title.x = element_text(margin=margin(t=6)))
      print(gg)
      dev.off()
    }
  }
}  
