plot_mixed_by_HC <- function(ddf,toalign,toprocess,totest,behavmodel,model_iter){
  
  ## Check plots
  if (strcmp(toalign,"clock")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC')
  } else if (strcmp(toalign,"feedback")){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC')
  }
  
  message("\nPlotting streams decoding")
  library(viridis)
  
  epoch_label = paste("Time relative to trial onset [s]")
  ddq <- ddf
  ddf <- as_tibble(ddf$coef_df_reml)
  if (strcmp(totest,'Explore-')){
    ddf$bin_num <- ddf$atlas_value
  }
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
  
  
  
  # if (strcmp(toprocess,'axis')){
  #   ddf <- ddf %>% mutate(bin_num1 = case_when(
  #     bin_num==1 ~ 12,
  #     bin_num==2 ~ 11,
  #     bin_num==3 ~ 10,
  #     bin_num==4 ~ 9,
  #     bin_num==5 ~ 8,
  #     bin_num==6 ~ 7,
  #     bin_num==7 ~ 6,
  #     bin_num==8 ~ 5,
  #     bin_num==9 ~ 4,
  #     bin_num==10 ~ 3,
  #     bin_num==11 ~ 2,
  #     bin_num==12 ~ 1
  #   ))
  # }
  if (strcmp(toprocess,'axis')){
    pal = wes_palette("Zissou1", 12, type = "continuous")
  } else if (strcmp(toprocess,'region') | strcmp(toprocess,'region-by-MF')){
    pal = palette()
    pal[2] <- '#C9D9F9'
    pal[1] <- '#818589'
  }
  library(wesanderson)
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed' & t < 8)
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    
    if (strcmp(toprocess,'region')){
      if (all(edf$p_level_fdr=='p < .001')){
        # fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        # pdf(fname, width = 9, height = 3.5)
        # gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
        #   # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        #   scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        # print(gg)
        # dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
        #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        #pal1 <- palette(c('#F21A00','#3B9AB2'))
        gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
          geom_vline(xintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1)+ xlab(epoch_label) +
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0, color="black") +
          geom_line(size = 1) + 
          ylab('Response [AU]') +
          scale_x_continuous(breaks = c(-4,-2,0,2,4)) + 
          #scale_color_gradientn(colors = pal, guide = 'none') + 
          scale_color_manual(values = pal) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          theme(legend.title = element_blank(),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
        print(gg)
        dev.off()
      } else {
        # fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        # pdf(fname, width = 9, height = 3.5)
        # gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
        #   # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        #   scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        # print(gg)
        # dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
        #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        #pal1 <- palette(c('#F21A00','#3B9AB2'))
        gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
          geom_hline(yintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) +
          geom_vline(xintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) + 
          geom_point(position=position_dodge(width=0.33),aes(size=p_level_fdr, alpha = p_level_fdr)) +
          geom_errorbar(position=position_dodge(width=0.33),aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0, color="black") +
          geom_line(position=position_dodge(width=0.33),size = 1) + 
          xlab(epoch_label) + ylab('Response [AU]') +
          scale_x_continuous(breaks = c(-4,-2,0,2,4)) + 
          #scale_color_gradientn(colors = pal, guide = 'none') + 
          scale_color_manual(values = pal) +  
          xlab(epoch_label) + ylab('Response [AU]') +
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          theme(legend.title = element_blank(),
                axis.title.y = element_text(margin=margin(r=6),size=22),
                axis.title.x = element_text(margin=margin(t=6),size=22),
                legend.text = element_text(size=22),
                axis.text.x = element_text(size=22),
                axis.text.y = element_text(size=22)
          )
        print(gg)
        dev.off()
      }
      
      
      
    } else if (strcmp(toprocess,'region-by-MF')){
      if (all(edf$p_level_fdr=='p < .001')){
        # fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        # pdf(fname, width = 9, height = 3.5)
        # gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
        #   # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        #   scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        # print(gg)
        # dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
        #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        #pal1 <- palette(c('#F21A00','#3B9AB2'))
        gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
          geom_vline(xintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1)+ xlab(epoch_label) +
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0, color="black") +
          geom_line(size = 1) + 
          ylab('Response [AU]') +
          scale_x_continuous(breaks = c(-4,-2,0,2,4)) + 
          facet_grid(~female) + 
          #scale_color_gradientn(colors = pal, guide = 'none') + 
          scale_color_manual(values = pal) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          theme(legend.title = element_blank(),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
        print(gg)
        dev.off()
      } else {
        # fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        # pdf(fname, width = 9, height = 3.5)
        # gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
        #   # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        #   scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        # print(gg)
        # dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
        #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        #pal1 <- palette(c('#F21A00','#3B9AB2'))
        gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
          geom_hline(yintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) +
          geom_vline(xintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) + 
          geom_point(position=position_dodge(width=0.33),aes(size=p_level_fdr, alpha = p_level_fdr)) +
          geom_errorbar(position=position_dodge(width=0.33),aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0, color="black") +
          geom_line(position=position_dodge(width=0.33),size = 1) + 
          xlab(epoch_label) + ylab('Response [AU]') +
          scale_x_continuous(breaks = c(-4,-2,0,2,4)) + 
          facet_grid(~female) + 
          #scale_color_gradientn(colors = pal, guide = 'none') + 
          scale_color_manual(values = pal) +  
          xlab(epoch_label) + ylab('Response [AU]') +
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          theme(legend.title = element_blank(),
                axis.title.y = element_text(margin=margin(r=6),size=22),
                axis.title.x = element_text(margin=margin(t=6),size=22),
                legend.text = element_text(size=22),
                axis.text.x = element_text(size=22),
                axis.text.y = element_text(size=22)
          )
        print(gg)
        dev.off()
      }
      
      
    } else if (strcmp(toprocess,'region-by-scanner')){
      if (all(edf$p_level_fdr=='p < .001')){
        # fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        # pdf(fname, width = 9, height = 3.5)
        # gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
        #   # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        #   scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        # print(gg)
        # dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
        #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        #pal1 <- palette(c('#F21A00','#3B9AB2'))
        gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
          geom_vline(xintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1)+ xlab(epoch_label) +
          geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_errorbar(aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0, color="black") +
          geom_line(size = 1) + 
          ylab('Response [AU]') +
          scale_x_continuous(breaks = c(-4,-2,0,2,4)) + 
          facet_grid(~scan_which) + 
          #scale_color_gradientn(colors = pal, guide = 'none') + 
          scale_color_manual(values = pal) + 
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          theme(legend.title = element_blank(),
                axis.title.y = element_text(margin=margin(r=6)),
                axis.title.x = element_text(margin=margin(t=6)))
        print(gg)
        dev.off()
      } else {
        # fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        # pdf(fname, width = 9, height = 3.5)
        # gg<- ggplot(edf, aes(t, HC_region)) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
        #   # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
        #   scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        # print(gg)
        # dev.off()
        
        fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
        pdf(fname, width = 9, height = 3.5)
        # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
        #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
        #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
        #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
        #pal1 <- palette(c('#F21A00','#3B9AB2'))
        gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(HC_region), group = as.factor(HC_region))) + 
          geom_hline(yintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) +
          geom_vline(xintercept = 0, lty = 'dashed', color = '#A9A9A9', size = 1) + 
          geom_point(position=position_dodge(width=0.33),aes(size=p_level_fdr, alpha = p_level_fdr)) +
          geom_errorbar(position=position_dodge(width=0.33),aes(ymin=estimate-std.error,ymax=estimate+std.error),width=0, color="black") +
          geom_line(position=position_dodge(width=0.33),size = 1) + 
          xlab(epoch_label) + ylab('Response [AU]') +
          scale_x_continuous(breaks = c(-4,-2,0,2,4)) + 
          facet_grid(~scan_which) + 
          #scale_color_gradientn(colors = pal, guide = 'none') + 
          scale_color_manual(values = pal) +  
          xlab(epoch_label) + ylab('Response [AU]') +
          #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
          theme_bw(base_size=13) +
          theme(legend.title = element_blank(),
                axis.title.y = element_text(margin=margin(r=6),size=22),
                axis.title.x = element_text(margin=margin(t=6),size=22),
                legend.text = element_text(size=22),
                axis.text.x = element_text(size=22),
                axis.text.y = element_text(size=22)
          )
        print(gg)
        dev.off()
      }
      
      
      
      
      } else if (strcmp(toprocess,'axis')){
      fname = paste(behavmodel,'-',totest,"_",toalign, "_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 5)
      print(ggplot(edf, aes(t, as.factor(bin_num))) + geom_tile(aes(fill = estimate, alpha = pfdr), size = 1) +
              # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab(""))
      dev.off()
      
      fname = paste(behavmodel,'-',totest,"_",toalign, "_line_", toprocess, "_", termstr,'-',model_iter, ".pdf", sep = "")
      pdf(fname, width = 9, height = 5)
      # gg <- ggplot(edf, aes(x=t, y=estimate,color=as.factor(bin_num))) +
      #   geom_line(size = 1) + geom_point(aes(size=pfdr), fill="red") +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #facet_wrap(~visuomotor_grad) +
      #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + xlab(epoch_label) +
      #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      gg<-ggplot(edf, aes(x=t, y=estimate,color = as.factor(bin_num), group = as.factor(bin_num))) + 
        geom_point(aes(size=p_level_fdr, alpha = p_level_fdr)) +
        # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
        geom_line(size = 1) + 
        geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('')
      #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
      print(gg)
      dev.off()
    }
  }
}  
