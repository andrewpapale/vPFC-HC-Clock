plot_mixed_by_simple <- function(ddf){
  
  epoch_label <- 'time peri feedback (s)'
  
  message("\nPlotting streams decoding")
  library(wesanderson)
  library(pracma)
  library(ggplot2)
  library(tidyverse)
  library(viridis)
  library(grid)
  ddf <- ddf$coef_df_reml
  
  ddf$atlas_value <- as.factor(ddf$atlas_value)
  ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
    padj_fdr_term > .05 ~ '1',
    padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
    padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
    padj_fdr_term <.001 ~ '4')))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  
  ddf$t <- ddf$evt_time
  ddf$t <- as.numeric(ddf$t)
  for (fe in terms) {
    # fe <- terms[1] # test only
    edf <- ddf %>% filter(term == paste(fe) & ddf$effect=='fixed')
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste0(termstr,'-',".pdf")
    pdf(fname, width = 9, height = 3.5)
    # gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=network1, size=`p, FDR-corrected`)) +
    #   geom_line(size = 1) + geom_point() +
    #   geom_errorbar() +
    #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~HC_region) +
    #   geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
    #   scale_color_manual(labels=c('DMN','CTR','LIM'),values=c('red','green','blue')) + xlab(epoch_label) +
    #   labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
    if (all(edf$`p, FDR-corrected`=='p < .001')){
      gg<-ggplot(edf, aes(x=t, y=estimate)) + 
        geom_point(aes(size=`p, FDR-corrected`, alpha = `p, FDR-corrected`,color=atlas_value)) + 
        # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
        geom_line(size = 1,aes(group=atlas_value)) + theme(legend.position = "none") + scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
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
      gg<-ggplot(edf, aes(x=t, y=estimate,color=atlas_value)) + 
        geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab(epoch_label) + ylab('') +
        geom_hline(yintercept = 0, lty = 'dashed',color = 'gray', size=1) + 
        geom_point(aes(size=`p, FDR-corrected`, alpha = `p, FDR-corrected`)) +
        # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
        geom_line(size = 1,aes(group=atlas_value)) + theme(legend.position = "none") +
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
  }
} 
