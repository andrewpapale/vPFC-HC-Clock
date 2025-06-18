# 2023-05-30 AndyP
# Plot t-stats for HCwithin, see if DMN > CTR > LIM order holds

library(wesanderson)
library(ggplot2)
library(grid)

nT <- 3
plot_dir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction'
paper_version <- 'v11'
doExplore = FALSE

for (iT in 1:3){
  if (doExplore){
    load(paste0('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-26-Explore-vPFC-HC-network-clock-HConly-',iT,'.Rdata'))
  } else {
    load(paste0('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-05-25-vmPFC-HC-network-clock-',iT,'.Rdata'))
  }
  ddf <- ddf$coef_df_reml
  hc_w <- ddf %>% filter(term=='HCwithin',effect=='fixed')
  # plot HCwithin
  hc_w <- hc_w %>% mutate(p_level_fdr = as.factor(case_when(
    padj_fdr_term > .05 ~ '1',
    padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
    padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
    padj_fdr_term <.001 ~ '4')))
  hc_w$p_level_fdr <- factor(hc_w$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < 0.05", "p < 0.01", "p < 0.001"))
  
  pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
  pal = palette()
  pal[1] = pal3[2]
  pal[2] = pal3[1]
  pal[3] = pal3[3]
  
  pal1 = palette()
  pal1[1] = '#ff484d'
  pal1[2] = '#2a52ff'
  
  setwd(plot_dir)
  
  if (doExplore){
    title_str <- paste0(paper_version,'-Explore-clock-network-by-HC-by-plot_t_stats-HCwithin-',iT,'.pdf')
  } else {
    title_str <- paste0(paper_version,'-MMClock-clock-network-by-HC-by-plot_t_stats-HCwithin-',iT,'.pdf')
  }
  pdf(title_str,height=8,width=16)
  
  if (all(hc_w$p_level_fdr=='p < 0.001')){
    gg1 <- ggplot(hc_w,aes(x=evt_time,y=statistic)) + 
      geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
      geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
      facet_grid(~HC_region) + 
      geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1) + 
      scale_color_manual(values = pal) +
      scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
      theme_bw(base_size=13) +
      theme(legend.title = element_blank(),
            panel.grid.major = element_line(colour = "grey45"), 
            panel.grid.minor = element_line(colour = "grey45"), 
            panel.background = element_rect(fill = 'grey40'),
            axis.title.y = element_text(margin=margin(r=6),size=22),
            axis.title.x = element_text(margin=margin(t=6),size=22),
            legend.text = element_text(size=22),
            axis.text.x = element_text(size=22),
            axis.text.y = element_text(size=22))
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
  } else {
    gg1 <- ggplot(hc_w,aes(x=evt_time,y=statistic)) + 
      geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
      geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
      facet_grid(~HC_region) + 
      geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1) +
      scale_color_manual(values = pal) +
      theme_bw(base_size=13) +
      theme(legend.title = element_blank(),
            panel.grid.major = element_line(colour = "grey45"), 
            panel.grid.minor = element_line(colour = "grey45"), 
            panel.background = element_rect(fill = 'grey40'),
            axis.title.y = element_text(margin=margin(r=6),size=22),
            axis.title.x = element_text(margin=margin(t=6),size=22),
            legend.text = element_text(size=22),
            axis.text.x = element_text(size=22),
            axis.text.y = element_text(size=22))
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
  
  
}