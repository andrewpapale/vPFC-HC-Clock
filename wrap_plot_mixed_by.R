# 2022-04-06 AndyP
# Wrapper script to plot vmPFC data

library(pracma)
library(tidyverse)

do_vPFC_fb = FALSE
do_vPFC_clock = FALSE
do_HC_fb = FALSE
do_HC_clock = FALSE
do_HC2vPFC_fb = FALSE
do_HC2vPFC_clock = TRUE
do_anat_fb = FALSE
do_anat_clock = FALSE
do_symmetry = FALSE
do_network = TRUE
do_Explore = TRUE

if (do_vPFC_fb){
  if (do_network){
    source('~/vmPFC/plot_mixed_by_vmPFC.R')
    source('~/vmPFC/plot_emmeans_vmPFC.R')
    for (i in 1:2){
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      model_str <- paste0('-vmPFC-network-feedback-',i,'.Rdata')
      model_str <- Sys.glob(paste0('*',model_str))
      load(model_str)
      model_iter <- i
      totest <- 'final-model'
      toprocess <- 'network'
      toalign <- 'feedback'
      behavmodel <- 'compressed'
      plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
      #plot_emmeans_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
    }
  }
  # source('~/vmPFC/plot_mixed_by_vmPFC.R')
  # for (i in 1){
  #   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  #   model_str <- paste0('-vmPFC-HC-feedback-network-outcome-',i,'.Rdata')
  #   model_str <- Sys.glob(paste0('*',model_str))
  #   load(model_str)
  #   model_iter <- i
  #   totest <- 'outcome'
  #   toprocess <- 'network'
  #   toalign <- 'feedback'
  #   behavmodel <- 'compressed'
  #   plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  # }
  if (do_symmetry){
    source('~/vmPFC/plot_mixed_by_vmPFC.R')
    for (i in 1:2){
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      model_str <- paste0('-vmPFC-symmetry-feedback-',i,'.Rdata')
      model_str <- Sys.glob(paste0('*',model_str))
      load(model_str)
      model_iter <- i
      totest <- 'final-model'
      toprocess <- 'symmetry_group'
      toalign <- 'feedback'
      behavmodel <- 'compressed'
      plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
    }
  }
}
if (do_vPFC_clock){
  if (do_network){
    source('~/vmPFC/plot_mixed_by_vmPFC.R')
    source('~/vmPFC/plot_emmeans_vmPFC.R')
    for (i in 1:2){
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      model_str <- paste0('-vmPFC-network-clock-',i,'.Rdata')
      model_str <- Sys.glob(paste0('*',model_str))
      load(model_str)
      model_iter <- i
      totest <- 'final-model'
      toprocess <- 'network'
      toalign <- 'clock'
      behavmodel <- 'compressed'
      plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
      #plot_emmeans_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
    }
  }
  if (do_symmetry){
    source('~/vmPFC/plot_mixed_by_vmPFC.R')
    for (i in 1:2){
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      model_str <- paste0('-vmPFC-symmetry-clock-',i,'.Rdata')
      model_str <- Sys.glob(paste0('*',model_str))
      load(model_str)
      model_iter <- i
      totest <- 'final-model'
      toprocess <- 'symmetry_group'
      toalign <- 'clock'
      behavmodel <- 'compressed'
      plot_mixed_by_vmPFC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
    }
  }
}
if (do_HC_fb){
  source('~/vmPFC/plot_mixed_by_HC.R')
  for (i in 1:2){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    if (do_Explore){
      model_str <- paste0('-Explore-HC-axis-feedback-',i,'.Rdata')
      totest <- 'Explore-'
    } else {
      model_str <- paste0('-HC-axis-feedback-',i,'.Rdata')
      totest <- 'testing-'
    }
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    model_iter <- i
    toprocess <- 'axis'
    toalign <- 'feedback'
    behavmodel <- 'compressed'
    plot_mixed_by_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  }
}
if (do_HC_clock){
  source('~/vmPFC/plot_mixed_by_HC.R')
  for (i in 1:2){
    setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    if (do_Explore){
      model_str <- paste0('-Explore-HC-region-clock-',i,'.Rdata')
      totest <- 'Explore-'
    } else {
      model_str <- paste0('-HC-axis-clock-',i,'.Rdata')
      totest <- 'testing-'
    }
    model_str <- Sys.glob(paste0('*',model_str))
    load(model_str)
    model_iter <- i
    toprocess <- 'region'
    toalign <- 'clock'
    behavmodel <- 'compressed'
    plot_mixed_by_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter)
  }
}
if (do_HC2vPFC_fb){
  if (do_network){
    source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
    source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
    for (i in 1:2){
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
      model_str <- paste0('-vmPFC-HC-network-feedback-',i,'.Rdata')
      model_str <- Sys.glob(paste0('*',model_str))
      load(model_str)
      model_iter <- i
      totest <- 'final-model'
      #toprocess <- 'symmetry-by-HC'
      toprocess <- 'network-by-HC'
      toalign <- 'feedback'
      behavmodel <- 'compressed'
      hc_LorR <- 'LR'
      plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
      #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    }
  }
  if (do_symmetry){
    source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
    source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
    for (i in 1:2){
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      model_str <- paste0('-vmPFC-HC-symmetry-feedback-',i,'.Rdata')
      #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
      model_str <- Sys.glob(paste0('*',model_str))
      load(model_str)
      model_iter <- i
      totest <- 'final-model'
      toprocess <- 'symmetry-by-HC'
      #toprocess <- 'network-by-HC'
      toalign <- 'feedback'
      behavmodel <- 'compressed'
      hc_LorR <- 'LR'
      plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
      #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    }
  }
}
if (do_HC2vPFC_clock){
  if (do_network){
    source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
    source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
    source('~/vmPFC/plot_emmeans_vmPFC_HC.R')
    source('~/vmPFC/plot_emmeans_vmPFC_HC2.R')
    for (i in 1){
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      #model_str <- paste0('-vmPFC-HC-full-symmetry-',i,'.Rdata')
      if (do_Explore){
        model_str <- paste0('-Explore-vPFC-HC-network-clock-rewFunc-',i,'.Rdata')
        totest <- 'Explore-'
      } else {
        model_str <- paste0('-vmPFC-HC-network-clock-rewFunc-',i,'.Rdata')
        totest <- 'testing-'
      }
      model_str <- Sys.glob(paste0('*',model_str))
      load(model_str)
      model_iter <- i
      #toprocess <- 'symmetry-by-HC'
      toprocess <- 'network-by-HC-by-rewFunc'
      toalign <- 'clock'
      behavmodel <- 'rewFunc'
      hc_LorR <- 'LR'
      plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
      #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
      #plot_emmeans_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
      #plot_emmeans_vmPFC_HC2(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    }
  }
  if (do_symmetry){
    source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
    source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
    for (i in 1:2){
      setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      model_str <- paste0('-vmPFC-HC-symmetry-clock-',i,'.Rdata')
      #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
      model_str <- Sys.glob(paste0('*',model_str))
      load(model_str)
      model_iter <- i
      totest <- 'final-model'
      toprocess <- 'symmetry-by-HC'
      #toprocess <- 'network-by-HC'
      toalign <- 'clock'
      behavmodel <- 'compressed'
      hc_LorR <- 'LR'
      plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
      #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
    }
  }
}


# nice HCwithin plot
setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
if (do_Explore){
  model_str <- paste0('-Explore-vPFC-HC-network-clock-',1,'.Rdata')
} else {
  model_str <- paste0('-vmPFC-HC-network-clock-',1,'.Rdata')
}
model_str <- Sys.glob(paste0('*',model_str))
load(model_str)

qdf <- ddf$coef_df_reml
qdf <- qdf %>% filter(effect=='fixed' & term=='HCwithin')
if (do_Explore){
  qdf <- qdf %>% mutate(network1 = network)
} else {
  qdf <- qdf %>% mutate(network1 = case_when(network=='C'~'CTR', network=='D'~'DMN',network=='L'~'LIM'))
}
pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal = palette()
pal[1] = pal3[2]
pal[2] = pal3[1]
pal[3] = pal3[3]
if (do_Explore){
  pdf('Nice-HCwihin-plot-Explore.pdf',height=3.5,width=6)
} else {
  pdf('Nice-HCwihin-plot.pdf',height=3.5,width=6)
}
  gg1 <- ggplot(qdf,aes(x=HC_region,y=estimate,color=network1)) + geom_violin() + 
  scale_color_manual(values = pal) + facet_wrap(~network1) +  
  ylab('Network Response') +
  xlab('') + guides(color='none') + 
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
print(gg1)
dev.off()


library(ggpattern)
library(grid)

emt <- ddf$emtrends_list$H_HC
emt <- emt %>% filter(network!='L')
emt <- emt %>% mutate(entropy_levels = case_when(v_entropy_sc=='-1.5' ~ 'low entropy', v_entropy_sc=='1.5' ~ 'high entropy'))
emt <- emt %>% mutate(network1 = case_when(network=='C'~'CTR', network=='D'~'DMN',network=='L'~'LIM'))
emt$entropy_levels <- factor(emt$entropy_levels,levels=c('high entropy','low entropy'))
pdf('Nice-HC-vPFC-Entropy-emt.pdf',height=4,width=7)
gg1 <- ggplot(emt, aes(x=HC_region,y=HCwithin.trend,lty=entropy_levels,color=network1)) + 
  geom_violin_pattern(width=1.25,aes(color=network1,pattern_spacing=entropy_levels)) + scale_color_manual(values = pal) +
  facet_wrap(~network1) + guides(color='none') + scale_pattern_spacing_discrete(range=c(0.05,1)) + 
  ylab('HC-vPFC coupling') +
  xlab('Hippocampus region') + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        strip.text.x = element_text(size=30),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=18),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22)
  )
gg2 <- ggplot_gtable(ggplot_build(gg1))
stripr <- which(grepl('strip-t', gg2$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal[k]
  k <- k+1
}
grid.draw(gg2)
dev.off()

library(wesanderson)

pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
pal2 = palette()
pal2[1] = pal3[2]
pal2[2] = pal3[1]
pal2[3] = pal3[3]

emt <- ddf$emtrends_list$V_HC
#emt <- emt %>% filter(network!='D')
emt <- emt %>% filter(v_max_wi!=5.0)
emt <- emt %>% mutate(value_levels = case_when(v_max_wi==-1.5 ~ 'low', v_max_wi==1.0 ~ 'high'))
emt <- emt %>% mutate(network1 = case_when(network=='C'~'CTR', network=='D'~'DMN',network=='L'~'LIM'))
emt$value_levels <- factor(emt$value_levels,levels=c('high','low'))
emt <- emt %>% filter(HC_region=='AH')
pdf('Nice-HC-vPFC-Value-emt',height=4,width=12)
gg1 <-ggplot(emt, aes(x=value_levels,y=HCwithin.trend,color=network1)) + 
  geom_violin(width=1.25,aes(color=network1)) + scale_color_manual(values = pal2) +
  facet_wrap(~network1) + guides(color='none') +
  ylab('HC-vPFC coupling') +
  xlab('Value Level') + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        strip.text.x = element_text(size=30),
        axis.title.y = element_text(margin=margin(r=6),size=22),
        axis.title.x = element_text(margin=margin(t=6),size=22),
        legend.text = element_text(size=18),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22)
  )
gg2 <- ggplot_gtable(ggplot_build(gg1))
stripr <- which(grepl('strip-t', gg2$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal2[k]
  k <- k+1
}
grid.draw(gg2)
dev.off()

#model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')

# source('~/vmPFC/plot_mixed_by_vmPFC_HC_entropy_split.R')
# source('~/vmPFC/plot_emtrends_vmPFC_HC_entropy_split.R')
# for (i in 1){
#   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   model_str <- paste0('-vmPFC-HC-network-clock-Hsplit-',i,'.Rdata')
#   #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
#   model_str <- Sys.glob(paste0('*',model_str))
#   load(model_str)
#   model_iter <- i
#   totest <- 'final-model'
#   toprocess <- 'network-by-HC'
#   #toprocess <- 'network-by-HC'
#   toalign <- 'clock'
#   behavmodel <- 'compressed'
#   hc_LorR <- 'LR'
#   plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
#   #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
# }

# # vmPFC-HC-network-testing
# source('~/vmPFC/plot_mixed_by_vmPFC_HC.R')
# source('~/vmPFC/plot_emtrends_vmPFC_HC.R')
# for (i in 1:2){
#   setwd('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
#   model_str <- paste0('-vmPFC-HC-network-clock-ranslopes-',i,'.Rdata')
#   #model_str <- paste0('-vmPFC-HC-full-network-',i,'.Rdata')
#   model_str <- Sys.glob(paste0('*',model_str))
#   load(model_str)
#   model_iter <- i
#   totest <- 'final-model'
#   toprocess <- 'network-by-HC'
#   #toprocess <- 'network-by-HC'
#   toalign <- 'clock'
#   behavmodel <- 'compressed'
#   hc_LorR <- 'LR'
#   plot_mixed_by_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
#   #plot_emtrends_vmPFC_HC(ddf,toalign,toprocess,totest,behavmodel,model_iter,hc_LorR)
# }

