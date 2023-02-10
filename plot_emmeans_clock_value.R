# 2023-02-06 AndyP
# Plot emmeans of Entropy for DMN 

library(tidyverse)
library(wesanderson)

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-01-07-vmPFC-network-clock-2.Rdata')

emm <- ddf$emmeans_list$V
emm <- emm %>% filter(network=='D')
emm <- emm %>% mutate(value = case_when(v_max_wi=='-1.5'~'low',
                                          v_max_wi=='1.5'~'high'))

pal = wes_palette("FantasticFox1", 3, type = "discrete")
pal1 <- palette()
pal1[1] = pal[1]
pal1[2] = '#f1d1a9'
  
setwd('/Volumes/Users/Andrew/Fig3_Value_vPFC')

pdf('DMN-value-MMClock.pdf',width=6,height=3.5)
gg1 <- ggplot(emm, aes(x=evt_time,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error,color=value)) + 
  geom_point(size=5) + 
  geom_errorbar(size=1.5) + 
  geom_line(size=1.5) +
  xlab('time relative to trial onset [s]') + 
  scale_color_manual(values = pal1) + 
  theme_bw(base_size=13) +
  ylab('DMN Response') +
  geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1) +
  theme(
    panel.grid.major = element_line(colour = "grey45"), 
    panel.grid.minor = element_line(colour = "grey45"), 
    panel.background = element_rect(fill = 'grey40'),
    axis.title.y = element_text(margin=margin(r=6),size=22),
    axis.title.x = element_text(margin=margin(t=6),size=22),
    legend.text = element_text(size=22),
    legend.title = element_text(size=22),
    axis.text.x = element_text(size=22),
    axis.text.y = element_text(size=22)
  )
print(gg1)
dev.off()




load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-01-07-Explore-vmPFC-network-clock-2.Rdata')

emm <- ddf$emmeans_list$H
emm <- emm %>% filter(network=='DMN')
emm <- emm %>% filter(v_max_wi=='-2' | v_max_wi=='2')
emm <- emm %>% mutate(value = case_when(v_max_wi=='-2'~'low',
                                          v_max_wi=='2'~'high'))

pal = wes_palette("FantasticFox1", 3, type = "discrete")
pal1 <- palette()
pal1[1] = pal[1]
pal1[2] = '#f1d1a9'
  
setwd('/Volumes/Users/Andrew/Fig3_Value_vPFC')

pdf('DMN-value-Explore.pdf',width=6,height=3.5)
gg1 <- ggplot(emm, aes(x=evt_time,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error,color=value)) + 
  geom_point(size=5) + 
  geom_errorbar(size=1) + 
  geom_line(size=1.5) +
  xlab('time relative to trial onset [s]') + 
  scale_color_manual(values = pal1) + 
  theme_bw(base_size=13) +
  ylab('DMN Response') +
  geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1) +
  theme(
    panel.grid.major = element_line(colour = "grey45"), 
    panel.grid.minor = element_line(colour = "grey45"), 
    panel.background = element_rect(fill = 'grey40'),
    axis.title.y = element_text(margin=margin(r=6),size=22),
    axis.title.x = element_text(margin=margin(t=6),size=22),
    legend.text = element_text(size=22),
    legend.title = element_text(size=22),
    axis.text.x = element_text(size=22),
    axis.text.y = element_text(size=22)
  )
print(gg1)
dev.off()
