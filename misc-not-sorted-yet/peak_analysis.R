# 2023-02-13 AndyP
# N peaks and peak prominence entropy relationship

library(tidyverse)

source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory = '~/clock_analysis',dataset='mmclock_fmri')

dP <- read_csv(file=file.path('~/vmPFC/','peaks1.csv'),col_names = FALSE)
dP <- dP %>% rename(peaks=X1,meanpp=X2,maxpp=X3,run=X4,run_trial=X5,id=X6)
dP <- inner_join(dP,df,by=c('id','run','run_trial'))
dP$peaks <- as.numeric(dP$peaks)

setwd('/Volumes/Users/Andrew/Fig_S#_Entropy_Peaks')
pdf('Entropy_peaks.pdf',height=6,width=9)
gg1<-ggplot(dP,aes(x=peaks,y=v_entropy_wi,group=peaks)) +
  theme_bw(base_size=13) +
  geom_hline(yintercept = 0, lty = 'dashed',color = 'gray', size=1) + 
  geom_violin() + geom_boxplot(width=0.1) +
  ylab('Scaled Entropy') +
  theme(legend.title = element_blank(),
      axis.title.y = element_text(margin=margin(r=6),size=48),
      axis.title.x = element_text(margin=margin(t=6),size=48),
      legend.text = element_text(size=24),
      axis.text.x = element_text(size=24),
      axis.text.y = element_text(size=24)
  )
print(gg1)
dev.off()

S <- aov(v_entropy_wi ~ peaks, data=dP)
anova(S)


