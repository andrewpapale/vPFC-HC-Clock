# 2022-06-07 AndyP
# plot Boxplot PH-AH

setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
load('2022-05-16-vmPFC-HC-network-clock-1.Rdata')
ddq <- ddf$coef_df_reml
edf <- ddq %>% filter(term=='HCwithin' & effect=='fixed')
edf <- edf %>% group_by(evt_time,network) %>% complete(HC_region=c('AH','PH')) %>% mutate(dH = estimate[HC_region=='PH'] - estimate[HC_region=='AH'])
edf <- edf %>% filter(HC_region=='AH')
edf <- edf %>% mutate(network1 = case_when(network=='C' ~ 'CTR', network=='D' ~ 'DMN', network=='L' ~ 'LIM'))
edf$network1 <- factor(edf$network1,levels=c('DMN','CTR','LIM'),ordered=TRUE)
pdf('test-boxplot.pdf',height=9,width=12)
gg1 <- ggplot(edf, aes(x=network1,y=dH,color=network1)) + 
  geom_boxplot(notch=FALSE,outlier.color='black',lwd=4) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
  ylab('PH - AH') + theme(legend.position="none") + xlab('network') + theme(text = element_text(size=40))
print(gg1)
dev.off()