# plot_beta_emt
# 2022-08-01 AndyP
session = "meg"
model = "slo"
clock_folder <- "~/clock_analysis" 
# source('~/code/Rhelpers/')
fmri_dir <- '~/vmPFC/MEDUSA Schaefer Analysis/parcel_maps_l2/'
setwd(fmri_dir)

library(grid)
library(readxl)
library(ggrepel)

if (session == "meg") {
  temp <- readRDS("/Volumes/Users/Andrew/parcel_maps_l2/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2_v_entropy_wi_meg_mixed_by.rds")  
} else if (session == "fmri") {
  temp <- readRDS("/Volumes/Users/Andrew/parcel_maps_l2/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2_v_entropy_wi_fmri_mixed_by.rds")
}

library(wesanderson)
pal = wes_palette("FantasticFox1", 3, type = "discrete")
pal_vPFC = palette()
pal_vPFC[2] = '#dd8d29'
pal_vPFC[3] = '#e2d200'
pal_vPFC[1] = '#46acc8'


df <- temp$coef_df_reml

if (model=='slo'){
  emt <- temp$emtrends_list$rt_lag_slo
  df <- df %>% filter(model_name=='slo' & l1_cope_name=='EV_entropy_wiz_clock' & mask_value %in% c(67,171,65,170,66,89,194,88,192,84,191,86,161,55,160,159,56))
} else if (model=='int'){
  emt <- temp$emtrends_list$rt_lag_int
  df <- df %>% filter(model_name=='int' & l1_cope_name=='EV_entropy_wiz_clock' & mask_value %in% c(67,171,65,170,66,89,194,88,192,84,191,86,161,55,160,159,56))
}
emt <- emt %>% filter(mask_value %in% c(67,171,65,170,66,89,194,88,192,84,191,86,161,55,160,159,56))
label_vPFC <- read_excel('~/vmPFC/PFC_region_gradients.xlsx')
label_vPFC <- label_vPFC %>% select(atlas_value, region,network) %>% mutate(mask_value=as.factor(atlas_value))
emt$mask_value <- as.factor(emt$mask_value)
emt <- inner_join(emt,label_vPFC, by='mask_value')
emt <- emt %>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
emt <- emt %>% filter(l1_cope_name=='EV_entropy_wiz_clock' & fmri_beta!=0)

df0 <- df %>% filter(term=='rt_lag:fmri_beta' | term=='rt_lag:fmri_beta:last_outcomeReward')
df0 <- df0 %>% mutate(last_outcome = case_when(term=='rt_lag:fmri_beta' ~ 'Omission', 
                                               term=='rt_lag:fmri_beta:last_outcomeReward' ~ 'Reward'))

df0$mask_value <- as.factor(df0$mask_value)
df0 <- inner_join(df0,label_vPFC, by='mask_value')

Q <- inner_join(emt,df0,by=c('atlas_value','last_outcome'))

Q <- Q  %>% group_by(network1) %>% mutate(padj_BY_term = p.adjust(p.value.y, method = 'bonferroni')) %>% ungroup() 

Q <- Q %>% 
  mutate(p_level_fdr = as.factor(case_when(
    padj_BY_term >= .05 ~ '1',
    padj_BY_term < .05 ~ '2' # & padj_BY_term > .01 ~ '2',
    # padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
    # padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
    # padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
    # padj_BY_term <.00001 ~ '6'
  )))
#Q$p_level_fdr <- factor(Q$p_level_fdr, levels = c('1', '2', '3', '4', '5', '6'), labels = c("NS","p < .05", "p < .01", "p < .001", "p < .0001", "p < .00001")) 
Q$p_level_fdr <- factor(Q$p_level_fdr,levels = c('1','2'),labels=c("NS","p < 0.05"))


pdf(paste0(model,'-',session,'-','Exploration-b2b.pdf'),height=12,width=12)
gg1 <- ggplot(Q, aes(x=fmri_beta,y=rt_lag.trend,color=statistic.y,alpha=p_level_fdr)) + 
  geom_violin(aes(fmri_beta, rt_lag.trend,group=fmri_beta), alpha = .2) + 
  geom_jitter(width=0.1,height=0,aes(size=p_level_fdr)) + facet_grid(last_outcome~network1,scales='free_y') + 
  geom_line(aes(group=region.x)) +
  xlab('entropy beta') +
  ylab('RT Convergence RT(N-1) to RT(N)') +
  theme(axis.text = element_text(size=15), axis.title=element_text(size=30), 
        legend.text=element_text(size=15),legend.title=element_text(size=30),
        strip.text.x = element_text(size=30),strip.text.y=element_text(size=30)) +
  geom_text_repel(aes(fmri_beta, rt_lag.trend, alpha = p_level_fdr, label=region.x), 
                  point.padding = 8, force = 3, color="#4FC3F7", size = 7,max.overlaps=10); 

gg2 <- ggplot_gtable(ggplot_build(gg1))
stripr <- which(grepl('strip-t', gg2$layout$name))
k <- 3
for (i in stripr) {
  j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_vPFC[k]
  k <- k-1
}
grid.draw(gg2)

dev.off()


if (model=='slo'){
  emt <- temp$emtrends_list$rt_vmax_slo
  df <- df %>% filter(model_name=='slo' & l1_cope_name=='EV_entropy_wiz_clock' & mask_value %in% c(67,171,65,170,66,89,194,88,192,84,191,86,161,55,160,159,56))
} else if (model=='int'){
  emt <- temp$emtrends_list$rt_vmax_int
  df <- df %>% filter(model_name=='int' & l1_cope_name=='EV_entropy_wiz_clock' & mask_value %in% c(67,171,65,170,66,89,194,88,192,84,191,86,161,55,160,159,56))
}
emt <- emt %>% filter(mask_value %in% c(67,171,65,170,66,89,194,88,192,84,191,86,161,55,160,159,56) & last_outcome=='Reward')
label_vPFC <- read_excel('~/vmPFC/PFC_region_gradients.xlsx')
label_vPFC <- label_vPFC %>% select(atlas_value, region,network) %>% mutate(mask_value=as.factor(atlas_value))
emt$mask_value <- as.factor(emt$mask_value)
emt <- inner_join(emt,label_vPFC, by='mask_value')
emt <- emt %>% mutate(network1 = case_when(network=='D'~'DMN', network=='C'~'CTR',network=='L'~'LIM'))
emt <- emt %>% filter(l1_cope_name=='EV_entropy_wiz_clock' & fmri_beta!=0)

df0 <- df %>% filter(term=="fmri_beta:rt_vmax_lag")

df0$mask_value <- as.factor(df0$mask_value)
df0 <- inner_join(df0,label_vPFC, by='mask_value')

Q <- inner_join(emt,df0,by=c('atlas_value'))

Q <- Q  %>% group_by(network1) %>% mutate(padj_BY_term = p.adjust(p.value.y, method = 'bonferroni')) %>% ungroup() %>% 
  mutate(p_level_fdr = as.factor(case_when(
    padj_BY_term > .05 ~ '1',
    padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
    padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
    padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
    padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
    padj_BY_term <.00001 ~ '6'
  )))
Q$p_level_fdr <- factor(Q$p_level_fdr, levels = c('1', '2', '3', '4', '5', '6'), labels = c("NS","p < .05", "p < .01", "p < .001", "p < .0001", "p < .00001")) 



pdf(paste0(model,'-',session,'-','Exploitation-b2b.pdf'),height=12,width=12)
gg1 <- ggplot(Q, aes(x=fmri_beta,y=rt_vmax_lag.trend,color=statistic.y,alpha=p_level_fdr)) + 
  geom_violin(aes(fmri_beta, rt_vmax_lag.trend,group=fmri_beta), alpha = .2) + 
  geom_jitter(width=0.1,height=0,aes(size=p_level_fdr)) + facet_wrap(~network1,scales='free_y') + 
  geom_line(aes(group=region.x)) +
  xlab('entropy beta') +
  ylab('Convergence on RT(Vmax)') +
  theme(axis.text = element_text(size=15), axis.title=element_text(size=30), 
        legend.text=element_text(size=15),legend.title=element_text(size=30),
        strip.text.x = element_text(size=30),strip.text.y=element_text(size=30)) +
  geom_text_repel(aes(fmri_beta, rt_vmax_lag.trend, alpha = p_level_fdr, label=region.x), 
                  point.padding = 8, force = 3, color="#4FC3F7", size = 7,max.overlaps=10); 

gg2 <- ggplot_gtable(ggplot_build(gg1))
stripr <- which(grepl('strip-t', gg2$layout$name))
k <- 3
for (i in stripr) {
  j <- which(grepl('rect', gg2$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_vPFC[k]
  k <- k-1
}
grid.draw(gg2)
dev.off()
