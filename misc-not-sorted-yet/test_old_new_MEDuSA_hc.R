# 2023-04-14 AndyP
# Check MEDuSA hippocampus old vs new decon

library(tidyverse)

#old_hc_l <- read_csv('/Volumes/Users/Andrew/MMClock_MEDuSA_old/HC_left_feedback_long_decon_locked.csv.gz') %>% mutate(side = 'L')
#old_hc_r <- read_csv('/Volumes/Users/Andrew/MMClock_MEDuSA_old/HC_right_feedback_long_decon_locked.csv.gz') %>% mutate(side = 'R')
#old_hc <- rbind(old_hc_l,old_hc_r) %>% arrange(id,run,trial,atlas_value,evt_time,side)
#old_hc <- old_hc %>% select(!decon_median & !decon_sd & !atlas)
#old_hc <- old_hc %>% rename(bin_num = atlas_value)
#rm(old_hc_l,old_hc_r)
#gc()

old_hc <- read_csv('/Volumes/Users/Andrew/MMClock_MEDuSA_old/long_axis_l_2.3mm_feedback_long_decon_locked.csv')
uB <- sort(unique(old_hc$axis_bin))
old_hc <- old_hc %>% mutate(bin_num = case_when(axis_bin==uB[1]~1,
                                                axis_bin==uB[2]~2,
                                                axis_bin==uB[3]~3,
                                                axis_bin==uB[4]~4,
                                                axis_bin==uB[5]~5,
                                                axis_bin==uB[6]~6,
                                                axis_bin==uB[7]~7,
                                                axis_bin==uB[8]~8,
                                                axis_bin==uB[9]~9,
                                                axis_bin==uB[10]~10,
                                                axis_bin==uB[11]~11,
                                                axis_bin==uB[12]~12))
old_hc <- old_hc %>% mutate(side='L')

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_fb.Rdata') # loads hc
hc <- hc %>% arrange(id,run,trial,bin_num,evt_time,side) %>% 
  select(!atlas_value & !HC_region & !atlas_value0) %>% 
  mutate(side1 = case_when(side=='l'~'L',side=='r'~'R')) %>% select(!side) %>% rename(side=side1)

old_hc <- old_hc %>% filter(evt_time > -6 & evt_time < 10) %>% mutate(dataset = 'old')
hc <- hc %>% filter(evt_time > -6 & evt_time < 10) %>% mutate(dataset = 'new')

old_hc <- old_hc %>% rename(old_decon_mean = decon_interp)
hc <- hc %>% rename(new_decon_mean = decon_mean)
gc()

hc_merged <- inner_join(hc,old_hc,by=c('id','run','run_trial','side','evt_time','bin_num'))
rm(hc, old_hc)
gc()

hc_merged <- hc_merged %>% group_by(id,run,trial,side,evt_time,bin_num) %>% mutate(diff_new_old = new_decon_mean-old_decon_mean) %>% ungroup()
hist(hc_merged$diff_new_old,100)
pdf('test_diff_new_old_hc.pdf',height=12,width=30)
gg1 <- ggplot(hc_merged, aes(x=diff_new_old)) + geom_histogram(bins=100) + facet_grid(evt_time~bin_num)
print(gg1)
dev.off()

corr_val <- hc_merged %>% group_by(evt_time,bin_num) %>% summarize(corr_val = cor(new_decon_mean,old_decon_mean,use='complete.obs')) %>% ungroup()
corr_val <- corr_val %>% pivot_wider(names_from = evt_time, values_from = corr_val)
corr_val <- corr_val %>% select(!bin_num)


raw <- hc_merged %>% group_by(bin_num,side,evt_time) %>% 
  summarize(mD_old = mean(old_decon_mean,na.rm=TRUE),mD_new = mean(new_decon_mean,na.rm=TRUE)) %>% ungroup()
raw <- raw %>% pivot_longer(cols=c(mD_old,mD_new),names_to='dataset')
