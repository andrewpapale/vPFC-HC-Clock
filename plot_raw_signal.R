# 2022-06-29 AndyP
# Plot Raw deconvolved bold vPFC

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)

# start with vmPFC simple, add in term by term, eventually add HC interaction
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'

#######################
### vPFC - feedback ###
#######################

message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'feedback_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- fb_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(fb_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>% 
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
         ))
df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 15 ~ 'Early',
  run_trial > 15 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% group_by(id) %>% mutate(entropy_split = case_when(
  v_entropy < mean(v_entropy,na.rm=TRUE) ~ 'low',
  v_entropy > mean(v_entropy,na.rm=TRUE) ~ 'high'
)) %>% ungroup()
df <- df %>% select(entropy_split,v_max_above_median,id,run,trial_bin,rewFunc,rt_bin,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)

dQ <- Q %>% group_by(id,evt_time,network,entropy_split,trial_bin) %>% summarize(vP = mean(vmPFC_decon,na.rm=TRUE),dP = sd(vmPFC_decon,na.rm=TRUE), N = length(vmPFC_decon))
dQ <- dQ %>% ungroup() %>% group_by(evt_time,network,entropy_split,trial_bin) %>% summarize(vP1 = mean(vP,na.rm=TRUE),dP1 = sd(dP,na.rm=TRUE), N1 = length(N))
dQ <- dQ %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR',network=='L'~'LIM'))
dQ$network1 <- factor(dQ$network1,levels=c('DMN','CTR','LIM'))
dQ <- dQ %>% filter(!is.na(entropy_split))
dQ$trial_bin <- factor(dQ$trial_bin,levels=c('Early','Middle','Late'))
pal = wes_palette("FantasticFox1", 3, type = "discrete")

dQ$entropy <- dQ$entropy_split
dQ <- dQ %>% ungroup() %>% select(!network) %>% rename(network=network1)

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('plot_raw_signal_feedback_ventropy.pdf',height=12,width=24)
gg1 <- ggplot(dQ,aes(x=evt_time,y=vP1,color=network,linetype=entropy)) + 
  geom_point(size=5) + geom_line(size=2) + 
  geom_errorbar(aes(ymin=vP1-(dP1/sqrt(N1)),ymax=vP1+(dP1/sqrt(N1)))) + 
  geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 3) + 
  facet_wrap(~trial_bin) + 
  scale_color_manual(values=pal,labels=c("DMN","CTR","LIM")) +
  xlab('Time relative to feedback (s)') +
  ylab('Deconvolved Signal (AU)') +
  theme(legend.text = element_text(size=25)) +
  theme(legend.title = element_text(size=30)) +
  theme(axis.title = element_text(size=40)) +
  theme(axis.text = element_text(size=25)) +
  theme(strip.text.x = element_text(size=40))
print(gg1)
dev.off()

rm(dQ)
dQ <- Q %>% group_by(id,evt_time,network,v_max_above_median,trial_bin) %>% summarize(vP = mean(vmPFC_decon,na.rm=TRUE),dP = sd(vmPFC_decon,na.rm=TRUE), N = length(vmPFC_decon))
dQ <- dQ %>% ungroup() %>% group_by(evt_time,network,v_max_above_median,trial_bin) %>% summarize(vP1 = mean(vP,na.rm=TRUE),dP1 = sd(dP,na.rm=TRUE), N1 = length(N))
dQ <- dQ %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR',network=='L'~'LIM'))
dQ$network1 <- factor(dQ$network1,levels=c('DMN','CTR','LIM'))
dQ <- dQ %>% filter(!is.na(v_max_above_median))
dQ$trial_bin <- factor(dQ$trial_bin,levels=c('Early','Middle','Late'))
pal = wes_palette("FantasticFox1", 3, type = "discrete")

dQ <- dQ %>% mutate(v_max = case_when(v_max_above_median == TRUE ~ 'high', v_max_above_median == FALSE ~ 'low'))
dQ <- dQ %>% ungroup() %>% select(!network) %>% rename(network=network1)

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('plot_raw_signal_feedback_vmax.pdf',height=12,width=24)
gg1 <- ggplot(dQ,aes(x=evt_time,y=vP1,color=network,linetype=v_max)) + 
  geom_point(size=5) + geom_line(size=2) + 
  geom_errorbar(aes(ymin=vP1-(dP1/sqrt(N1)),ymax=vP1+(dP1/sqrt(N1)))) + 
  geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 3) + 
  facet_wrap(~trial_bin) + 
  scale_color_manual(values=pal,labels=c("DMN","CTR","LIM")) +
  xlab('Time relative to feedback (s)') +
  ylab('Deconvolved Signal (AU)') +
  theme(legend.text = element_text(size=25)) +
  theme(legend.title = element_text(size=30)) +
  theme(axis.title = element_text(size=40)) +
  theme(axis.text = element_text(size=25)) +
  theme(strip.text.x = element_text(size=40))
print(gg1)
dev.off()


####################
### vPFC - clock ###
####################
rm(Q)
rm(dQ)
message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
vmPFC <- clock_comb
vmPFC <- vmPFC %>% filter(evt_time > -6 & evt_time < 6)
rm(clock_comb)
vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
source('~/vPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
df <- df %>%
  group_by(id, run) %>% 
  mutate(iti_lag = lag(iti_ideal), rt_sec = rt_csv/1000) %>% ungroup() %>%
  mutate(v_chosen_sc = scale(v_chosen),
         abs_pe_max_sc = scale(abs(pe_max)),
         pe_max_sc = scale(pe_max),
         pe_max_lag_sc = scale(lag(pe_max)),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_vmax_sc = scale(rt_vmax),
         v_entropy_sc = scale(v_entropy),
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld3 = case_when(
           kld3 ==0 ~ NA_real_,
           kld3 >0 ~ log10(kld3)
         )) %>% mutate(log10kld3_lag = case_when(
           kld3_lag==0 ~NA_real_,
           kld3_lag>0 ~ log10(kld3_lag)
         ))
df <- df %>% group_by(id,run) %>% mutate(expl_longer =(case_when(
  rt_csv - rt_lag > 1 ~ 'Longer',
  rt_csv - rt_lag < -1 ~ '0',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run) %>% mutate(expl_shorter =(case_when(
  rt_csv - rt_lag > 1 ~ '0',
  rt_csv - rt_lag < -1 ~ 'Shorter',
  rt_csv - rt_lag < 1 & rt_csv - rt_lag > -1 ~ '0'
)))
df <- df %>% group_by(id,run)  %>% mutate(rt_bin = (case_when(
  rt_csv_sc <= -1 ~ '-1',
  rt_csv_sc > -1 & rt_csv_sc <= 0 ~ '-0.5',
  rt_csv_sc > 0 & rt_csv_sc <= 1 ~ '0.5',
  rt_csv_sc > 1 ~ '1'
)))
df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
  run_trial <= 15 ~ 'Early',
  run_trial > 15 & run_trial < 30 ~ 'Middle',
  run_trial >=30 ~ 'Late',
)))
df <- df %>% group_by(id) %>% mutate(entropy_split = case_when(
  v_entropy < mean(v_entropy,na.rm=TRUE) ~ 'low',
  v_entropy > mean(v_entropy,na.rm=TRUE) ~ 'high'
)) %>% ungroup()
df <- df %>% select(id,run,trial_bin,rewFunc,rt_bin,entropy_split,v_max_above_median,v_max_wi,expl_longer,expl_shorter,rt_csv_sc,v_entropy_wi, v_entropy_wi_change,run_trial,trial_neg_inv_sc,rt_vmax_change,kld3,abs_pe_max_sc,abs_pe_max_lag_sc,pe_max,pe_max_lag,pe_max_sc,pe_max_lag_sc,v_entropy_wi_change_lag)
Q <- merge(df, vmPFC, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
# test age & sex
demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
demo <- demo %>% rename(id=lunaid)
demo <- demo %>% select(!adult & !scandate)
Q <- inner_join(Q,demo,by=c('id'))
Q$female <- relevel(as.factor(Q$female),ref='0')
Q$age <- scale(Q$age)

dQ <- Q %>% group_by(id,evt_time,network,entropy_split,trial_bin) %>% summarize(vP = mean(vmPFC_decon,na.rm=TRUE),dP = sd(vmPFC_decon,na.rm=TRUE), N = length(vmPFC_decon))
dQ <- dQ %>% ungroup() %>% group_by(evt_time,network,entropy_split,trial_bin) %>% summarize(vP1 = mean(vP,na.rm=TRUE),dP1 = sd(dP,na.rm=TRUE), N1 = length(N))
dQ <- dQ %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR',network=='L'~'LIM'))
dQ$network1 <- factor(dQ$network1,levels=c('DMN','CTR','LIM'))
dQ <- dQ %>% filter(!is.na(entropy_split))
dQ$trial_bin <- factor(dQ$trial_bin,levels=c('Early','Middle','Late'))
pal = wes_palette("FantasticFox1", 3, type = "discrete")

dQ$entropy <- dQ$entropy_split
dQ <- dQ %>% ungroup() %>% select(!network) %>% rename(network=network1)

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('plot_raw_signal_clock_ventropy.pdf',height=12,width=24)
gg1 <- ggplot(dQ,aes(x=evt_time,y=vP1,color=network,linetype=entropy)) + 
  geom_point(size=5) + geom_line(size=2) + 
  geom_errorbar(aes(ymin=vP1-(dP1/sqrt(N1)),ymax=vP1+(dP1/sqrt(N1)))) + 
  geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 3) + 
  facet_wrap(~trial_bin) + 
  scale_color_manual(values=pal,labels=c("DMN","CTR","LIM")) +
  xlab('Time relative to clock (s)') +
  ylab('Deconvolved Signal (AU)') +
  theme(legend.text = element_text(size=25)) +
  theme(legend.title = element_text(size=30)) +
  theme(axis.title = element_text(size=40)) +
  theme(axis.text = element_text(size=25)) +
  theme(strip.text.x = element_text(size=40))
print(gg1)
dev.off()

rm(dQ)
dQ <- Q %>% group_by(id,evt_time,network,v_max_above_median,trial_bin) %>% summarize(vP = mean(vmPFC_decon,na.rm=TRUE),dP = sd(vmPFC_decon,na.rm=TRUE), N = length(vmPFC_decon))
dQ <- dQ %>% ungroup() %>% group_by(evt_time,network,v_max_above_median,trial_bin) %>% summarize(vP1 = mean(vP,na.rm=TRUE),dP1 = sd(dP,na.rm=TRUE), N1 = length(N))
dQ <- dQ %>% mutate(network1 = case_when(network=='D'~'DMN',network=='C'~'CTR',network=='L'~'LIM'))
dQ$network1 <- factor(dQ$network1,levels=c('DMN','CTR','LIM'))
dQ <- dQ %>% filter(!is.na(v_max_above_median))
dQ$trial_bin <- factor(dQ$trial_bin,levels=c('Early','Middle','Late'))
pal = wes_palette("FantasticFox1", 3, type = "discrete")

dQ <- dQ %>% mutate(v_max = case_when(v_max_above_median == TRUE ~ 'high', v_max_above_median == FALSE ~ 'low'))
dQ <- dQ %>% ungroup() %>% select(!network) %>% rename(network=network1)

setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
pdf('plot_raw_signal_clock_vmax.pdf',height=12,width=24)
gg1 <- ggplot(dQ,aes(x=evt_time,y=vP1,color=network,linetype=v_max)) + 
  geom_point(size=5) + geom_line(size=2) + 
  geom_errorbar(aes(ymin=vP1-(dP1/sqrt(N1)),ymax=vP1+(dP1/sqrt(N1)))) + 
  geom_vline(xintercept = 0, lty = "dashed", color = "#808080", size = 3) + 
  facet_wrap(~trial_bin) + 
  scale_color_manual(values=pal,labels=c("DMN","CTR","LIM")) +
  xlab('Time relative to clock (s)') +
  ylab('Deconvolved Signal (AU)') +
  theme(legend.text = element_text(size=25)) +
  theme(legend.title = element_text(size=30)) +
  theme(axis.title = element_text(size=40)) +
  theme(axis.text = element_text(size=25)) +
  theme(strip.text.x = element_text(size=40))
print(gg1)
dev.off()



