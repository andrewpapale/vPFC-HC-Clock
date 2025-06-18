# 2023-02-28 AndyP
# Lead/Lag HC-vPFC analysis

vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
ncores = 26

##################################
### MMClock - clock hc -> vPFC ###
##################################

load('/Volumes/Users/Andrew/2023-03-02-MMClock-vPFC-clock.Rdata')
ddf <- ddf %>% group_by(id,run,trial,network) %>% mutate(vPFC_scaled = scale(vPFC_mean), vPFC_mean = mean(vPFC_mean,na.rm=TRUE)) %>% ungroup()
ddf <- ddf %>% mutate(run_trial = case_when(run==1 ~ trial,
                                            run==2 ~ trial-50,
                                            run==3 ~ trial-100,
                                            run==4 ~ trial-150,
                                            run==5 ~ trial-200,
                                            run==6 ~ trial-250,
                                            run==7 ~ trial-300,
                                            run==8 ~ trial-350
                                            ))
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock.Rdata')
hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
hc <- hc %>% group_by(id,run,HC_region) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()

Q <- merge(ddf,hc,by=c("id","run","run_trial","evt_time"))
Q <- Q %>% select(!decon1)

Q <- Q %>% rename(evt_time0 = evt_time)

### add lead/lag evt_times for hc/vPFC
Q <- Q %>% group_by(id,run,run_trial,network,HC_region) %>% mutate(evt_time_lag1 = lag(evt_time0, 1),
                                                 evt_time_lag2 = lag(evt_time0, 2),
                                                 evt_time_lag3 = lag(evt_time0, 3),
                                                 evt_time_lag4 = lag(evt_time0, 4),
                                                 evt_time_lag5 = lag(evt_time0, 5),
                                                 evt_time_lag6 = lag(evt_time0, 6),
                                                 evt_time_lag7 = lag(evt_time0, 7),
                                                 evt_time_lag8 = lag(evt_time0, 8),
                                                 evt_time_lead1 = lead(evt_time0,1),
                                                 evt_time_lead2 = lead(evt_time0,2),
                                                 evt_time_lead3 = lead(evt_time0,3),
                                                 evt_time_lead4 = lead(evt_time0,4),
                                                 evt_time_lead5 = lead(evt_time0,5)) %>% ungroup()

source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory = '~/clock_analysis',dataset='mmclock_fmri')
df <- df %>% select(id,run,run_trial,iti_prev,rt_csv,iti_ideal)

Q <- inner_join(Q,df,by=c('id','run','run_trial'))

Q$vPFC_scaled[Q$evt_time0 < -Q$iti_prev] = NA
Q$HCwithin[Q$evt_time0 < -Q$iti_prev] = NA
Q$vPFC_scaled[Q$evt_time0 > Q$rt_csv + Q$iti_ideal] = NA
Q$HCwithin[Q$evt_time0 > Q$rt_csv + Q$iti_ideal] = NA

save(Q,file='2023-03-02-HC-vPFC-intermediate.Rdata')

substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

iC = 1
col1 <- c('evt_time0','evt_time_lag1','evt_time_lag2','evt_time_lag3','evt_time_lag4','evt_time_lag5',
          'evt_time_lag6','evt_time_lag7','evt_time_lag8','evt_time_lead1','evt_time_lead2','evt_time_lead3',
          'evt_time_lead4','evt_time_lead5')
col2 <- c('evt_time0','evt_time_lag1','evt_time_lag2','evt_time_lag3','evt_time_lag4','evt_time_lag5',
          'evt_time_lag6','evt_time_lag7','evt_time_lag8','evt_time_lead1','evt_time_lead2','evt_time_lead3',
          'evt_time_lead4','evt_time_lead5')
Q2 <- NULL
for (iT1 in 1:14){
  for (iT2 in 1:14){
    message(col1[iT1], col2[iT2])
    iC <- iC + 1
    message(iC)
    Q0 <- Q %>% select(id,run,run_trial,HCwithin,HCbetween,HC_region,vPFC_scaled,vPFC_mean,network,col1[iT1],col2[iT2])
    vPFC <- Q0 %>% filter(!!as.name(col1[iT1]) > -3 & !!as.name(col1[iT1]) < 3) %>% select(!HCwithin & !HCbetween & !col1[iT1] & !col2[iT2] & !HC_region) # get window around zero
    HC <- Q0 %>% filter(!!as.name(col2[iT2]) > -3 & !!as.name(col2[iT2]) < 3) %>% select(!vPFC_scaled & !vPFC_mean & !col1[iT1] & !col2[iT2] & !network)
    vPFC <- vPFC[!duplicated(vPFC),]
    HC <- HC[!duplicated(HC),]
    
    temp_str1 <- substrRight(col1[iT1],2)
    temp_str2 <- substrRight(col2[iT2],2)
    
    if (grepl('^e',temp_str1)) { # case when we are at evt_time0
      temp_str1 <- substrRight(col1[iT1],1)
    } else if (grepl('^g',temp_str1)) {  # case where we are negative
      temp_str1 <- substrRight(col1[iT1],1)
      temp_str1 <- strcat('-',temp_str1)
    } else if (grepl('^d',temp_str1)) {  # case where we are positive
      temp_str1 <- substrRight(col1[iT1],1)
    }
    if (grepl('^e',temp_str2)) { # case when we are at evt_time0
      temp_str2 <- substrRight(col2[iT2],1)
    } else if (grepl('^g',temp_str2)) {  # case where we are negative
      temp_str2 <- substrRight(col2[iT2],1)
      temp_str2 <- strcat('-',temp_str2)
    } else if (grepl('^d',temp_str2)) {  # case where we are positive
      temp_str2 <- substrRight(col2[iT2],1)
    }
    vPFC <- vPFC %>% mutate(lead_lag1 = temp_str1)
    HC <- HC %>% mutate(lead_lag2 = temp_str2)
    
    Q1 <- inner_join(vPFC,HC,by=c('id','run','run_trial'))
    
    Q2 <- rbind(Q2,Q1)
    
    # setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
    # save(Q1,file='Lead_lag_HC_to_vPFC_intermediate.Rdata')   
    
    if (iT1==5 && iT2==14){
       setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
       save(Q2,file='Lead_lag_HC_to_vPFC_intermediate_1_of_3.Rdata')
       rm(Q2)
       gc()
       Q1 <- NULL
       Q2 <- NULL
     } else if (iT1==10 && iT2==14) {
       setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
       save(Q2,file='Lead_lag_HC_to_vPFC_intermediate_2_of_3.Rdata')
       rm(Q2)
       gc()
       Q1 <- NULL
       Q2 <- NULL
     } else if (iT1==14 && iT2==14) {
       setwd('~/vmPFC/MEDUSA Schaefer Analysis/')
       save(Q2,file='Lead_lag_HC_to_vPFC_intermediate_3_of_3.Rdata')
       rm(Q2)
       gc()
     }
  }
} 
  
  splits <- c('lead_lag','network','HC_region')

  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Lead_lag_HC_to_vPFC_intermediate_1_of_3.Rdata')
  Q2$lead_lag <- paste(Q2$lead_lag1,Q2$lead_lag2)
  
  decode_formula <- formula(~ HCwithin + (1|id/run/run_trial))
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  source("~/fmri.pipeline/R/mixed_by.R")
  hc2vPFC1 <- mixed_by(Q2, outcomes = "vPFC_scaled", rhs_model_formulae = decode_formula , split_on = splits,
                  padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                  tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  
  save(hc2vPFC1,file='hc2vPFC1.Rdata')  
  
  decode_formula <- formula(~ vPFC_scaled + (1|id/run/run_trial))
  
  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  source("~/fmri.pipeline/R/mixed_by.R")
  vPFC2hc1 <- mixed_by(Q2, outcomes = "HCwithin", rhs_model_formulae = decode_formula , split_on = splits,
                       padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                       tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  
  save(vPFC2hc1,file='vPFC2hc1.Rdata') 
  
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Lead_lag_HC_to_vPFC_intermediate_2_of_3.Rdata')
  Q2$lead_lag <- paste(Q2$lead_lag1,Q2$lead_lag2)

  decode_formula <- formula(~ HCwithin + (1|id/run/run_trial))

  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  source("~/fmri.pipeline/R/mixed_by.R")
  hc2vPFC2 <- mixed_by(Q2, outcomes = "vPFC_scaled", rhs_model_formulae = decode_formula , split_on = splits,
                       padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                       tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  
  save(hc2vPFC2,file='hc2vPFC2.Rdata') 

  decode_formula <- formula(~ vPFC_scaled + (1|id/run/run_trial))

  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  source("~/fmri.pipeline/R/mixed_by.R")
  vPFC2hc2 <- mixed_by(Q2, outcomes = "HCwithin", rhs_model_formulae = decode_formula , split_on = splits,
                       padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                       tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))

  save(vPFC2hc2,file='vPFC2hc2.Rdata')
  
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Lead_lag_HC_to_vPFC_intermediate_3_of_3.Rdata')
  Q2$lead_lag <- paste(Q2$lead_lag1,Q2$lead_lag2)

  decode_formula <- formula(~ HCwithin + (1|id/run/run_trial))

  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  source("~/fmri.pipeline/R/mixed_by.R")
  hc2vPFC3 <- mixed_by(Q2, outcomes = "vPFC_scaled", rhs_model_formulae = decode_formula , split_on = splits,
                       padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                       tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))
  
  save(hc2vPFC3,file='hc2vPFC3.Rdata')

  decode_formula <- formula(~ vPFC_scaled + (1|id/run/run_trial))

  setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
  source("~/fmri.pipeline/R/mixed_by.R")
  vPFC2hc3 <- mixed_by(Q2, outcomes = "HCwithin", rhs_model_formulae = decode_formula , split_on = splits,
                       padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                       tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE))

  save(vPFC2hc3,file='vPFC2hc3.Rdata')  

  rm(Q2)
  gc()
  
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/vPFC2hc1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/vPFC2hc2.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/vPFC2hc3.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/hc2vPFC1.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/hc2vPFC2.Rdata')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/hc2vPFC3.Rdata')

df1 <- vPFC2hc1$coef_df_reml
df2 <- vPFC2hc2$coef_df_reml
df3 <- vPFC2hc3$coef_df_reml
df1 <- df1 %>% filter(effect=='fixed' & term=='vPFC_scaled')
df2 <- df2 %>% filter(effect=='fixed' & term=='vPFC_scaled')
df3 <- df3 %>% filter(effect=='fixed' & term=='vPFC_scaled')
df <- rbind(df1,df2,df3)
df <- df %>% mutate(direction = 'vPFC to HC')


dq1 <- hc2vPFC1$coef_df_reml
dq2 <- hc2vPFC2$coef_df_reml
dq3 <- hc2vPFC3$coef_df_reml
dq1 <- dq1 %>% filter(effect=='fixed' & term=='HCwithin')
dq2 <- dq2 %>% filter(effect=='fixed' & term=='HCwithin')
dq3 <- dq3 %>% filter(effect=='fixed' & term=='HCwithin')
dq <- rbind(dq1,dq2,dq3)
dq <- dq %>% mutate(direction = 'HC to vPFC')

rm(hc2vPFC1,hc2vPFC2,hc2vPFC3,vPFC2hc1,vPFC2hc2,vPFC2hc3)
gc()

df <- df %>% filter(!grepl('NA',lead_lag))
df <- df %>% mutate(hc_time = as.numeric(case_when(grepl('^-8',lead_lag)~-8,
                                       grepl('^-7',lead_lag)~-7,
                                       grepl('^-6',lead_lag)~-6,
                                       grepl('^-5',lead_lag)~-5,
                                       grepl('^-4',lead_lag)~-4,
                                       grepl('^-3',lead_lag)~-3,
                                       grepl('^-2',lead_lag)~-2,
                                       grepl('^-1',lead_lag)~-1,
                                       grepl('^0',lead_lag)~0,
                                       grepl('^1',lead_lag)~1,
                                       grepl('^2',lead_lag)~2,
                                       grepl('^3',lead_lag)~3,
                                       grepl('^4',lead_lag)~4,
                                       grepl('^5',lead_lag)~5,
                                       grepl('^6',lead_lag)~6,
                                       grepl('^7',lead_lag)~7,
                                       grepl('^8',lead_lag)~8)),
                    vPFC_time = as.numeric(case_when(grepl('-8$',lead_lag)~-8,
                                                   grepl('-7$',lead_lag)~-7,
                                                   grepl('-6$',lead_lag)~-6,
                                                   grepl('-5$',lead_lag)~-5,
                                                   grepl('-4$',lead_lag)~-4,
                                                   grepl('-3$',lead_lag)~-3,
                                                   grepl('-2$',lead_lag)~-2,
                                                   grepl('-1$',lead_lag)~-1,
                                                   grepl('0$',lead_lag)~0,
                                                   grepl('1$',lead_lag)~1,
                                                   grepl('2$',lead_lag)~2,
                                                   grepl('3$',lead_lag)~3,
                                                   grepl('4$',lead_lag)~4,
                                                   grepl('5$',lead_lag)~5,
                                                   grepl('6$',lead_lag)~6,
                                                   grepl('7$',lead_lag)~7,
                                                   grepl('8$',lead_lag)~8)))
dq <- dq %>% filter(!grepl('NA',lead_lag))
dq <- dq %>% mutate(hc_time = as.numeric(case_when(grepl('^-8',lead_lag)~-8,
                                                   grepl('^-7',lead_lag)~-7,
                                                   grepl('^-6',lead_lag)~-6,
                                                   grepl('^-5',lead_lag)~-5,
                                                   grepl('^-4',lead_lag)~-4,
                                                   grepl('^-3',lead_lag)~-3,
                                                   grepl('^-2',lead_lag)~-2,
                                                   grepl('^-1',lead_lag)~-1,
                                                   grepl('^0',lead_lag)~0,
                                                   grepl('^1',lead_lag)~1,
                                                   grepl('^2',lead_lag)~2,
                                                   grepl('^3',lead_lag)~3,
                                                   grepl('^4',lead_lag)~4,
                                                   grepl('^5',lead_lag)~5,
                                                   grepl('^6',lead_lag)~6,
                                                   grepl('^7',lead_lag)~7,
                                                   grepl('^8',lead_lag)~8)),
                    vPFC_time = as.numeric(case_when(grepl('-8$',lead_lag)~-8,
                                                     grepl('-7$',lead_lag)~-7,
                                                     grepl('-6$',lead_lag)~-6,
                                                     grepl('-5$',lead_lag)~-5,
                                                     grepl('-4$',lead_lag)~-4,
                                                     grepl('-3$',lead_lag)~-3,
                                                     grepl('-2$',lead_lag)~-2,
                                                     grepl('-1$',lead_lag)~-1,
                                                     grepl('0$',lead_lag)~0,
                                                     grepl('1$',lead_lag)~1,
                                                     grepl('2$',lead_lag)~2,
                                                     grepl('3$',lead_lag)~3,
                                                     grepl('4$',lead_lag)~4,
                                                     grepl('5$',lead_lag)~5,
                                                     grepl('6$',lead_lag)~6,
                                                     grepl('7$',lead_lag)~7,
                                                     grepl('8$',lead_lag)~8)))


# 
ddf <- rbind(dq,df)

ddf <- ddf %>% arrange(hc_time,vPFC_time)
ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term, 
                       p_level_fdr = as.factor(case_when(
                         # p_fdr > .1 ~ '0',
                         # p_fdr < .1 & p_fdr > .05 ~ '1',
                         p_fdr > .05 ~ '1',
                         p_fdr < .05 & p_fdr > .01 ~ '2',
                         p_fdr < .01 & p_fdr > .001 ~ '3',
                         p_fdr <.001 ~ '4')))
# ddf <- ddf %>% group_by(network,HC_region,direction) %>% mutate(estimate_scaled = scale(estimate,center=FALSE)) %>% ungroup()

setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction/')
pdf('Lead_Lag_HC_vPFC_Analysis-HC-to-vPFC.pdf',height=24,width=12)
ggplot() + 
  geom_errorbar(data=ddf,aes(x=hc_time,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error,group=network)) + 
  geom_point(data=ddf,aes(x=hc_time,y=estimate,color=network,group=network,alpha=p_level_fdr,size=3)) + 
  geom_line(data=ddf,aes(x=hc_time,y=estimate,group=network,color=network)) +  
  facet_grid(vPFC_time~HC_region)
dev.off()
pdf('Lead_Lag_HC_vPFC_Analysis-vPFC-to-HC.pdf',height=24,width=12)
ggplot() + 
  geom_errorbar(data=ddq,aes(x=vPFC_time,y=estimate,ymin=estimate-std.error,ymax=estimate+std.error,group=network)) + 
  geom_point(data=ddq,aes(x=vPFC_time,y=estimate,color=network,group=network,alpha=p_level_fdr,size=3)) + 
  geom_line(data=ddq,aes(x=vPFC_time,y=estimate,group=network,color=network)) + 
  facet_grid(hc_time~HC_region)
dev.off()

