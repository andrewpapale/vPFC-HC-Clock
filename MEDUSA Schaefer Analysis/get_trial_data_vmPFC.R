# 2024-07-14 AndyP
# test 3 way of entropy*HC*lastoutcome

if (do_HC2vPFC_clock){
  rm(Q)
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMclock_clock_Aug2023.Rdata')
  vmPFC <- clock_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(clock_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock_Aug2023.Rdata')
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
  
  Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
  Q <- Q %>% select(!decon1)
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory=repo_directory,dataset='mmclock_fmri')
  df <- df %>% 
    group_by(id, run) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           iti_lag_sc = scale(iti_prev),
           score_sc = scale(score_csv),
           score_lag_sc = scale(lag(score_csv)),
           iti_sc = scale(iti_ideal),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           rt_vmax_sc = scale(rt_vmax),
           ev_sc = scale(ev),
           ev_lag_sc = scale(lag(ev)),
           v_entropy_sc = scale(v_entropy),
           v_entropy_lag_sc = scale(lag(v_entropy)),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
           abs_rt_vmax_change = abs(rt_vmax_change),
           rt_vmax_change_bin = case_when(
             abs_rt_vmax_change < 4/24 ~ 'No Change',
             abs_rt_vmax_change >= 4/24 ~ 'Change',
           ),
           v_entropy_wi_change_lag_bin = case_when(
             v_entropy_wi_change_lag < -0.5 ~ 'Decrease',
             v_entropy_wi_change_lag > 0.5  ~ 'Increase',
             v_entropy_wi_change_lag >= -0.5 & v_entropy_wi_change_lag <= 0.5 ~ 'No Change',
           ),
           rt_vmax_change_sc = scale(rt_vmax_change)) %>% arrange(id, run, run_trial) %>% mutate(log10kld4 = case_when(
             kld4 ==0 ~ NA_real_,
             kld4 >0 ~ log10(kld4)
           )) %>% mutate(log10kld4_lag = case_when(
             kld4_lag==0 ~NA_real_,
             kld4_lag>0 ~ log10(kld4_lag)
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
  df <- df %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                       total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))
  
  #df <- df %>% filter(!is.na(rt_vmax_change_bin) | !is.na(v_entropy_wi_change_lag_bin))
  df <- df %>% select(id,run,trial,run_trial,rt_lag_sc,total_earnings_split,outcome,iti_ideal,iti_prev,rt_csv,abs_pe_max_lag_sc,v_entropy_wi,rt_vmax_change_sc,trial_bin,rewFunc,trial_neg_inv_sc,rt_csv_sc,v_entropy_sc,expl_longer,expl_shorter,rt_bin,trial_bin,last_outcome,v_max_wi,v_entropy_wi_change_lag,score_lag_sc,iti_sc,iti_lag_sc,ev_lag_sc)
  Q <- inner_join(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
  Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;
  Q$HCbetween[Q$evt_time < -(Q$iti_prev)] = NA;
  Q$expl_longer <- relevel(as.factor(Q$expl_longer),ref='0')
  Q$expl_shorter <- relevel(as.factor(Q$expl_shorter),ref='0')
  Q$rt_bin <- relevel(as.factor(Q$rt_bin),ref='-0.5')
  Q$trial_bin <- relevel(as.factor(Q$trial_bin),ref='Middle')
  #Q$v_entropy_wi_change_lag_bin <- relevel(as.factor(Q$v_entropy_wi_change_lag_bin),ref='No Change')
  #Q$rt_vmax_change_bin <- relevel(as.factor(Q$rt_vmax_change_bin),ref='No Change')
  
  # test age & sex
  demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
  demo <- demo %>% rename(id=lunaid)
  demo <- demo %>% select(!adult & !scandate)
  Q <- inner_join(Q,demo,by=c('id'))
  Q$female <- relevel(as.factor(Q$female),ref='0')
  Q$age <- scale(Q$age)
  if (remove_1to10){
    Q <- Q %>% filter(trial > 10)
  }
  rm(decode_formula)
  decode_formula <- NULL
  #decode_formula[[1]] = formula(~age * HCwithin + female * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + v_max_wi * HCwithin + v_entropy_wi_change_lag * HCwithin + rt_csv_sc * HCwithin + iti_lag_sc * HCwithin + iti_sc * HCwithin + last_outcome * HCwithin + rt_vmax_change_sc * HCwithin + HCwithin * HCbetween + (1 | id/run)) 
  #decode_formula[[3]] = formula(~ age + female + abs_pe_max_lag_sc*HCwithin + trial_neg_inv_sc*HCwithin  + rt_csv_sc*HCwithin  + iti_lag_sc*HCwithin + iti_sc + HCwithin*HCbetween + (1| id/run))
  #decode_formula[[1]] = formula(~ v_max_wi*HCwithin + HCbetween + (1|id/run))
  #decode_formula[[2]] = formula(~ v_entropy_wi*HCwithin + HCbetween + (1|id/run))
  #decode_formula[[1]] = formula(~v_max_wi + (1|id/run))
  #decode_formula[[1]] = formula(~v_max_wi*HCwithin + v_entropy_wi*HCwithin + HCbetween + (1|id/run))
  decode_formula[[1]] = formula(~age * HCwithin + female * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run)) 
  decode_formula[[2]] = formula(~age * HCwithin + female * HCwithin + trial_neg_inv_sc * HCwithin + v_max_wi * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run)) 
  #decode_formula[[3]] = formula(~age * HCwithin + female * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + v_max_wi * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run)) 
  #decode_formula[[4]] = formula(~age * HCwithin + female * HCwithin + abs_pe_max_lag_sc * HCwithin + trial_neg_inv_sc * HCwithin + rt_lag_sc*HCwithin + HCbetween + (1|id/run))
  
  qT <- c(-0.7,0.43)
  
  if (do_network){
    
    splits = c('evt_time','network','HC_region')
    source("~/fmri.pipeline/R/mixed_by.R")
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)#,
                      # emmeans_spec = list(
                      #   H = list(outcome='vmPFC_decon', model_name='model1', 
                      #            specs=formula(~v_entropy_sc:HCwithin), at = list(v_entropy_sc=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                      #   V = list(outcome='vmPFC_decon', model_name='model1',
                      #            specs=formula(~v_max_wi:HCwithin), at=list(v_max_wi=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                      #   Tr = list(outcome='vmPFC_decon',model_name='model1',
                      #             specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5)))
                      # ),
                      # emtrends_spec = list(
                      #   H_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                      #               specs=formula(~v_entropy_sc:HCwithin), at=list(HCwithin=c(-1.5,1.5),v_entropy_sc=c(-1.5,1.5))),
                      #   T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                      #               specs=formula(~trial_neg_inv_sc:HCwithin),at=list(trial_neg_inv_sc=qT)),
                      #   V_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                      #               specs=formula(~v_max_wi:HCwithin),at=list(v_max_wi=c(-1.5,1,5))),
                      #   O_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                      #               specs=formula(~last_outcome:HCwithin)),
                      #   HCbw = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                      #               specs=formula(~HCbetween:HCwithin),at=list(HCbetween=c(-1.5,1,5))),
                      #   RT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                      #                specs=formula(~rt_csv_sc:HCwithin),at=list(rt_csv_sc=c(-1.5,1.5)))
                      # )
      )
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      if (i==4){
        save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-pe_max_lag_sc-',i,'.Rdata'))
      } else {
        save(ddf,file=paste0(curr_date,'-vmPFC-HC-network-clock-',i,'.Rdata'))
      }
    }
    
    if (do_vif){
      nE = unique(Q2$evt_time)
      nHC = unique(Q2$HC_region)
      nN = unique(Q2$network)
      vdf <- NULL
      network <- NULL
      HC_region <- NULL
      evt_time <- NULL
      for (iE in 1:length(nE)){
        for (iHC in 1:length(nHC)){
          for (iN in 1:length(nN)){
            temp <- Q %>% filter(evt_time == nE[iE])
            temp <- temp %>% filter(HC_region == nHC[iHC])
            temp <- temp %>% filter(network == nN[iN])
            m1 <- lmerTest::lmer(data=temp, vmPFC_decon ~ age * HCwithin + female * HCwithin + v_entropy_wi * HCwithin + trial_neg_inv_sc * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run))
            v0 <- car::vif(m1)
            vdf <- rbind(vdf,v0)
            network <- rbind(network,nN[iN])
            HC_region <- rbind(HC_region,nHC[iHC])
            evt_time <- rbind(evt_time,nE[iE])
          }
        }
      }
      vdf1 <- data.frame(vdf,network=network,HC_region=HC_region,evt_time=evt_time)
    }
    
    
  }
  
  if (do_symmetry){
    
    splits = c('evt_time','symmetry_group','HC_region')
    for (i in 1:length(decode_formula)){
      setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
      df0 <- decode_formula[[i]]
      print(df0)
      ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                      padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                      tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE),
                      emmeans_spec = list(
                        H = list(outcome='vmPFC_decon', model_name='model1', 
                                 specs=formula(~v_entropy_wi:HCwithin), at = list(v_entropy_wi=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        V = list(outcome='vmPFC_decon', model_name='model1',
                                 specs=formula(~v_max_wi:HCwithin), at=list(v_max_wi=c(-1.5,1.5),HCwithin=c(-1.5,1.5))),
                        Tr = list(outcome='vmPFC_decon',model_name='model1',
                                  specs=formula(~trial_neg_inv_sc:HCwithin), at=list(trial_neg_inv_sc=qT,HCwithin=c(-1.5,1.5)))
                      ),
                      emtrends_spec = list(
                        H_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_entropy_wi:HCwithin), at=list(v_entropy_wi=c(-1.5,1.5))),
                        T_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~trial_neg_inv_sc:HCwithin),at=list(trial_neg_inv_sc=qT)),
                        V_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~v_max_wi:HCwithin),at=list(v_max_wi=c(-1.5,1,5))),
                        A_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~age:HCwithin),at=list(age=c(-1.5,1.5))),
                        O_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~last_outcome:HCwithin)),
                        dH_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                     specs=formula(~v_entropy_wi_change_lag:HCwithin)),
                        dRT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                      specs=formula(~rt_vmax_change_sc:HCwithin),at=list(HCwithin=c(-1.5,1.5))),
                        HCbw = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                    specs=formula(~HCbetween:HCwithin),at=list(HCwithin=c(-1.5,1,5),HCbetween=c(-1.5,1,5))),
                        RT_HC = list(outcome='vmPFC_decon',model_name='model1',var='HCwithin',
                                     specs=formula(~rt_csv_sc:HCwithin),at=list(rt_csv_sc=c(-1.5,1.5),HCwithin=c(-1.5,1.5)))
                      )
      )
      curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
      save(ddf,file=paste0(curr_date,'-vmPFC-HC-symmetry-clock-',i,'.Rdata'))
    }
  }
}