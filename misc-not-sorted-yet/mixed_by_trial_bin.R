# 2023-05-25 AndyP
# Run Trial Test
# bin run_trial and do mixed_by, look at time course of activation across trial bins
# c.f. Fig 4c-d in HC paper

library(tidyverse)
library(wesanderson)

repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
repo_str <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/'
plot_dir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction'
ncores <- 26
paper_version <- 'v11'
do_MMClock <- FALSE
do_plot_MMClock <- TRUE
do_Explore = FALSE
do_plot_Explore = TRUE
source("~/fmri.pipeline/R/mixed_by.R")

if (do_MMClock){
  if (exists('Q')){rm(Q)}
  message("Loading vmPFC medusa data from cache: ", vmPFC_cache_dir)
  load(file.path(vmPFC_cache_dir,  'clock_vmPFC_Schaefer_tall_ts_1.Rdata'))
  vmPFC <- clock_comb
  vmPFC <- vmPFC %>% filter(evt_time > -5 & evt_time < 5)
  rm(clock_comb)
  vmPFC <- vmPFC %>% select(id,run,run_trial,decon_mean,atlas_value,evt_time,region,symmetry_group,network)
  vmPFC <- vmPFC %>% rename(vmPFC_decon = decon_mean)
  load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock.Rdata')
  hc <- hc %>% filter(evt_time > -5 & evt_time < 5)
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon1),HCbetween=mean(decon1,na.rm=TRUE)) %>% ungroup()
  
  
  
  Q <- merge(vmPFC,hc,by=c("id","run","run_trial","evt_time"))
  Q <- Q %>% select(!decon1)
  source('~/vmPFC/get_trial_data_vmPFC.R')
  df <- get_trial_data_vmPFC(repo_directory=repo_directory,dataset='mmclock_fmri')
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
  df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
    run_trial <= 10 ~ '1-10',
    run_trial > 10 & run_trial <= 20 ~ '10-20',
    run_trial > 20 & run_trial <= 30 ~ '20-30',
    run_trial > 30 & run_trial <= 40 ~ '30-40',
    run_trial > 40 & run_trial <= 50 ~ '40-50'
  )))
  df <- df %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                       total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))
  
  #df <- df %>% filter(!is.na(rt_vmax_change_bin) | !is.na(v_entropy_wi_change_lag_bin))
  df <- df %>% select(id,run,run_trial,rt_lag_sc,total_earnings_split,outcome,iti_ideal,iti_prev,rt_csv,abs_pe_max_lag_sc,v_entropy_wi,rt_vmax_change_sc,trial_bin,rewFunc,trial_neg_inv_sc,rt_csv_sc,v_entropy_sc,trial_bin,last_outcome,v_max_wi,v_entropy_wi_change_lag,score_lag_sc,iti_sc,iti_lag_sc,ev_lag_sc)
  Q <- inner_join(df, Q, by = c("id", "run", "run_trial")) %>% arrange("id","run","run_trial","evt_time")
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA;
  Q$vmPFC_decon[Q$evt_time < -(Q$iti_prev)] = NA;
  Q$HCwithin[Q$evt_time < -(Q$iti_prev)] = NA;
  Q$HCbetween[Q$evt_time < -(Q$iti_prev)] = NA;
  
  Q <- Q %>% mutate(online = case_when(evt_time < 0 ~'offline',
                                       evt_time >=0 ~ 'online'))
  
  # test age & sex
  demo <- read.table(file=file.path(repo_directory, 'fmri/data/mmy3_demographics.tsv'),sep='\t',header=TRUE)
  demo <- demo %>% rename(id=lunaid)
  demo <- demo %>% select(!adult & !scandate)
  Q <- inner_join(Q,demo,by=c('id'))
  Q$female <- relevel(as.factor(Q$female),ref='0')
  Q$age <- scale(Q$age)
  
  if (exists('decode_formula')){rm(decode_formula)}
  decode_formula <- NULL
  decode_formula[[1]] = formula(~v_entropy_wi * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run)) 
  decode_formula[[2]] = formula(~v_max_wi * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run)) 
  decode_formula[[3]] = formula(~v_entropy_wi * HCwithin + v_max_wi * HCwithin + rt_lag_sc*HCwithin + iti_lag_sc * HCwithin + last_outcome * HCwithin + HCbetween + (1 | id/run)) 
  
  
  qT <- c(-0.7,0.43)
  
  
  splits = c('online','network','HC_region','trial_bin')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)#,
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'-vPFC-HC-network-clock-trial_bin-',i,'.Rdata'))
  }
}

##################
### PLOT #########
##################
if (do_plot_MMClock){
  for (iM in 1:3){
    load(paste0(repo_str,'/','2023-05-25-vPFC-HC-network-clock-trial_bin-',iM,'.Rdata')) # ddf
    
    # x-axis = trial bin, y-axis = estimate, color = network, facet wrap = HC_region / online, size/alpha = p_level_fdr
    ddf <- ddf$coef_df_reml
    
    # plot HCwithin
    hc_w <- ddf %>% filter(effect=='fixed' & term=='HCwithin')
    
    hc_w <- hc_w %>% mutate(p_level_fdr = as.factor(case_when(
      padj_fdr_term > .05 ~ '1',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
      padj_fdr_term <.001 ~ '4')))
    hc_w$p_level_fdr <- factor(hc_w$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < 0.05", "p < 0.01", "p < 0.001"))
    
    
    # plot HCwithin*v_entropy_wi
    hc_en <- ddf %>% filter(effect=='fixed' & (term=='v_entropy_wi:HCwithin' | term=='HCwithin:v_entropy_wi'))
    
    hc_en <- hc_en %>% mutate(p_level_fdr = as.factor(case_when(
      padj_fdr_term > .05 ~ '1',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
      padj_fdr_term <.001 ~ '4')))
    hc_en$p_level_fdr <- factor(hc_en$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < 0.05", "p < 0.01", "p < 0.001"))
    
    # plot HCwithin*v_max_wi
    hc_vmax <- ddf %>% filter(effect=='fixed' & (term=='HCwithin:v_max_wi' | term=='v_max_wi:HCwithin'))
    
    hc_vmax <- hc_vmax %>% mutate(p_level_fdr = as.factor(case_when(
      padj_fdr_term > .05 ~ '1',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
      padj_fdr_term <.001 ~ '4')))
    hc_vmax$p_level_fdr <- factor(hc_vmax$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < 0.05", "p < 0.01", "p < 0.001"))
    
    pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
    pal = palette()
    pal[1] = pal3[2]
    pal[2] = pal3[1]
    pal[3] = pal3[3]
    
    pal1 = palette()
    pal1[1] = '#ff484d'
    pal1[2] = '#2a52ff'
    
    setwd(plot_dir)
    
    title_str <- paste0(paper_version,'-MMClock-clock-network-by-HC-by-trial_bin-HCwithin-',iM,'.pdf')
    pdf(title_str,height=8,width=16)
    
    if (all(hc_w$p_level_fdr=='p < 0.001')){
      gg1 <- ggplot(hc_w,aes(x=trial_bin,y=estimate)) + 
        geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
        geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
        geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
        facet_grid(online~HC_region) + 
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
      gg1 <- ggplot(hc_w,aes(x=trial_bin,y=estimate)) + 
        geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
        geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
        geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
        facet_grid(online~HC_region) + 
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
    
    if (nrow(hc_en)>0){
      title_str <- paste0(paper_version,'-MMClock-clock-network-by-HC-by-trial_bin-HCwithin-v_entropy_wi-',iM,'.pdf')
      pdf(title_str,height=8,width=16)
      
      if (all(hc_en$p_level_fdr=='p < 0.001')){
        gg1 <- ggplot(hc_en,aes(x=trial_bin,y=estimate)) + 
          geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
          geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
          geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
          facet_grid(online~HC_region) + 
          scale_color_manual(values = pal) +
          scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
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
        gg1 <- ggplot(hc_en,aes(x=trial_bin,y=estimate)) + 
          geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
          geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
          geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
          facet_grid(online~HC_region) + 
          scale_color_manual(values = pal) +
          geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
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
    if (nrow(hc_vmax)>0){
      title_str <- paste0(paper_version,'-MMClock-clock-network-by-HC-by-trial_bin-HCwithin-v_max_wi-',iM,'.pdf')
      pdf(title_str,height=8,width=16)
      
      if (all(hc_vmax$p_level_fdr=='p < 0.001')){
        gg1 <- ggplot(hc_vmax,aes(x=trial_bin,y=estimate)) + 
          geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
          geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
          geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
          facet_grid(online~HC_region) + 
          scale_color_manual(values = pal) +
          scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
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
        gg1 <- ggplot(hc_vmax,aes(x=trial_bin,y=estimate)) + 
          geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
          geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
          geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
          facet_grid(online~HC_region) + 
          scale_color_manual(values = pal) +
          geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
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
  }
}


if (do_Explore){
  if (exists('Q')){rm(Q)}
  load('/Volumes/Users/Andrew/MEDuSA_data_Explore/clock-vPFC.Rdata')
  md <- md %>% filter(evt_time > -4 & evt_time < 4)
  # load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock.Rdata')
  # hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  # hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
  # hc <- hc %>% mutate(block = case_when(trial <= 40 ~ 1, 
  #                                     trial > 40 & trial <= 80 ~ 2,
  #                                     trial > 80 & trial <=120 ~ 3, 
  #                                     trial > 120 & trial <=160 ~ 4,
  #                                     trial > 160 & trial <=200 ~ 5,
  #                                     trial > 200 & trial <=240 ~ 6))
  # hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup() # 12 -> 2  hc <- hc %>% rename(decon_mean=decon_mean1)
  #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_clock.Rdata')
  hc <- read_csv('/Volumes/Users/Bea/StriatumHippThalamus/clock_aligned_striatum_hipp_thalamus.csv.gz')
  hc <- hc %>% mutate(run1 = as.integer(str_sub(run,4,4))) %>% select(!run) %>% rename(run=run1)
  hc <- hc %>% filter(atlas_value %in% c(223,224,225,226,227,228,229,230))
  hc <- hc %>% select(!decon_median & !decon_sd)
  #hc #load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/Explore_HC_fb.Rdata')
  #hc <- hc %>% select(id,run,trial,run_trial,decon_mean,evt_time,side,HC_region,atlas_value)
  hc <- hc %>% filter(evt_time > -4 & evt_time < 4)
  hc <- hc %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                        trial > 40 & trial <= 80 ~ 2,
                                        trial > 80 & trial <=120 ~ 3, 
                                        trial > 120 & trial <=160 ~ 4,
                                        trial > 160 & trial <=200 ~ 5,
                                        trial > 200 & trial <=240 ~ 6))
  hc <- hc %>% mutate(HC_region = case_when(atlas_value==223 ~ 'PH',
                                            atlas_value==224 ~ 'PH',
                                            atlas_value==225 ~ 'AH',
                                            atlas_value==226 ~ 'AH',
                                            atlas_value==227 ~ 'PH',
                                            atlas_value==228 ~ 'PH',
                                            atlas_value==229 ~ 'AH',
                                            atlas_value==230 ~ 'AH'))
  hc <- hc %>% mutate(run_trial = case_when(trial <= 40 ~ trial,
                                            trial > 40 & trial <= 80 ~ trial - 40,
                                            trial > 80 & trial <=120 ~ trial - 80,
                                            trial > 120 & trial <= 160  ~ trial - 120,
                                            trial > 160 & trial <=200 ~ trial - 160,
                                            trial > 200 & trial <=240 ~ trial - 200))
  hc <- hc %>% select(!atlas_value)
  hc <- hc %>% group_by(id,run,run_trial,evt_time,HC_region) %>% summarize(decon1 = mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  hc <- hc %>% rename(decon_mean=decon1)
  hc <- hc %>% group_by(id,run) %>% mutate(HCwithin = scale(decon_mean),HCbetween=mean(decon_mean,na.rm=TRUE)) %>% ungroup()
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory='/Volumes/Users/Andrew/MEDuSA_data_Explore',dataset='explore')
  df <- df %>%
    group_by(id,run_number) %>%
    arrange(id, run_number, trial) %>%
    mutate(condition_trial = run_trial-floor(run_trial/40.5)*40,
           condition_trial_neg_inv = -1000 / condition_trial,
           condition_trial_neg_inv_sc = as.vector(scale(condition_trial_neg_inv)),
    ) %>% ungroup()
  df <- df %>%
    group_by(id) %>% 
    mutate(v_entropy_sc_r = scale(v_entropy)) %>% ungroup() %>%
    group_by(id, run) %>% 
    mutate(v_chosen_sc = scale(v_chosen),
           abs_pe_max_sc = scale(abs(pe_max)),
           score_sc = scale(score_csv),
           iti_sc = scale(iti_ideal),
           iti_lag_sc = lag(iti_sc),
           v_chosen_sc = scale(v_chosen),
           pe_max_sc = scale(pe_max),
           pe_max_lag_sc = scale(lag(pe_max)),
           abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
           rt_vmax_sc = scale(rt_vmax),
           ev_sc = scale(ev),
           v_entropy_sc = scale(v_entropy),
           v_max_lag_sc = scale(lag(v_max)),
           v_entropy_wi_change_lag = lag(v_entropy_wi_change),
           abs_rt_vmax_change = abs(rt_vmax_change),
           rt_vmax_change_bin = case_when(
             abs_rt_vmax_change < 4/24 ~ 'No Change',
             abs_rt_vmax_change >= 4/24 ~ 'Change'
           ),
           v_entropy_wi_change_bin = case_when(
             v_entropy_wi_change < -0.5 ~ 'Decrease',
             v_entropy_wi_change > 0.5  ~ 'Increase',
             v_entropy_wi_change >= -0.5 & v_entropy_wi_change <= 0.5 ~ 'No Change'
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
    run_trial <= 13 ~ 'Early',
    run_trial > 13 & run_trial < 26 ~ 'Middle',
    run_trial >=26 ~ 'Late',
  )))
  df <- df %>% group_by(id,run) %>% mutate(trial_bin = (case_when(
    run_trial <= 10 ~ '1-10',
    run_trial > 10 & run_trial <= 20 ~ '10-20',
    run_trial > 20 & run_trial <= 30 ~ '20-30',
    run_trial > 30 & run_trial <= 40 ~ '30-40',
    run_trial > 40 & run_trial <= 50 ~ '40-50'
  )))
  df <- df %>% mutate(total_earnings_split = case_when(total_earnings >= median(df$total_earnings,na.rm=TRUE)~'richer',
                                                       total_earnings < median(df$total_earnings,na.rm=TRUE)~'poorer'))
  
  df <- df %>% select(total_earnings_split,condition_trial_neg_inv_sc,iti_ideal,rt_lag_sc,iti_lag_sc,iti_prev,iti_sc,v_entropy_wi_change_lag,outcome,ev_sc,v_chosen_sc,last_outcome,rt_csv,rt_bin,rt_vmax_change_sc,trial_bin,rt_csv_sc,run_trial,id,run,v_entropy_wi,v_max_wi,trial_neg_inv_sc,trial,rewFunc)
  df$id <- as.character(df$id)
  Q <- full_join(md,df,by=c('id','run','trial'))
  Q <- Q %>% rename(vmPFC_decon = decon_mean) %>% select(!decon_median & !decon_sd)
  rm(md)
  hc$id <- as.character(hc$id)
  Q <- Q %>% mutate(block = case_when(trial <= 40 ~ 1, 
                                      trial > 40 & trial <= 80 ~ 2,
                                      trial > 80 & trial <=120 ~ 3, 
                                      trial > 120 & trial <=160 ~ 4,
                                      trial > 160 & trial <=200 ~ 5,
                                      trial > 200 & trial <=240 ~ 6))
  Q <- inner_join(Q,hc,by=c('id','run','run_trial','evt_time'))
  rm(hc)
  Q$HCwithin[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$vmPFC_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCbetween[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA
  Q$HCwithin[Q$evt_time < -Q$rt_csv] = NA
  Q$vmPFC_decon[Q$evt_time < -Q$rt_csv] = NA
  Q <- Q %>% mutate(network = case_when(
    atlas_value==67 | atlas_value==171 | atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 'CTR',
    atlas_value==89 | atlas_value==194 | atlas_value==88 | atlas_value==192 | atlas_value==84 | atlas_value==191 | atlas_value==86 | atlas_value==161 ~ 'DMN',
    atlas_value==55 | atlas_value==160 | atlas_value==56 | atlas_value==159 ~ 'LIM'
  )) %>% mutate(symmetry_group=case_when(
    atlas_value==67 | atlas_value==171 ~ 6,
    atlas_value==65 | atlas_value==66 | atlas_value==170 ~ 1,
    atlas_value==89 | atlas_value==194 ~ 7,
    atlas_value==88 | atlas_value==192 ~ 5,
    atlas_value==84 | atlas_value==191 ~ 4,
    atlas_value==86 | atlas_value==161 ~ 2,
    atlas_value==55 | atlas_value==160 ~ 8,
    atlas_value==56 | atlas_value==159 ~ 3
  ))
  Q <- Q %>% arrange(id,run,trial,evt_time)
  Q <- Q %>% filter(evt_time > -4 & evt_time < 4)
  
  demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
  demo <- demo %>% select(registration_redcapid,age,gender,registration_group,wtar,education_yrs)
  demo <- demo %>% rename(id=registration_redcapid,group=registration_group)
  demo$gender <- relevel(as.factor(demo$gender),ref='M')
  demo$age <- scale(demo$age)
  demo$wtar <- scale(demo$wtar)
  demo$education_yrs <- scale(demo$education_yrs)
  
  Q <- merge(demo,Q,by='id')
  #Q <- Q %>% filter(total_earnings_split=='richer')
  Q$group <- relevel(factor(Q$group),ref='HC')
  Q <- Q %>% filter(group!='ATT')
  Q <- Q %>% filter(group=='HC')
  Q <- Q %>% filter(!is.na(rewFunc))
  Q <- Q %>% filter(trial > 10)
  
  Q <- Q %>% mutate(online = case_when(evt_time < 0 ~'offline',
                                       evt_time >=0 ~ 'online'))
  
  rm(decode_formula)
  decode_formula <- NULL
  decode_formula[[1]] = formula(~age*HCwithin + gender*HCwithin + wtar*HCwithin + v_max_wi*HCwithin + rt_lag_sc*HCwithin + HCbetween + last_outcome*HCwithin+ (1|id))
  decode_formula[[2]] = formula(~age*HCwithin + gender*HCwithin + wtar*HCwithin + v_entropy_wi*HCwithin + rt_lag_sc*HCwithin + HCbetween + last_outcome*HCwithin+ (1|id))
  decode_formula[[3]] = formula(~age*HCwithin + gender*HCwithin + wtar*HCwithin + v_max_wi*HCwithin + v_entropy_wi*HCwithin + rt_lag_sc*HCwithin + HCbetween + last_outcome*HCwithin + (1|id))
  decode_formula[[4]] = formula(~v_max_wi*HCwithin + HCbetween + (1|id))
  
  qT <- c(-0.7,0.43)
  
  
  splits = c('online','network','HC_region','trial_bin')
  source("~/fmri.pipeline/R/mixed_by.R")
  for (i in 1:length(decode_formula)){
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection')
    df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "vmPFC_decon", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3,
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE)#,
    )
    curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
    save(ddf,file=paste0(curr_date,'vPFC-HC-Explore-network-clock-trial_bin-',i,'.Rdata'))
  }
}

##################
### PLOT #########
##################
if (do_plot_Explore){
  for (iM in 1:3){
    load(paste0(repo_str,'/','2023-05-25-vPFC-HC-Explore-network-clock-trial_bin-',iM,'.Rdata')) # ddf
    
    # x-axis = trial bin, y-axis = estimate, color = network, facet wrap = HC_region / online, size/alpha = p_level_fdr
    ddf <- ddf$coef_df_reml
    
    # plot HCwithin
    hc_w <- ddf %>% filter(effect=='fixed' & term=='HCwithin')
    
    hc_w <- hc_w %>% mutate(p_level_fdr = as.factor(case_when(
      padj_fdr_term > .05 ~ '1',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
      padj_fdr_term <.001 ~ '4')))
    hc_w$p_level_fdr <- factor(hc_w$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < 0.05", "p < 0.01", "p < 0.001"))
    
    
    # plot HCwithin*v_entropy_wi
    hc_en <- ddf %>% filter(effect=='fixed' & (term=='v_entropy_wi:HCwithin' | term=='HCwithin:v_entropy_wi'))
    
    hc_en <- hc_en %>% mutate(p_level_fdr = as.factor(case_when(
      padj_fdr_term > .05 ~ '1',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
      padj_fdr_term <.001 ~ '4')))
    hc_en$p_level_fdr <- factor(hc_en$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < 0.05", "p < 0.01", "p < 0.001"))
    
    # plot HCwithin*v_max_wi
    hc_vmax <- ddf %>% filter(effect=='fixed' & (term=='HCwithin:v_max_wi' | term=='v_max_wi:HCwithin'))
    
    hc_vmax <- hc_vmax %>% mutate(p_level_fdr = as.factor(case_when(
      padj_fdr_term > .05 ~ '1',
      padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
      padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
      padj_fdr_term <.001 ~ '4')))
    hc_vmax$p_level_fdr <- factor(hc_vmax$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < 0.05", "p < 0.01", "p < 0.001"))
    
    pal3 = wes_palette("FantasticFox1", 3, type = "discrete")
    pal = palette()
    pal[1] = pal3[2]
    pal[2] = pal3[1]
    pal[3] = pal3[3]
    
    pal1 = palette()
    pal1[1] = '#ff484d'
    pal1[2] = '#2a52ff'
    
    setwd(plot_dir)
    
    title_str <- paste0(paper_version,'-Explore-clock-network-by-HC-by-trial_bin-HCwithin-',iM,'.pdf')
    pdf(title_str,height=8,width=16)
    
    if (all(hc_w$p_level_fdr=='p < 0.001')){
      gg1 <- ggplot(hc_w,aes(x=trial_bin,y=estimate)) + 
        geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
        geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
        geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
        facet_grid(online~HC_region) + 
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
      gg1 <- ggplot(hc_w,aes(x=trial_bin,y=estimate)) + 
        geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
        geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
        geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
        facet_grid(online~HC_region) + 
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
    
    if (nrow(hc_en)>0){
      title_str <- paste0(paper_version,'-Explore-clock-network-by-HC-by-trial_bin-HCwithin-v_entropy_wi-',iM,'.pdf')
      pdf(title_str,height=8,width=16)
      
      if (all(hc_en$p_level_fdr=='p < 0.001')){
        gg1 <- ggplot(hc_en,aes(x=trial_bin,y=estimate)) + 
          geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
          geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
          geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
          facet_grid(online~HC_region) + 
          scale_color_manual(values = pal) +
          scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
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
        gg1 <- ggplot(hc_en,aes(x=trial_bin,y=estimate)) + 
          geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
          geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
          geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
          facet_grid(online~HC_region) + 
          scale_color_manual(values = pal) +
          geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
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
    if (nrow(hc_vmax)>0){
      title_str <- paste0(paper_version,'-Explore-clock-network-by-HC-by-trial_bin-HCwithin-v_max_wi-',iM,'.pdf')
      pdf(title_str,height=8,width=16)
      
      if (all(hc_vmax$p_level_fdr=='p < 0.001')){
        gg1 <- ggplot(hc_vmax,aes(x=trial_bin,y=estimate)) + 
          geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
          geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
          geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
          facet_grid(online~HC_region) + 
          scale_color_manual(values = pal) +
          scale_alpha_discrete(range=c(1,1)) + scale_size_manual(values=c(6)) +
          geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
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
        gg1 <- ggplot(hc_vmax,aes(x=trial_bin,y=estimate)) + 
          geom_line(position=position_dodge(width=0.2),aes(color=network,group=network)) + 
          geom_point(position=position_dodge(width=0.2),aes(color=network,group=network,size=p_level_fdr,alpha=p_level_fdr)) + 
          geom_errorbar(position=position_dodge(width=0.2),width=0,color='white',aes(ymin=estimate-std.error,ymax=estimate+std.error,group=network)) +
          facet_grid(online~HC_region) + 
          scale_color_manual(values = pal) +
          geom_hline(yintercept = 0, lty = 'dashed', color = 'white', size = 1) +
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
  }
}



