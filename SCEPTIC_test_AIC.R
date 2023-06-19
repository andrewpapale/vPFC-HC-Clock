# 2023-06-13 AndyP
# AIC model comparison compressed/uncompressed

do_mmc_vPFC = TRUE
do_mmc_HC = TRUE
do_mmc_vPFC_HC = TRUE
do_entropy = FALSE
do_value = TRUE
if (do_mmc_vPFC){
  
  if (do_entropy){
    # test entropy AIC vPFC clock
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str_u <- '-vmPFC-uncompressed-network-clock-1.Rdata'
    model_str_c <- '-vmPFC-compressed-network-clock-1.Rdata'
    model_str_u <- Sys.glob(paste0('*',model_str_u))
    model_str_c <- Sys.glob(paste0('*',model_str_c))
    load(model_str_u)
    uncomp <- ddf;
    load(model_str_c)
    comp <- ddf;
    rm(ddf)
    gc()
    
    uncomp <- uncomp$fit_df;
    comp <- comp$fit_df;
    
    uncomp <- uncomp %>% mutate(RL_model = 'uncompressed')
    comp <- comp %>% mutate(RL_model = 'compressed')
    
    ddq <- rbind(uncomp,comp)
    ddq <- ddq %>% select(evt_time,network,AIC,RL_model)
    ddq <- ddq %>% pivot_wider(names_from=RL_model,values_from = AIC)
    
    ggplot(ddq,aes(x=evt_time,y=uncompressed-compressed,color=network,group=network)) + geom_point() + geom_line()
  }
  
  if (do_value){
    # test v_max AIC vPFC clock
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str_u <- '-vmPFC-uncompressed-network-clock-2.Rdata'
    model_str_c <- '-vmPFC-compressed-network-clock-2.Rdata'
    model_str_u <- Sys.glob(paste0('*',model_str_u))
    model_str_c <- Sys.glob(paste0('*',model_str_c))
    load(model_str_u)
    uncomp <- ddf;
    model_str_c <- model_str_c[1]
    load(model_str_c)
    comp <- ddf;
    rm(ddf)
    gc()
    
    uncomp <- uncomp$fit_df;
    comp <- comp$fit_df;
    
    uncomp <- uncomp %>% mutate(RL_model = 'uncompressed')
    comp <- comp %>% mutate(RL_model = 'compressed')
    
    ddq <- rbind(uncomp,comp)
    ddq <- ddq %>% select(evt_time,network,AIC,RL_model)
    ddq <- ddq %>% pivot_wider(names_from=RL_model,values_from = AIC)
    
    ggplot(ddq,aes(x=evt_time,y=uncompressed-compressed,color=network,group=network)) + geom_point() + geom_line()
  }
}

if (do_mmc_HC){
  if (do_entropy){
    # test entropy AIC vPFC clock
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str_u <- '-HC-uncompressed-region-clock-1.Rdata'
    model_str_c <- '-HC-compressed-region-clock-1.Rdata'
    model_str_u <- Sys.glob(paste0('*',model_str_u))
    model_str_c <- Sys.glob(paste0('*',model_str_c))
    load(model_str_u)
    uncomp <- ddf;
    load(model_str_c)
    comp <- ddf;
    rm(ddf)
    gc()
    
    uncomp <- uncomp$fit_df;
    comp <- comp$fit_df;
    
    uncomp <- uncomp %>% mutate(RL_model = 'uncompressed')
    comp <- comp %>% mutate(RL_model = 'compressed')
    
    ddq <- rbind(uncomp,comp)
    ddq <- ddq %>% select(evt_time,HC_region,AIC,RL_model)
    ddq <- ddq %>% pivot_wider(names_from=RL_model,values_from = AIC)
    
    ggplot(ddq,aes(x=evt_time,y=uncompressed-compressed,color=HC_region,group=HC_region)) + geom_point() + geom_line()
  }
  if (do_value){
    # test v_max AIC vPFC clock
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str_u <- '-HC-uncompressed-region-clock-2.Rdata'
    model_str_c <- '-HC-compressed-region-clock-2.Rdata'
    model_str_u <- Sys.glob(paste0('*',model_str_u))
    model_str_c <- Sys.glob(paste0('*',model_str_c))
    load(model_str_u)
    uncomp <- ddf;
    model_str_c <- model_str_c[1]
    load(model_str_c)
    comp <- ddf;
    rm(ddf)
    gc()
    
    uncomp <- uncomp$fit_df;
    comp <- comp$fit_df;
    
    uncomp <- uncomp %>% mutate(RL_model = 'uncompressed')
    comp <- comp %>% mutate(RL_model = 'compressed')
    
    ddq <- rbind(uncomp,comp)
    ddq <- ddq %>% select(evt_time,HC_region,AIC,RL_model)
    ddq <- ddq %>% pivot_wider(names_from=RL_model,values_from = AIC)
    
    ggplot(ddq,aes(x=evt_time,y=uncompressed-compressed,color=HC_region,group=HC_region)) + geom_point() + geom_line()
  }
}

if (do_mmc_vPFC_HC){
  if (do_entropy){
    # test entropy AIC vPFC clock
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str_u <- '-vmPFC-HC-uncompressed-network-clock-1.Rdata'
    model_str_c <- '-vmPFC-HC-compressed-network-clock-1.Rdata'
    model_str_u <- Sys.glob(paste0('*',model_str_u))
    model_str_c <- Sys.glob(paste0('*',model_str_c))
    load(model_str_u)
    uncomp <- ddf;
    load(model_str_c)
    comp <- ddf;
    rm(ddf)
    gc()
    
    uncomp <- uncomp$fit_df;
    comp <- comp$fit_df;
    
    uncomp <- uncomp %>% mutate(RL_model = 'uncompressed')
    comp <- comp %>% mutate(RL_model = 'compressed')
    
    ddq <- rbind(uncomp,comp)
    ddq <- ddq %>% select(evt_time,HC_region,network,AIC,RL_model)
    ddq <- ddq %>% pivot_wider(names_from=RL_model,values_from = AIC)
    
    ggplot(ddq,aes(x=evt_time,y=uncompressed-compressed,color=network,group=network)) + geom_point() + geom_line() + facet_wrap(~HC_region)
  }
  if (do_value){
    # test v_max AIC vPFC clock
    setwd('~/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/')
    model_str_u <- '-vmPFC-HC-uncompressed-network-clock-2.Rdata'
    model_str_c <- '-vmPFC-HC-compressed-network-clock-2.Rdata'
    model_str_u <- Sys.glob(paste0('*',model_str_u))
    model_str_c <- Sys.glob(paste0('*',model_str_c))
    load(model_str_u)
    uncomp <- ddf;
    model_str_c <- model_str_c[1]
    load(model_str_c)
    comp <- ddf;
    rm(ddf)
    gc()
    
    uncomp <- uncomp$fit_df;
    comp <- comp$fit_df;
    
    uncomp <- uncomp %>% mutate(RL_model = 'uncompressed')
    comp <- comp %>% mutate(RL_model = 'compressed')
    
    ddq <- rbind(uncomp,comp)
    ddq <- ddq %>% select(evt_time,HC_region,network,AIC,RL_model)
    ddq <- ddq %>% pivot_wider(names_from=RL_model,values_from = AIC)
    
    ggplot(ddq,aes(x=evt_time,y=uncompressed-compressed,color=network,group=network)) + geom_point() + geom_line() + facet_wrap(~HC_region)
  }
}