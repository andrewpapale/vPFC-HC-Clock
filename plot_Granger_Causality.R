plot_Granger_Causality <- function(Xuni,Xbi,Yuni,Ybi,toalign,toprocess,totest){

  
Xuni <- Xuni$residuals %>% 
  filter(!is.na(ddf$residuals$V1)) %>% 
  group_by(evt_time,network,HC_region) %>% 
  summarize(resid = var(V1,na.rm=TRUE)) %>% ungroup()    

Yuni <- Yuni$residuals %>% 
  filter(!is.na(ddf$residuals$V1)) %>% 
  group_by(evt_time,network,HC_region) %>% 
  summarize(resid = var(V1,na.rm=TRUE)) %>% ungroup()    

Xbi <- Xbi$residuals %>% 
  filter(!is.na(ddf$residuals$V1)) %>% 
  group_by(evt_time,network,HC_region) %>% 
  summarize(resid = var(V1,na.rm=TRUE)) %>% ungroup()    

Ybi <- Ybi$residuals %>% 
  filter(!is.na(ddf$residuals$V1)) %>% 
  group_by(evt_time,network,HC_region) %>% 
  summarize(resid = var(V1,na.rm=TRUE)) %>% ungroup()    
  
  
}