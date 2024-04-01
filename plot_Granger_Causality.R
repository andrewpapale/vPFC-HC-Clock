#plot_Granger_Causality <- function(Xuni,Xbi,Yuni,Ybi,toalign,toprocess,totest){

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-04-01-MMClock-Granger-Causality-pred-vmPFC-network-clock-10.Rdata')
Xuni <- ddf;
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-04-01-MMClock-Granger-Causality-pred-vmPFC-network-clock-4.Rdata')
Xbi <- ddf;
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-04-01-MMClock-Granger-Causality-pred-HC-network-clock-7.Rdata')
Yuni <- ddf;
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2024-04-01-MMClock-Granger-Causality-pred-HC-network-clock-1.Rdata')
Ybi <- ddf;
  
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
  
Xuni <- Xuni %>% rename(Xuni_resid = resid)
Xbi <- Xbi %>% rename(Xbi_resid = resid)
Yuni <- Yuni %>% rename(Yuni_resid = resid)
Ybi <- Ybi %>% rename(Ybi_resid = resid)  

Fyx <- inner_join(Xuni,Xbi,by=c('evt_time','network','HC_region'))
Fxy <- inner_join(Yuni,Ybi,by=c('evt_time','network','HC_region'))

Fyx <- Fyx %>% mutate(Fyx = log(Xuni_resid/Xbi_resid)) %>% select(!Xuni_resid & !Xbi_resid)
Fxy <- Fxy %>% mutate(Fxy = log(Yuni_resid/Ybi_resid)) %>% select(!Yuni_resid & !Ybi_resid)

dF <- inner_join(Fyx,Fxy,by=c('evt_time','network','HC_region'))
dF <- dF %>% mutate(dF = Fyx-Fxy)

ggplot(dF,aes(x=evt_time,y=dF,color=network,group=network)) + geom_point() + geom_line() + facet_wrap(~HC_region)
#}