

Q <- Q %>% group_by(id,run,trial,network, HC_region) %>% 
  summarize(vmPFC_decon = mean(vmPFC_decon,na.rm=TRUE),HCwithin = mean(HCwithin,na.rm=TRUE), HCbetween = mean(HCbetween,na.rm=TRUE)) %>% 
  ungroup()
Q <- Q %>% group_by(id,run,trial) %>% 
  pivot_wider(values_from=c(network,HC_region,vmPFC_decon,HCwithin,HCbetween),names_from = 'network') %>% 
  ungroup()