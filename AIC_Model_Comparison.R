# 2022-06-03 AndyP
# AIC Model Comparison


if (exists('D0')){ # start fresh
  rm(D0)
}
for (i in 1:length(decode_formula)){
  curr_model <- paste0('-vmPFC-HC-network-testing-trxH',i,'.Rdata')
  curr_model <- Sys.glob(paste0('*',curr_model))
  load(curr_model)
  if (i==1){
    D0 <- ddf$fit_df
    new_AIC <- paste0('AIC','-',i)
    D0 <- D0 %>% select(evt_time,network,AIC,HC_region) %>% rename(!! new_AIC := AIC)
  } else {
    D1 <- ddf$fit_df
    new_AIC <- paste0('AIC','-',i)
    D1 <- D1 %>% select(evt_time,network,AIC,HC_region) %>% rename(!! new_AIC := AIC)
    D0 <- inner_join(D0,D1,by=c('evt_time','network','HC_region'))
    rm(D1)
  }
  rm(ddf)
  disp(i)
}

D0 <- D0 %>% group_by(evt_time,network,HC_region) %>% mutate(AIC21 = `AIC-2` - `AIC-1`,
                                                   AIC31 = `AIC-3` - `AIC-1`,
                                                   AIC32 = `AIC-3` - `AIC-2`,
                                                   AIC41 = `AIC-4` - `AIC-1`,
                                                   AIC42 = `AIC-4` - `AIC-2`,
                                                   AIC43 = `AIC-4` - `AIC-3`
                                                   # AIC51 = `AIC-5` - `AIC-1`,
                                                   # AIC52 = `AIC-5` - `AIC-2`,
                                                   # AIC53 = `AIC-5` - `AIC-3`,
                                                   # AIC54 = `AIC-5` - `AIC-4`
                                                   ) %>% ungroup()

D0 <- D0 %>% select(! `AIC-1` & !`AIC-2` & !`AIC-3` & !`AIC-4`)# & !`AIC-5`)
D0 <- pivot_longer(D0,cols=starts_with('AIC'),names_to='model',names_prefix="AIC",values_to="AIC")

D1 <- D0 %>% filter(model==41 | model==42 | model==43)# | model==51)

epoch_label = paste("Time relative to",'feedback', "[s]")
pdf(paste0('AIC-Model-Comparison','.pdf'),width=9,height=6)
gg <- ggplot(D1, group=model) + facet_wrap(HC_region~network) +
  geom_line(size=0.1, aes(x=evt_time,y=(AIC),group=model)) + 
  geom_point(size=2,aes(x=evt_time,y=(AIC),color=model)) +
  xlab(epoch_label) + ylab('AIC') 
print(gg)
dev.off()
