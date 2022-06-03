# 2022-06-03 AndyP
# AIC Model Comparison

if (exists('D0')){ # start fresh
  rm(D0)
}
for (i in 1:length(decode_formula)){
  curr_model <- paste0('-',i,'.Rdata')
  curr_model <- Sys.glob(paste0('*',curr_model))
  load(curr_model)
  if (i==1){
    D0 <- ddf$fit_df
    new_AIC <- paste0('AIC','-',i)
    D0 <- D0 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
  } else {
    D1 <- ddf$fit_df
    new_AIC <- paste0('AIC','-',i)
    D1 <- D1 %>% select(evt_time,network,AIC) %>% rename(!! new_AIC := AIC)
    D0 <- inner_join(D0,D1,by=c('evt_time','network'))
    rm(D1)
  }
  rm(ddf)
  disp(i)
}