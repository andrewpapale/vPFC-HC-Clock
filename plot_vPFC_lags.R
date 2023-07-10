# 2023-07-03 AndyP
# plot vPFC_lags

library(tidyverse)
source("~/vmPFC/plot_mixed_by_vmPFC_HC.R", echo=TRUE)
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-07-03-lagged-models-vmPFC-HC-network-clock-1.Rdata')
plot_mixed_by_vmPFC_HC(ddf,'clock','network-by-HC','lag0','v11',1,'LR')
lag0 <- ddf$fit_df %>% mutate(model = 'lag0')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-07-03-lagged-models-vmPFC-HC-network-clock-2.Rdata')
plot_mixed_by_vmPFC_HC(ddf,'clock','network-by-HC','lag1','v11',1,'LR')
lag1 <- ddf$fit_df %>% mutate(model = 'lag1')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-07-03-lagged-models-vmPFC-HC-network-clock-3.Rdata')
plot_mixed_by_vmPFC_HC(ddf,'clock','network-by-HC','lag2','v11',1,'LR')
lag2 <- ddf$fit_df %>% mutate(model = 'lag2')
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-07-03-lagged-models-vmPFC-HC-network-clock-4.Rdata')
lag12 <- ddf$fit_df %>% mutate(model = 'lag12')
plot_mixed_by_vmPFC_HC(ddf,'clock','network-by-HC','lag12','v11',1,'LR')

df <- rbind(lag0,lag1,lag2,lag12)
df <- df %>% select(evt_time,network,HC_region,AIC,model)
df <- df %>% pivot_wider(names_from = model, values_from = AIC)
df <- df %>% mutate(lag1_m_lag0 = lag1-lag0,
                    lag2_m_lag0 = lag2-lag0,
                    lag2_m_lag1 = lag2-lag1,
                    lag12_m_lag1 = lag12-lag1,
                    lag12_m_lag2 = lag12-lag2,
                    lag12_m_lag0 = lag12-lag0)


ggplot(df, aes(x=evt_time,y=lag1_m_lag0,color=network,group=network)) + geom_point(size=10) + geom_line() + facet_wrap(~HC_region)
ggplot(df, aes(x=evt_time,y=log10(-lag1)-log10(-lag0),color=network,group=network)) + geom_point(size=10) + geom_line() + facet_wrap(~HC_region)
