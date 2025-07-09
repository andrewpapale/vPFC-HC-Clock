# 2025-07-09 AndyP
# Plot AIC models

plot_aic_comparison <- function(model_list) {
  
  # Extract AIC values into a tidy dataframe
  all_df <- imap_dfr(model_list, function(model, name) {
    data.frame(
      Index = seq_len(nrow(model$fit_df)),
      AIC = model$fit_df$AIC,
      evt_time = model$fit_df$evt_time,
      network = model$fit_df$network,
      HC_region = model$fit_df$HC_region,
      name = name
    )
  })
  return(all_df)
}

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection/MMClock-model-comparison-sex.Rdata')

all_df <- plot_aic_comparison(models)
all_df0 <- all_df %>% filter(evt_time==0)

ggplot(data=all_df0,aes(x=name,y=AIC,color=HC_region,group=HC_region)) + 
  geom_point(size=2) + geom_line() + facet_wrap(~network,scales='free_y') + ggtitle('MMClock sex')


load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Age_vPFC_HC_model_selection/Bsocial-model-comparison-sex.Rdata')

all_df <- plot_aic_comparison(models)
all_df0 <- all_df %>% filter(evt_time==0)

ggplot(data=all_df0,aes(x=name,y=AIC,color=HC_region,group=HC_region)) + 
  geom_point(size=2) + geom_line() + facet_wrap(~network,scales='free_y') + ggtitle('BSocial sex')
