# 2025-03-11 AndyP
# subject level random slopes for Trust

library(tidyverse)
source('/Users/dnplserv/vmPFC/mixed_by_clmm.R')

df <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/Tim_trust_trialdf_clean.csv')
df <- df %>% group_by(id,trustee) %>% mutate(tdec_prev = lag(t_decides), pdec_prev = lag(p_decides)) %>% ungroup()
df <- df %>% mutate(pdec_prev_char = case_when(pdec_prev == 0 ~ 'keep', pdec_prev == 1 ~ 'share'))
df$pdec_prev_char <- relevel(as.factor(df$pdec_prev_char),ref='keep')
df <- df %>% mutate(dataset_dummy = 'dataset')
df$block <- as.factor(df$block)
df$block <- relevel(df$block,ref='50')
df <- df %>% arrange(id,trial) %>% group_by(id,trustee,block) %>% mutate(block_trial = seq(1,length(unique(trial)),length.out=length(unique(trial)))) %>% ungroup()
df$trustee <- relevel(as.factor(df$trustee),ref='neutral')

df <- df %>% mutate(dataset_dummy = 'dataset')
df0 <- df %>% mutate(dataset_dummy = 'dataset_copy')
df1 <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/Tim_trust_trialdf_clean.csv')
df1 <- rbind(df,df0)
splits = 'dataset_dummy'

ddf <- mixed_by_clmm(data=df1,doclmm=TRUE,outcomes='p_decides', 
                     rhs_model_formula = formula(~pdec_prev_char*tdec_prev + scale(block_trial) + PE_mult_z + V_mult_z*trustee + tdec_prev + (1|id)), 
                     padjust_by = "term", ncores = 26, refit_on_nonconvergence = 3, split_on = splits,
                     return_models=TRUE,tidy_args = list(effects=c("fixed","ran_vals")),
                     emmeans_spec = list( 
                       VT = list(outcome='p_decides', model_name='model1', specs=formula(~V_mult_z:trustee), at = list(V_mult_z = c(-1.5,0.1,5)))
                     )
)
                     