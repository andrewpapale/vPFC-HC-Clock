# 2023-02-10 AndyP
# Get stats for vPFC-HC paper

library(tidyverse)

#####################
### Entropy ########
####################


# network - MMClock
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-01-07-vmPFC-network-clock-1.Rdata')

ddf <- ddf$coef_df_reml
#if (strcmp(toprocess,"network")){
ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

h_term <- ddf %>% filter(term=='v_entropy_wi' & effect=='fixed')
D <- h_term %>% filter(network=='D')
C <- h_term %>% filter(network=='C')

summary(D)
summary(C)

# network - EXP
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-01-07-Explore-vmPFC-network-clock-1.Rdata')

ddf <- ddf$coef_df_reml
#if (strcmp(toprocess,"network")){
ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

h_term <- ddf %>% filter(term=='v_entropy_wi' & effect=='fixed')
D <- h_term %>% filter(network=='DMN')
L <- h_term %>% filter(network=='LIM')
C <- h_term %>% filter(network=='CTR')
summary(D)
summary(L)
summary(C)


# symmetry - MMClock
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-01-07-vmPFC-symmetry-clock-1.Rdata')

ddf <- ddf$coef_df_reml
#if (strcmp(toprocess,"network")){
ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

h_term <- ddf %>% filter(term=='v_entropy_wi' & effect=='fixed')
D <- h_term %>% filter(symmetry_group %in% c(2,4,5,7))
L <- h_term %>% filter(symmetry_group == 3)
rl9_10 <- h_term %>% filter(symmetry_group==6)
r11_47 <- h_term %>% filter(symmetry_group==1)
m142532 <- h_term %>% filter(symmetry_group == 4)
summary(D)
summary(L)
summary(rl9_10)
summary(r11_47)
summary(m142532)

############################
###### Value ###############
############################

# network - MMClock
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-01-07-vmPFC-network-clock-2.Rdata')

ddf <- ddf$coef_df_reml
#if (strcmp(toprocess,"network")){
ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

v_term <- ddf %>% filter(term=='v_max_wi' & effect=='fixed')
D <- v_term %>% filter(network=='D')
C <- v_term %>% filter(network=='C')
L <- v_term %>% filter(network=='L')

summary(D)
summary(C)
summary(L)
# network - EXP
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-01-07-Explore-vmPFC-network-clock-2.Rdata')

ddf <- ddf$coef_df_reml
#if (strcmp(toprocess,"network")){
ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

v_term <- ddf %>% filter(term=='v_max_wi' & effect=='fixed')
D <- v_term %>% filter(network=='DMN')
L <- v_term %>% filter(network=='LIM')
C <- v_term %>% filter(network=='CTR')
summary(D)
summary(L)
summary(C)


# symmetry - MMClock
load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-01-07-vmPFC-symmetry-clock-2.Rdata')

ddf <- ddf$coef_df_reml
#if (strcmp(toprocess,"network")){
ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

v_term <- ddf %>% filter(term=='v_max_wi' & effect=='fixed')
D <- v_term %>% filter(symmetry_group %in% c(4,5,7))
fp10 <- v_term %>% filter(symmetry_group == 2)
r11_13 <- v_term %>% filter(symmetry_group == 8)
rl9_10 <- v_term %>% filter(symmetry_group==6)
r11_47 <- v_term %>% filter(symmetry_group==1)
m142532 <- v_term %>% filter(symmetry_group == 4)
summary(D)
summary(fp10)
summary(r11_13)
summary(rl9_10)
summary(r11_47)
summary(m142532)

##############################
##### B2B - Entropy ##########
##############################

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-02-09-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-both-1.Rdata')

ddq <- ddq$coef_df_reml
#if (strcmp(toprocess,"network")){
ddq <- ddq  %>% mutate(p_fdr = padj_fdr_term)
ddq <- ddq %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

rt_term <- ddq %>% filter(term=='rt_lag:subj_level_rand_slope')
D <- rt_term %>% filter(network=='D')
C <- rt_term %>% filter(network=='C')
summary(D)
summary(C)

rt_vmax_term <- ddq %>% filter(term=='subj_level_rand_slope:rt_vmax_lag')
D <- rt_vmax_term %>% filter(network=='D')
summary(D)


load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-02-10-vmPFC-network-ranslopes-replication-clock-pred-rt_csv_sc-both-1.Rdata')

ddq <- ddq$coef_df_reml
#if (strcmp(toprocess,"network")){
ddq <- ddq  %>% mutate(p_fdr = padj_fdr_term)
ddq <- ddq %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

rt_term <- ddq %>% filter(term=='rt_lag:subj_level_rand_slope')
D <- rt_term %>% filter(network=='D')
C <- rt_term %>% filter(network=='C')
summary(D)
summary(C)

rt_vmax_term <- ddq %>% filter(term=='subj_level_rand_slope:rt_vmax_lag')
D <- rt_vmax_term %>% filter(network=='D')
L <- rt_vmax_term %>% filter(network=='L')
summary(D)
summary(L)

##############################
##### B2B - Value ##########
##############################

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-02-09-vmPFC-network-ranslopes-clock-pred-rt_csv_sc-both-2.Rdata')

ddq <- ddq$coef_df_reml
#if (strcmp(toprocess,"network")){
ddq <- ddq  %>% mutate(p_fdr = padj_fdr_term)
ddq <- ddq %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

rt_term <- ddq %>% filter(term=='rt_lag:subj_level_rand_slope')
D <- rt_term %>% filter(network=='D')
summary(D)

C <- rt_term %>% filter(network=='C')
summary(C)

L <- rt_term %>% filter(network=='L')
summary(L)

rt_vmax_term <- ddq %>% filter(term=='subj_level_rand_slope:rt_vmax_lag')
D <- rt_vmax_term %>% filter(network=='D')
summary(D)


load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2023-02-10-vmPFC-network-ranslopes-replication-clock-pred-rt_csv_sc-both-2.Rdata')

ddq <- ddq$coef_df_reml
#if (strcmp(toprocess,"network")){
ddq <- ddq  %>% mutate(p_fdr = padj_fdr_term)
ddq <- ddq %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddq$p_level_fdr <- factor(ddq$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

rt_term <- ddq %>% filter(term=='rt_lag:subj_level_rand_slope')
D <- rt_term %>% filter(network=='D')
summary(D)

C <- rt_term %>% filter(network=='C')
summary(C)

L <- rt_term %>% filter(network=='L')
summary(L)

rt_vmax_term <- ddq %>% filter(term=='subj_level_rand_slope:rt_vmax_lag')
D <- rt_vmax_term %>% filter(network=='D')
summary(D)


############################
#### HC ###################
###########################

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2022-11-28-HC-axis-clock-2.Rdata')

ddf <- ddf$coef_df_reml
#if (strcmp(toprocess,"network")){
ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

v_term <- ddf %>% filter(term=='v_max_wi' & effect=='fixed')
AH <- v_term %>% filter(HC_region=='AH')
PH <- v_term %>% filter(HC_region=='PH')
summary(AH)
summary(PH)

load('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2022-11-28-vmPFC-HC-network-clock-1.Rdata')
ddf <- ddf$coef_df_reml
ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term)
ddf <- ddf %>% mutate(p_level_fdr = as.factor(case_when(
  # p_fdr > .1 ~ '0',
  # p_fdr < .1 & p_fdr > .05 ~ '1',
  p_fdr > .05 ~ '1',
  p_fdr < .05 & p_fdr > .01 ~ '2',
  p_fdr < .01 & p_fdr > .001 ~ '3',
  p_fdr <.001 ~ '4')))
ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

a_term <- ddf %>% filter(term=='HCwithin' & effect=='fixed')
S <- aov(estimate ~ HC_region*network, data=a_term)
anova(S)


h_term <- ddf %>% filter(term=='HCwithin:v_entropy_sc' & effect=='fixed')
D = h_term %>% filter(network=='D')
summary(D)

C = h_term %>% filter(network=='C')
summary(C)


v_term <- ddf %>% filter(term=='HCwithin:v_max_wi' & effect=='fixed')
D = v_term %>% filter(network=='D')
summary(D)

C = v_term %>% filter(network=='C')
summary(C)

L = v_term %>% filter(network=='L')
summary(L)
