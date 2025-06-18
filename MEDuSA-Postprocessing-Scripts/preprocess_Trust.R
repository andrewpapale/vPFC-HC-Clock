# 2025-03-03 AndyP
# preprocess Trust MEDuSA 

library(tidyverse)

#hc_r <- read_csv('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/R/Explore_HC_r_clock.csv.gz')
#hc_l <- read_csv('/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Explore_HC/L/Explore_HC_l_clock.csv.gz')

hc_r <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/outcome_aligned_trust-transformed-hc-r.csv.gz')
hc_l <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/outcome_aligned_trust-transformed-hc-l.csv.gz')

df <- read_csv('/Volumes/Users/Andrew/MEDuSA_data_Trust/Tim_trust_trialdf_clean.csv')

hc_r <- hc_r %>% mutate(side='r')
hc_l <- hc_l %>% mutate(side='l')

hc_l_m <- oro.nifti::readNIfTI('/Volumes/Users/Andrew/MEDuSA_data_Trust/trust-transformed-hc-l.nii.gz',reorient=FALSE)
hc_r_m <- oro.nifti::readNIfTI('/Volumes/Users/Andrew/MEDuSA_data_Trust/trust-transformed-hc-r.nii.gz',reorient=FALSE)
mi_l <- which(hc_l_m > 0, arr.ind=TRUE)
mi_r <- which(hc_r_m > 0, arr.ind=TRUE)
bin_cuts_l <- seq(min(hc_l_m[mi_l])-5e-3,max(hc_l_m[mi_l])+5e-3,length.out=12+1)
bin_cuts_r <- seq(min(hc_r_m[mi_r])-1e-5,max(hc_r_m[mi_r])+1e-5,length.out=12+1)
hc_l <- hc_l %>% mutate(atlas_value0 = cut(atlas_value, bin_cuts_l))
hc_r <- hc_r %>% mutate(atlas_value0 = cut(atlas_value, bin_cuts_r))

uav <- sort(unique(hc_l$atlas_value0))

hc_l <- hc_l %>% mutate(atlas_value=case_when(
  atlas_value0==uav[1] ~ 1,
  atlas_value0==uav[2] ~ 2,
  atlas_value0==uav[3] ~ 3,
  atlas_value0==uav[4] ~ 4,
  atlas_value0==uav[5] ~ 5,
  atlas_value0==uav[6] ~ 6,
  atlas_value0==uav[7] ~ 7,
  atlas_value0==uav[8] ~ 8,
  atlas_value0==uav[9] ~ 9,
  atlas_value0==uav[10]~ 10,
  atlas_value0==uav[11]~ 11,
  atlas_value0==uav[12]~ 12,
))

uav <- sort(unique(hc_r$atlas_value0))

hc_r <- hc_r %>% mutate(atlas_value=case_when(
  atlas_value0==uav[1] ~ 1,
  atlas_value0==uav[2] ~ 2,
  atlas_value0==uav[3] ~ 3,
  atlas_value0==uav[4] ~ 4,
  atlas_value0==uav[5] ~ 5,
  atlas_value0==uav[6] ~ 6,
  atlas_value0==uav[7] ~ 7,
  atlas_value0==uav[8] ~ 8,
  atlas_value0==uav[9] ~ 9,
  atlas_value0==uav[10]~ 10,
  atlas_value0==uav[11]~ 11,
  atlas_value0==uav[12]~ 12,
))



hc <- rbind(hc_r,hc_l)
rm(hc_r,hc_l)
gc()
hc <- hc %>% mutate(run1 = case_when(run=='run1'~1)) %>% select(!run) %>% rename(run=run1)

df <- df %>% select(id,trial)

hc$id <- as.double(hc$id)
hc <- inner_join(hc,df,by=c('id','trial'))

hc <- hc %>% select(!decon_median & !decon_sd)
hc <- hc %>% mutate(HC_region = case_when(atlas_value < 5 ~ 'PH', atlas_value >=5 ~ 'AH'))

save(hc,file='/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/Trust/Trust_HC_outcome_TRdiv2.Rdata')



