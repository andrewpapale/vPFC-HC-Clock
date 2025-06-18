# script to preprocess HC data from BSOC
# adapted from script by Andrew Papale for Explore data
# Feb 2025

library(tidyverse)
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
repo_directory <- file.path('/Volumes/Users/Andrew/MEDuSA_data_BSOC')

hc_r <- read_csv(file.path(repo_directory,'clock_aligned_bsocial_hc_r.csv.gz'))
hc_l <- read_csv(file.path(repo_directory,'clock_aligned_bsocial_hc_l.csv.gz'))

hc_r <- hc_r %>% mutate(side='r')
hc_l <- hc_l %>% mutate(side='l')

hc_l_m <- oro.nifti::readNIfTI(file.path(repo_directory,'bsocial_hc_l.nii.gz'),reorient=FALSE)
hc_r_m <- oro.nifti::readNIfTI(file.path(repo_directory,'bsocial_hc_r.nii.gz'),reorient=FALSE)
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


hc$id <- as.character(hc$id)
hc <- hc %>% group_by(id, run) %>% mutate(run_trial = trial - min(trial) + 1) %>% ungroup()
hc <- hc %>% mutate(run1 = case_when(run=='run1'~1,run=='run2'~2)) %>% select(!run) %>% rename(run=run1)

hc <- hc %>% select(!decon_median & !decon_sd)
hc <- hc %>% mutate(HC_region = case_when(atlas_value < 5 ~ 'PH', atlas_value >=5 ~ 'AH'))

save(hc,file=file.path(rootdir,'BSOC_HC_clock_TRdiv2.Rdata'))
