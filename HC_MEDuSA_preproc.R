# 2023-02-21 AndyP
# HC-MMClock-MEDuSA

hc_r <- read_csv('/Volumes/Users/Andrew/mmc_clock_aligned_HC_R.csv.gz')
hc_l <- read_csv('/Volumes/Users/Andrew/mmc_clock_aligned_HC_L.csv.gz')

source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory='~/clock_analysis',dataset='mmclock_fmri')


hc_r <- hc_r %>% mutate(side='r')
hc_l <- hc_l %>% mutate(side='l')

hc_l_m <- oro.nifti::readNIfTI('/Volumes/Users/Andrew/long_axis_l_cobra_2.3mm.nii.gz',reorient=FALSE)
hc_r_m <- oro.nifti::readNIfTI('/Volumes/Users/Andrew/long_axis_r_cobra_2.3mm.nii.gz',reorient=FALSE)
mi_l <- which(hc_l_m > 0, arr.ind=TRUE)
mi_r <- which(hc_r_m > 0, arr.ind=TRUE)
bin_cuts_l <- seq(min(hc_l_m[mi_l])-5e-3,max(hc_l_m[mi_l])+5e-3,length.out=12+1)
bin_cuts_r <- seq(min(hc_r_m[mi_r])-5e-3,max(hc_r_m[mi_r])+5e-3,length.out=12+1)
hc_l <- hc_l %>% mutate(atlas_value0 = cut(atlas_value, bin_cuts_l))
hc_r <- hc_r %>% mutate(atlas_value0 = cut(atlas_value, bin_cuts_r))

uav <- sort(unique(hc_l$atlas_value0))

hc_l <- hc_l %>% mutate(bin_num=case_when(
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

hc_r <- hc_r %>% mutate(bin_num=case_when(
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

hc <- hc %>% mutate(run1 = case_when(run=='run1'~1,
                                     run=='run2'~2,
                                     run=='run3'~3,
                                     run=='run4'~4,
                                     run=='run5'~5,
                                     run=='run6'~6,
                                     run=='run7'~7,
                                     run=='run8'~8)) %>% select(!run) %>% rename(run=run1)

df <- df %>% select(id,run,trial,run_trial)

hc <- inner_join(hc,df,by=c('id','run','trial'))

hc <- hc %>% select(!decon_median & !decon_sd)
hc <- hc %>% mutate(HC_region = case_when(bin_num < 5 ~ 'PH', bin_num >=5 ~ 'AH'))

save(hc,file='/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_clock.Rdata')


################
## feedback ###
###############

hc_r <- read_csv('/Volumes/Users/Andrew/mmc_rt_aligned_HC_R.csv.gz')
hc_l <- read_csv('/Volumes/Users/Andrew/mmc_rt_aligned_HC_L.csv.gz')

source('~/vmPFC/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory='~/clock_analysis',dataset='mmclock_fmri')


hc_r <- hc_r %>% mutate(side='r')
hc_l <- hc_l %>% mutate(side='l')

hc_l_m <- oro.nifti::readNIfTI('/Volumes/Users/Andrew/long_axis_l_cobra_2.3mm.nii.gz',reorient=FALSE)
hc_r_m <- oro.nifti::readNIfTI('/Volumes/Users/Andrew/long_axis_r_cobra_2.3mm.nii.gz',reorient=FALSE)
mi_l <- which(hc_l_m > 0, arr.ind=TRUE)
mi_r <- which(hc_r_m > 0, arr.ind=TRUE)
bin_cuts_l <- seq(min(hc_l_m[mi_l])-5e-3,max(hc_l_m[mi_l])+5e-3,length.out=12+1)
bin_cuts_r <- seq(min(hc_r_m[mi_r])-5e-3,max(hc_r_m[mi_r])+5e-3,length.out=12+1)
hc_l <- hc_l %>% mutate(atlas_value0 = cut(atlas_value, bin_cuts_l))
hc_r <- hc_r %>% mutate(atlas_value0 = cut(atlas_value, bin_cuts_r))

uav <- sort(unique(hc_l$atlas_value0))

hc_l <- hc_l %>% mutate(bin_num=case_when(
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

hc_r <- hc_r %>% mutate(bin_num=case_when(
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

hc <- hc %>% mutate(run1 = case_when(run=='run1'~1,
                                     run=='run2'~2,
                                     run=='run3'~3,
                                     run=='run4'~4,
                                     run=='run5'~5,
                                     run=='run6'~6,
                                     run=='run7'~7,
                                     run=='run8'~8)) %>% select(!run) %>% rename(run=run1)

df <- df %>% select(id,run,trial,run_trial)

hc <- inner_join(hc,df,by=c('id','run','trial'))

hc <- hc %>% select(!decon_median & !decon_sd)
hc <- hc %>% mutate(HC_region = case_when(bin_num < 5 ~ 'PH', bin_num >=5 ~ 'AH'))

save(hc,file='/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/HC_fb.Rdata')


