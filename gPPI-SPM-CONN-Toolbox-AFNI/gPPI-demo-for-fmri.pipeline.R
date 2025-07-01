
repo_directory_mmclock <- "~/clock_analysis"
source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
dfmmc <- get_trial_data(repo_directory=repo_directory_mmclock,dataset='mmclock_fmri')

trial_df <- dfmmc

files <- list.files('/Volumes/bierka_root/datamesh/PROC/MMClock/2025-06-27-gPPI-voxelwisedecon',recursive=TRUE,full.names=TRUE)

decon_df <- NULL;
for (iF in 1:length(files)){
  currF <- read_csv(files[iF])
  temp_str <- str_split(files[iF],'/')[[1]][10]
  id <- substring(temp_str,4,8)
  run <- substring(temp_str,10,13)
  if (grepl('hippocampus_l',files[iF])){
    currF <- currF %>% mutate(side = 'l',id = id, run=run)
  } else {
    currF <- currF %>% mutate(side = 'r',id = id, run=run)
  }
  decon_df <- rbind(decon_df,currF)
  print(iF)
}


hc_l_m <- oro.nifti::readNIfTI('/Volumes/Users/Andrew/long_axis_l_cobra_2.3mm.nii.gz',reorient=FALSE)
mi_l <- which(hc_l_m > 0, arr.ind=TRUE)
bin_cuts_l <- seq(min(hc_l_m[mi_l])-5e-3,max(hc_l_m[mi_l])+5e-3,length.out=12+1)
hc_r_m <- oro.nifti::readNIfTI('/Volumes/Users/Andrew/long_axis_r_cobra_2.3mm.nii.gz',reorient=FALSE)
mi_r <- which(hc_r_m > 0, arr.ind=TRUE)
bin_cuts_r <- seq(min(hc_r_m[mi_l])-5e-3,max(hc_r_m[mi_l])+5e-3,length.out=12+1)
decon_df <- decon_df %>% group_by(side) %>% mutate(atlas_value0 = case_when(side == 'l' ~ cut(atlas_value, bin_cuts_l),
                                                                                          side == 'r' ~ cut(atlas_value, bin_cuts_r))
) %>% ungroup()

uav <- sort(unique(decon_df$atlas_value0))

decon_df <- decon_df %>% mutate(bin_num=case_when(
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

decon_df$atlas_value0 <- as.numeric(decon_df$atlas_value0)

decon_df <- decon_df %>% mutate(HC_region = case_when(atlas_value0 < 5 ~ 'PH', atlas_value0 >=5 ~ 'AH'))

decon_df <- decon_df %>% group_by(id, run, time, HC_region) %>% summarize(decon = mean(decon,na.rm=TRUE)) %>% ungroup()

setwd("/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R")
write.csv(decon_df, file = 'MMClock-gPPI.csv')
