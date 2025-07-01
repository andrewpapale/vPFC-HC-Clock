
repo_directory_mmclock <- "~/clock_analysis"
dfmmc <- get_trial_data(repo_directory=repo_directory_mmclock,dataset='mmclock_fmri')

trial_df <- dfmmc %>% filter(id== c('10637','10997','11279'))

files <- list.files('/Volumes/bierka_root/datamesh/PROC/MMClock/2025-06-27-gPPI-voxelwisedecon',recursive=TRUE,full.names=TRUE)

gPPI_decon_demo <- NULL;
for (iF in 1:length(files)){
  currF <- read_csv(files[iF])
  temp_str <- str_split(files[iF],'/')[[1]][9]
  id <- substring(temp_str,4,8)
  run <- substring(temp_str,10,13)
  if (grepl('hippocampus_l',files[iF])){
    currF <- currF %>% mutate(side = 'l',id = id, run=run)
  } else {
    currF <- currF %>% mutate(side = 'r',id = id, run=run)
  }
  gPPI_decon_demo <- rbind(gPPI_decon_demo,currF)
  print(iF)
}


hc_l_m <- oro.nifti::readNIfTI('/Volumes/Users/Andrew/long_axis_l_cobra_2.3mm.nii.gz',reorient=FALSE)
mi_l <- which(hc_l_m > 0, arr.ind=TRUE)
bin_cuts_l <- seq(min(hc_l_m[mi_l])-5e-3,max(hc_l_m[mi_l])+5e-3,length.out=12+1)
gPPI_decon_demo <- gPPI_decon_demo %>% mutate(atlas_value0 = cut(atlas_value, bin_cuts_l))

uav <- sort(unique(gPPI_decon_demo$atlas_value0))

gPPI_decon_demo <- gPPI_decon_demo %>% mutate(bin_num=case_when(
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

gPPI_decon_demo$atlas_value0 <- as.numeric(gPPI_decon_demo$atlas_value0)

gPPI_decon_demo <- gPPI_decon_demo %>% mutate(HC_region = case_when(atlas_value0 < 5 ~ 'PH', atlas_value0 >=5 ~ 'AH'))
