# 2025-08-01 AndyP
# Create MMClock-Light Directory for gPPI

library(stringr)

# this file enumeration takes about 30 min, so saved output.
# elig_files <- list.files("/Volumes/bierka_root/datamesh/RAW/MMClock/MR_Proc", pattern = "nfaswuktm_clock[1-8]_5_drop2_trunc[0-9]+\\.nii\\.gz$",full.names = TRUE, recursive=TRUE)
#conf_files <- list.files("/Volumes/bierka_root/datamesh/RAW/MMClock/MR_Proc", pattern = "nuisance_regressors.txt",full.names = TRUE, recursive=TRUE)
setwd('/Volumes/Users/Andrew/2025-05-29-MMC-gPPI-Prototype-forMNH/')
load('2025-08-04-elig_files.Rdata')
path_out0 <- '/Users/dnplserv/MMClock_lite'

dirs <- list.dirs(path = '/Volumes/bierka_root/datamesh/RAW/MMClock/MR_Proc',full.names=FALSE,recursive=FALSE)
nD <- length(dirs)

for (iD in 1:nD){
  dir.create(file.path(path_out0, dirs[iD]))
}
for (iD in 1:nD){
  dir.create(file.path(path_out0, dirs[iD],'mni_5mm_aroma'))
}

dirs_new <- list.dirs(path = path_out0,full.names = FALSE, recursive = FALSE)
dirs_new_full <- list.dirs(path = path_out0,full.names = TRUE, recursive = FALSE)
# testing for 8 runs / subject 578 / 8 != an integer
nF = length(elig_files)
for (iF in 1:nF){
  id <- str_split(elig_files[iF],'/')[[1]][8]
  run <- str_split(elig_files[iF],'/')[[1]][10]
  file <- str_split(elig_files[iF],'/')[[1]][11]
  for (iD in 1:length(dirs_new)){
    if (id==dirs_new[iD]){
      dir.create(paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run))
      if (!file.exists(paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run,'/',file))){
        file.copy(elig_files[iF],paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run))
      }
      break;
    }
  }
  print(id)
}

nF = length(conf_files)
conf_files0 <- NULL
iR <- 1;
for (iF in 1:nF){
  tempfold <- str_split(conf_files[iF],'/')[[1]][9]
  if (tempfold == 'mni_5mm_aroma'){
      conf_files0[iR] <- conf_files[iF]
      iR <- iR + 1;
  }
}

nF = length(conf_files0)
for (iF in 1:nF){
  id <- str_split(conf_files0[iF],'/')[[1]][8]
  run <- str_split(conf_files0[iF],'/')[[1]][10]
  file <- str_split(conf_files0[iF],'/')[[1]][11]
  for (iD in 1:length(dirs_new)){
    if (id==dirs_new[iD]){
      dir.create(paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run))
      if (!file.exists(paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run,'/',file))){
        file.copy(conf_files0[iF],paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run))
      }
      break;
    }
  }
  print(id)
}


motion_files <- list.files("/Volumes/bierka_root/datamesh/RAW/MMClock/MR_Proc", pattern = "motion.par",full.names = TRUE, recursive=TRUE)

nF = length(motion_files)
motion_files0 <- NULL
iR <- 1;
for (iF in 1:nF){
  tempfold <- str_split(motion_files[iF],'/')[[1]][9]
  if (tempfold == 'mni_5mm_aroma'){
    motion_files0[iR] <- motion_files[iF]
    iR <- iR + 1;
  }
}

nF = length(motion_files0)
for (iF in 1:nF){
  id <- str_split(motion_files0[iF],'/')[[1]][8]
  run <- str_split(motion_files0[iF],'/')[[1]][10]
  file <- str_split(motion_files0[iF],'/')[[1]][11]
  for (iD in 1:length(dirs_new)){
    if (id==dirs_new[iD]){
      dir.create(paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run))
      if (!file.exists(paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run,'/',file))){
        file.copy(motion_files0[iF],paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run))
      }
      break;
    }
  }
  print(id)
}

nF <- length(motion_files0)
nD <- length(elig_files)
ix2rem_conf <- NULL
ix2rem_motion <- NULL
for (iF in 1:nF){
  motion_temp_id <- str_split(motion_files0[iF],'/')[[1]][8]
  motion_temp_run <- str_split(motion_files0[iF],'/')[[1]][10]
  conf_temp_id <- str_split(conf_files0[iF],'/')[[1]][8]
  conf_temp_run <- str_split(conf_files0[iF],'/')[[1]][10]
  in_the_set <- FALSE
  for (iD in 1:nD){
    mri_temp_id <- str_split(elig_files[iD],'/')[[1]][8]
    mri_temp_run <- str_split(elig_files[iD],'/')[[1]][10]
    if (conf_temp_id==mri_temp_id && conf_temp_run==mri_temp_run){
      in_the_set <- TRUE
    }
  }
  if (in_the_set==FALSE){
    ix2rem_conf <- rbind(ix2rem_conf,iF)
  }
  in_the_set <- FALSE
  for (iD in 1:nD){
    mri_temp_id <- str_split(elig_files[iD],'/')[[1]][8]
    mri_temp_run <- str_split(elig_files[iD],'/')[[1]][10]
    if (motion_temp_id==mri_temp_id && motion_temp_run==mri_temp_run){
      in_the_set <- TRUE
    }
  }
  if (in_the_set==FALSE){
    ix2rem_motion <- rbind(ix2rem_motion,iF)
  }  
}

save(ix2rem_motion, file='ix2rem_motion.Rdata')
save(ix2rem_conf, file='ix2rem_conf.Rdata')

conf_files0 <- conf_files0[-ix2rem_conf]
motion_files0 <- motion_files0[-ix2rem_motion]

save(motion_files0,file='motion_files.Rdata')
save(conf_files0,file='confound_files.Rdata')
