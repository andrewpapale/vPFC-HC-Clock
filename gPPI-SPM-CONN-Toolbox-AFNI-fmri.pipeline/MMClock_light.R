# 2025-08-01 AndyP
# Create MMClock-Light Directory for gPPI

library(stringr)

# this file enumeration takes about 30 min, so saved output.
# elig_files <- list.files("/Volumes/bierka_root/datamesh/RAW/MMClock/MR_Proc", pattern = "nfaswuktm_clock[1-8]_5_drop2_trunc[0-9]+\\.nii\\.gz$",full.names = TRUE, recursive=TRUE)
setwd('/Volumes/Users/Andrew/2025-05-29-MMC-gPPI-Prototype-forMNH/')
load('2025-08-04-elig_files.Rdata')
path_out0 <- '/Volumes/bierka_root/datamesh/PROC/MMClock/MMClock_lite'

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
  for (iD in 1:length(dirs_new)){
    if (id==dirs_new[iD]){
      dir.create(paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run))
      file.copy(elig_files[iF],paste0(dirs_new_full[iD],'/mni_5mm_aroma/',run))
      break;
    }
  }
}

