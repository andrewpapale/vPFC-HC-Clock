# 2025-04-21 AndyP
# Extract L2 betas from the ABCD MID dataset.  These files are CIFTIs and need to be converted to GIFTIS
# The conversion takes place using tools from the connectomics workbench: https://www.humanconnectome.org/software/connectome-workbench
# Then, use AFNI to extract L2 betas in a Schaefer mask
library(tidyverse)

Sys.setenv(PATH = paste(Sys.getenv("PATH"),"/Users/dnplserv/abin/", sep = ":"))
datadir <- '/Volumes/Users/Andrew/ABCC_MID_L2/fmriresults01/derivatives/abcd-hcp-pipeline'

# I did two passes on the data, one with the fsLR_L.gii atlas and one with the resampled atlas (but stopped partway through as I was getting mostly warnings about dimension mismatches). 
# It is possible different subjects were computed with the resampled atlas...I resampled with 3dresample, but Jo Etzel advises against this.  The original contrast 10 .nii CIFTI file was in a different dimension.
# Apparently there are some "standard" CIFTI file dimensions that we should be able to find for Schaefer masks.
#atlas <- '/Volumes/Users/Andrew/ABCC_MID_L2/Schaefer2018_400Parcels_7Networks_order_fsLR_resampled.gii'
atlas <- '/Volumes/Users/Andrew/ABCC_MID_L2/Schaefer2018_400Parcels_7Networks_order_fsLR_L.label.gii'

subj <- list.files('/Volumes/Users/Andrew/ABCC_MID_L2/fmriresults01/derivatives/abcd-hcp-pipeline',full.names = TRUE)
nSubj <- length(subj)
for (iS in 1:nSubj){
  print(subj[iS])
  setwd(subj[iS])
  ses0 <- list.files(subj[iS],full.names=TRUE)
  setwd(ses0)
  func0 <- list.files(ses0,full.names=TRUE)
  setwd(func0)
  file <- list.files(getwd(),pattern=glob2rx("*29*.nii"),full.names=FALSE)
  if (length(file)>0){
    file_name <- strsplit(file, "\\.")[[1]][1]
    system(paste0('wb_command -cifti-separate ', file, ' COLUMN -metric CORTEX_LEFT ', file_name,'_cortex_left.func.gii'))
    system(paste0('wb_command -cifti-separate ', file, ' COLUMN -metric CORTEX_RIGHT ', file_name,'_cortex_right.func.gii'))
    system(paste0('3dROIstats -mask ', atlas, ' ', file_name,'_cortex_right.func.gii > ', file_name,'_cortex_right.txt'))
    system(paste0('3dROIstats -mask ', atlas, ' ', file_name,'_cortex_left.func.gii > ', file_name,'_cortex_left.txt'))
  }
}

L2roifiles <- list.files('/Volumes/Users/Andrew/ABCC_MID_L2/fmriresults01/derivatives/abcd-hcp-pipeline', pattern = glob2rx("*29*.txt"), recursive = TRUE, full.names = TRUE)

df <- NULL
nSubj <- length(L2roifiles)
for (iS in 1:nSubj){
  print(L2roifiles[iS])
  curr_file <- read.table(L2roifiles[iS], sep="\t", header=TRUE)
  df0 <- curr_file %>% pivot_longer(cols = starts_with("Mean_")) %>% 
    mutate(atlas_value = as.numeric(str_extract(name, "(?<=_)[0-9]+")), id = str_extract(File, "[^_]+")) %>%
    mutate(id = str_extract(id, "(?<=-)[A-za-z0-9]+")) %>%
    select(!File & !name & !Sub.brick)
  
  if (grepl('_right',L2roifiles[iS])){
    df0 <- df0 %>% mutate(atlas_value = atlas_value + 200)
  }
  
  df0 <- df0[, c(3,2,1)]
  df <- rbind(df, df0)
}

setwd('/Volumes/Users/Andrew/ABCC_MID_L2/fmriresults01/derivatives/abcd-hcp-pipeline')
write.csv(df, file='2025-04-25-ABCC_MID_L2_contrast_29.csv')
