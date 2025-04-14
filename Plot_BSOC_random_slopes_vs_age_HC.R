# 2025-04-14 AndyP
# Plotting random slopes vs. age and AH/PH for devlopmental project

# libraries we'll need
library(tidyverse)
library(fmri.pipeline)
library(MplusAutomation)
# set root directory
rootdir1 <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/MMClock_BSOC_MPlus'
rootdir <- '/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/BSOCIAL'
repo_directory <- file.path('/Volumes/Users/Andrew/MEDuSA_data_BSOC')
repo_directory1 <- '~/clock_analysis'
ncores <- 26
# load mixed_by function for analyses

##################################
##### Load in and format data ####
#####     Clock-aligned     ######
##################################

load("/Users/dnplserv/vmPFC/MEDUSA Schaefer Analysis/vmPFC_HC_model_selection/2025-03-05-vmPFC-HC-network-clock-ranslopes-nofixedeffect-noHCbetween-3.Rdata")
Q <- ddf$coef_df_reml %>% filter(effect=='ran_vals' & term=='HCwithin')
Q <- Q %>% rename(id=level)
Q$id <- as.character(Q$id)

Q1ah <- Q %>% filter(HC_region == "AH")
Q1ph <- Q %>% filter(HC_region == "PH")