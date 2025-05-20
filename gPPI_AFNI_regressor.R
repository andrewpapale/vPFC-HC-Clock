# 2025-05-09 AndyP
# Generate AFNI formatted time files for a parametric regressor of entropy

library(tidyverse)

repo_directory <- "~/clock_analysis"


demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
HC <- demo %>% filter(registration_group == 'HC')
ids <- HC$registration_redcapid

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='explore')

for (iD in ids){
  for (run0 in 1:2){
    df0 <- df %>% filter(id == iD & run == run0) %>% select(trial,clock_onset, rt_csv,v_entropy)
    
    if (nrow(df0) > 0){
      
      str <- NULL; for (iD in 1:nrow(df0)){str <- paste0(str,paste0(df0$clock_onset[iD],":",df0$rt_csv[iD],' '))}
      
      
      setwd(paste0('/Users/dnplserv/gPPI/Explore_HC_only/sub-',iD,'/func'))
      
      fileConn <- file("clock_onset_regressor.1D")
      writeLines(str,fileConn)
      close(fileConn)
      
      fileConn <- file(paste0(iD,'-run-',run0,'-clock_onset_regressor-script.tcsh'))
      
      str <- paste0(
        '
#!/bin/tcsh \\ 

3dDeconvolve                                                                 \\
    -input           ',dir2,'  \\
    -polort          -1                                                      \\
    -num_stimts      1                                                       \\
    -GOFORIT         120                                                     \\
    -stim_times_AM1  1 clock_onset_regressor.1D "dmUBLOCK(1)"                                    \\
    -x1D             clock_onset_regressor_dmU.1D                                                
\\
3dDeconvolve                                                                 \\
    -input           ',dir2,'  \\
    -polort          -1                                                      \\
    -num_stimts      1                                                       \\
    -GOFORIT         120                                                     \\
    -stim_times_AM1  1 clock_onset_regressor.1D "dmUBLOCK(-1)"                                   \\
    -x1D             clock_onset_regressor_dmUn.1D                                               

\\
# simplify formatting: output number-only files \\
 \\
1dcat clock_onset_regressor_dmU.1D  > r_dmU.1D \\
1dcat clock_onset_regressor_dmUn.1D > r_dmUn.1D \\'
      )      
      writeLines(str,fileConn)
      close(fileConn)
      
      system(paste0('tcsh ',paste0(iD,'-run-',run0,'-clock_onset_regressor_script.tcsh')))
    }
    
    
    
    df0 <- df %>% filter(id == iD & run == run0) %>% select(trial,clock_onset,v_entropy)
    
    if (nrow(df0) > 0){
      
      str <- NULL; for (iD in 1:nrow(df0)){str <- paste0(str,paste0(df0$clock_onset[iD],":",0.7,' '))}
      
      
      setwd(paste0('/Users/dnplserv/gPPI/Explore_HC_only/sub-',iD,'/func'))
      dir0 <- list.files(pattern = 'nfas')
      
      dir1 <- NULL
      iC <- 1
      for (iD0 in 1:length(dir0)){
        if (grepl('.nii.gz',dir0[iD0])){
          dir1[iC] <- dir0[iD0]
          iC <- iC + 1
        }
      }
      
      dir2 <- NULL
      iC <- 1
      for (iD0 in 1:length(dir1)){
        if (grepl(paste0('run-',run0),dir1[iD0])){
          dir2[iC] <- dir0[iD0]
          iC <- iC + 1
        }
      }
      
      checkmate::assert(length(dir2)==1)
      
      fileConn <- file("feedback_regressor.1D")
      writeLines(str,fileConn)
      close(fileConn)
      
      fileConn <- file(paste0(iD,'-run-',run0,'-feedback_regressor-script.tcsh'))
      
      str <- paste0(
        '
#!/bin/tcsh \\ 

3dDeconvolve                                                                 \\
    -input           ',dir2,'  \\
    -polort          -1                                                      \\
    -num_stimts      1                                                       \\
    -GOFORIT         120                                                     \\
    -stim_times_AM1  1 feedback_regressor.1D "dmUBLOCK(1)"                                    \\
    -x1D             feedback_regressor_dmU.1D                                                
\\
3dDeconvolve                                                                 \\
    -input           ',dir2,'  \\
    -polort          -1                                                      \\
    -num_stimts      1                                                       \\
    -GOFORIT         120                                                     \\
    -stim_times_AM1  1 feedback_regressor.1D "dmUBLOCK(-1)"                                   \\
    -x1D             feedback_regressor_dmUn.1D                                               

\\
# simplify formatting: output number-only files \\
 \\
1dcat feedback_regressor_dmU.1D  > r_dmU.1D \\
1dcat feedback_regressor_dmUn.1D > r_dmUn.1D \\'
      )      
      writeLines(str,fileConn)
      close(fileConn)
      
      system(paste0('tcsh ',paste0(iD,'-run-',run0,'-feedback_regressor_script.tcsh')))
    }
    
    
    
  }
}