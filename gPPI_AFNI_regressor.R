# 2025-05-09 AndyP
# Generate AFNI formatted time files for a parametric regressor of entropy

library(tidyverse)

repo_directory <- "~/clock_analysis"

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
df <- get_trial_data(repo_directory=repo_directory,dataset='explore')
df0 <- df %>% filter(id == '202200' & run == 1) %>% select(trial,clock_onset, rt_csv,v_entropy)

str <- NULL; for (iD in 1:nrow(df0)){str <- paste0(str,paste0(df0$clock_onset[iD],":",df0$rt_csv[iD],' '))}


setwd('~/gPPI')

fileConn <- file("p.1D")
writeLines(str,fileConn)
close(fileConn)

fileConn <- file(paste0('202200','-','regressor-script.tcsh'))

str <-
'
#!/bin/tcsh \\ 

3dDeconvolve                                                                 \\
    -input           nfas-sub-202200_task-clockRev_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii  \\
    -polort          -1                                                      \\
    -num_stimts      1                                                       \\
    -GOFORIT         120                                                     \\
    -stim_times_AM1  1 p.1D "dmUBLOCK(1)"                                    \\
    -x1D             q_dmU.1D                                                
\\
3dDeconvolve                                                                 \\
    -input           nfas-sub-202200_task-clockRev_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii  \\
    -polort          -1                                                      \\
    -num_stimts      1                                                       \\
    -GOFORIT         120                                                     \\
    -stim_times_AM1  1 p.1D "dmUBLOCK(-1)"                                   \\
    -x1D             q_dmUn.1D                                               

\\
# simplify formatting: output number-only files \\
 \\
1dcat q_dmU.1D  > r_dmU.1D \\
1dcat q_dmUn.1D > r_dmUn.1D \\'

writeLines(str,fileConn)
close(fileConn)

system(paste0('tcsh ',paste0('202200','-','regressor-script.tcsh')))
