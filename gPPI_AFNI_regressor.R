# 2025-05-09 AndyP
# Generate AFNI formatted time files for a parametric regressor of entropy

library(tidyverse)

repo_directory <- "~/clock_analysis"
thread = 26
dofeedback = FALSE
doclock = TRUE
dodf = FALSE

demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
HC <- demo %>% filter(registration_group == 'HC')
ids <- HC$registration_redcapid

par_cl <- parallel::makeCluster(spec = thread,type = "FORK")
parallel::parLapply(par_cl,ids,function(iD){
  
  source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
  df <- get_trial_data(repo_directory=repo_directory,dataset='explore')
  
  for (run0 in 1:2){
    
    df0 <- df %>% filter(id == iD & run == run0) %>% select(trial,clock_onset, rt_csv)
    
    if (nrow(df0) > 0){
      
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
          dir2[iC] <- dir1[iD0]
          iC <- iC + 1
        }
      }
      
      checkmate::assert(length(dir2)==1)
      
      if (doclock){
        
        str <- NULL; for (iQ in 1:nrow(df0)){str <- paste0(str,paste0(df0$clock_onset[iQ],":",df0$rt_csv[iQ],' '))}
        
        fileConn <- file("clock_onset_regressor.1D")
        writeLines(str,fileConn)
        close(fileConn)
        
        fileConn <- file(paste0(iD,'-run-',run0,'-clock_onset_regressor_script.tcsh'))
        
        str <- paste0(
          '
#!/bin/tcsh

3dDeconvolve                                                                 \\
    -input          ',dir2,'  \\
    -polort          -1                                                      \\
    -num_stimts      1                                                       \\
    -GOFORIT         120                                                     \\
    -stim_times_AM1  1 clock_onset_regressor.1D "dmUBLOCK(1)"                                    \\
    -x1D             clock_onset_regressor_dmU.1D                                                
\\
3dDeconvolve                                                                 \\
    -input          ',dir2,'  \\
    -polort          -1                                                      \\
    -num_stimts      1                                                       \\
    -GOFORIT         120                                                     \\
    -stim_times_AM1  1 clock_onset_regressor.1D "dmUBLOCK(-1)"                                   \\
    -x1D             clock_onset_regressor_dmUn.1D                                               

\\
# simplify formatting: output number-only files

1dcat clock_onset_regressor_dmU.1D  > r_clock_onset_regressor_dmU_',run0,'.1D

1dcat clock_onset_regressor_dmUn.1D > r_clock_onset_regressor_dmUn_',run0,'.1D'
        )      
        writeLines(str,fileConn)
        close(fileConn)
        
        system(paste0('tcsh ',paste0(iD,'-run-',run0,'-clock_onset_regressor_script.tcsh')))
      }
      
    }
    
    if (dofeedback){
      
      df0 <- df %>% filter(id == iD & run == run0) %>% select(trial,feedback_onset,v_entropy)
      
      if (nrow(df0) > 0){
        
        
        checkmate::assert(length(dir2)==1)
        
        str <- NULL; for (iQ in 1:nrow(df0)){str <- paste0(str,paste0(df0$feedback_onset[iQ],":",0.7,' '))}
        fileConn <- file("feedback_regressor.1D")
        writeLines(str,fileConn)
        close(fileConn)
        
        fileConn <- file(paste0(iD,'-run-',run0,'-feedback_regressor_script.tcsh'))
        
        str <- paste0(
          '
 #!/bin/tcsh 
 
 3dDeconvolve                                                                 \\
     -input          ',dir2,'  \\
     -polort          -1                                                      \\
     -num_stimts      1                                                       \\
     -GOFORIT         120                                                     \\
     -stim_times_AM1  1 feedback_regressor.1D "dmUBLOCK(1)"                                    \\
     -x1D             feedback_regressor_dmU.1D                                                
 \\
 3dDeconvolve                                                                 \\
     -input          ',dir2,'  \\
     -polort          -1                                                      \\
     -num_stimts      1                                                       \\
     -GOFORIT         120                                                     \\
     -stim_times_AM1  1 feedback_regressor.1D "dmUBLOCK(-1)"                                   \\
     -x1D             feedback_regressor_dmUn.1D                                               
 
 \\
 # simplify formatting: output number-only files \\
 
 1dcat feedback_regressor_dmU.1D  > r_feedback_regressor_dmU_',run0,'.1D

 1dcat feedback_regressor_dmUn.1D > r_feedback_regressor_dmUn_',run0,'.1D'
        )      
        writeLines(str,fileConn)
        close(fileConn)
        
        system(paste0('tcsh ',paste0(iD,'-run-',run0,'-feedback_regressor_script.tcsh')))
      }
    }
    
    if (dodf){
      
      # write entropy (will construct into parametric modulator using MATLAB/SPM)
      
      df0 <- df %>% filter(id == iD & run == run0) %>% select(trial,clock_onset,v_entropy)
      write.table(df0,file=paste0(iD,'-run-',run0,'-entropy-PM.csv'),col.names=FALSE)
    }
  }
})
parallel::stopCluster(par_cl)