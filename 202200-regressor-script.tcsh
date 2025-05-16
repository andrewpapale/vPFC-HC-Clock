
#!/bin/tcsh 

3dDeconvolve                                                                 \
    -input           nfas-sub-202200_task-clockRev_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii  \
    -polort          -1                                                      \
    -num_stimts      1                                                       \
    -GOFORIT         120                                                     \
    -stim_times_AM1  1 p.1D "dmUBLOCK(1)"                                    \
    -x1D             q_dmU.1D                                                

3dDeconvolve                                                                 \
    -input           nfas-sub-202200_task-clockRev_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii  \
    -polort          -1                                                      \
    -num_stimts      1                                                       \
    -GOFORIT         120                                                     \
    -stim_times_AM1  1 p.1D "dmUBLOCK(-1)"                                   \
    -x1D             q_dmUn.1D                                               


# simplify formatting: output number-only files

1dcat q_dmU.1D  > r_dmU.1D
1dcat q_dmUn.1D > r_dmUn.1D
