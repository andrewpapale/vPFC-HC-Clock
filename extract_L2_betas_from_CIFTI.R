# 2025-04-21 AndyP
# Extract L2 betas from the ABCD MID dataset.  These filese are CIFTIs and need to be converted to GIFTIS
# The conversion takes place using tools from the connectomics workbench: https://www.humanconnectome.org/software/connectome-workbench
# Then, use AFNI to extract L2 betas in a Schaefer mask
Sys.setenv(PATH = paste(Sys.getenv("PATH"),"/Users/dnplserv/abin/", sep = ":"))
datadir <- '/Volumes/Users/Andrew/ABCD_MID_Task_L2s/abcd-hcp-pipeline'
atlas <- '/Volumes/Users/Andrew/ABCD_MID_Task_L2s/Schaefer2018_400Parcels_7Networks_order_fsLR_L.label.gii'

subj <- list.files('/Volumes/Users/Andrew/ABCD_MID_Task_L2s/abcd-hcp-pipeline',full.names = TRUE)
nSubj <- length(subj)
for (iS in 1:nSubj){
    disp(subj[iS])
    setwd(subj[iS])
    ses0 <- list.files(subj[iS],full.names=TRUE)
    setwd(ses0)
    func0 <- list.files(ses0,full.names=TRUE)
    setwd(func0)
    file <- list.files(func0,pattern="\\.nii",full.names=FALSE)
    file_name <- strsplit(file, "\\.")[[1]][1]
    system(paste0('wb_command -cifti-separate ', file, ' COLUMN -metric CORTEX_LEFT ', file_name,'_cortex_left.func.gii'))
    system(paste0('wb_command -cifti-separate ', file, ' COLUMN -metric CORTEX_RIGHT ', file_name,'_cortex_right.func.gii'))
    system(paste0('3dROIstats -mask ', atlas, ' ', file_name,'_cortex_right.func.gii > ', file_name,'_cortex_right.txt'))
    system(paste0('3dROIstats -mask ', atlas, ' ', file_name,'_cortex_left.func.gii > ', file_name,'_cortex_left.txt'))
}
