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

# Get a list of full directory paths for all subjects in the specified base pipeline folder
subj <- list.files('/Volumes/Users/Andrew/ABCC_MID_L2/fmriresults01/derivatives/abcd-hcp-pipeline',full.names = TRUE)

# Calculate the total number of subjects to establish the upper limit for the loop
nSubj <- length(subj)

# Initialize a loop to iterate through each subject by their numerical index
for (iS in 1:nSubj){
  
  # Output the current subject's directory path to the console to track the script's progress
  print(subj[iS])
  
  # Change the active working directory to the current subject's root folder
  setwd(subj[iS])
  
  # List the full paths of the contents inside the subject directory (expected to be the session folder)
  ses0 <- list.files(subj[iS],full.names=TRUE)
  
  # Change the working directory to navigate inside that specific session folder
  setwd(ses0)
  
  # List the full paths of the contents inside the session folder (expected to be the functional data folder)
  func0 <- list.files(ses0,full.names=TRUE)
  
  # Change the working directory to navigate inside the functional data folder
  setwd(func0)
  
  # Search the current working directory for files matching the pattern "*29*.nii" 
  # Returning full.names=FALSE means it only grabs the filename, not the path
  file <- list.files(getwd(),pattern=glob2rx("*29*.nii"),full.names=FALSE)
  
  # Proceed only if at least one matching file was found in the directory
  if (length(file)>0){
    
    # Strip the file extension to get the base name by splitting the string at the first period
    file_name <- strsplit(file, "\\.")[[1]][1]
    
    # Run a shell command: Use Connectome Workbench to separate the LEFT cortex data from the CIFTI file into a .func.gii metric file
    system(paste0('wb_command -cifti-separate ', file, ' COLUMN -metric CORTEX_LEFT ', file_name,'_cortex_left.func.gii'))
    
    # Run a shell command: Use Connectome Workbench to separate the RIGHT cortex data from the CIFTI file into a .func.gii metric file
    system(paste0('wb_command -cifti-separate ', file, ' COLUMN -metric CORTEX_RIGHT ', file_name,'_cortex_right.func.gii'))
    
    # Run a shell command: Use AFNI's 3dROIstats to extract summary statistics for the RIGHT hemisphere using an atlas mask, saving output to a .txt file
    system(paste0('3dROIstats -mask ', atlas, ' ', file_name,'_cortex_right.func.gii > ', file_name,'_cortex_right.txt'))
    
    # Run a shell command: Use AFNI's 3dROIstats to extract summary statistics for the LEFT hemisphere using an atlas mask, saving output to a .txt file
    system(paste0('3dROIstats -mask ', atlas, ' ', file_name,'_cortex_left.func.gii > ', file_name,'_cortex_left.txt'))
  }
}

# Recursively search the pipeline directory for all .txt files containing "29" in the name
# This gathers all the ROI statistic files generated in the previous step into one list
L2roifiles <- list.files('/Volumes/Users/Andrew/ABCC_MID_L2/fmriresults01/derivatives/abcd-hcp-pipeline', 
                         pattern = glob2rx("*29*.txt"), recursive = TRUE, full.names = TRUE)

# Initialize an empty object to store the final aggregated data
df <- NULL

# Determine the total number of text files found to set the loop limit
nSubj <- length(L2roifiles)

# Loop through each individual text file in the list
for (iS in 1:nSubj){
  
  # Print the current file path to the console for progress tracking
  print(L2roifiles[iS])
  
  # Read the tab-delimited text file into a temporary data frame
  curr_file <- read.table(L2roifiles[iS], sep="\t", header=TRUE)
  
  # Process the data using tidyverse functions:
  # 1. pivot_longer: Convert columns starting with "Mean_" (the ROIs) into a long format (name and value)
  # 2. mutate (atlas_value): Extract the numeric ROI ID from the column name (e.g., "Mean_12" becomes 12)
  # 3. mutate (id): Extract the subject ID from the "File" column using regex (looking for characters after the "-")
  # 4. select: Remove the original 'File', 'name', and 'Sub.brick' columns to keep only the clean data
  df0 <- curr_file %>% pivot_longer(cols = starts_with("Mean_")) %>% 
    mutate(atlas_value = as.numeric(str_extract(name, "(?<=_)[0-9]+")), id = str_extract(File, "[^_]+")) %>%
    mutate(id = str_extract(id, "(?<=-)[A-za-z0-9]+")) %>%
    select(!File & !name & !Sub.brick)
  
  # Check if the current filename contains the string '_right'
  if (grepl('_right', L2roifiles[iS])){
    # If it is the right hemisphere, add 200 to the atlas value to differentiate it from the left hemisphere IDs
    df0 <- df0 %>% mutate(atlas_value = atlas_value + 200)
  }
  
  # Reorder the columns: moves 'id' to the first column, followed by 'atlas_value' and the data 'value'
  df0 <- df0[, c(3,2,1)]
  
  # Append the current subject's processed data to the master data frame (df)
  df <- rbind(df, df0)
}

setwd('/Volumes/Users/Andrew/ABCC_MID_L2/fmriresults01/derivatives/abcd-hcp-pipeline')
write.csv(df, file='2025-04-25-ABCC_MID_L2_contrast_29.csv')
