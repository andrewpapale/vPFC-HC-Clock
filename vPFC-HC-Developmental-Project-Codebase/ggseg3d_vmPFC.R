# 2026-10-19 AndyP
# plot 3d rendered images of 17 Network Yeo Schaefer 400 parcellated cortex
# and Choi 2012 Striatum

# https://github.com/ggsegverse/ggsegSchaefer
#pak::pak("ggsegverse/ggsegSchaefer")
#install.packages("ggseg3d")

# 2026-11-20 This script was used to generate Figure 2a, Figure 2b Figure 5b, in the original submission for the CLIMATE Corticostriatal Paper
# Papale, et al (circa ~ 2026-2027) Corticostriatal and Cortical Network Reorganization Associated with Clozapine Response in Treatment-Resistant Schizophrenia
# I may have had to manually rearrange some colors in the masks (in Inkscape) to get everything to match


# Render an interactive 3D mesh of the 17-network (400 parcels variant)
library(ggseg3d)
library(ggsegSchaefer)
library(ggbrain)
# Use the 3d data function
plot(schaefer17_400())

library(ggseg3d)
library(tidyverse)

# Standard Yeo/Schaefer 17-Network Hex Palette 
# (Ordered Network 1 to Network 17)
# yeo17_colors <- c(
#   #"NONE"          = "#000000",
#   "VisCent"       = "#781286",
#   "VisPeri"       = "#FF0000",
#   "SomMotA"       = "#4682B4",
#   "SomMotB"       = "#2ACCA4",
#   "DorsAttnA"     = "#4A9B3C",
#   "DorsAttnB"     = "#00760E",
#   "SalVentAttnA"  = "#C43AFA", 
#   "SalVentAttnB"  = "#FF98D5",
#   "LimbicB"       = "#7A8732", 
#   "LimbicA"       = "#DCF8A4", 
#   "ContA"         = "#E69422", 
#   "ContB"         = "#87324A", 
#   "ContC"         = "#778CB0",
#   "DefaultA"      = "#FFFF00",
#   "DefaultB"      = "#CD3E4E",
#   "DefaultC"      = "#000082",
#   "TempPar"       = "#0C30FF" 
# )

library(wesanderson)
pal = wes_palette("FantasticFox1", 3, type = "discrete")

names(yeo17_colors) <- as.character(1:17)


# started with the FSL whole brain template, original idea was to plot a translucent whole brain with the Choi 2012 17 Network in opaque colors matching
# the original Yeo 17 network colors.  We decided to eliminate the translucent whole brian.


# 2. Establish clean file paths
underlay_path <- "~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Data/masks/FSL_MNI152_FreeSurferConformed_1mm.nii.gz"
vmPFC_path <- '~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Data/masks/vmPFC_mask.nii.gz'

library(oro.nifti)
library(neurobase)
library(rgl)
library(misc3d)
library(RNifti)

# 1. Read your NIFTI files
# Ensure 'cortex.nii.gz' and 'striatum.nii.gz' are aligned in the same spatial space (e.g., MNI152)
cortex_nifti   <- RNifti::readNifti(underlay_path)
vmPFC_nifti <- RNifti::readNifti(vmPFC_path)

# 2. Extract raw arrays for contouring
# RNifti images act like standard 3D arrays directly compatible with misc3d
cortex_arr   <- as.array(cortex_nifti)
vmpfc_arr <- array(as.numeric(vmPFC_nifti), dim = dim(vmPFC_nifti))
vmpfc_arr <- drop(vmpfc_arr)


# 3. Open interactive 3D RGL Canvas
open3d()
bg3d("white")

# 4. Render Transparent Cortex Surface
contour3d(
  cortex_arr,
  level = 0.25,           # Threshold boundary (adjust if probability map)
  color = "peachpuff",   # Natural outer shell tone
  alpha = 0.05,          # High transparency to view interior
  smooth = FALSE,
  add = TRUE
)

iC <- 1
# 4. Loop through each Yeo network label (1 through 17)
# These are the 9 networks in the Choi 17 network striatal mask
map_roi2col <- c(3,4,7,8,9,11,12,14,15)
for (roi_val in c(3,4,7,8,10,12,13,16,17)) {
  
  # Create a temporary binary mask for the current network label
  roi_mask <- array(0, dim = dim(vmpfc_arr))
  roi_mask[vmpfc_arr == roi_val] <- 1
  
  # Only attempt to plot if the network label exists in your striatum volume
  if (max(roi_mask) > 0) {
    print(roi_val)
    print(yeo17_colors[map_roi2col[iC]])
    contour3d(
      roi_mask,
      level = 0.5,
      color = yeo17_colors[map_roi2col[iC]],
      alpha = 1,        # Fully opaque for solid colors
      smooth = TRUE,
      add = TRUE,
      specular = "black",
      shininess = 0,
    )
  }
  iC <- iC + 1
}


library(misc3d)
library(rgl)

setwd('~/Desktop')

par3d(windowRect = c(0, 0, 1200, 1200))

# 2. Capture the exact window output into a PDF structure, this doesn't always work for some reason.  I think for the final figures I just took screenshots.
snapshot3d("vmPFC.pdf", fmt = "pdf", webshot = TRUE)



# # 4. Loop through each Yeo network label (1 through 17)
iC <- 1
# 4. Loop through each Yeo network label (1 through 17)
# Create a temporary binary mask for the current network label
#roi_mask <- array(0, dim = dim(striatum_arr))
#roi_mask[striatum_arr == roi_val] <- 1

#for (roi_val in c(3,4,7,8,10,12,13,16,17)) {

# Create a temporary binary mask for the current network label
roi_mask <- array(0, dim = dim(striatum_arr))
roi_mask[striatum_arr > 0] <- 1


# Only attempt to plot if the network label exists in your striatum volume
if (max(roi_mask) > 0) {
  #print(roi_val)
  #print(yeo17_colors[map_roi2col[iC]])
  contour3d(
    roi_mask,
    level = 0.3,
    color = "grey60",
    alpha = 0.1,        # Fully opaque for solid colors
    smooth = 30,
    add = TRUE,
    specular = "black",
    shininess = 0
  )
}
iC <- iC + 1

#}

# # Create a temporary binary mask for the current network label
# mask2 <- array(0, dim = dim(lseed_arr))
# mask2[lseed_arr == 1] <- 1
# 
# # Create a temporary binary mask for the current network label
# mask3 <- array(0, dim = dim(rseed_arr))
# mask3[rseed_arr == 1] <- 1
# 
# contour3d(
#   mask3,
#   level = 0.5,
#   color = 'black',
#   alpha = 1,        # Fully opaque for solid colors
#   smooth = TRUE,
#   add = TRUE
# )  
# 
# contour3d(
#   mask2,
#   level = 0.5,
#   color = 'black',
#   alpha = 1,        # Fully opaque for solid colors
#   smooth = TRUE,
#   add = TRUE
# )  

# # 4. Loop through each Yeo network label (1 through 17)
iC <- 1
# 4. Loop through each Yeo network label (1 through 17)
for (roi_val in c(3,4,7,8,10,12,13,16,17)) {
  
  # Create a temporary binary mask for the current network label
  roi_mask <- array(0, dim = dim(ass_str_arr))
  roi_mask[ass_str_arr == roi_val] <- 1
  
  
  # Only attempt to plot if the network label exists in your striatum volume
  if (max(roi_mask) > 0) {
    print(roi_val)
    print(yeo17_colors[map_roi2col[iC]])
    contour3d(
      roi_mask,
      level = 0.5,
      color = yeo17_colors[map_roi2col[iC]],
      alpha = 1,        # Fully opaque for solid colors
      smooth = 10,
      add = TRUE,
      specular = "black",
      shininess = 0
    )
  }
  iC <- iC + 1
  
}


###########################
# dorsal caudate ##########
###########################

# # 4. Loop through each Yeo network label (1 through 17)
iC <- 1
# 4. Loop through each Yeo network label (1 through 17)
#for (roi_val in c(3,4,7,8,10,12,13,16,17)) {

# Create a temporary binary mask for the current network label
roi_mask <- array(0, dim = dim(lseed_arr))
roi_mask[lseed_arr > 0] <- 1


# Only attempt to plot if the network label exists in your striatum volume
if (max(roi_mask) > 0) {
  #print(roi_val)
  #print(yeo17_colors[map_roi2col[iC]])
  contour3d(
    roi_mask,
    level = 0.5,
    color = "black",
    alpha = 1,        # Fully opaque for solid colors
    smooth = 5,
    add = TRUE,
    specular = "black",
    shininess = 0
  )
}
iC <- iC + 1

#}

# # 4. Loop through each Yeo network label (1 through 17)
iC <- 1
# 4. Loop through each Yeo network label (1 through 17)
#for (roi_val in c(3,4,7,8,10,12,13,16,17)) {

# Create a temporary binary mask for the current network label
roi_mask <- array(0, dim = dim(rseed_arr))
roi_mask[rseed_arr > 0] <- 1


# Only attempt to plot if the network label exists in your striatum volume
if (max(roi_mask) > 0) {
  #print(roi_val)
  #print(yeo17_colors[map_roi2col[iC]])
  contour3d(
    roi_mask,
    level = 0.5,
    color = "black",
    alpha = 1,        # Fully opaque for solid colors
    smooth = TRUE,
    add = TRUE,
    specular = "black",
    shininess = 0
  )
}
iC <- iC + 1

#}




