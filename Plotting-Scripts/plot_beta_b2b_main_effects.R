# plot_beta_emt
# 2022-08-01 AndyP
session = "fmri"
model = "slo"
clock_folder <- "~/clock_analysis" 
# source('~/code/Rhelpers/')
fmri_dir <- '~/vmPFC/MEDUSA Schaefer Analysis/parcel_maps_l2/'
setwd(fmri_dir)

library(grid)
library(readxl)
library(ggrepel)

if (session == "meg") {
  temp <- readRDS("/Volumes/Users/Andrew/parcel_maps_l2/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2_v_entropy_wi_meg_mixed_by.rds")  
} else if (session == "fmri") {
  temp <- readRDS("/Volumes/Users/Andrew/parcel_maps_l2/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2_v_entropy_wi_fmri_mixed_by.rds")
}

library(wesanderson)
pal = wes_palette("FantasticFox1", 3, type = "discrete")
pal_vPFC = palette()
pal_vPFC[2] = '#dd8d29'
pal_vPFC[3] = '#e2d200'
pal_vPFC[1] = '#46acc8'


df <- temp$coef_df_reml

if (model=='slo'){
  df <- df %>% filter(model_name=='slo' & l1_cope_name=='EV_entropy_wiz_clock' & mask_value %in% c(67,171,65,170,66,89,194,88,192,84,191,86,161,55,160,159,56))
} else if (model=='int'){
  df <- df %>% filter(model_name=='int' & l1_cope_name=='EV_entropy_wiz_clock' & mask_value %in% c(67,171,65,170,66,89,194,88,192,84,191,86,161,55,160,159,56))
}
label_vPFC <- read_excel('~/vmPFC/PFC_region_gradients.xlsx')
label_vPFC <- label_vPFC %>% select(atlas_value, region,network) %>% mutate(mask_value=as.factor(atlas_value))
df$mask_value <- as.factor(df$mask_value)
df <- inner_join(df,label_vPFC, by='mask_value')








