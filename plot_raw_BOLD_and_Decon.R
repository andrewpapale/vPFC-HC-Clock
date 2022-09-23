# 2022-09-20 AndyP
# Clock script (for now)
# Plot raw BOLD and raw decon on same graph.  Useful for checking alignment.

library(data.table)
library(tidyverse)
library(oro.nifti)
library(abind)
library(pracma)

medusa_dir <- '/Users/andypapale/vmPFC/Schaefer_deconvolved/Schaefer_Default' # should be a directory with a list of raw (unaligned) decon files 1 / subject
bold_dir <- '/Users/andypapale/vmPFC/MR_Proc/' # should be a directory with a NIFTI of post-preprocessed BOLD from all subjects
atlas_value1 <- 89 # what region or atlas_value do we want to look at?  Using Schaefer parcellation for now.
mask_dir <- '/Users/andypapale/vmPFC/Schaefer2018/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants.nii.gz'

subj <- '10637' # pick a random subject
run0 <- 1
align_evt <- 'feedback' # what do we want to align to

med_files <- list.files(path=medusa_dir)
bold_files <- list.files(path=bold_dir)

# /todo code to make sure subj is in med_files and bold_files
med_files <- med_files[1]

if (any(subj %like% bold_files)){
  bold_files <- list.files(paste0(bold_dir,subj,'/','func'),pattern='nfas')
}

# extract decon atlas_values to truncate bold_dir, average
decon_df <- fread(paste0(medusa_dir,'/',med_files))
decon_df <- decon_df %>% filter(atlas_value==atlas_value1)
decon_df <- decon_df %>% group_by(time) %>% summarize(decon_mean = mean(decon), decon_stderr = sd(decon)/sqrt(length(decon))) %>% ungroup()
# extract NIFTI voxels that match atlas_value, average

bold_df <- oro.nifti::readNIfTI(paste0(bold_dir,subj,'/func','/',bold_files))

mask <- readNIfTI(mask_dir, reorient=FALSE)
#ensure that mask is binary
mask[mask != atlas_value1] <- 0.0
mask[mask == atlas_value1] <- 1.0
aindices <- which(mask==1, arr.ind=TRUE)

v <- apply(bold_df, 4, function(x) {
  vox_in_mask <- x[aindices]
  c(mean=mean(vox_in_mask), sd=sd(vox_in_mask), N=sum(mask==1))
})

# to assign (not required here)
# mm <- repmat(aindices)
# my4d[mm] <- results

# another approach
# dlist <- lapply(1:dim(bold_df)[4], function(time) {
#   data <- bold_df[,,,time]
# })

# this is slow
# for (iT in 1:dim(bold_df)[4]){ 
#   aimg <- abind(aimg,bold_df[1:dim(bold_df)[1],1:dim(bold_df)[2],1:dim(bold_df)[3],iT]*mask,along=4)
# }
# a_indices <- which(aimg > zero_thresh, arr.ind=TRUE)

# plot on same axis (range is large enough to catch decon signal and BOLD signal (later))
merged_df1 <- data.frame(time=decon_df$time,decon_mean=decon_df$decon_mean,decon_stderr=decon_df$decon_stderr,bold_mean=v[1,],bold_stderr=v[2,]/sqrt(v[3,]))

source('~/vPFC-HC-Clock/get_trial_data_vmPFC.R')
df <- get_trial_data_vmPFC(repo_directory = '~/clock_analysis',dataset="mmclock_fmri")
df <- df %>% filter(id==subj & run==run0)

if (strcmp(align_evt,'feedback')){
  align_to <- round(df$feedback_onset)
}

bold <- matrix(as.double(NA),nrow=15+15+1,ncol=nrow(df))
decon <- matrix(as.double(NA),nrow=15+15+1,ncol=nrow(df))
for (iT in 1:nrow(df)){
    idx <- merged_df1$time >= align_to[iT]-15 & merged_df1$time <= align_to[iT]+15
    if (sum(idx)==15+15+1){
      bold[,iT] <- merged_df1$bold_mean[idx]
      decon[,iT] <- merged_df1$decon_mean[idx]
    } else if (min(merged_df1$time) > align_to[iT]-15){
      bold[,iT] <- c(rep(NA,min(merged_df1$time)-(align_to[iT]-15)),merged_df1$bold_mean[idx])
      decon[,iT] <- c(rep(NA,min(merged_df1$time)-(align_to[iT]-15)),merged_df1$decon_mean[idx])
    } else if (max(merged_df$time) < align_to[iT]+15){
      bold[,iT] <- c(merged_df1$bold_mean[idx],rep(NA,align_to[iT]+15-max(merged_df1$time)))
      decon[,iT] <- c(merged_df1$decon_mean[idx],rep(NA,align_to[iT]+15-max(merged_df1$time)))
    }
}

bold <- data.frame(bold)
bold <- bold %>% mutate(evt_time = -15:1:15)
bold <- bold %>% pivot_longer(cols = starts_with("X"))
bold <- bold %>% mutate(run_trial = rep(1:nrow(df),15+15+1))
bold <- bold %>% rename(bold = value)
bold$bold <- scale(bold$bold)
decon <- data.frame(decon)
decon <- decon %>% mutate(evt_time = -15:1:15)
decon <- decon %>% pivot_longer(cols = starts_with("X"))
decon <- decon %>% mutate(run_trial = rep(1:nrow(df),15+15+1))
decon <- decon %>% rename(decon = value)
decon$decon <- scale(decon$decon)
merged_df <- merge(decon,bold,by=c('evt_time','run_trial'))
merged_df <- merged_df %>% select(decon,bold,evt_time,run_trial) %>% arrange(run_trial,evt_time)

merged_df <- merged_df %>% pivot_longer(cols= decon | bold,names_to='data_type',values_to='data_values')

sm_df <- merged_df %>% group_by(evt_time,data_type) %>% summarize(mean0 = mean(data_values,na.rm=TRUE), sd0 = sd(data_values,na.rm=TRUE), N0 = sum(!is.na(data_values)))%>% ungroup()
c0 <- ccf(sm_df$mean0[sm_df$data_type=='decon'],sm_df$mean0[sm_df$data_type=='bold'],lag=14,pl=TRUE)
c0 <- data.frame(lag=c0$lag,acf=c0$acf)
setwd('~/Documents')
pdf(paste0(subj,'-','BOLD-decon.pdf'),height=6,width=12)
gg1 <- ggplot() + 
  geom_line(data=sm_df,aes(x=evt_time,y=mean0,color=data_type,group=data_type)) + 
              geom_smooth(data=sm_df,aes(x=evt_time,y=mean0,color=data_type,group=data_type)) +
  geom_line(data=c0,aes(x=lag,y=acf/2.5)) +
  scale_y_continuous(sec.axis=sec_axis(~.*2.5,name='acf'))
print(gg1)
dev.off()


