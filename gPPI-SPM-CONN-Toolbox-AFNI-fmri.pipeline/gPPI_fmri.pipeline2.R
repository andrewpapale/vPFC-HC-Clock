# 2025-06-03 AndyP and Angela I
# script to run gPPI in fmri.pipeline

# PPI prototype

#Sys.setenv(PATH = paste("/usr/lib64/", Sys.getenv("PATH"), sep = ":"))
#Sys.setenv(LD_LIBRARY_PATH = paste("/usr/lib64/libjpeg.so", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
library(data.table)
library(tidyr)
library(dplyr)
library(fmri.pipeline)
library(readr)
library(stringr)

create_new_gpa_object = TRUE

#Sys.setenv(PATH = paste("/ihome/crc/install/fsl/6.0.7/fsl/bin", Sys.getenv("PATH"), sep = ":"))
#Sys.setenv(FSLDIR = paste("/ihome/crc/install/fsl/6.0.7/fsl", Sys.getenv("FSLDIR")))
#Sys.setenv(PATH = paste("/ihome/crc/install/fsl/6.0.7/fsl/bin", Sys.getenv("PATH"), sep = ":"))

run_nifti <- list.files('/ix1/adombrovski/DNPL_DataMesh/Data/MMClock/MMClock_lite',pattern = "nfaswuktm_clock[1-8]_5_drop2_trunc[0-9]+\\.nii\\.gz$",full.names=TRUE,recursive=TRUE)

#physio_data <- data.table::fread("/proj/mnhallqlab/projects/mmy3_gppi_prototype/2025-05-29-HC_MEDuSA-gPPI-Prototype.csv")
#physio_data <- data.table::fread("/ix1/adombrovski/papalea/2025-07-07-gPPI/2025-06-27-gPPI-demo.csv")
#str(physio_data)
ppi_data <- data.table::fread('/ix1/adombrovski/papalea/2025-07-07-gPPI/2025-08-21-MMClock-decon.csv')
ppi_data$run<- as.character((ppi_data$run))

# for simplicity, pare this down and reshape it -- not quire sure what the intent is here since these are not deconvolved time series
# ppi_data <- physio_data %>%
#     select(id, run)

#elig_files <- list.files("/Volumes/bierka_root/datamesh/RAW/MMClock/MR_Proc", pattern = "nfaswuktm_clock[1-8]_5_drop2_trunc[0-9]+\\.nii\\.gz$",
#                         #elig_files <- list.files("/proj/mnhallqlab/projects/mmy3_gppi_prototype/MMClock/MR_Proc", pattern = "nfaswuktm_clock.*", 
#                         full.names = TRUE, recursive=TRUE
#)

# build_info <- function(input) {
#     require(RNifti)

#     id <- as.integer(sub(".*/MR_Proc/([^_]+)_.*", "\\1", input))
#     run_number <- as.integer(sub(".*/MR_Proc/[^/]+/mni_5mm_aroma/clock([1-8]+)/.*", "\\1", input))
#     # f <- readNifti(input, internal = TRUE)
#     # volumes <- dim(f)[4]
#     volumes <- lookup_run_volumes(input)
#     return(list(id=id, run_number=run_number, volumes=volumes))
# }


# info_list <- do.call(rbind.data.frame, lapply(elig_files, function(e) build_info(e)))

# # fake ppi data for now
# ppi_data <- info_list %>%
#     rowwise() %>%
#     reframe(id = id, run_number = run_number, volume = 1:volumes, hippo_ph = rnorm(volumes), hippo_ah = rnorm(volumes)) %>%
#     ungroup()

# real data
#ppi_data <- fread("/proj/mnhallqlab/projects/mmy3_gppi_prototype/MMClock_HippoDecon_27Jun2025.csv.gz")


#elig_files <- read.csv('/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/2025-07-03-MMClock-MRProc.csv')

ppi_data <- read.csv('/ix1/adombrovski/papalea/2025-07-07-gPPI/2025-08-21-MMClock-decon.csv')
ppi_data$volume <- ppi_data$time + 1
ppi_data$time <- NULL
ppi_data <- ppi_data %>%
  rename(run_number = run) %>%
  pivot_wider(names_from="HC_region", values_from="decon")

source('/ix1/adombrovski/DNPL_DataMesh/Data/bea_demo/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
trial_df <- get_trial_data(repo_directory = "/ix1/adombrovski/DNPL_DataMesh/Data/bea_demo/clock_analysis/", dataset = "mmclock_fmri", groupfixed = TRUE) %>%
  dplyr::rename(run_number = "run") %>% mutate(feedback_duration = iti_onset - feedback_onset) %>% mutate(isi_duration = trial_df$feedback_onset - trial_df$isi_onset)


temp_folder <- data.frame(mr_dir = list.dirs('/ix1/adombrovski/DNPL_DataMesh/Data/MMClock/MMClock_lite',recursive = FALSE))

id <- NULL
for (iD in 1:nrow(temp_folder)){
  id[iD] <- str_split(str_split(temp_folder$mr_dir[iD],'/')[[1]][8],'_')[[1]][1]
}

temp_folder <- temp_folder %>% mutate(id = id)

subject_df <- temp_folder %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) #make sure we're looking in the right folder



# subject_df0 <- list.dirs("/ix1/adombrovski/papalea/2025-07-07-gPPI/MR_Proc",recursive = FALSE, full.names=TRUE)  
# mr_dir <- NULL
# id <- NULL
# run <- NULL
# for (iD in 1:length(subject_df0)){
#   mr_dir <- rbind(paste0(subject_df0[[iD]][1],'/mni_5mm_aroma'),mr_dir)
#   id <- rbind(id,str_split(str_split(subject_df0[[iD]][1],'/')[[1]][7],'_')[[1]][1])
# }
# subject_df <- data.frame(mr_dir = mr_dir, id=id) %>% filter(id %in% c(10637, 10997, 11279))

# run_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_run_data.rds")

if (create_new_gpa_object==TRUE){
  gpa <- setup_glm_pipeline(
    analysis_name = "ppi_demo2", scheduler = "slurm",
    output_directory = "/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output8",
    subject_data = subject_df, trial_data = trial_df, run_data = NULL, ppi_data = ppi_data,
    tr = 1.0, drop_volumes = 0,
    n_expected_runs = 8,
    l1_models = NULL, l2_models = NULL, l3_models = NULL,
    fmri_file_regex = "nfaswuktm_clock[1-8]_5_drop2_trunc[0-9]+\\.nii\\.gz",
    fmri_path_regex = "clock[1-9]",
    run_number_regex = ".*clock([0-9]+)_5.*",
    confound_settings = list(
      motion_params_file = "motion.par",
      confound_input_file = "nuisance_regressors.txt",
      confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
      l1_confound_regressors = c("csf", "dcsf", "wm", "dwm"),
      exclude_run = "max(framewise_displacement) > 5 | sum(framewise_displacement > .9)/length(framewise_displacement) > .10", # this must evaluate to a scalar per run
      exclude_subject = "n_good_runs < 4",
      # truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
      spike_volumes = NULL
    ),
    parallel = list(
      fsl = list(
        l1_feat_alljobs_time = "48:00:00"
      )
    )
  )
  
  gpa <- setup_compute_environment(gpa) # need to fill in the right module load commandgpa$output_locations
  
  # gpa$run_data <- gpa$run_data %>% mutate(run_nifti = run_nifti, exclude_run = FALSE)
  save(gpa, file="/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output8/MMClock-2025-08-22-ppi.Rdata")
  gpa <- build_l1_models(gpa)
  
  #gpa$parallel$l1_setup_cores <- 1 # fallback to serial for debugging
  
  #gpa <- readRDS('/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output4/ppi_demo2/ppi_MMClock1.rds')
  
  gpa <- build_l2_models(gpa)
  gpa$run_data <- gpa$run_data %>% mutate(exclude_run = FALSE, exclude_subject = FALSE)
  
  gpa <- build_l3_models(gpa)
  
  gpa$subject_data <- gpa$subject_data  %>% mutate(exclude_subject = FALSE)
  
  gpa$run_data <- NULL
  
  id <- NULL
  run_number <- NULL
  for (iD in 1:length(run_nifti)){
    id[iD] <- str_split(str_split(run_nifti[iD],'/')[[1]][8],'_')[[1]][1]
    run0 <- str_split(str_split(run_nifti[iD],'/')[[1]][11],'_')[[1]][2]
    run_number[iD] <- as.numeric(gsub('clock','',run0)) 
  }
  
} else {
  
  load(gpa, file="/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output8/MMClock-2025-08-22-ppi.Rdata")
  
}

id <- NULL
run_number <- NULL
for (iD in 1:length(run_nifti)){
  id[iD] <- str_split(str_split(run_nifti[iD],'/')[[1]][8],'_')[[1]][1]
  run0 <- str_split(str_split(run_nifti[iD],'/')[[1]][11],'_')[[1]][2]
  run_number[iD] <- as.numeric(gsub('clock','',run0)) 
}
gpa$run_data <- data.frame(mr_dir = run_nifti,id = as.numeric(id), run_number = run_number, session = 1)
gpa$subject_data$id <- as.numeric(gpa$subject_data$id)
gpa$run_data$id <- as.numeric(gpa$run_data$id)
gpa$run_data$tr <- 1
#gpa <- run_glm_pipeline(gpa)

load('/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output4/motion_files.Rdata')
load('/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output4/confound_files.Rdata')

gpa$run_data$confound_input_file <- conf_files0
gpa$run_data$motion_params_file <- motion_files0

### MNH checks

#gpa <- readRDS("/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output2/ppi_demo2/ppi_demo2.rds")
# gpa$output_directory <- "/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output8"
gpa$scheduler <- "slurm"
gpa$sqlite_con <- NULL # seems to be holding an invalid database connection
# gpa$output_locations <- lapply(gpa$output_locations, function(x) {
#   if (is.character(x)) sub("gPPI-output4/", "gPPI-output4/", x, fixed=TRUE)
#   else x
# })

# has files, not directories
gpa$run_data$mr_dir <- dirname(gpa$run_data$mr_dir)

#gpa <- finalize_pipeline_configuration(gpa)

save(gpa, file="/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output8/MMClock-2025-08-22-ppi.Rdata")
export_glm_config(gpa, "/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output8/model_config.yaml")

#run_glm_pipeline(gpa)
#pa <- readRDS("/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output8/ppi_MMClock1.rds")


#debug script
setwd('/ix1/adombrovski/papalea/2025-07-07-gPPI/')

# source('get_output_directory.R')
# source('sqlite_functions.R')
# source('get_l1_confounds.R')
# source('glm_helper_functions.R')
# source('lookup_nifti_inputs.R')
# source('finalize_pipeline_configuration.R')

library(tidyverse)
library(glue)
library(stringr)

setwd('/ix1/adombrovski/papalea/2025-07-07-gPPI/gPPI-output4/')

# gpa <- readRDS('MMClock-gPPI-08-06.rds')
# 
# lg <- lgr::get_logger("glm_pipeline/setup_glm_pipeline")
# gpa <- finalize_confound_settings(gpa, lg)

elig_files <- list.files("/ix1/adombrovski/DNPL_DataMesh/Data/MMClock/MMClock_lite/", pattern = "nfaswuktm_clock[1-8]_5_drop2_trunc[0-9]+\\.nii\\.gz$",full.names = TRUE, recursive=TRUE)


conf_files <- list.files("/ix1/adombrovski/DNPL_DataMesh/Data/MMClock/MMClock_lite/", pattern = "nuisance_regressors.txt",full.names = TRUE, recursive=TRUE)


nF = length(conf_files)
conf_files0 <- NULL
iR <- 1;
for (iF in 1:nF){
  tempfold <- str_split(conf_files[iF],'/')[[1]][10]
  if (tempfold == 'mni_5mm_aroma'){
    conf_files0[iR] <- conf_files[iF]
    iR <- iR + 1;
  }
}




motion_files <- list.files("/ix1/adombrovski/DNPL_DataMesh/Data/MMClock/MMClock_lite/", pattern = "motion.par",full.names = TRUE, recursive=TRUE)

nF = length(motion_files)
motion_files0 <- NULL
iR <- 1;
for (iF in 1:nF){
  tempfold <- str_split(motion_files[iF],'/')[[1]][10]
  if (tempfold == 'mni_5mm_aroma'){
    motion_files0[iR] <- motion_files[iF]
    iR <- iR + 1;
  }
}


nF <- length(motion_files0)
nD <- length(elig_files)
ix2rem_conf <- NULL
ix2rem_motion <- NULL
for (iF in 1:nF){
  motion_temp_id <- str_split(motion_files0[iF],'/')[[1]][9]
  motion_temp_run <- str_split(motion_files0[iF],'/')[[1]][11]
  conf_temp_id <- str_split(conf_files0[iF],'/')[[1]][9]
  conf_temp_run <- str_split(conf_files0[iF],'/')[[1]][11]
  in_the_set <- FALSE
  for (iD in 1:nD){
    mri_temp_id <- str_split(elig_files[iD],'/')[[1]][9]
    mri_temp_run <- str_split(elig_files[iD],'/')[[1]][11]
    if (conf_temp_id==mri_temp_id && conf_temp_run==mri_temp_run){
      in_the_set <- TRUE
    }
  }
  if (in_the_set==FALSE){
    ix2rem_conf <- rbind(ix2rem_conf,iF)
  }
  in_the_set <- FALSE
  for (iD in 1:nD){
    mri_temp_id <- str_split(elig_files[iD],'/')[[1]][9]
    mri_temp_run <- str_split(elig_files[iD],'/')[[1]][11]
    if (motion_temp_id==mri_temp_id && motion_temp_run==mri_temp_run){
      in_the_set <- TRUE
    }
  }
  if (in_the_set==FALSE){
    ix2rem_motion <- rbind(ix2rem_motion,iF)
  }  
}

conf_files <- conf_files[-ix2rem_conf]
motion_files <- motion_files[-ix2rem_motion]

