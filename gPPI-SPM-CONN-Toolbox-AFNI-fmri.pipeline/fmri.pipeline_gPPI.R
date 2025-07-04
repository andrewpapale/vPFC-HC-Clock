# 2025-06-03 AndyP and Angela I
# script to run gPPI in fmri.pipeline

# PPI prototype

library(data.table)
library(tidyr)
library(dplyr)
library(fmri.pipeline)
library(readr)

#physio_data <- data.table::fread("/proj/mnhallqlab/projects/mmy3_gppi_prototype/2025-05-29-HC_MEDuSA-gPPI-Prototype.csv")
physio_data <- data.table::fread("/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/2025-06-27-gPPI-demo.csv")
str(physio_data)

ppi_data <- physio_data %>% group_by(id, run, time, HC_region) %>% summarize(decon=mean(decon)) %>% ungroup()

# for simplicity, pare this down and reshape it -- not quire sure what the intent is here since these are not deconvolved time series
# ppi_data <- physio_data %>%
#     select(id, run)

elig_files <- list.files("/Volumes/bierka_root/datamesh/RAW/MMClock/MR_Proc", pattern = "nfaswuktm_clock[1-8]_5_drop2_trunc[0-9]+\\.nii\\.gz$",
                         #elig_files <- list.files("/proj/mnhallqlab/projects/mmy3_gppi_prototype/MMClock/MR_Proc", pattern = "nfaswuktm_clock.*", 
                         full.names = TRUE, recursive=TRUE
)

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
ppi_data$volume <- ppi_data$time + 1
ppi_data$time <- NULL
ppi_data <- ppi_data %>%
  rename(run_number = run) %>%
  mutate(run_number = readr::parse_number(run_number)) %>%
  pivot_wider(names_from="HC_region", values_from="decon")

source('/Users/dnplserv/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R')
trial_df <- get_trial_data(repo_directory = "/Users/dnplserv/clock_analysis/", dataset = "mmclock_fmri", groupfixed = TRUE) %>%
  dplyr::rename(run_number = "run")

subject_df <- readRDS("/Volumes/Users/Andrew/v19-2025-05-27-JNeuro-postR2R/mmclock_subject_data.rds") %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) %>% #make sure we're looking in the right folder
  filter(id %in% c(10637, 10997, 11279))

# run_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_run_data.rds")

gpa <- setup_glm_pipeline(
  analysis_name = "ppi_demo", scheduler = "local",
  output_directory = "/Volumes/Users/Andrew/2025-05-29-MMC-gPPI-Prototype-forMNH/MMClock/gPPI-output",
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

gpa <- build_l1_models(gpa)

gpa <- finalize_pipeline_configuration(gpa)
gpa$parallel$l1_setup_cores <- 1 # fallback to serial for debugging
gpa <- setup_l1_models(gpa)