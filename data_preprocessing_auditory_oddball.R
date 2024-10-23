# Script purpose: READ AND ANALYZE DATA - SEGA PROJECT - AUDITORY ODDBALL
# Author: Nico Bast
# Date Created: `r paste(Sys.Date())`
# Copyright (c) Nico Bast, `r paste(format(Sys.Date(), "%Y"))`
# Email: nico.bast@kgu.de

###notes on

# - 027 without task data
# - 037 task aborted
# - 073 without ET data
# - 056 program failed
# - 070 no data due to subject-related issues
# - 073 problems with starting the task
# - 077 complete data but script fails
# - 139 no data due to subject-related issues
# - 122 no eeg data due to eeg pc problem
# - 114 no eeg data due to eeg pc problem
# - 103 only .psydat, no csv file
# - 164 no eeg due to subject's aversion
# - 161 no egg because wearing wig


## SETUP ####

sessionInfo()

# REQUIRED PACKAGES
library(rhdf5, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE) # efficient due to parallelization
library(zoo, warn.conflicts = FALSE) # used for na.approx
library(pbapply, warn.conflicts = FALSE) # progress bar for apply functions
library(lme4, warn.conflicts = FALSE) # linear-mixed-effects models
library(lmerTest, warn.conflicts = FALSE) # linear-mixed-effects models
library(emmeans, warn.conflicts = FALSE) # estimated marginal means (EMMs)
library(ggplot2, warn.conflicts = FALSE) # creating graphs
library(hexbin, warn.conflicts = FALSE) # binning + plotting functions
library(gridExtra, warn.conflicts = FALSE) # multiple plots arrangement
library(dplyr, warn.conflicts = FALSE) # for %>% operator
library(cowplot, warn.conflicts = FALSE) # get only legend from plot
library(ggpubr, warn.conflicts = FALSE) # save legend from plot
library(stringr)
library(tidyverse)

# PATHS
if (Sys.info()["sysname"] == "Linux") {
  home_path <- "~"
  project_path <- "/PowerFolders/project_sega"
  data_path <- "/PowerFolders/project_sega/data/AuditoryOddball"
  data_path_eeg <- "/PowerFolders/project_sega/data/AuditoryOddball_EEG"
  datapath <- paste0(home_path, data_path) # .csv + .hdf5 input files
  datapath_eeg <- paste0(home_path, data_path_eeg) # .txt input files (eeg)
  # List all .hdf and .csv files
  data_files <- list.files(path = datapath, full.names = TRUE)
}

if (Sys.info()["sysname"] == "Windows") {
  home_path <- "C:/Users/Nico"
  project_path <- "/PowerFolders/project_sega"
  data_path <- "/PowerFolders/project_sega/data/et auditory oddball"
  data_path_task<-"/PowerFolders/project_sega/data/task auditory oddball"
  data_path_eeg <- "/PowerFolders/project_sega/data/tests/eeg preprocessed auditory oddball"
  datapath <- paste0(home_path, data_path) # hdf5 input files
  datapath_task <- paste0(home_path, data_path_task) # hdf5 input files
  datapath_eeg <- paste0(home_path, data_path_eeg) # .txt input files (eeg)
  # List all .hdf and .csv files
  data_files <- c(list.files(path = datapath, full.names = TRUE),
                  list.files(path = datapath_task, full.names = TRUE))
}

if (Sys.info()["sysname"] == "Darwin") {
  home_path <- "~"
  project_path <- "/code"
  data_path <- "/code/input/AuditoryOddball/eyetracking"
  data_path_task <- "/code/input/AuditoryOddball/taskdata"
  data_path_eeg <- "/code/input/AuditoryOddball/eeg"
  datapath <- paste0(home_path, data_path) # .hdf5 input files
  datapath_task <- paste0(home_path, data_path_task) # .csv input files
  datapath_eeg <- paste0(home_path, data_path_eeg) # .txt input files (eeg)
  data_path_single_trial_eeg <- "/code/input/AuditoryOddball/eeg_single_trial"
  datapath_single_eeg <- paste0(home_path, data_path_single_trial_eeg)
  # List all .hdf and .csv files
  data_files <- c(
    list.files(path = datapath, full.names = TRUE),
    list.files(path = datapath_task, full.names = TRUE))
}

# DATA IMPORT AND RESHAPING ####

# Get eye tracking data and store them in a list of df (one per subject)
data_files_et <- data_files[grepl(".hdf5", data_files)]

#exclude files that cannot be read
exclude_unreadable<-c("C:/Users/Nico/PowerFolders/project_sega/data/et auditory oddball/auditory_103_2024-04-28-1236.hdf5",
                      "C:/Users/Nico/PowerFolders/project_sega/data/et auditory oddball/auditory_2022-10-04-0911.hdf5",
                      "C:/Users/Nico/PowerFolders/project_sega/data/et auditory oddball/auditory_2022-10-05-1440.hdf5",
                      "C:/Users/Nico/PowerFolders/project_sega/data/et auditory oddball/auditory_2022-10-19-1220.hdf5")


data_files_et<-data_files_et[!data_files_et %in% exclude_unreadable]

list_et_data <- list(0)
for (i in 1:length(data_files_et)) {
  print(paste0("now reading: ", data_files_et[i]))
  list_et_data[[i]] <- h5read(
    file = data_files_et[i],
    name = "data_collection/events/eyetracker/BinocularEyeSampleEvent")
  #print(paste0("completed reading: ", data_files_et[i]))
 }
h5closeAll()

# List names for each subject are unique including date and time of recording.
#NICO: removed magic number here

id_names <- substr(
  data_files_et,
  nchar(datapath) + 2,
  nchar(data_files_et))

names(list_et_data) <- id_names

# These eye tracker variables are being dropped, keeping variables
# left_pupil_measure1, right_pupil_measure1, logged_time,
# left_gaze_x and left_gaze_y.
constant_variables <- c(
  "experiment_id",
  "status",
  "session_id",
  "device_id",
  "type",
  "device_time",
  "time",
  "delay",
  "confidence_interval",
  "filter_id",
  "left_gaze_z", "right_gaze_z",
  "left_angle_x", "right_angle_x",
  "left_angle_y", "right_angle_y",
  "left_raw_x", "right_raw_x",
  "left_raw_y", "right_raw_y",
  "left_pupil_measure1_type", "right_pupil_measure1_type",
  "left_pupil_measure2_type", "right_pupil_measure2_type",
  "left_ppd_x", "right_ppd_x",
  "left_ppd_y", "right_ppd_y",
  "left_velocity_x", "right_velocity_x",
  "left_velocity_y", "right_velocity_y",
  "left_velocity_xy", "right_velocity_xy",
  "left_pupil_measure2", "right_pupil_measure2",
  "left_eye_cam_x", "right_eye_cam_x",
  "left_eye_cam_y", "right_eye_cam_y",
  "left_eye_cam_z", "right_eye_cam_z"
  )

list_et_data <- lapply(
  list_et_data, function(x) {
    x[!(names(x) %in% constant_variables)]})


# Get trial data and store them in a list of df (one per subject)
data_files_trial <- data_files[grepl(".csv", data_files)]
list_trial_data <- list(0)
trial_variables <- c(
  ".thisRepN",
  ".thisTrialN",
  "phase",
  "block_counter",
  "stimulus_duration",
  "baseline_trial_counter",
  "grip_strength",
  "trial",
  "timestamp_exp",
  "oddball_trial_counter",
  "ISI_duration",
  "manipulation_trial_counter",
  "id",
  "group",
  "oddball_frequency",
  "standard_frequency"
  )

for (i in 1:length(data_files_trial)) {
  list_trial_data[[i]] <- fread(data_files_trial[i], select = trial_variables)
  print(paste0("read TRIAL data file: ", i))
  }
list_trial_data <- lapply(list_trial_data, data.frame)

# List names for each subject are unique including date and time of recording.
id_names <- substr(
  data_files_trial,
  nchar(datapath) + 4,
  nchar(data_files_trial))
names(list_trial_data) <- id_names

#check: which ET data and task data per participant do not match
#task_data<-substr(data_files_trial,nchar(datapath_task)+2,nchar(data_files_trial)-4)
#et_data<-substr(data_files_et,nchar(datapath)+2,nchar(data_files_et)-5)
#unmatched_et<-et_data[!(et_data %in% task_data)]
#unmatched_task<-task_data[!(task_data %in% et_data)]

# Check for .csv- + .hdf5- file matching (path-independent):
unmatched_et <- tools::file_path_sans_ext(basename(data_files_et))[
  !(tools::file_path_sans_ext(basename(data_files_et)) %in% tools::file_path_sans_ext(basename(data_files_trial)))]
unmatched_task <- tools::file_path_sans_ext(basename(data_files_trial))[
  !(tools::file_path_sans_ext(basename(data_files_trial)) %in% tools::file_path_sans_ext(basename(data_files_et)))]

# print unmatching files in console
cat(unmatched_task, "do not have a matching et file", sep = "\n") 
cat(unmatched_et, "do not have a matching trial file", sep = "\n")

#reduce data files for participant with ET and task data
list_et_data<-list_et_data[!(names(list_et_data) %in% paste0(unmatched_et,'.hdf5'))]
list_trial_data<-list_trial_data[!(names(list_trial_data) %in% paste0(unmatched_task,'.csv'))]
  
  #TODO: function below fails with this ID
  list_et_data<-list_et_data[names(list_et_data)!="auditory_77_2023-04-17-1831.hdf5"]
  list_trial_data<-list_trial_data[names(list_trial_data)!="auditory_77_2023-04-17-1831.csv"]

# Grip strength variable
list_trial_data <- lapply(list_trial_data, function(x) {
  mean_grip_strength <- mean(x$grip_strength, na.rm = TRUE)
  x$mean_grip_strength <- mean_grip_strength
  return(x)
  })
  
# subject data frames are row-wise combined
df_trial <- plyr::rbind.fill(list_trial_data)

# Eye tracking data (logged_time) are assigned to trials (timestamp_exp).
# Before it needs to be checked that trial data matches to ET data
fun_merge_all_ids <- function(et_data, trial_data) {
# Time variables: eye tracking (logged_time) + trial data (timestamp_exp)
  start_ts <- trial_data$timestamp_exp # trial start
  end_ts <- c(trial_data$timestamp_exp[-1], NA) # trial end
  et_ts <- et_data$logged_time
  split_trial_data <- split(trial_data, seq(nrow(trial_data)))
  
  fun_merge_data <- function(ts_1, ts_2, trial_data_splitted) {
    matched_time <- which(et_ts >= ts_1 & et_ts < ts_2)
    if (trial_data_splitted$baseline_trial_counter == "6" & !is.na(trial_data_splitted$baseline_trial_counter)) {
      matched_time <- which(et_ts >= ts_1)
    } 
    selected_et_data <- et_data[matched_time, ] # et data for trial duration
    # trial data: 1 row == 1 trial -> is repeated for each eye tracking event
    repeated_trial_data <- data.frame(
      sapply(trial_data_splitted, function(x) {
        rep(x, length(matched_time))}, simplify = FALSE))
        merged_data <- data.frame(repeated_trial_data, selected_et_data)
  }
  
  print(paste0("merge: ", unique(trial_data$id))) #debugging print
  
  df_one_id <- mapply(
    fun_merge_data,
    ts_1 = start_ts,
    ts_2 = end_ts,
    trial_data_splitted = split_trial_data,
    SIMPLIFY = FALSE)
    df_one_id <- dplyr::bind_rows(df_one_id) # faster than rbind.fill
}

# Calling function
df_list <- pbmapply(
  fun_merge_all_ids,
  et_data = list_et_data,
  trial_data = list_trial_data, SIMPLIFY = FALSE)

#drop empty list elements (unmatched)
df_list<-df_list[sapply(df_list,function(x){length(x)!=0})]

#create trial_index and trial_number variables
df_list <- pblapply(df_list, function(x) {
  # New variable trial_index
  x$trial_index <- with(x, droplevels(interaction(
    .thisTrialN, .thisRepN, block_counter)))
    # New variable: trial_number
    x$trial_number <- with(x, rep(seq_along(table(trial_index)),
    times = table(trial_index)))
    return(x)})

# Delete 3 preceding standard trials at the beginning of each of the 4 blocks
standards_index <- c("-1.0.3", "-1.0.5", "-1.0.10", "-1.0.12")
df_list <- pblapply(df_list, function(x) {
  x <- subset(x, ! trial_index %in% standards_index)
  return(x)
})

# Split trial data per trial -> list of all trials of all subjects.
# Each list element contains a df of trial data
list_split_trial <- pblapply(df_list, function(x) {
  split(x, x$trial_number)})
list_split_trial <- unlist(list_split_trial, recursive = FALSE)

# New variable ts_trial: Timestamp for each et event within a trial
list_split_trial <- pblapply(list_split_trial, function(x) {
  x$ts_trial <- x$logged_time - x$timestamp_exp
  return(x)
  })

# DATA PREPROCESSING ####
# Blinks are defined as consecutive missing et data for 75–250 ms
fun_blink_cor <- function(
  signal, lower_threshold = 23, upper_threshold = 75,
  samples_before = 8, samples_after = 8) {
  # Replace Na with 999
  findna <- ifelse(is.na(signal), 999, signal)
  repets <- rle(findna) # gives number of repetitions
  # Repeat number of repetition as often the value is
  repets <- rep(repets[["lengths"]], times = repets[["lengths"]])
  # 75 ms / 3,33 sampling interval = 23 samples (rows of et data)
  # 250 ms / 3,33 ms sampling interval = 75 samples (rows of eye tracking data)
  # If value is consecutively repeted >= 23 and <= 75, coding is "1", else "0"
  repets <- ifelse(repets >= lower_threshold & repets <= upper_threshold, 1, 0)
  # Repeated values other than Na are set to "0"
  repets[findna != 999 & repets == 1] <- 0
  # Differences between consecutive values indicate blink artefact bounderies
  changes <- c(diff(repets), 0)
  change_start <- which(changes == 1)
  # Blink sequence includes 8 samples before after blink, repectively.
  start_seq <- unlist(lapply(change_start, function(x) {
    seq(max(x - (samples_before - 1), 1), x)
    }
    ))
  repets[start_seq] <- 1
  changes_end <- which(changes == -1) + 1
  end_seq <- unlist(lapply(changes_end, function(x) {
    seq(x, min(x + (samples_before - 1), length(repets)))
    }
    ))
  repets[end_seq] <- 1
  # Data in blink interval is replaced with Na.
  signal[repets == 1] <- NA
  return(signal)
}

func_pd_preprocess <- function(x) {
  left_diameter <- x$left_pupil_measure1
  right_diameter <- x$right_pupil_measure1
  remote_time <- x$ts_trial * 1000 # *1000 to convert s -> ms format
  # Pupil diameter outliers (< 2 mm or > 8 mm) are replaced with Na.
  pl <- ifelse((left_diameter < 2 | left_diameter > 8), NA, left_diameter)
  pr <- ifelse((right_diameter < 2 | right_diameter > 8), NA, right_diameter)
  # Dilation speed outliers: > constant * median change values are excluded
  constant <- 3
  # speed defined as movement / time
  # Dilatation speed for left eye
  pl_speed1 <- diff(pl) / diff(remote_time) # compared to previous et event
  pl_speed2 <- diff(rev(pl)) / diff(rev(remote_time)) # compared to next event
  pl_speed1 <- c(NA, pl_speed1)
  pl_speed2 <- c(rev(pl_speed2), NA)
  pl_speed <- pmax(pl_speed1, pl_speed2, na.rm = TRUE)
  rm(pl_speed1, pl_speed2)
  # Dilatation speed for right eye
  pr_speed1 <- diff(pr) / diff(remote_time) # compared to previous et event
  pr_speed2 <- diff(rev(pr)) / diff(rev(remote_time)) # compared to next event
  pr_speed1 <- c(NA, pr_speed1)
  pr_speed2 <- c(rev(pr_speed2), NA)
  pr_speed <- pmax(pr_speed1, pr_speed2, na.rm = TRUE)
  rm(pr_speed1, pr_speed2)
  # Threshold (in mm/ms): dilation speed median + 3 * median absolute deviation
  # Left eye
  pl_speed_med <- median(pl_speed, na.rm = TRUE)
  pl_mad <- median(abs(pl_speed - pl_speed_med), na.rm = TRUE)
  pl_treshold_speed <- pl_speed_med + constant * pl_mad
  # Right eye
  pr_speed_med <- median(pr_speed, na.rm = TRUE)
  pr_mad <- median(abs(pr_speed - pr_speed_med), na.rm = TRUE)
  pr_treshold_speed <- pr_speed_med + constant * pr_mad
  # Replace pupil data higher than threshold with Na
  pl <- ifelse(abs(pl_speed) > pl_treshold_speed, NA, pl)
  pr <- ifelse(abs(pr_speed) > pr_treshold_speed, NA, pr)
  # Calling function for blink correction
  pl <- fun_blink_cor(pl)
  pr <- fun_blink_cor(pr)
  # Two pass approach. 1st pass: Exclude deviation from trend
  # line derived from all samples. 2nd pass: Exclude deviation from trend
  # line derived from samples passing. Reintroduction of sample that might
  # have been falsely excluded due to outliers estimate smooth size based
  # on sampling rate
  smooth_length <- 150 # in ms
  # take sampling rate into account (300 vs. 120):
  smooth_size <- round(
    smooth_length / median(diff(remote_time), # remote_time is ts_trial in ms
    na.rm = TRUE))
  is_even <- function(x) {
    x %% 2 == 0
    }
  smooth_size <- ifelse(
    is_even(smooth_size) == TRUE,
    smooth_size + 1, smooth_size) # odd values for runmed()-function
  # for left and right eye:
  # giving the smooth function Na would raise an error
  pl_smooth <- na.approx(pl, na.rm = FALSE, rule = 2)
  # Robust Scatter Plot Smoothing
  if (sum(!is.na(pl_smooth)) != 0) {
    pl_smooth <- runmed(pl_smooth, k = smooth_size)
    }
  pl_mad <- median(abs(pl - pl_smooth), na.rm = TRUE)
  # Giving the smooth function Na would raise an error
  pr_smooth <- na.approx(pr, na.rm = FALSE, rule = 2)
  # Robust Scatter Plot Smoothing
  if (sum(!is.na(pr_smooth)) != 0) {
    pr_smooth <- runmed(pr_smooth, k = smooth_size)
    }
  pr_mad <- median(abs(pr - pr_smooth), na.rm = TRUE)
  # correct pupil dilation for size outliers - 1st pass
  pl_pass1 <- ifelse(
    (pl > pl_smooth + constant * pl_mad) | (pl < pl_smooth - constant * pl_mad),
    NA, pl)
  pr_pass1 <- ifelse(
    (pr > pr_smooth + constant * pr_mad) | (pr < pr_smooth - constant * pr_mad),
    NA, pr)
  # for left and right eye:
  # giving the smooth function Na would raise an error
  pl_smooth <- na.approx(pl_pass1, na.rm = FALSE, rule = 2)
  # Robust Scatter Plot Smoothing
  if (sum(!is.na(pl_smooth)) != 0) {
    pl_smooth <- runmed(pl_smooth, k = smooth_size)
    }
  pl_mad <- median(abs(pl - pl_smooth), na.rm = TRUE)
  # Giving the smooth function Na would raise an error
  pr_smooth <- na.approx(pr_pass1, na.rm = FALSE, rule = 2)
  # Robust Scatter Plot Smoothing
  if (sum(!is.na(pr_smooth)) != 0) {
    pr_smooth <- runmed(pr_smooth, k = smooth_size)
    }
  pr_mad <- median(abs(pr - pr_smooth), na.rm = TRUE)
  # correct pupil dilation for size outliers - 2nd pass
  pl_pass2 <- ifelse(
    (pl > pl_smooth + constant * pl_mad) | (pl < pl_smooth - constant * pl_mad),
    NA, pl)
  pr_pass2 <- ifelse(
    (pr > pr_smooth + constant * pr_mad) | (pr < pr_smooth - constant * pr_mad),
    NA, pr)
  pl <- pl_pass2
  pr <- pr_pass2
  # Fill Na with offset value
  pd_offset <- pl - pr
  pd_offset <- na.approx(pd_offset, rule = 2)
  pl <- ifelse(is.na(pl) == FALSE, pl, pr + pd_offset)
  pr <- ifelse(is.na(pr) == FALSE, pr, pl - pd_offset)
  # Interpolation of missing values < 300 ms
  pl <- na.approx(pl, na.rm = FALSE, maxgap = 90, rule = 2)
  pr <- na.approx(pr, na.rm = FALSE, maxgap = 90, rule = 2)
  # mean pupil dilation across both eyes
  pd <- (pl + pr) / 2
  x[, "pd"] <- pd
  return(x)
}

list_split_trial <- pblapply(
  list_split_trial, func_pd_preprocess)

# Bind trials together to a df
df <- dplyr::bind_rows(list_split_trial)

# Create df with 6 baseline means
baseline_means <- data.frame()
subjects <- unique(df$id)
baseline_means_subj <- NULL

for (subject in subjects){
  subject_df <- df[df$id == subject, ]
  baseline_trial_counter <- c(1, 2, 3, 4, 5, 6)
  for (i in baseline_trial_counter){
    baseline_df <- subject_df[subject_df$trial == "baseline" & subject_df$baseline_trial_counter == i, ]
    # block_baseline_mean is the new variable
    block_baseline_mean <- mean(baseline_df$pd, na.rm = TRUE)
    # id for merging with df_trial
    id <- subject
    # baseline_trial_counter for plot
    baseline_trial_counter <- i
    # block_counter of following oddball-Block for merging
    if (i == 1){
      block_counter <- 3}
    if (i == 2){
      block_counter <- 5}    
    if (i == 3){
      block_counter <- NA}
    if (i == 4){
      block_counter <- 10}
    if (i == 5){
      block_counter <- 12}
    if (i == 6){
      block_counter <- NA}
    baseline_means_subj <- cbind(id, block_counter, block_baseline_mean, baseline_trial_counter)
    baseline_means <- rbind(baseline_means, baseline_means_subj)
  }}

# baseline_trial_counter as factor for boxplot
baseline_means$baseline_trial_counter <- as.factor(baseline_means$baseline_trial_counter)

# Pupil size for baselines
ggplot(
  baseline_means,
  aes(x = baseline_trial_counter,
      y = block_baseline_mean)) + 
  geom_boxplot(fill = "steelblue")

# split by block and id
list_split_blocks <- split(df,
droplevels(interaction(df$block_counter, df$id)))

list_split_blocks <- pblapply(list_split_blocks, function(x) {
  # name to identify individual trials within BLOCK
  x$trial_index_in_block <- with(x,
  droplevels(interaction(.thisTrialN, .thisRepN)))
  # all trials including empty trials
  x$trial_number_in_block <- with(x,
  rep(seq_along(table(trial_index_in_block)),
  times = table(trial_index_in_block)))
  return(x)})

df <- data.table::rbindlist(list_split_blocks)

# melt to data.frame
df$id <- as.character(df$id)

# split by trial (for further processing)
list_split_trial <- split(df, droplevels(interaction(df$id, df$trial_number)))

  # Number of et event per phase
  table(df$phase)
  # Number of et events per trial
  table(df$trial)
  # Frequency of pupil diameter (mm) in 5 size categories
  hist(df$pd)
  # preprocessed pupil diameter for each participant
  # ggplot(
  #   df[is.finite(df$pd), ],
  #   aes(x = pd)) + geom_histogram(bins = 100) + facet_wrap(~id)
  
  
# reduce ET data to per trial data (and merge with df_trial)
# pupil diameter late in trial duration
# rpd_high <- sapply(list_split_trial, function(x) {
#   mean(x$pd[x$ts_trial > 0.875 & x$ts_trial < 1.125])})

require(DescTools) #AUC function  
rpd_auc <- sapply(list_split_trial, function(x) {
  AUC(x$ts_trial,x$pd,na.rm=T)})
rpd_auc<-ifelse(rpd_auc>15,NA,rpd_auc)

rpd_high <- sapply(list_split_trial, function(x) {
  mean(x$pd[x$ts_trial > 0.500 & x$ts_trial < 1.5])})

rpd_low <- sapply(list_split_trial, function(x) {
  mean(x$pd[x$ts_trial < 0.250])})

# news variables to index specific trials
trial_number <- sapply(list_split_trial, function(x) {
  unique(x$trial_number)})
trial_number_in_block <- sapply(list_split_trial, function(x) {
  unique(x$trial_number_in_block)[1]})

# merge data
merger_id <- sapply(list_split_trial, function(x) {
  interaction(x$id, x$trial_index)[1]})

df_et_trial <- data.frame(
  merger_id,
  rpd_auc,
  rpd_high,
  rpd_low,
  trial_number,
  trial_number_in_block)

##save trialwise structure
df_trial_backup<-df_trial

# name to identify individual trials
df_trial$trial_index <- with(df_trial,
interaction(.thisTrialN, .thisRepN, block_counter))
df_trial$merger_id <- interaction(df_trial$id, df_trial$trial_index)
df_trial <- merge(df_trial, df_et_trial, id = "merger_id")

# rpd = scaled trial-baseline corrected pupil diameter
df_trial$rpd <- df_trial$rpd_high - df_trial$rpd_low

# Add baseline means to df_trial
# Before, get rid of column "baseline_trial_counter" from df baseline_means. Was necessary for the plot but
# now prevent to have this column twice due to merging.
baseline_means <- subset(baseline_means, select = -c(baseline_trial_counter))
df_trial <- merge(x = baseline_means, y = df_trial, by = c("id", "block_counter"), all.y = TRUE)

# new variable rpd_block is block baseline corrected rpd_high
df_trial$rpd_block <- df_trial$rpd_high - df_trial$block_baseline_mean

# new variables manipulation indicates whether before or after manipulation
df_trial$manipulation <- factor(
  ifelse(df_trial$block_counter < 8, "before",
  ifelse(df_trial$block_counter > 8, "after",
  ifelse(df_trial$block_counter == 8, "during", "manipulation"))),
  levels = c("before", "after", "during"))

# distinguish pitch: oddball or standard frequency
df_trial$pitch <- ifelse(df_trial$trial == "oddball",
df_trial$oddball_frequency, df_trial$standard_frequency)
# distinguish between (forward) oddball blocks and reverse oddball blocks
df_trial$block <- ifelse(
  df_trial$phase == "oddball_block_rev", "reverse", "forward")
# distinguish between oddball and standards trials
df_trial$trial <- as.factor(ifelse(grepl(
  "oddball",
  df_trial$trial),
  "oddball", "standard"))


### VISUALIZATION ####

##manipulation check
sampled_rows<-sample(1:nrow(df),nrow(df)/1000) #randomly select rows
ggplot(df[sampled_rows & df$phase=='oddball_block' & df$ts_trial<2,],aes(x=ts_trial,y=pd,group=trial,color=trial))+geom_smooth()+facet_wrap(~block_counter)+theme_bw()  
ggplot(df[sampled_rows & df$phase=='oddball_block_rev' & df$ts_trial<2,],aes(x=ts_trial,y=pd,group=trial,color=trial))+geom_smooth()+facet_wrap(~block_counter)+theme_bw()  
###--> select rpd high based on data inspection


###--> show auditory oddball response in first block before manipulation
df$rpd<-with(df,pd-rpd_low)
  ###--> also shows a pupillary response in reverse trials

ggplot(df[sampled_rows & df$ts_trial<1 & df$phase=='oddball_block',],aes(x=ts_trial,y=rpd,group=trial,color=trial))+geom_smooth()+
  facet_wrap(~block_counter)+
  theme_bw()  

ggplot(df[sampled_rows & df$ts_trial<1 & df$phase=='oddball_block_rev',],aes(x=ts_trial,y=rpd,group=trial,color=trial))+geom_smooth()+
  facet_wrap(~block_counter)+
  theme_bw()  


require(wesanderson) #custom color palettes
custom_condition_colors <- wes_palette('Darjeeling1',2,type='discrete') #reverse custom colors to match color coding in other figures

ggplot(df[sampled_rows & df$phase=='oddball_block' & df$ts_trial<0.7 & df$block_counter<8,],
       aes(x=ts_trial,y=rpd,group=trial,color=trial))+geom_smooth()+
       labs(x='trial duration (s)',y='pupillary response (mm)')+
       scale_color_manual(values = custom_condition_colors)+
       theme_bw()  

# Association between pupil variables
hist(df_trial$rpd_block, xlim = c(-2, 2))
hist(df_trial$rpd, xlim = c(-1, 1))

### add experimental group data ####

###--> where is this file?

# Read experimental group data from .csv file
exp_groups <- read.csv("exp_groups.csv", header = TRUE)
exp_groups$SEGA_ID <- as.character(exp_groups$SEGA_ID)
# Include experimental group (ASD, CON, MHC) from .csv file in df_trial
for (row in 1:length(df_trial$id)) {
  SEGA_ID_df <- df_trial[row, "id"]
  row_number <- which(exp_groups$SEGA_ID == SEGA_ID_df)
  group <- exp_groups[row_number, "group"]
  df_trial[row, "group"] <- group}

# Read handdynamometer z-values from .csv files
for (row in 1:length(df_trial$id)) {
  SEGA_ID_df <- df_trial[row, "id"]
  row_number <- which(exp_groups$SEGA_ID == SEGA_ID_df)
  z_handdynamometer <- exp_groups[row_number, "z_handdynamometer"]
  df_trial[row, "z_handdynamometer"] <- z_handdynamometer}

##--> save preprocess df_trial ####
#saveRDS(df_trial,file=paste0(home_path,project_path,'/data/preprocessed_auditory_ETdata.rds'))
saveRDS(df_trial,file=paste0(home_path,project_path,'/data/preprocessed_auditory_ETdata_24092024.rds'))

# # Can be used to skip preprocessing and directly read proprocessed data from .rds file:
# df_trial <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ETdata.rds'))
# 
# #changed random intercept to a factor
# df_trial$id<-as.factor(df_trial$id) #change ID to factor
# #df_trial$trial<-as.factor(df_trial$trial)
# 
# require(performance)

### Data plausibility check ####

table(df_trial$block_counter, df_trial$phase)
table(df_trial$trial)
# Number of trials for forward and reverse oddball blocks
table(df_trial$block, df_trial$trial)
hist(df_trial$rpd, 50)
with(df_trial, by(rpd, trial, mean, na.rm = TRUE))

# Read EEG single trial data
# Single trial EEG data is stored in 2 files (before + after manipulation) per subject
# List all files separately for before + after
files_single_eeg_before <- list.files(
  path = datapath_single_eeg, full.names = TRUE, pattern = "before")
files_single_eeg_after <- list.files(
  path = datapath_single_eeg, full.names = TRUE, pattern = "after")

# For removing 3 preceeding standards
standards <- c(1, 2, 3, 104, 105, 106, 207, 208, 209, 310, 311, 312)

# Read files before manipulation
MMN_df_trial_before <- data.frame()
P3a_df_trial_before <- data.frame()

for (file in files_single_eeg_before) {
  if (grepl("MMN", file)) {
    eeg_data_before <- read.table(file, header = TRUE, fill = TRUE)
    colnames(eeg_data_before)[colnames(eeg_data_before) == "Filename"] <- "SEGA_ID"
    eeg_data_before$SEGA_ID <- as.factor(substr(eeg_data_before$SEGA_ID, 6, 8))
    eeg_data_before$SEGA_ID <- sub("^0+", "", eeg_data_before$SEGA_ID)
    colnames(eeg_data_before)[colnames(eeg_data_before) == "MinMMN.L"] <- "MMN_latency"
    colnames(eeg_data_before)[colnames(eeg_data_before) == "MinMMN.V"] <- "MMN_amplitude"
    colnames(eeg_data_before)[colnames(eeg_data_before) == "Segment"] <- "oddball_trial_counter"
    eeg_data_before <- subset(eeg_data_before, !(oddball_trial_counter %in% standards))
    MMN_df_trial_before <- rbind(MMN_df_trial_before, eeg_data_before)
  }
  if (grepl("P3a", file)) {
    eeg_data_before <- read.table(file, header = TRUE, fill = TRUE)
    colnames(eeg_data_before)[colnames(eeg_data_before) == "Filename"] <- "SEGA_ID"
    eeg_data_before$SEGA_ID <- as.factor(substr(eeg_data_before$SEGA_ID, 6, 8))
    eeg_data_before$SEGA_ID <- sub("^0+", "", eeg_data_before$SEGA_ID)
    colnames(eeg_data_before)[colnames(eeg_data_before) == "MaxP3a.L"] <- "P3a_latency"
    colnames(eeg_data_before)[colnames(eeg_data_before) == "MaxP3a.V"] <- "P3a_amplitude"
    colnames(eeg_data_before)[colnames(eeg_data_before) == "Segment"] <- "oddball_trial_counter"
    eeg_data_before <- subset(eeg_data_before, !(oddball_trial_counter %in% standards))
    P3a_df_trial_before <- rbind(P3a_df_trial_before, eeg_data_before)}
}
# P3a + MMN before manipulation in one df
single_trial_eeg_before <- merge(MMN_df_trial_before, P3a_df_trial_before, by = c("SEGA_ID", "oddball_trial_counter"))

# Read files after manipulation
MMN_df_trial_after <- data.frame()
P3a_df_trial_after <- data.frame()

for (file in files_single_eeg_after) {
  if (grepl("MMN", file)){
    eeg_data_after <- read.table(file, header = TRUE, fill = TRUE)
    colnames(eeg_data_after)[colnames(eeg_data_after) == "Filename"] <- "SEGA_ID"
    eeg_data_after$SEGA_ID <- as.factor(substr(eeg_data_after$SEGA_ID, 6, 8))
    eeg_data_after$SEGA_ID <- sub("^0+", "", eeg_data_after$SEGA_ID)
    colnames(eeg_data_after)[colnames(eeg_data_after) == "MinMMN.L"] <- "MMN_latency"
    colnames(eeg_data_after)[colnames(eeg_data_after) == "MinMMN.V"] <- "MMN_amplitude"
    colnames(eeg_data_after)[colnames(eeg_data_after) == "Segment"] <- "oddball_trial_counter"
    eeg_data_after$oddball_trial_counter <- seq(from = 207, length.out = 206) # trials are continuously numbered throughout the experiment
    eeg_data_after <- subset(eeg_data_after, !(oddball_trial_counter %in% standards)) # after renumbering!
    MMN_df_trial_after <- rbind(MMN_df_trial_after, eeg_data_after)
  }
  if (grepl("P3a", file)){
    eeg_data_after <- read.table(file, header = TRUE, fill = TRUE)
    colnames(eeg_data_after)[colnames(eeg_data_after) == "Filename"] <- "SEGA_ID"
    eeg_data_after$SEGA_ID <- as.factor(substr(eeg_data_after$SEGA_ID, 6, 8))
    eeg_data_after$SEGA_ID <- sub("^0+", "", eeg_data_after$SEGA_ID)
    colnames(eeg_data_after)[colnames(eeg_data_after) == "MaxP3a.L"] <- "P3a_latency"
    colnames(eeg_data_after)[colnames(eeg_data_after) == "MaxP3a.V"] <- "P3a_amplitude"
    colnames(eeg_data_after)[colnames(eeg_data_after) == "Segment"] <- "oddball_trial_counter"
    eeg_data_after$oddball_trial_counter <- seq(from = 207, length.out = 206) # trials are continuously numbered throughout the experiment
    eeg_data_after <- subset(eeg_data_after, !(oddball_trial_counter %in% standards)) # after renumbering!
    P3a_df_trial_after <- rbind(P3a_df_trial_after, eeg_data_after)}
}
# P3a + MMN after manipulation in one df
single_trial_eeg_after <- merge(MMN_df_trial_after, P3a_df_trial_after, by = c("SEGA_ID", "oddball_trial_counter"))

# Bind ERP data before + after manipulation together
# EEG_df_trial contains EEG data on trial level for oddball_blocks + reverse_blocks
EEG_df_trial <- rbind(single_trial_eeg_before, single_trial_eeg_after)

# Merge ET and ERP data on trial level
# ET_df_trial = subset of df_trial with ET data only from oddball_block(rev)
ET_df_trial <- df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev"), ]
colnames(ET_df_trial)[colnames(ET_df_trial) == "id"] <- "SEGA_ID" #colname SEGA_ID for further matching
# left join because experimental infos (pitch, phase, ect.) come with ET data.
ET_ERP_trial <- merge(ET_df_trial, EEG_df_trial, by = c("SEGA_ID", "oddball_trial_counter"), all.x = T)

# Correct data types in ET_ERP_trial data frame
ET_ERP_trial$phase <- as.factor(ET_ERP_trial$phase)
ET_ERP_trial$pitch <- as.factor(ET_ERP_trial$pitch)
ET_ERP_trial$block <- as.factor(ET_ERP_trial$block)
ET_ERP_trial$MMN_latency <- as.numeric(ET_ERP_trial$MMN_latency)
ET_ERP_trial$MMN_amplitude <- as.numeric(ET_ERP_trial$MMN_amplitude)
ET_ERP_trial$P3a_amplitude <- as.numeric(ET_ERP_trial$P3a_amplitude)
ET_ERP_trial$P3a_latency <- as.numeric(ET_ERP_trial$P3a_latency)
ET_ERP_trial$group <- as.factor(ET_ERP_trial$group)
ET_ERP_trial$SEGA_ID <- as.numeric(ET_ERP_trial$SEGA_ID)
str(ET_ERP_trial)

# Checkliste.csv contains age + gender
checkliste_import <- fread("checkliste.csv", select = c("ID_Bado", "ID_Studie", "Geburt_Index", "Geschlecht_Index"))
colnames(checkliste_import)[colnames(checkliste_import) == "ID_Studie"] <- "SEGA_ID" # same column name as in exp_groups for matching
checkliste_import$SEGA_ID <- sub(".*SEGA_", "", checkliste_import$SEGA_ID)
checkliste_import$SEGA_ID <- str_replace(checkliste_import$SEGA_ID, "^0+", "")

# SCQ
SCQ_import <- read.csv("FSK.csv", header = T, sep = ";", dec = ",", fill = T)
SCQ <- SCQ_import[c("ID_Bado", "FSK_ges")]
# SRS
SRS_import <- read.csv("SRS.csv", header = T, sep = ";", dec = ",", fill = T)
SRS <- SRS_import[c("ID_Bado", "Gesamtwert_N_k_RW")]
colnames(SRS)[colnames(SRS) == "Gesamtwert_N_k_RW"] <- "SRS_Gesamtwert_N_k_RW" # add questionnaire name to column names
# CBCL
CBCL_import <- read.csv("CBCL_4_18.csv", header = T, sep = ";", dec = ",", fill = T)
CBCL <- CBCL_import[c("ID_Bado", "CBCL_T_INT", "CBCL_T_EXT", "CBCL_T_GES")]
# YSR
YSR_import <- read.csv("YSR.csv", header = T, sep = ";", dec = ",", fill = T)
YSR <- YSR_import[c("ID_Bado", "YSR_T_GES")]
# SP2
SP2_import <- read.csv("SP2_Eltern.csv", header = T, sep = ";", dec = ",", fill = T)
SP2 <- SP2_import[c("ID_Bado", "SP2_Eltern_RW_Auditiv")]

# Merge sample characteristics 
questionnaires <- list(checkliste_import, SCQ, SRS, CBCL, YSR, SP2) %>%
  reduce(full_join, by = "ID_Bado")
sample_characteristics <- merge(questionnaires, exp_groups, by = "SEGA_ID")
# Calculated age in separate column
sample_characteristics$age <- as.numeric(
  difftime(sample_characteristics$date_data_collection, sample_characteristics$Geburt_Index,
           units = "weeks"))/52.25

# Add sample characteristics to ET_ERP_trial
for (row in 1:length(ET_ERP_trial$SEGA_ID)) {
  SEGA_ID_df <- ET_ERP_trial[row, "SEGA_ID"]
  row_number <- which(sample_characteristics$SEGA_ID == SEGA_ID_df)
  gender <- sample_characteristics[row_number, "Geschlecht_Index"]
  age <- sample_characteristics[row_number, "age"]
  SCQ <- sample_characteristics[row_number, "FSK_ges"]
  SRS <- sample_characteristics[row_number, "SRS_Gesamtwert_N_k_RW"]
  CBCL <- sample_characteristics[row_number, "CBCL_T_GES"]
  YSR <- sample_characteristics[row_number, "YSR_T_GES"]
  SP2 <- sample_characteristics[row_number, "SP2_Eltern_RW_Auditiv"]
  ET_ERP_trial[row, "gender"] <- gender
  ET_ERP_trial[row, "age"] <- age
  ET_ERP_trial[row, "SCQ"] <- SCQ
  ET_ERP_trial[row, "SRS"] <- SRS
  ET_ERP_trial[row, "CBCL"] <- CBCL
  ET_ERP_trial[row, "YSR"] <- YSR
  ET_ERP_trial[row, "SP2"] <- SP2
}

# Correct data types of ET_ERP_trial for analyses
ET_ERP_trial$SEGA_ID <- as.factor(ET_ERP_trial$SEGA_ID)
ET_ERP_trial$gender <- as.factor(ET_ERP_trial$gender)

# Remove unused factor level
ET_ERP_trial$manipulation <- droplevels(ET_ERP_trial$manipulation)

# Scaled values for correlation plots of dependend variables
ET_ERP_trial$z_rpd <- as.numeric(scale(ET_ERP_trial$rpd))
ET_ERP_trial$z_rpd_low <- as.numeric(scale(ET_ERP_trial$rpd_low))
ET_ERP_trial$z_MMN_amplitude <- as.numeric(scale(ET_ERP_trial$MMN_amplitude))
ET_ERP_trial$z_MMN_latency <- as.numeric(scale(ET_ERP_trial$MMN_latency))
ET_ERP_trial$z_P3a_amplitude <- as.numeric(scale(ET_ERP_trial$P3a_amplitude))
ET_ERP_trial$z_P3a_latency <- as.numeric(scale(ET_ERP_trial$P3a_latency))

# Save .Rds file for analysis on trial level
saveRDS(ET_ERP_trial,file=paste0(home_path,project_path,'/data/preprocessed_auditory_ET_ERP_trial.rds'))

## CONDITION LEVEL
# List and read EEG files
data_files_eeg <- list.files(
  path = datapath_eeg, full.names = TRUE)

# store diff-files in separate variable because of different data structure
data_files_diff_MMN <- data_files_eeg[grepl("*_diff", data_files_eeg)]
data_files_eeg <- data_files_eeg[!grepl('diff', data_files_eeg)] # remove the diff-files

list_eeg_data <- lapply(
  data_files_eeg,
  read.table,
  header = TRUE,
  fill = TRUE,
  skip = 2)

df_MMN_500oddball <- list_eeg_data[[1]]
df_MMN_750oddball <- list_eeg_data[[2]]
df_P3a_500oddball <- list_eeg_data[[3]]
df_P3a_750oddball <- list_eeg_data[[4]]
df_P3b_500oddball <- list_eeg_data[[5]]
df_P3b_750oddball <- list_eeg_data[[6]]

# Change column names for 500Hz-MMN data frame
names(df_MMN_500oddball) <- c(
  "id",
  "L-Peak_MMN_Oddball_500_Hz_before",
  "L-Peak_MMN_Standard_750_Hz_before",
  "L-Peak_MMN_Oddball_rev_750_Hz_before",
  "L-Peak_MMN_Standard_rev_500_Hz_before",
  "L-Peak_MMN_Oddball_500_Hz_after",
  "L-Peak_MMN_Standard_750_Hz_after",
  "L-Peak_MMN_Oddball_rev_750_Hz_after",
  "L-Peak_MMN_Standard_rev_500_Hz_after",
  "Peak_MMN_Oddball_500_Hz_before",
  "Peak_MMN_Standard_750_Hz_before",
  "Peak_MMN_Oddball_rev_750_Hz_before",
  "Peak_MMN_Standard_rev_500_Hz_before",
  "Peak_MMN_Oddball_500_Hz_after",
  "Peak_MMN_Standard_750_Hz_after",
  "Peak_MMN_Oddball_rev_750_Hz_after",
  "Peak_MMN_Standard_rev_500_Hz_after"
)

# Change column names for 750Hz-MMN data frame
names(df_MMN_750oddball) <- c(
  "id",
  "L-Peak_MMN_Oddball_750_Hz_before",
  "L-Peak_MMN_Standard_500_Hz_before",
  "L-Peak_MMN_Oddball_rev_500_Hz_before",
  "L-Peak_MMN_Standard_rev_750_Hz_before",
  "L-Peak_MMN_Oddball_750_Hz_after",
  "L-Peak_MMN_Standard_500_Hz_after",
  "L-Peak_MMN_Oddball_rev_500_Hz_after",
  "L-Peak_MMN_Standard_rev_750_Hz_after",
  "Peak_MMN_Oddball_750_Hz_before",
  "Peak_MMN_Standard_500_Hz_before",
  "Peak_MMN_Oddball_rev_500_Hz_before",
  "Peak_MMN_Standard_rev_750_Hz_before",
  "Peak_MMN_Oddball_750_Hz_after",
  "Peak_MMN_Standard_500_Hz_after",
  "Peak_MMN_Oddball_rev_500_Hz_after",
  "Peak_MMN_Standard_rev_750_Hz_after"
)

# Change column names for 500Hz-P3a data frame
names(df_P3a_500oddball) <- c(
  "id",
  "L-Peak_P3a_Oddbal_500_Hz_before",
  "L-Peak_P3a_Standard_750_Hz_before",
  "L-Peak_P3a_Oddball_rev_750_Hz_before",
  "L-Peak_P3a_Standard_rev_500_Hz_before",
  "L-Peak_P3a_Oddball_500_Hz_after",
  "L-Peak_P3a_Standard_750_Hz_after",
  "L-Peak_P3a_Oddball_rev_750_after",
  "L-Peak_P3a_Standard_rev_500_Hz_after",
  "Peak_P3a_Oddbal_500_Hz_before",
  "Peak_P3a_Standard_750_Hz_before",
  "Peak_P3a_Oddball_rev_750_Hz_before",
  "Peak_P3a_Standard_rev_500_Hz_before",
  "Peak_P3a_Oddball_500_Hz_after",
  "Peak_P3a_Standard_750_Hz_after",
  "Peak_P3a_Oddball_rev_750_after",
  "Peak_P3a_Standard_rev_500_Hz_after"
)

# Change column names for 750Hz-P3a data frame
names(df_P3a_750oddball) <- c(
  "id",
  "L-Peak_P3a_Oddball_750_Hz_before",
  "L-Peak_P3a_Standard_500_Hz_before",
  "L-Peak_P3a_Oddball_rev_500_Hz_before",
  "L-Peak_P3a_Standard_rev_750_Hz_before",
  "L-Peak_P3a_Oddball_750_Hz_after",
  "L-Peak_P3a_Standard_500_Hz_after",
  "L-Peak_P3a_Oddball_rev_500_Hz_after",
  "L-Peak_P3a_Standard_rev_750_Hz_after",
  "Peak_P3a_Oddball_750_Hz_before",
  "Peak_P3a_Standard_500_Hz_before",
  "Peak_P3a_Oddball_rev_500_Hz_before",
  "Peak_P3a_Standard_rev_750_Hz_before",
  "Peak_P3a_Oddball_750_Hz_after",
  "Peak_P3a_Standard_500_Hz_after",
  "Peak_P3a_Oddball_rev_500_Hz_after",
  "Peak_P3a_Standard_rev_750_Hz_after"
)

# Change column names for 500Hz-P3b data frame
names(df_P3b_500oddball) <- c(
  "id",
  "L-Peak_P3b_Oddball_500_Hz_before",
  "L-Peak_P3b_Standard_750_Hz_before",
  "L-Peak_P3b_Oddball_rev_750_Hz_before",
  "L-Peak_P3b_Standard_rev_500_Hz_before",
  "L-Peak_P3b_Oddball_500_Hz_after",
  "L-Peak_P3b_Standard_750_Hz_after",
  "L-Peak_P3b_Oddball_rev_750_Hz_after",
  "L-Peak_P3b_Standard_rev_500_Hz_after",
  "Peak_P3b_Oddball_500_Hz_before",
  "Peak_P3b_Standard_750_Hz_before",
  "Peak_P3b_Oddball_rev_750_Hz_before",
  "Peak_P3b_Standard_rev_500_Hz_before",
  "Peak_P3b_Oddball_500_Hz_after",
  "Peak_P3b_Standard_750_Hz_after",
  "Peak_P3b_Oddball_rev_750_Hz_after", 
  "Peak_P3b_Standard_rev_500_Hz_after"
)

# Change column names for 750Hz-P3b data frame
names(df_P3b_750oddball) <- c(
  "id",
  "L-Peak_P3b_Oddball_750_Hz_before",
  "L-Peak_P3b_Standard_500_Hz_before",
  "L-Peak_P3b_Oddball_rev_500_Hz_before",
  "L-Peak_P3b_Standard_rev_750_Hz_before",
  "L-Peak_P3b_Oddball_750_Hz_after",
  "L-Peak_P3b_Standard_500_Hz_after",
  "L-Peak_P3b_Oddball_rev_500_Hz_after",
  "L-Peak_P3b_Standard_rev_750_Hz_after",
  "Peak_P3b_Oddball_750_Hz_before",
  "Peak_P3b_Standard_500_Hz_before",
  "Peak_P3b_Oddball_rev_500_Hz_before",
  "Peak_P3b_Standard_rev_750_Hz_before",
  "Peak_P3b_Oddball_750_Hz_after",
  "Peak_P3b_Standard_500_Hz_after",
  "Peak_P3b_Oddball_rev_500Hz_after",
  "Peak_P3b_Standard_rev_750_Hz_after"
)

# Remove empty columns + rows
df_MMN_500oddball <- Filter(function(x)!all(is.na(x)), df_MMN_500oddball)
df_MMN_500oddball <- df_MMN_500oddball[!(df_MMN_500oddball$`L-Peak_MMN_Oddball_500_Hz_before` == "???"), ]
df_MMN_750oddball <- Filter(function(x)!all(is.na(x)), df_MMN_750oddball)
df_MMN_750oddball <- df_MMN_750oddball[!(df_MMN_750oddball$`L-Peak_MMN_Oddball_750_Hz_before` == "???"), ]

df_P3a_500oddball <- Filter(function(x)!all(is.na(x)), df_P3a_500oddball)
df_P3a_500oddball <- df_P3a_500oddball[!(df_P3a_500oddball$`L-Peak_P3a_Oddbal_500_Hz_before` == "???"), ]
df_P3a_750oddball <- Filter(function(x)!all(is.na(x)), df_P3a_750oddball)
df_P3a_750oddball <- df_P3a_750oddball[!(df_P3a_750oddball$`L-Peak_P3a_Oddball_750_Hz_before` == "???"), ]

df_P3b_500oddball <- Filter(function(x)!all(is.na(x)), df_P3b_500oddball)
df_P3b_500oddball <- df_P3b_500oddball[!(df_P3b_500oddball$`L-Peak_P3b_Oddball_500_Hz_before` == "???"), ]
df_P3b_750oddball <- Filter(function(x)!all(is.na(x)), df_P3b_750oddball)
df_P3b_750oddball <- df_P3b_750oddball[!(df_P3b_750oddball$`L-Peak_P3b_Oddball_750_Hz_before` == "???"), ]

# Replace commas with point decimals
for (i in 2:ncol(df_MMN_500oddball)) {
  df_MMN_500oddball[, i] <- as.numeric(gsub(",", ".", df_MMN_500oddball[, i]))
}
for (i in 2:ncol(df_MMN_750oddball)) {
  df_MMN_750oddball[, i] <- as.numeric(gsub(",", ".", df_MMN_750oddball[, i]))
}
for (i in 2:ncol(df_P3a_500oddball)) {
  df_P3a_500oddball[, i] <- as.numeric(gsub(",", ".", df_P3a_500oddball[, i]))
}
for (i in 2:ncol(df_P3a_750oddball)) {
  df_P3a_750oddball[, i] <- as.numeric(gsub(",", ".", df_P3a_750oddball[, i]))
}
for (i in 2:ncol(df_P3b_500oddball)) {
  df_P3b_500oddball[, i] <- as.numeric(gsub(",", ".", df_P3b_500oddball[, i]))
}
for (i in 2:ncol(df_P3b_750oddball)) {
  df_P3b_750oddball[, i] <- as.numeric(gsub(",", ".", df_P3b_750oddball[, i]))
}

# Split the 4 data frames based on amplitude/latency 
# MMN AMPLITUDE
MMN_amplitude_500 <- df_MMN_500oddball[c(
  "id",
  "Peak_MMN_Oddball_500_Hz_before",
  "Peak_MMN_Standard_750_Hz_before",
  "Peak_MMN_Oddball_rev_750_Hz_before",
  "Peak_MMN_Standard_rev_500_Hz_before",
  "Peak_MMN_Oddball_500_Hz_after",
  "Peak_MMN_Standard_750_Hz_after",
  "Peak_MMN_Oddball_rev_750_Hz_after",
  "Peak_MMN_Standard_rev_500_Hz_after"
)]
reshaped_MMN_amplitude_500 <- reshape::melt(MMN_amplitude_500, id = "id")

# Additional variables and correct 3 SEGA-IDs
reshaped_MMN_amplitude_500$block <- as.factor(ifelse(grepl("rev", reshaped_MMN_amplitude_500$variable), "reverse", "forward"))
reshaped_MMN_amplitude_500$trial <- as.factor(str_extract(reshaped_MMN_amplitude_500$variable, "Oddball|Oddbal|Standard"))
reshaped_MMN_amplitude_500$manipulation <- as.factor(str_extract(reshaped_MMN_amplitude_500$variable, "before|after"))
reshaped_MMN_amplitude_500$pitch <- as.factor(str_extract(reshaped_MMN_amplitude_500$variable, "500|750"))
reshaped_MMN_amplitude_500$id <- str_replace(
  reshaped_MMN_amplitude_500$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_MMN_amplitude_500$id <- str_replace(
  reshaped_MMN_amplitude_500$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_MMN_amplitude_500$id <- str_replace(
  reshaped_MMN_amplitude_500$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_MMN_amplitude_500$SEGA_ID <- as.factor(substr(reshaped_MMN_amplitude_500$id, 6, 8))
reshaped_MMN_amplitude_500$SEGA_ID <- sub("^0+", "", reshaped_MMN_amplitude_500$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_MMN_amplitude_500 <- subset(reshaped_MMN_amplitude_500, select = -c(variable, id))
# Rename column with MMN amplitude values
colnames(reshaped_MMN_amplitude_500)[colnames(reshaped_MMN_amplitude_500) == "value"] <- "MMN_amplitude"

MMN_amplitude_750 <- df_MMN_750oddball[c(
  "id",
  "Peak_MMN_Oddball_750_Hz_before",
  "Peak_MMN_Standard_500_Hz_before",
  "Peak_MMN_Oddball_rev_500_Hz_before",
  "Peak_MMN_Standard_rev_750_Hz_before",
  "Peak_MMN_Oddball_750_Hz_after",
  "Peak_MMN_Standard_500_Hz_after",
  "Peak_MMN_Oddball_rev_500_Hz_after",
  "Peak_MMN_Standard_rev_750_Hz_after"
)]
reshaped_MMN_amplitude_750 <- reshape::melt(MMN_amplitude_750, id = "id")

# Additional variables an correct 3 SEGA-IDs
reshaped_MMN_amplitude_750$block <- as.factor(ifelse(grepl("rev", reshaped_MMN_amplitude_750$variable), "reverse", "forward"))
reshaped_MMN_amplitude_750$trial <- as.factor(str_extract(reshaped_MMN_amplitude_750$variable, "Oddball|Oddbal|Standard"))
reshaped_MMN_amplitude_750$manipulation <- as.factor(str_extract(reshaped_MMN_amplitude_750$variable, "before|after"))
reshaped_MMN_amplitude_750$pitch <- as.factor(str_extract(reshaped_MMN_amplitude_750$variable, "500|750"))
reshaped_MMN_amplitude_750$id <- str_replace(
  reshaped_MMN_amplitude_750$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_MMN_amplitude_750$id <- str_replace(
  reshaped_MMN_amplitude_750$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_MMN_amplitude_750$id <- str_replace(
  reshaped_MMN_amplitude_750$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_MMN_amplitude_750$SEGA_ID <- as.factor(substr(reshaped_MMN_amplitude_750$id, 6, 8))
reshaped_MMN_amplitude_750$SEGA_ID <- sub("^0+", "", reshaped_MMN_amplitude_750$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_MMN_amplitude_750 <- subset(reshaped_MMN_amplitude_750, select = -c(variable, id))

# Rename column with MMN amplitude values
colnames(reshaped_MMN_amplitude_750)[colnames(reshaped_MMN_amplitude_750) == "value"] <- "MMN_amplitude"

# combine MMN Amplitude in MMN_amplitude_df
MMN_amplitude_df <- rbind(reshaped_MMN_amplitude_500, reshaped_MMN_amplitude_750)

# MMM LATENCY
MMN_latency_500 <- df_MMN_500oddball[c(
  "id",
  "L-Peak_MMN_Oddball_500_Hz_before",
  "L-Peak_MMN_Standard_750_Hz_before",
  "L-Peak_MMN_Oddball_rev_750_Hz_before",
  "L-Peak_MMN_Standard_rev_500_Hz_before",
  "L-Peak_MMN_Oddball_500_Hz_after",
  "L-Peak_MMN_Standard_750_Hz_after",
  "L-Peak_MMN_Oddball_rev_750_Hz_after",
  "L-Peak_MMN_Standard_rev_500_Hz_after"
)]
reshaped_MMN_latency_500 <- reshape::melt(MMN_latency_500, id = "id")

# Additional variables an correct 3 SEGA-IDs
reshaped_MMN_latency_500$block <- as.factor(ifelse(grepl("rev", reshaped_MMN_latency_500$variable), "reverse", "forward"))
reshaped_MMN_latency_500$trial <- as.factor(str_extract(reshaped_MMN_latency_500$variable, "Oddball|Oddbal|Standard"))
reshaped_MMN_latency_500$manipulation <- as.factor(str_extract(reshaped_MMN_latency_500$variable, "before|after"))
reshaped_MMN_latency_500$pitch <- as.factor(str_extract(reshaped_MMN_latency_500$variable, "500|750"))
reshaped_MMN_latency_500$id <- str_replace(
  reshaped_MMN_latency_500$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_MMN_latency_500$id <- str_replace(
  reshaped_MMN_latency_500$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_MMN_latency_500$id <- str_replace(
  reshaped_MMN_latency_500$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_MMN_latency_500$SEGA_ID <- as.factor(substr(reshaped_MMN_latency_500$id, 6, 8))
reshaped_MMN_latency_500$SEGA_ID <- sub("^0+", "", reshaped_MMN_latency_500$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_MMN_latency_500 <- subset(reshaped_MMN_latency_500, select = -c(variable, id))
# Rename column with MMN latency values
colnames(reshaped_MMN_latency_500)[colnames(reshaped_MMN_latency_500) == "value"] <- "MMN_latency"

MMN_latency_750 <- df_MMN_750oddball[c(
  "id",
  "L-Peak_MMN_Oddball_750_Hz_before",
  "L-Peak_MMN_Standard_500_Hz_before",
  "L-Peak_MMN_Oddball_rev_500_Hz_before",
  "L-Peak_MMN_Standard_rev_750_Hz_before",
  "L-Peak_MMN_Oddball_750_Hz_after",
  "L-Peak_MMN_Standard_500_Hz_after",
  "L-Peak_MMN_Oddball_rev_500_Hz_after",
  "L-Peak_MMN_Standard_rev_750_Hz_after"
)]
reshaped_MMN_latency_750 <- reshape::melt(MMN_latency_750, id = "id")

# Additional variables an correct 3 SEGA-IDs
reshaped_MMN_latency_750$block <- as.factor(ifelse(grepl("rev", reshaped_MMN_latency_750$variable), "reverse", "forward"))
reshaped_MMN_latency_750$trial <- as.factor(str_extract(reshaped_MMN_latency_750$variable, "Oddball|Oddbal|Standard"))
reshaped_MMN_latency_750$manipulation <- as.factor(str_extract(reshaped_MMN_latency_750$variable, "before|after"))
reshaped_MMN_latency_750$pitch <- as.factor(str_extract(reshaped_MMN_latency_750$variable, "500|750"))
reshaped_MMN_latency_750$id <- str_replace(
  reshaped_MMN_latency_750$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_MMN_latency_750$id <- str_replace(
  reshaped_MMN_latency_750$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_MMN_latency_750$id <- str_replace(
  reshaped_MMN_latency_750$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_MMN_latency_750$SEGA_ID <- as.factor(substr(reshaped_MMN_latency_750$id, 6, 8))
reshaped_MMN_latency_750$SEGA_ID <- sub("^0+", "", reshaped_MMN_latency_750$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_MMN_latency_750 <- subset(reshaped_MMN_latency_750, select = -c(variable, id))
# Rename column with MMN amplitude values
colnames(reshaped_MMN_latency_750)[colnames(reshaped_MMN_latency_750) == "value"] <- "MMN_latency"

MMN_latency_df <- rbind(reshaped_MMN_latency_500, reshaped_MMN_latency_750)

MMN_df <- merge(MMN_amplitude_df, MMN_latency_df, by = c("SEGA_ID", "block", "manipulation", "trial", "pitch"))

# P3a Amplitude
P3a_amplitude_500 <- df_P3a_500oddball[c(
  "id",
  "Peak_P3a_Oddbal_500_Hz_before",
  "Peak_P3a_Standard_750_Hz_before",
  "Peak_P3a_Oddball_rev_750_Hz_before",
  "Peak_P3a_Standard_rev_500_Hz_before",
  "Peak_P3a_Oddball_500_Hz_after",
  "Peak_P3a_Standard_750_Hz_after",
  "Peak_P3a_Oddball_rev_750_after",
  "Peak_P3a_Standard_rev_500_Hz_after"
)]
reshaped_P3a_amplitude_500 <- reshape::melt(P3a_amplitude_500, id = "id")

# Additional variables an correct 3 SEGA-IDs
reshaped_P3a_amplitude_500$block <- as.factor(ifelse(grepl("rev", reshaped_P3a_amplitude_500$variable), "reverse", "forward"))
reshaped_P3a_amplitude_500$trial <- as.factor(str_extract(reshaped_P3a_amplitude_500$variable, "Oddball|Oddbal|Standard"))
reshaped_P3a_amplitude_500$trial <- as.factor(str_replace(reshaped_P3a_amplitude_500$trial, "Oddbal(?!l)", "Oddball")) # correct typo
reshaped_P3a_amplitude_500$manipulation <- as.factor(str_extract(reshaped_P3a_amplitude_500$variable, "before|after"))
reshaped_P3a_amplitude_500$pitch <- as.factor(str_extract(reshaped_P3a_amplitude_500$variable, "500|750"))
reshaped_P3a_amplitude_500$id <- str_replace(
  reshaped_P3a_amplitude_500$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_P3a_amplitude_500$id <- str_replace(
  reshaped_P3a_amplitude_500$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_P3a_amplitude_500$id <- str_replace(
  reshaped_P3a_amplitude_500$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_P3a_amplitude_500$SEGA_ID <- as.factor(substr(reshaped_P3a_amplitude_500$id, 6, 8))
reshaped_P3a_amplitude_500$SEGA_ID <- sub("^0+", "", reshaped_P3a_amplitude_500$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_P3a_amplitude_500 <- subset(reshaped_P3a_amplitude_500, select = -c(variable, id))
# Rename column with MMN amplitude values
colnames(reshaped_P3a_amplitude_500)[colnames(reshaped_P3a_amplitude_500) == "value"] <- "P3a_amplitude"

P3a_amplitude_750 <- df_P3a_750oddball[c(
  "id",
  "Peak_P3a_Oddball_750_Hz_before",
  "Peak_P3a_Standard_500_Hz_before",
  "Peak_P3a_Oddball_rev_500_Hz_before",
  "Peak_P3a_Standard_rev_750_Hz_before",
  "Peak_P3a_Oddball_750_Hz_after",
  "Peak_P3a_Standard_500_Hz_after",
  "Peak_P3a_Oddball_rev_500_Hz_after",
  "Peak_P3a_Standard_rev_750_Hz_after"
)]
reshaped_P3a_amplitude_750 <- reshape::melt(P3a_amplitude_750, id = "id")

reshaped_P3a_amplitude_750$block <- as.factor(ifelse(grepl("rev", reshaped_P3a_amplitude_750$variable), "reverse", "forward"))
reshaped_P3a_amplitude_750$trial <- as.factor(str_extract(reshaped_P3a_amplitude_750$variable, "Oddball|Standard"))
reshaped_P3a_amplitude_750$manipulation <- as.factor(str_extract(reshaped_P3a_amplitude_750$variable, "before|after"))
reshaped_P3a_amplitude_750$pitch <- as.factor(str_extract(reshaped_P3a_amplitude_750$variable, "500|750"))
reshaped_P3a_amplitude_750$id <- str_replace(
  reshaped_P3a_amplitude_750$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_P3a_amplitude_750$id <- str_replace(
  reshaped_P3a_amplitude_750$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_P3a_amplitude_750$id <- str_replace(
  reshaped_P3a_amplitude_750$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_P3a_amplitude_750$SEGA_ID <- as.factor(substr(reshaped_P3a_amplitude_750$id, 6, 8))
reshaped_P3a_amplitude_750$SEGA_ID <- sub("^0+", "", reshaped_P3a_amplitude_750$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_P3a_amplitude_750 <- subset(reshaped_P3a_amplitude_750, select = -c(variable, id))
# Rename column with MMN amplitude values
colnames(reshaped_P3a_amplitude_750)[colnames(reshaped_P3a_amplitude_750) == "value"] <- "P3a_amplitude"

# Combine amplitude + latency in P3a_df
P3_amplitude_df <- rbind(reshaped_P3a_amplitude_500, reshaped_P3a_amplitude_750)

# P3a LATENCY
P3a_latency_500 <- df_P3a_500oddball[c(
  "id",
  "L-Peak_P3a_Oddbal_500_Hz_before",
  "L-Peak_P3a_Standard_750_Hz_before",
  "L-Peak_P3a_Oddball_rev_750_Hz_before",
  "L-Peak_P3a_Standard_rev_500_Hz_before",
  "L-Peak_P3a_Oddball_500_Hz_after",
  "L-Peak_P3a_Standard_750_Hz_after",
  "L-Peak_P3a_Oddball_rev_750_after",
  "L-Peak_P3a_Standard_rev_500_Hz_after"
)]
reshaped_P3a_latency_500 <- reshape::melt(P3a_latency_500, id = "id")

# Additional variables an correct 3 SEGA-IDs
reshaped_P3a_latency_500$block <- as.factor(ifelse(grepl("rev", reshaped_P3a_latency_500$variable), "reverse", "forward"))
reshaped_P3a_latency_500$trial <- as.factor(str_extract(reshaped_P3a_latency_500$variable, "Oddball|Oddbal|Standard"))
reshaped_P3a_latency_500$trial <- as.factor(str_replace(reshaped_P3a_latency_500$trial, "Oddbal(?!l)", "Oddball")) # correct typo
reshaped_P3a_latency_500$manipulation <- as.factor(str_extract(reshaped_P3a_latency_500$variable, "before|after"))
reshaped_P3a_latency_500$pitch <- as.factor(str_extract(reshaped_P3a_latency_500$variable, "500|750"))
reshaped_P3a_latency_500$id <- str_replace(
  reshaped_P3a_latency_500$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_P3a_latency_500$id <- str_replace(
  reshaped_P3a_latency_500$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_P3a_latency_500$id <- str_replace(
  reshaped_P3a_latency_500$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_P3a_latency_500$SEGA_ID <- as.factor(substr(reshaped_P3a_latency_500$id, 6, 8))
reshaped_P3a_latency_500$SEGA_ID <- sub("^0+", "", reshaped_P3a_latency_500$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_P3a_latency_500 <- subset(reshaped_P3a_latency_500, select = -c(variable, id))
# Rename column with MMN amplitude values
colnames(reshaped_P3a_latency_500)[colnames(reshaped_P3a_latency_500) == "value"] <- "P3a_latency"

P3a_latency_750 <- df_P3a_750oddball [c(
  "id",
  "L-Peak_P3a_Oddball_750_Hz_before",
  "L-Peak_P3a_Standard_500_Hz_before",
  "L-Peak_P3a_Oddball_rev_500_Hz_before",
  "L-Peak_P3a_Standard_rev_750_Hz_before",
  "L-Peak_P3a_Oddball_750_Hz_after",
  "L-Peak_P3a_Standard_500_Hz_after",
  "L-Peak_P3a_Oddball_rev_500_Hz_after",
  "L-Peak_P3a_Standard_rev_750_Hz_after"
)]
reshaped_P3a_latency_750 <- reshape::melt(P3a_latency_750, id = "id")

# Additional variables an correct 3 SEGA-IDs
reshaped_P3a_latency_750$block <- as.factor(ifelse(grepl("rev", reshaped_P3a_latency_750$variable), "reverse", "forward"))
reshaped_P3a_latency_750$trial <- as.factor(str_extract(reshaped_P3a_latency_750$variable, "Oddball|Oddbal|Standard"))
reshaped_P3a_latency_750$trial <- as.factor(str_replace(reshaped_P3a_latency_750$trial, "Oddbal(?!l)", "Oddball")) # correct typo
reshaped_P3a_latency_750$manipulation <- as.factor(str_extract(reshaped_P3a_latency_750$variable, "before|after"))
reshaped_P3a_latency_750$pitch <- as.factor(str_extract(reshaped_P3a_latency_750$variable, "500|750"))
reshaped_P3a_latency_750$id <- str_replace(
  reshaped_P3a_latency_750$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_P3a_latency_750$id <- str_replace(
  reshaped_P3a_latency_750$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_P3a_latency_750$id <- str_replace(
  reshaped_P3a_latency_750$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_P3a_latency_750$SEGA_ID <- as.factor(substr(reshaped_P3a_latency_750$id, 6, 8))
reshaped_P3a_latency_750$SEGA_ID <- sub("^0+", "", reshaped_P3a_latency_750$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_P3a_latency_750 <- subset(reshaped_P3a_latency_750, select = -c(variable, id))
# Rename column with MMN amplitude values
colnames(reshaped_P3a_latency_750)[colnames(reshaped_P3a_latency_750) == "value"] <- "P3a_latency"

P3a_latency_df <- rbind(reshaped_P3a_latency_500, reshaped_P3a_latency_750)

# Combine amplitude + latency in P3a_df
P3a_df <- merge(P3_amplitude_df, P3a_latency_df, by = c("SEGA_ID", "block", "manipulation", "trial", "pitch"))

# P3b AMPLITUDE
P3b_amplitude_500 <- df_P3b_500oddball[c(
  "id",
  "Peak_P3b_Oddball_500_Hz_before",
  "Peak_P3b_Standard_750_Hz_before",
  "Peak_P3b_Oddball_rev_750_Hz_before",
  "Peak_P3b_Standard_rev_500_Hz_before",
  "Peak_P3b_Oddball_500_Hz_after",
  "Peak_P3b_Standard_750_Hz_after",
  "Peak_P3b_Oddball_rev_750_Hz_after",
  "Peak_P3b_Standard_rev_500_Hz_after"
)]
reshaped_P3b_amplitude_500 <- reshape::melt(P3b_amplitude_500, id = "id")

# Additional variables and correct 3 SEGA-IDs
reshaped_P3b_amplitude_500$block <- as.factor(ifelse(grepl("rev", reshaped_P3b_amplitude_500$variable), "reverse", "forward"))
reshaped_P3b_amplitude_500$trial <- as.factor(str_extract(reshaped_P3b_amplitude_500$variable, "Oddball|Oddbal|Standard"))
reshaped_P3b_amplitude_500$manipulation <- as.factor(str_extract(reshaped_P3b_amplitude_500$variable, "before|after"))
reshaped_P3b_amplitude_500$pitch <- as.factor(str_extract(reshaped_P3b_amplitude_500$variable, "500|750"))
reshaped_P3b_amplitude_500$id <- str_replace(
  reshaped_P3b_amplitude_500$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_P3b_amplitude_500$id <- str_replace(
  reshaped_P3b_amplitude_500$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_P3b_amplitude_500$id <- str_replace(
  reshaped_P3b_amplitude_500$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_P3b_amplitude_500$SEGA_ID <- as.factor(substr(reshaped_P3b_amplitude_500$id, 6, 8))
reshaped_P3b_amplitude_500$SEGA_ID <- sub("^0+", "", reshaped_P3b_amplitude_500$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_P3b_amplitude_500 <- subset(reshaped_P3b_amplitude_500, select = -c(variable, id))
# Rename column with MMN amplitude values
colnames(reshaped_P3b_amplitude_500)[colnames(reshaped_P3b_amplitude_500) == "value"] <- "P3b_amplitude"

P3b_amplitude_750 <- df_P3b_750oddball[c(
  "id",
  "Peak_P3b_Oddball_750_Hz_before",
  "Peak_P3b_Standard_500_Hz_before",
  "Peak_P3b_Oddball_rev_500_Hz_before",
  "Peak_P3b_Standard_rev_750_Hz_before",
  "Peak_P3b_Oddball_750_Hz_after",
  "Peak_P3b_Standard_500_Hz_after",
  "Peak_P3b_Oddball_rev_500Hz_after",
  "Peak_P3b_Standard_rev_750_Hz_after"
)]
reshaped_P3b_amplitude_750 <- reshape::melt(P3b_amplitude_750, id = "id")

# Additional variables an correct 3 SEGA-IDs
reshaped_P3b_amplitude_750$block <- as.factor(ifelse(grepl("rev", reshaped_P3b_amplitude_750$variable), "reverse", "forward"))
reshaped_P3b_amplitude_750$trial <- as.factor(str_extract(reshaped_P3b_amplitude_750$variable, "Oddball|Oddbal|Standard"))
reshaped_P3b_amplitude_750$manipulation <- as.factor(str_extract(reshaped_P3b_amplitude_750$variable, "before|after"))
reshaped_P3b_amplitude_750$pitch <- as.factor(str_extract(reshaped_P3b_amplitude_750$variable, "500|750"))
reshaped_P3b_amplitude_750$id <- str_replace(
  reshaped_P3b_amplitude_750$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_P3b_amplitude_750$id <- str_replace(
  reshaped_P3b_amplitude_750$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_P3b_amplitude_750$id <- str_replace(
  reshaped_P3b_amplitude_750$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_P3b_amplitude_750$SEGA_ID <- as.factor(substr(reshaped_P3b_amplitude_750$id, 6, 8))
reshaped_P3b_amplitude_750$SEGA_ID <- sub("^0+", "", reshaped_P3b_amplitude_750$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_P3b_amplitude_750 <- subset(reshaped_P3b_amplitude_750, select = -c(variable, id))

# Rename column with MMN amplitude values
colnames(reshaped_P3b_amplitude_750)[colnames(reshaped_P3b_amplitude_750) == "value"] <- "P3b_amplitude"

# combine MMN Amplitude in MMN_amplitude_df
P3b_amplitude_df <- rbind(reshaped_P3b_amplitude_500, reshaped_P3b_amplitude_750)

# P3b LATENCY
P3b_latency_500 <- df_P3b_500oddball[c(
  "id",
  "L-Peak_P3b_Oddball_500_Hz_before",
  "L-Peak_P3b_Standard_750_Hz_before",
  "L-Peak_P3b_Oddball_rev_750_Hz_before",
  "L-Peak_P3b_Standard_rev_500_Hz_before",
  "L-Peak_P3b_Oddball_500_Hz_after",
  "L-Peak_P3b_Standard_750_Hz_after",
  "L-Peak_P3b_Oddball_rev_750_Hz_after",
  "L-Peak_P3b_Standard_rev_500_Hz_after"
)]
reshaped_P3b_latency_500 <- reshape::melt(P3b_latency_500, id = "id")

# Additional variables an correct 3 SEGA-IDs
reshaped_P3b_latency_500$block <- as.factor(ifelse(grepl("rev", reshaped_P3b_latency_500$variable), "reverse", "forward"))
reshaped_P3b_latency_500$trial <- as.factor(str_extract(reshaped_P3b_latency_500$variable, "Oddball|Oddbal|Standard"))
reshaped_P3b_latency_500$manipulation <- as.factor(str_extract(reshaped_P3b_latency_500$variable, "before|after"))
reshaped_P3b_latency_500$pitch <- as.factor(str_extract(reshaped_P3b_latency_500$variable, "500|750"))
reshaped_P3b_latency_500$id <- str_replace(
  reshaped_P3b_latency_500$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_P3b_latency_500$id <- str_replace(
  reshaped_P3b_latency_500$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_P3b_latency_500$id <- str_replace(
  reshaped_P3b_latency_500$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_P3b_latency_500$SEGA_ID <- as.factor(substr(reshaped_P3b_latency_500$id, 6, 8))
reshaped_P3b_latency_500$SEGA_ID <- sub("^0+", "", reshaped_P3b_latency_500$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_P3b_latency_500 <- subset(reshaped_P3b_latency_500, select = -c(variable, id))
# Rename column with MMN latency values
colnames(reshaped_P3b_latency_500)[colnames(reshaped_P3b_latency_500) == "value"] <- "P3b_latency"

P3b_latency_750 <- df_P3b_750oddball[c(
  "id",
  "L-Peak_P3b_Oddball_750_Hz_before",
  "L-Peak_P3b_Standard_500_Hz_before",
  "L-Peak_P3b_Oddball_rev_500_Hz_before",
  "L-Peak_P3b_Standard_rev_750_Hz_before",
  "L-Peak_P3b_Oddball_750_Hz_after",
  "L-Peak_P3b_Standard_500_Hz_after",
  "L-Peak_P3b_Oddball_rev_500_Hz_after",
  "L-Peak_P3b_Standard_rev_750_Hz_after"
)]
reshaped_P3b_latency_750 <- reshape::melt(P3b_latency_750, id = "id")

# Additional variables an correct 3 SEGA-IDs
reshaped_P3b_latency_750$block <- as.factor(ifelse(grepl("rev", reshaped_P3b_latency_750$variable), "reverse", "forward"))
reshaped_P3b_latency_750$trial <- as.factor(str_extract(reshaped_P3b_latency_750$variable, "Oddball|Oddbal|Standard"))
reshaped_P3b_latency_750$manipulation <- as.factor(str_extract(reshaped_P3b_latency_750$variable, "before|after"))
reshaped_P3b_latency_750$pitch <- as.factor(str_extract(reshaped_P3b_latency_750$variable, "500|750"))
reshaped_P3b_latency_750$id <- str_replace(
  reshaped_P3b_latency_750$id, "SEGA_AuditoryOddball_022_24102022", "SEGA_022_AuditoryOddball_24102022")
reshaped_P3b_latency_750$id <- str_replace(
  reshaped_P3b_latency_750$id, "SEGA_AuditoryOddball_046_24102022", "SEGA_046_AuditoryOddball_24102022")
reshaped_P3b_latency_750$id <- str_replace(
  reshaped_P3b_latency_750$id, "SEGA_003_05.01.23_Auditory", "SEGA_003_AuditoryOddball_05012023")
reshaped_P3b_latency_750$SEGA_ID <- as.factor(substr(reshaped_P3b_latency_750$id, 6, 8))
reshaped_P3b_latency_750$SEGA_ID <- sub("^0+", "", reshaped_P3b_latency_750$SEGA_ID)

# Delete unnecessary 2 columns
reshaped_P3b_latency_750 <- subset(reshaped_P3b_latency_750, select = -c(variable, id))
# Rename column with MMN amplitude values
colnames(reshaped_P3b_latency_750)[colnames(reshaped_P3b_latency_750) == "value"] <- "P3b_latency"

P3b_latency_df <- rbind(reshaped_P3b_latency_500, reshaped_P3b_latency_750)

P3b_df <- merge(P3b_amplitude_df, P3b_latency_df, by = c("SEGA_ID", "block", "manipulation", "trial", "pitch"))

# Combine MMN + P3a in ERP_df (outer join because of subjects having only one ERP)
ERP_df <- merge(P3a_df, MMN_df, by = c("SEGA_ID", "block", "manipulation", "trial", "pitch"), all = TRUE)
# Add P3b to ERP_df
ERP_df <- merge(ERP_df, P3b_df, by = c("SEGA_ID", "block", "manipulation", "trial", "pitch"), all = TRUE)

# Include experimental group from .csv file in ERP_df
for (row in 1:length(ERP_df$SEGA_ID)) {
  SEGA_ID_df <- ERP_df[row, "SEGA_ID"]
  row_number <- which(exp_groups$SEGA_ID == SEGA_ID_df)
  group <- exp_groups[row_number, "group"]
  ERP_df[row, "group"] <- group
}
ERP_df$group <- as.factor(ERP_df$group)

# Include handdynamometer z-value from .csv file
for (row in 1:length(ERP_df$SEGA_ID)) {
  SEGA_ID_df <- ERP_df[row, "SEGA_ID"]
  row_number <- which(exp_groups$SEGA_ID == SEGA_ID_df)
  z_handdynamometer <- exp_groups[row_number, "z_handdynamometer"]
  ERP_df[row, "z_handdynamometer"] <- z_handdynamometer
}

# ET_df contains mean BPS + SEPR for 8 conditions per subject.
unique_ids <- unique(ET_df_trial$SEGA_ID)
ET_df <- NULL

for (i in unique_ids){
  i_subject_df <- ET_df_trial[ET_df_trial$SEGA_ID == i, ]
  manipulation_block_df <- i_subject_df[i_subject_df$manipulation == "before" & i_subject_df$block == "forward", ]
  stimulus_df <- manipulation_block_df[manipulation_block_df$trial == "oddball", ]
  
  SEGA_ID <- unique(as.numeric(i))
  block <- unique(stimulus_df$block)
  manipulation <- unique(stimulus_df$manipulation)
  trial <- unique(stimulus_df$trial)
  rpd <- mean(stimulus_df$rpd, na.rm = TRUE)
  rpd_block <- mean(stimulus_df$rpd_block, na.rm = TRUE)
  rpd_low <- mean(stimulus_df$rpd_low, na.rm = TRUE)
  print(i)
  ET_df <- rbind(ET_df, data.frame(SEGA_ID, block, manipulation, trial, rpd, rpd_block, rpd_low))
  manipulation_block_df <- i_subject_df[i_subject_df$manipulation == "before" & i_subject_df$block == "forward", ]
  stimulus_df <- manipulation_block_df[manipulation_block_df$trial == "standard", ]
  
  SEGA_ID <- unique(as.numeric(i))
  block <- unique(stimulus_df$block)
  manipulation <- unique(stimulus_df$manipulation)
  trial <- unique(stimulus_df$trial)
  rpd <- mean(stimulus_df$rpd, na.rm = TRUE)
  rpd_block <- mean(stimulus_df$rpd_block, na.rm = TRUE)
  rpd_low <- mean(stimulus_df$rpd_low, na.rm = TRUE)
  print(paste0("2nd", i))
  ET_df <- rbind(ET_df, data.frame(SEGA_ID, block, manipulation, trial, rpd, rpd_block, rpd_low))
  
  manipulation_block_df <- i_subject_df[i_subject_df$manipulation == "before" & i_subject_df$block == "reverse", ]
  stimulus_df <- manipulation_block_df[manipulation_block_df$trial == "oddball", ]
  
  SEGA_ID <- unique(as.numeric(i))
  block <- unique(stimulus_df$block)
  manipulation <- unique(stimulus_df$manipulation)
  trial <- unique(stimulus_df$trial)
  rpd <- mean(stimulus_df$rpd, na.rm = TRUE)
  rpd_block <- mean(stimulus_df$rpd_block, na.rm = TRUE)
  rpd_low <- mean(stimulus_df$rpd_low, na.rm = TRUE)
  print(paste0("3rd", i))
  ET_df <- rbind(ET_df, data.frame(SEGA_ID, block, manipulation, trial, rpd, rpd_block, rpd_low))
  
  manipulation_block_df <- i_subject_df[i_subject_df$manipulation == "before" & i_subject_df$block == "reverse", ]
  stimulus_df <- manipulation_block_df[manipulation_block_df$trial == "standard", ]
  
  SEGA_ID <- unique(as.numeric(i))
  block <- unique(stimulus_df$block)
  manipulation <- unique(stimulus_df$manipulation)
  trial <- unique(stimulus_df$trial)
  rpd <- mean(stimulus_df$rpd, na.rm = TRUE)
  rpd_block <- mean(stimulus_df$rpd_block, na.rm = TRUE)
  rpd_low <- mean(stimulus_df$rpd_low, na.rm = TRUE)
  print(paste0("4th", i))
  ET_df <- rbind(ET_df, data.frame(SEGA_ID, block, manipulation, trial, rpd, rpd_block, rpd_low))
  
  manipulation_block_df <- i_subject_df[i_subject_df$manipulation == "after" & i_subject_df$block == "forward", ]
  stimulus_df <- manipulation_block_df[manipulation_block_df$trial == "oddball", ]
  
  SEGA_ID <- unique(as.numeric(i))
  block <- unique(stimulus_df$block)
  manipulation <- unique(stimulus_df$manipulation)
  trial <- unique(stimulus_df$trial)
  rpd <- mean(stimulus_df$rpd, na.rm = TRUE)
  rpd_block <- mean(stimulus_df$rpd_block, na.rm = TRUE)
  rpd_low <- mean(stimulus_df$rpd_low, na.rm = TRUE)
  print(paste0("5th", i))
  ET_df <- rbind(ET_df, data.frame(SEGA_ID, block, manipulation, trial, rpd, rpd_block, rpd_low))
  
  manipulation_block_df <- i_subject_df[i_subject_df$manipulation == "after" & i_subject_df$block == "forward", ]
  stimulus_df <- manipulation_block_df[manipulation_block_df$trial == "standard", ]
  
  SEGA_ID <- unique(as.numeric(i))
  block <- unique(stimulus_df$block)
  manipulation <- unique(stimulus_df$manipulation)
  trial <- unique(stimulus_df$trial)
  rpd <- mean(stimulus_df$rpd, na.rm = TRUE)
  rpd_block <- mean(stimulus_df$rpd_block, na.rm = TRUE)
  rpd_low <- mean(stimulus_df$rpd_low, na.rm = TRUE)
  print(paste0("6th", i))
  ET_df <- rbind(ET_df, data.frame(SEGA_ID, block, manipulation, trial, rpd, rpd_block, rpd_low))
  
  manipulation_block_df <- i_subject_df[i_subject_df$manipulation == "after" & i_subject_df$block == "reverse", ]
  stimulus_df <- manipulation_block_df[manipulation_block_df$trial == "oddball", ]
  
  SEGA_ID <- unique( as.numeric(i))
  block <- unique(stimulus_df$block)
  manipulation <- unique(stimulus_df$manipulation)
  trial <- unique(stimulus_df$trial)
  rpd <- mean(stimulus_df$rpd, na.rm = TRUE)
  rpd_block <- mean(stimulus_df$rpd_block, na.rm = TRUE)
  rpd_low <- mean(stimulus_df$rpd_low, na.rm = TRUE)
  print(paste0("7th", i))
  ET_df <- rbind(ET_df, data.frame(SEGA_ID, block, manipulation, trial, rpd, rpd_block, rpd_low))
  
  manipulation_block_df <- i_subject_df[i_subject_df$manipulation == "after" & i_subject_df$block == "reverse", ]
  stimulus_df <- manipulation_block_df[manipulation_block_df$trial == "standard", ]
  
  SEGA_ID <- unique(as.numeric(i))
  block <- unique(stimulus_df$block)
  manipulation <- unique(stimulus_df$manipulation)
  trial <- unique(stimulus_df$trial)
  rpd <- mean(stimulus_df$rpd, na.rm = TRUE)
  rpd_block <- mean(stimulus_df$rpd_block, na.rm = TRUE)
  rpd_low <- mean(stimulus_df$rpd_low, na.rm = TRUE)
  print(paste0("8th", i))
  ET_df <- rbind(ET_df, data.frame(SEGA_ID, block, manipulation, trial, rpd, rpd_block, rpd_low))
}

# Trial variable with capital letters due to congruence with eeg format
ET_df$trial <- str_to_sentence(ET_df$trial)

# Include experimental group from .csv file in ET_df
for (row in 1:length(ET_df$SEGA_ID)) {
  SEGA_ID_df <- ET_df[row, "SEGA_ID"]
  row_number <- which(exp_groups$SEGA_ID == SEGA_ID_df)
  group <- exp_groups[row_number, "group"]
  ET_df[row, "group"] <- group
}
ET_df$group <- as.factor(ET_df$group)

# Include IQ from .csv file
for (row in 1:length(ET_df$SEGA_ID)) {
  SEGA_ID_df <- ET_df[row, "SEGA_ID"]
  print(paste("id", SEGA_ID_df))
  row_number <- which(exp_groups$SEGA_ID == SEGA_ID_df)
  MT <- exp_groups[row_number, "MT"]
  GF <- exp_groups[row_number, "GF"]
  MZ <- exp_groups[row_number, "MZ"]
  WT <- exp_groups[row_number, "WT"]
  GF_WT <- c(GF, WT)
  MT_MZ <- c(MT, MZ)
  verbal_IQ <- mean(GF_WT, na.rm = T)
  non_verbal_IQ <- mean(MT_MZ, na.rm = T)
  ET_df[row, "verbal_IQ"] <- verbal_IQ
  ET_df[row, "non_verbal_IQ"] <- non_verbal_IQ}

ET_ERP_subject <- merge(ET_df, ERP_df, by = c("SEGA_ID", "block", "trial", "manipulation", 'group'), all = T)

# Add sample characteristics to ET_ERP_subject
for (row in 1:length(ET_ERP_subject$SEGA_ID)) {
  SEGA_ID_df <- ET_ERP_subject[row, "SEGA_ID"]
  row_number <- which(sample_characteristics$SEGA_ID == SEGA_ID_df)
  gender <- sample_characteristics[row_number, "Geschlecht_Index"]
  age <- sample_characteristics[row_number, "age"]
  SCQ <- sample_characteristics[row_number, "FSK_ges"]
  SRS <- sample_characteristics[row_number, "SRS_Gesamtwert_N_k_RW"]
  CBCL <- sample_characteristics[row_number, "CBCL_T_GES"]
  YSR <- sample_characteristics[row_number, "YSR_T_GES"]
  SP2 <- sample_characteristics[row_number, "SP2_Eltern_RW_Auditiv"]
  ET_ERP_subject[row, "gender"] <- gender
  ET_ERP_subject[row, "age"] <- age
  ET_ERP_subject[row, "SCQ"] <- SCQ
  ET_ERP_subject[row, "SRS"] <- SRS
  ET_ERP_subject[row, "CBCL"] <- CBCL
  ET_ERP_subject[row, "YSR"] <- YSR
  ET_ERP_subject[row, "SP2"] <- SP2
}

# Correct data types for analysis
ET_ERP_subject$trial <- as.factor(ET_ERP_subject$trial)
ET_ERP_subject$gender <- as.factor(ET_ERP_subject$gender)
ET_ERP_subject$SEGA_ID <- as.factor(ET_ERP_subject$SEGA_ID)
ET_ERP_subject$block <- as.factor(ET_ERP_subject$block)

# Scaled values for correlation plots of dependent variables
ET_ERP_subject$z_rpd <- as.numeric(scale(ET_ERP_subject$rpd))
ET_ERP_subject$z_rpd_low <- as.numeric(scale(ET_ERP_subject$rpd_low))
ET_ERP_subject$z_MMN_amplitude <- as.numeric(scale(ET_ERP_subject$MMN_amplitude))
ET_ERP_subject$z_MMN_latency <- as.numeric(scale(ET_ERP_subject$MMN_latency))
ET_ERP_subject$z_P3a_amplitude <- as.numeric(scale(ET_ERP_subject$P3a_amplitude))
ET_ERP_subject$z_P3a_latency <- as.numeric(scale(ET_ERP_subject$P3a_latency))
ET_ERP_subject$z_P3b_amplitude <- as.numeric(scale(ET_ERP_subject$P3b_amplitude))
ET_ERP_subject$z_P3b_latency <- as.numeric(scale(ET_ERP_subject$P3b_latency))
ET_ERP_subject$z_rpd_block <- as.numeric(scale(ET_ERP_subject$rpd_block))

# Remove unused factor level
ET_ERP_subject$manipulation <- droplevels(ET_ERP_subject$manipulation)

# Save .Rds file for analysis on subject level
saveRDS(ET_ERP_subject,file=paste0(home_path,project_path,'/data/preprocessed_auditory_ET_ERP_subject.rds'))


saveRDS(df, file = paste0(home_path, project_path, '/data/preprocessed_auditory_df.rds'))
