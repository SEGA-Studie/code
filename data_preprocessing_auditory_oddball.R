# Script purpose: READ AND ANALYZE DATA - SEGA PROJECT - AUDITORY ODDBALL
# Author: Nico Bast
# Date Created: `r paste(Sys.Date())`
# Copyright (c) Nico Bast, `r paste(format(Sys.Date(), "%Y"))`
# Email: nico.bast@kgu.de

###notes on

# - 027 without task data
# - 037 task aborted
# - 056 program failed
# - 016 --> et daten fehlten bisher

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
  # List all .hdf and .csv files
  data_files <- c(
    list.files(path = datapath, full.names = TRUE),
    list.files(path = datapath_task, full.names = TRUE))
}

# DATA IMPORT AND RESHAPING ####

# Get eye tracking data and store them in a list of df (one per subject)
data_files_et <- data_files[grepl(".hdf5", data_files)]
list_et_data <- list(0)
for (i in 1:length(data_files_et)) {
  list_et_data[[i]] <- h5read(
    file = data_files_et[i],
    name = "data_collection/events/eyetracker/BinocularEyeSampleEvent")
  print(paste0("read ET data file: ", i))
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

# Split trial data per trial -> list of all trials of all subjects.
# Each list element contains a df of trial data
list_split_trial <- pblapply(df_list, function(x) {
  split(x, x$trial_number)})
list_split_trial <- unlist(list_split_trial, recursive = FALSE)

# New variable ts_trial: Timestamp for each et event within a trial
list_split_trial <- lapply(list_split_trial, function(x) {
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

# Sample size
n_asd <- 0
n_td <- 0

for (subject in df_list){
  if (subject$group[1] == 'TD')
    n_td <<- n_td + 1
  if (subject$group[1] == 'ASD')
    n_asd <<- n_asd + 1
}

# trial-baseline correction and number of missing trial baselines
counter <- 0
counter_td <- 0
counter_asd <- 0

list_split_trial <- lapply(list_split_trial, function(x) {
  rpd_low <- mean(x$pd[x$ts_trial < 0.250])
  if (is.na(unique(rpd_low))) {
    counter <<- counter + 1
    if (x$group[1] == 'TD')
      counter_td <<- counter_td + 1
    if (x$group[1] == 'ASD')
      counter_asd <<- counter_asd + 1
  }
  
  trial_corr_pd <- x$pd - rpd_low
  x[, "rpd_low"] <- rep(rpd_low, times = nrow(x))
  x[, "trial_corr_pd"] <- trial_corr_pd
  return(x)
})

# Percentage of missing trial baseline
percentage_trial <- 100 / (length(list_split_trial))
print(paste0(round(percentage_trial * counter, digits = 2),
             " % of all trials with trial baseline == NA: ",
             round(percentage_trial * counter_asd, digits = 2), "% (ASD), ",
             round(percentage_trial * counter_td, digits = 2), "% (TD)"
             ))

# Bind trials together to a df
df <- dplyr::bind_rows(list_split_trial)

# split by block and id
list_split_blocks <- split(df,
droplevels(interaction(df$block_counter, df$id)))

list_split_blocks <- lapply(list_split_blocks, function(x) {
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
  ggplot(
    df[is.finite(df$pd), ],
    aes(x = pd)) + geom_histogram(bins = 100) + facet_wrap(~id)
  
  
# reduce ET data to per trial data (and merge with df_trial)
# pupil diameter late in trial duration
# rpd_high <- sapply(list_split_trial, function(x) {
#   mean(x$pd[x$ts_trial > 0.875 & x$ts_trial < 1.125])})

require(DescTools) #AUC function  
rpd_auc <- sapply(list_split_trial, function(x) {
  AUC(x$ts_trial,x$pd,na.rm=T)})
rpd_auc<-ifelse(rpd_auc>15,NA,rpd_auc)

rpd_high <- sapply(list_split_trial, function(x) {
  mean(x$pd[x$ts_trial > 0.500 & x$ts_trial < 1])})

rpd_low <- sapply(list_split_trial, function(x) {
  mean(x$pd[x$ts_trial < 0.250])})

# news variables to index specific trials
trial_number <- sapply(list_split_trial, function(x) {
  unique(x$trial_number)})
trial_number_in_block <- sapply(list_split_trial, function(x) {
  unique(x$trial_number_in_block)})

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

# new variables manipulation indicates whether before or after manipulation
df_trial$manipulation <- factor(
  ifelse(df_trial$block_counter < 8, "before",
  ifelse(df_trial$block_counter > 8, "after", "manipulation")),
  levels = c("before", "after"))

# distinguish pitch: oddball or standard frequency
df_trial$pitch <- ifelse(df_trial$trial == "oddball",
df_trial$oddball_frequency, df_trial$standard_frequency)
# distinguish between (forward) oddball blocks and reverse oddball blocks
df_trial$reverse <- ifelse(
  df_trial$phase == "oddball_block_rev", "reverse", "forward")
# distinguish between oddball and standards trials
df_trial$oddball <- as.factor(ifelse(grepl(
  "oddball",
  df_trial$trial),
  "oddball", "standard"))


### VISUALIZATION ####

##manipulation check
sampled_rows<-sample(1:nrow(df),nrow(df)/100) #randomly select rows
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

length(unique(df$id))

##--> save preprocess df_trial ####
saveRDS(df_trial,file=paste0(home_path,project_path,'/data/preprocessed_auditory_ETdata_12092023.rds'))
# Can be used to skip preprocessing and directly read proprocessed data from .rds file:
df_trial <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ETdata_12092023.rds'))

#changed random intercept to a factor
df_trial$id<-as.factor(df_trial$id) #change ID to factor
#df_trial$trial<-as.factor(df_trial$trial)

require(performance)

### Data plausibility ####

table(df_trial$block_counter, df_trial$phase)
table(df_trial$trial)
# Number of trials for forward and reverse oddball blocks
table(df_trial$reverse, df_trial$oddball)
hist(df_trial$rpd, 50)
with(df_trial, by(rpd, trial, mean, na.rm = TRUE))

### DATA ANALYSIS  ####

lmm <- lmer(
  rpd ~ oddball * manipulation * reverse + trial_number_in_block + (1 | id),
  data = df_trial[
    df_trial$phase %in% c("oddball_block", "oddball_block_rev"), ])
summary(lmm)
contrast(emmeans(lmm, ~ oddball), "pairwise")

ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev") &
  is.finite(df_trial$rpd), ], aes(
    x = trial_number,
    y = rpd,
    group = oddball,
    color = oddball)) +
  geom_smooth() +
  facet_grid(rows = vars(reverse), cols = vars(manipulation))

ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev") &
  df_trial$trial_number_in_block >= 10 & is.finite(
    df_trial$rpd), ],
  aes(x = trial_number_in_block, y = rpd, group = oddball, color = oddball)) +
  geom_smooth() +
  facet_grid(rows = vars(reverse), cols = vars(manipulation)) + theme_bw()

ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev") &
  is.finite(df_trial$rpd), ],
  aes(x = as.factor(trial_number_in_block), y = rpd)) +
  geom_boxplot() +
  facet_grid(
    rows = vars(reverse),
    cols = vars(manipulation, oddball)) +
    theme_bw()

# DATA ANALYSIS: Linear mixed model - trial data ####


table(df_trial$phase)

#pupillary response - SEPR
lmm<-lmer(scale(rpd)~oddball*group*manipulation*reverse+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm)

  contrast(emmeans(lmm,~oddball),'pairwise') ##--> higher response to oddball
  contrast(emmeans(lmm,~reverse),'pairwise') ##--> higher response in forward trials
  contrast(emmeans(lmm,~oddball|reverse),'pairwise') ##--> oddball response is specific to forward blocks
  contrast(emmeans(lmm,~reverse|oddball),'pairwise') ###--> oddball larger in forward 
  
  # #these effects are not significant - as of Sept 2023
  # contrast(emmeans(lmm,~group),'pairwise') ##--> TD with higher response to all trials
  # contrast(emmeans(lmm,~group|manipulation),'pairwise') ##--> higher response in TD before manipulation
  # contrast(emmeans(lmm,~manipulation|group),'pairwise') ##--> manipulation has an  effect in ASD
  # contrast(emmeans(lmm,~group|oddball+manipulation+reverse),'pairwise')
  # ###--> TD compared to ASD with higher response to standards in foward trials before manipulation
  # ###--> TD comapred to ASD with higehr repsonse to oddballs in reverse trials before manipulation
  
#baseline  - BPS
lmm<-lmer(scale(rpd_low)~oddball*group*manipulation*reverse+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm)

    #--> does not differ by stimulus
    contrast(emmeans(lmm,~manipulation),'pairwise') ##--> higher pupil size after manipulation
    confint(contrast(emmeans(lmm,~manipulation|group),'pairwise')) ##--> this effect of manipulation is emphasized in ASD
    
    contrast(emmeans(lmm,~reverse),'pairwise') ##--> higher pupil size in foward trials
    confint(contrast(emmeans(lmm,~reverse|group),'pairwise')) ##--> this difference is emphasized in TD
    confint(contrast(emmeans(lmm,~reverse|manipulation),'pairwise')) ##--> this difference is emphasized before the manipulation
    
#area under the curve  
lmm<-lmer(rpd_auc~oddball*group*manipulation*reverse+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm)

    contrast(emmeans(lmm,~manipulation),'pairwise') ##--> higher response after manipulation
    confint(contrast(emmeans(lmm,~manipulation|group),'pairwise')) ##--> higher effect of manipulation in ASD
    
    confint(contrast(emmeans(lmm,~reverse|group),'pairwise')) 
    
# DATA ANALYSIS: Baseline phase ####
df_baseline <- df[df$phase %in% c("baseline", "baseline_calibration"), ]
df_baseline <- df_baseline[is.finite(df_baseline$pd), ]

# Number of pupil data per baseline trial (white, black and grey slide)
with(df_baseline[df_baseline$phase == "baseline_calibration", ], table(trial))
# Check: Block counter should be in line with aforementioned table
with(df_baseline, table(trial, block_counter))
with(df_baseline, by(
  timestamp_exp, interaction(baseline_trial_counter, phase, id),
  mean, na.rm = TRUE))

# subject baseline is defined as mean of baseline trial during calibration phase
df_meanpd <- with(
  df_baseline[df_baseline$phase == "baseline_calibration" &
  df_baseline$trial == "baseline" , ], by(as.numeric(pd), id, mean))
df_meanpd <- data.frame(names(df_meanpd), as.numeric(df_meanpd))
names(df_meanpd) <- c("id", "mean_pd")
df_baseline <- merge(df_baseline, df_meanpd, by = "id")

# Visualization: LAPR group comparison
# New df with baseline trial only from calibration phase
df_baseline_calibration <- df_baseline[
  df_baseline$phase == "baseline_calibration", ]

# Plot 1: LAPR_Subplots with trial_baseline corrected pd
ggplot(
# Baseline trial duration is theoretically 5 seconds.
df_baseline_calibration[df_baseline_calibration$ts_trial < 5.0, ],
aes(x = ts_trial,
y = trial_corr_pd,
group = group,
color = group)) +
geom_smooth() +
theme_bw() +
xlab("trial duration [s]") +
ylab("trial-corr. pupil dilation [mm]") +
labs(title = "Trial-corrected LAPR for group and trial") +
facet_wrap(~ factor(trial,
levels = c("baseline", "baseline_whiteslide", "baseline_blackslide")))

ggsave(
  "output/Plot_2_lapr_subplots_trial_corrected.tiff",
  device = "tiff",
  width = 6,
  height = 4,
  units = "in",
  compression = "lzw",
  dpi = 800
)

### plot with pd --> trial corrected effects are reversed
ggplot(
  # Baseline trial duration is theoretically 5 seconds.
  df_baseline_calibration[df_baseline_calibration$ts_trial < 5.0, ],
  aes(x = ts_trial,
      y = pd,
      group = group,
      color = group)) +
  geom_smooth() +
  theme_bw() +
  xlab("trial duration [s]") +
  ylab("trial-corr. pupil dilation [mm]") +
  labs(title = "LAPR for group and trial") +
  facet_wrap(~ factor(trial,
                      levels = c("baseline", "baseline_whiteslide", "baseline_blackslide")))


# DATA ANLAYSIS: Oddball phase ####
table(df$phase)

df_oddball <- df[df$phase %in% c("oddball_block", "oddball_block_rev"), ]
table(df_oddball$trial)

# merge with baseline PD
df_oddball <- merge(df_oddball, df_meanpd, by = "id")

# subject-baseline corrected pd
df_oddball$rpd <- df_oddball$pd - df_oddball$mean_pd
hist(df_oddball$rpd, 50)

# correct for non-finite values
df_oddball <- df_oddball[is.finite(df_oddball$rpd), ]

# trial type (standard versus oddball)
df_oddball$trial_type <- substr(df_oddball$trial, 1, 7)

table(df_oddball$block_counter, df_oddball$trial)
df_oddball$manipulation <- factor(ifelse(df_oddball$block_counter < 8, "before",
ifelse(df_oddball$block_counter > 8, "after", "manipulation")),
levels = c("before", "after"))

df_oddball$order <- ifelse(grepl("rev", df_oddball$trial), "reverse", "normal")
table(df_oddball$order, df_oddball$block_counter)

# define trials in oddball phase
df_oddball$trial_index_oddballphase <- with(
  # name to indentify individual trials
  df_oddball, interaction(block_counter, .thisRepN, .thisTrialN))
# all trials including empty trials
df_oddball$trial_number_oddballphase <- with(
  df_oddball, rep(seq_along(table(trial_index_oddballphase)),
  times = table(trial_index_oddballphase)))
df_oddball$trial_number_oddballphase <- with(
  df_oddball, rep(seq_along(table(trial_number_oddballphase)),
# rerun remove empty trials from trial number
times = table(trial_number_oddballphase)))

# Number of pupil data per trial number
table(df_oddball$trial_number)
hist(df_oddball$trial_number_oddballphase)

# pupillary response by trial
ggplot(
  df_oddball,
  aes(x = trial_number_oddballphase, y = rpd)) +
  geom_smooth(method = "lm") + facet_wrap(~trial + manipulation)

ggplot(df_oddball,
aes(rpd, fill = interaction(manipulation))) +
geom_density(alpha = 0.2) + facet_wrap(~trial)

# Visualization 1: Manipulation effect
tiff(file = paste0(
  home_path,
  project_path,
  "/output/figure_audio_effectofmanipulation_pupilresponse_testdata.tiff"),
  width = 8,
  height = 4,
  units = "in",
  res = 300,
  compression = "lzw")
# Only forward (normal) oddball blocks, not reverse

ggplot(
  df_oddball[df_oddball$ts_trial < 1.8 & df_oddball$order == "normal", ],
  aes(x = ts_trial,
      y = scale(rpd),
      group = trial_type,
      color = trial_type)) +
  geom_smooth() +
  theme_bw() +
  labs(x = "trial duration (s)",
       y = "standardized pupil response (z)",
  title = "effect of manipulation of pupil response") +
  facet_wrap(~manipulation)

dev.off()

# Visualization 2: Group effect
tiff(file = paste0(
  home_path,
  project_path,
  "/output/group_effect.tiff"),
  width = 8,
  height = 4,
  units = "in",
  res = 300,
  compression = "lzw")
# Only forward (normal) oddball blocks, not reverse
ggplot(
  df_oddball[df_oddball$ts_trial < 1.8 & df_oddball$order == "normal", ],
  aes(x = ts_trial,
      y = scale(rpd),
      group = interaction(group, trial_type),
      color = group,
      linetype = trial_type)) +
  geom_smooth() +
  theme_bw() +
  labs(x = "trial duration (s)",
       y = "standardized pupil response (z)",
       title = "effect of manipulation of pupil response")

dev.off()

# Visualization 3: Interaction effect
tiff(file = paste0(
  home_path,
  project_path,
  "/output/interaction_effect.tiff"),
  width = 8,
  height = 4,
  units = "in",
  res = 300,
  compression = "lzw")
# Only forward (normal) oddball blocks, not reverse
ggplot(
  df_oddball[df_oddball$ts_trial < 1.8 & df_oddball$order == "normal", ],
  aes(x = ts_trial,
      y = scale(rpd),
      group = interaction(group, trial_type),
      color = group,
      linetype = trial_type)) +
  geom_smooth() +
  theme_bw() +
  labs(x = "trial duration (s)",
       y = "standardized pupil response (z)",
       title = "effect of manipulation of pupil response") +
  facet_wrap(~manipulation)
dev.off()

# LMM: Interaction effect of manipulation, trial_type and group
# with random intercept id.
lmm_interaction <- lmerTest::lmer(
  rpd ~ trial_type * manipulation * group + (1|id),
  data = df_oddball)
anova(lmm_interaction) # p-value

# Get confidence interval
confint(contrast(emmeans::emmeans(lmm_interaction, ~ group + trial_type + manipulation), method = 'pairwise'))
confint(contrast(emmeans::emmeans(lmm_interaction, ~ manipulation + trial_type | group), method = 'pairwise'))

# neurophysiological habituation to standard trials within block
ggplot(
  df_oddball[df_oddball$ts_trial < 2 & is.finite(df_oddball$rpd) &
  df_oddball$trial_type == "standard", ],
  aes(x = .thisRepN,
  y = scale(rpd),
  group = interaction(manipulation, order),
  color = interaction(manipulation, order))) +
  geom_smooth() + theme_bw()

tiff(file = paste0(home_path,
project_path,
"/output/figure_audio_neurophysiological_habituation_standards_Jan2023.tiff"),
width = 8,
height = 4,
units = "in",
res = 300,
compression = "lzw")

ggplot(
  df_oddball[df_oddball$ts_trial < 2 & is.finite(df_oddball$rpd), ],
  aes(x = trial_number_in_block,
  y = scale(rpd),
  group = interaction(manipulation, order, trial_type),
  color = interaction(manipulation, order),
  linetype = trial_type)) +
  geom_smooth() + theme_bw() +
  xlab("trial within block") +
  ylab("relative pupil size (z)")

dev.off()

# DATA ANALYSIS: gaze data - raw pupil data ####
hist(df_oddball$left_pupil_measure1)
hist(df_oddball$right_pupil_measure1)
# raw gaze data
hist(
  df_oddball$left_gaze_x[df_oddball$left_gaze_x < 200 &
  df_oddball$left_gaze_x > -200])
hist(
  df_oddball$right_gaze_x[df_oddball$right_gaze_x < 200 &
  df_oddball$right_gaze_x > -200])
hist(
  df_oddball$left_gaze_y[df_oddball$left_gaze_y < 200 &
  df_oddball$left_gaze_y > -200])
hist(
  df_oddball$right_gaze_y[df_oddball$right_gaze_y < 200 &
  df_oddball$right_gaze_y > -200])

# refined heatmap
ggplot(
  df_oddball[df_oddball$ts_trial < 2 & is.finite(df_oddball$rpd) &
  is.finite(df_oddball$left_gaze_x) & is.finite(df_oddball$left_gaze_y), ],
  aes(x = left_gaze_x, y = left_gaze_y)) +
  stat_density_2d(aes(fill = after_stat(density)),
  n = 50,
  geom = "raster",
  contour = FALSE) +
  xlim(-200, 200) +
  ylim(-200, 200) +
  scale_fill_gradientn(colours = rev(rainbow(3)))

# DATA ANALYSIS: manipulation phase ####
df_manip <- df[df$phase %in% c("manipulation_block"), ]

table(df_manip$trial)
table(df_manip$manipulation_trial_counter)
with(df_manip, table(trial, manipulation_trial_counter))

# merge with baseline PD
df_manip <- merge(df_manip, df_meanpd, by = "id")

# baselinecorrected pd
df_manip$rpd <- df_manip$pd - df_manip$mean_pd
hist(df_manip$rpd)

# correct for non-finite values
df_manip <- df_manip[is.finite(df_manip$rpd), ]

# define different baseline scenarios
df_manip$trial_type <- with(
  df_manip, ifelse(
    manipulation_trial_counter %in% c(3, 7, 11, 15, 19),
    "baseline_after_squeeze", ifelse(
      manipulation_trial_counter %in% c(5, 9, 13, 17), "baseline_after_relax",
      ifelse(manipulation_trial_counter %in% c(2, 6, 10, 14, 18), "squeeze",
      ifelse(manipulation_trial_counter %in% c(4, 8, 12, 16, 20),
      "relax", "baseline_start")))))

# 18 is duration of squeeze phase
ggplot(
  df_manip[df_manip$ts_trial < 18 & df_manip$trial != "baseline", ],
  aes(x = ts_trial,
  y = rpd,
  group = trial,
  color = trial)) +
  geom_smooth() +
  theme_bw() +
  labs(x = "trial duration (s)",
  y = "standardized pupillary response (z)",
  title = "pupil response between manipulation phases")
ggplot(
  df_manip[df_manip$trial_type != "baseline_start", ],
  aes(x = trial_type,
  y = rpd,
  fill = trial_type)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "manipulation phase (s)",
  y = "standardized pupillary response (z)",
  title = "pupil response between manipulation phases")

ggplot(
  df_manip,
  aes(x = manipulation_trial_counter,
  y = rpd,
  group = manipulation_trial_counter,
  fill = trial_type)) +
  geom_boxplot() +
  theme_bw()

# Model
lmm <- lmer(
  rpd ~ trial_type + (1 | manipulation_trial_counter),
  data = df_manip[df_manip$trial_type != "baseline_start", ])
anova(lmm)
plot(contrast(emmeans(lmm, ~ trial_type), "pairwise"))

# DATA ANALYSIS: PRELIMINARY EEG analysis ####
# Reag eeg data
 data_files_eeg <- list.files(
   path = datapath_eeg, full.names = TRUE)
 list_eeg_data <- lapply(
   data_files_eeg,
   read.table,
   header = TRUE,
   fill = TRUE,
   skip = 2)

 df_MMN_500oddball <- list_eeg_data[[1]]
 df_MMN_750oddball <- list_eeg_data[[2]]
 df_P3_500oddball <- list_eeg_data[[3]]
 df_P3_750oddball <- list_eeg_data[[4]]
 
 
 df_MMN_500oddball <- reshape2::melt(df_MMN_500oddball, id.vars = "File")
 df_MMN_750oddball <- reshape2::melt(df_MMN_750oddball, id.vars = "File")
 df_P3_500oddball <- reshape2::melt(df_P3_500oddball, id.vars = "File")
 df_P3_750oddball <- reshape2::melt(df_P3_750oddball, id.vars = "File")
 
 df_erp <- rbind(
   df_MMN_500oddball,
   df_MMN_750oddball,
   df_P3_500oddball,
   df_P3_750oddball)
 
 # define variables
 df_erp$electrode <- substr(df_erp$variable, 1, 7)
 df_erp$manipulation <- substr(
   df_erp$variable,
   nchar(as.character(df_erp$variable)) - 5,
   nchar(as.character(df_erp$variable)))
 
 df_erp$pitch <- substr(
   df_erp$variable, nchar(as.character(df_erp$variable)) - 12,
   nchar(as.character(df_erp$variable)) - 9)
 df_erp$reverse <- as.logical(
   ifelse(
     nchar(as.character(df_erp$variable)) > 42, "T", "F"))
 
 df_erp$erp <- substr(df_erp$variable, 18, 20)
 df_erp$erp <- ifelse(df_erp$erp == "_MM", "MMN", df_erp$erp)
 df_erp$trial <- substr(df_erp$variable, 22, 26)
 df_erp$trial <- ifelse(grepl("db", df_erp$trial), "oddball", "standard")
 
 # remove characters
 df_erp$electrode <- gsub("_", "", df_erp$electrode)
 df_erp$manipulation <- gsub("_", "", df_erp$manipulation)
 df_erp$pitch <- gsub("_", "", df_erp$pitch)
 df_erp$erp <- gsub("_", "", df_erp$erp)
 
 # change decimal - ERP amplitude variable
 df_erp$value <- gsub(",", ".", df_erp$value)
 df_erp$value <- ifelse(df_erp$value == "???", NA, df_erp$value)
 df_erp$value <- as.numeric(df_erp$value)
 df_erp$erp_amplitude <- df_erp$value
 df_erp <- df_erp[, !(names(df_erp) == "value")]
 
 #extract id
 df_erp$id <- ifelse(
   grepl("SEGA_AuditoryOddball", df_erp$File),
   substr(df_erp$File, 21, 24), substr(df_erp$File, 6, 8))
 df_erp$id <- as.numeric(gsub("_", "", df_erp$id))
 
 # Which id has EEG and also ET data?
 names(table(df_erp$id)) %in% names(table(df$id)) # 3 # 4
 names(table(df$id)) %in% names(table(df_erp$id)) # 32 # 82
 table(df$id)
 
 # converge EEG with eye tracking oddball data on condition level
 # prepare ET data - allign variables
 df_oddball$manipulation <- ifelse(
   df_oddball$block_counter < 8, "before", "after")
 df_oddball$pitch <- ifelse(
   df_oddball$trial == "oddball",
   df_oddball$oddball_frequency,
   df_oddball$standard_frequency)
 df_oddball$reverse <- as.logical(ifelse(
   df_oddball$block_counter %in% c(5, 12), "T", "F"))
 
 # create aggregated varaibles per condition (id, manipulation, reverse, phase)
 df_oddball$merger_id <- with(
   df_oddball, interaction(id, manipulation, reverse, trial_type))
 oddball_per_condition_high <- with(
   df_oddball[df_oddball$ts_trial > 0.75 & df_oddball$ts_trial < 1.25, ],
   by(as.numeric(rpd), merger_id, mean, na.rm = TRUE))
 oddball_per_condition_low <- with(
   df_oddball[df_oddball$ts_trial > 0 & df_oddball$ts_trial < 0.25, ], by(
     as.numeric(rpd), merger_id, mean, na.rm = TRUE))
 oddball_per_condition <- as.numeric(
   oddball_per_condition_high - oddball_per_condition_low)
 
 id_per_condition <- as.character(
   with(df_oddball, by(id, merger_id, head, n = 1)))
 manipulation_per_condition <- as.character(
   with(df_oddball, by(manipulation, merger_id, head, n = 1)))
 reverse_per_condition <- as.logical(
   with(df_oddball, by(reverse, merger_id, head, n = 1)))
 trial_type_per_condition <- as.character(
   with(df_oddball, by(trial_type, merger_id, head, n = 1)))
 
 df_rpd_agg <- data.frame(
   id_per_condition,
   manipulation_per_condition,
   reverse_per_condition,
   oddball_per_condition,
   trial_type_per_condition)
 names(df_rpd_agg) <- c("id", "manipulation", "reverse", "rpd", "trial_type")
 df_rpd_agg$trial_type <- ifelse(
   df_rpd_agg$trial_type == "oddball", "oddball", "standard")
 
 # rather plausible
 # merge
 df_erp$merger_id <- with(
   df_erp, interaction(id, manipulation, trial, reverse))
 df_rpd_agg$merger_id <- with(
   df_rpd_agg, interaction(id, manipulation, trial_type, reverse))
 df_rpd_agg <- df_rpd_agg[, !(names(df_rpd_agg) %in% c(
   "id", "manipulation", "trial_type", "reverse"))]
 df_erp_rpd_agg <- merge(df_rpd_agg, df_erp, by = "merger_id")
 
 # pupillary response
 ggplot(
   df_erp_rpd_agg,
   aes(x = interaction(
     manipulation, trial),
     y = rpd,
     fill = interaction(manipulation, trial))) +
     geom_violin() +
     geom_boxplot(alpha = 0.4) +
     facet_wrap(~reverse)
 
 ggplot(
   df_erp_rpd_agg, aes(rpd, fill = trial)) +
   geom_density(alpha = 0.2,
   adjust = 2) +
   facet_wrap(~reverse + manipulation)
 
 # erp
 ggplot(
   df_erp_rpd_agg[df_erp_rpd_agg$electrode == "L.Peak", ],
   aes(x = interaction(manipulation, trial),
   y = erp_amplitude,
   fill = interaction(manipulation, trial))) +
   geom_violin() +
   geom_boxplot(alpha = 0.4) +
   facet_wrap(~reverse)
 
 # display
 ggplot(
   df_erp_rpd_agg[df_erp_rpd_agg$electrode == "L.Peak", ],
   aes(x = rpd,
   y = erp_amplitude)) +
   geom_point() +
   geom_smooth(method = "lm") +
   facet_wrap(~trial + reverse + manipulation)
 ggplot(
   df_erp_rpd_agg[df_erp_rpd_agg$electrode == "Pz.Peak", ],
   aes(x = rpd,
   y = erp_amplitude)) +
   geom_point() +
   geom_smooth(method = "lm") +
   facet_wrap(~trial + reverse + manipulation)
 ggplot(
   df_erp_rpd_agg[df_erp_rpd_agg$electrode == "L.Peak" &
   df_erp_rpd_agg$trial == "standard" &
   df_erp_rpd_agg$reverse == FALSE, ],
   aes(x = scale(rpd),
   y = scale(erp_amplitude))) +
   geom_point() +
   geom_smooth(method = "lm",
   color = "chocolate2",
   fill = "chocolate1") +
   xlab("pupillary response (z)") +
   ylab("P3 amplitude (z)") +
   theme_bw()
 ggplot(
   df_erp_rpd_agg[df_erp_rpd_agg$electrode == "L.Peak", ],
   aes(x = rpd,
   y = erp_amplitude)) +
   geom_point() +
   geom_smooth(method = "lm")

# Visualization
 g1 <- ggplot(
   df_baseline[df_baseline$phase == "baseline_calibration" &
   df_baseline$ts_trial < 6, ],
   aes(x = ts_trial,
   y = pd,
   group = trial,
   color = trial)) +
   geom_smooth() +
   theme_bw() +
   xlab("time (s)") +
   ylab("pupil size (mm)") +
   labs(title = "change of PD during initial baseline between IDs")
 g2 <- ggplot(
   df_manip[df_manip$ts_trial < 18 & df_manip$trial != "baseline", ],
   aes(x = ts_trial,
   y = rpd,
   group = trial,
   color = trial)) +
   geom_smooth() +
   theme_bw() +
   labs(x = "trial duration (s)",
   y = "standardized pupillary response (z)",
   title = "pupil response during manipulation")
 g3 <- ggplot(
   df_manip[df_manip$trial_type != "baseline_start", ],
   aes(x = trial_type,
   y = rpd,
   fill = trial_type)) +
   geom_boxplot() +
   theme_bw() +
   labs(x = "manipulation phase (s)",
   y = "standardized pupillary response (z)",
   title = "pupil response in manipulation phases")
 g4 <- ggplot(
   df_oddball[df_oddball$ts_trial < 1.8 & df_oddball$order == "normal", ],
   aes(x = ts_trial,
   y = scale(rpd),
   group = trial_type,
   color = trial_type)) +
   geom_smooth() +
   theme_bw() +
   labs(x = "trial duration (s)",
   y = "standardized pupil response (z)",
   title = "effect of manipulation on stimulus-evoked pupillary response
   (before vs. after manipulation)") +
   facet_wrap(~manipulation) +
   scale_color_brewer(palette = "Dark2")
 tiff(
   file = paste0(home_path, project_path,
   "/output/figure_audio_results_testdata.tiff"),
   width = 12,
   height = 16,
   units = "in",
   res = 300,
   compression = "lzw")
 grid.arrange(g1, g2, g3, g4,
 layout_matrix = rbind(c(1, 1), c(2, 3), c(4, 4)))
 dev.off()

print("The End")

# DATA ANALYSIS: PRELIMINARY CORTISOL analysis ####

### - load hair and saliva cortisol data from xlsx ####

require(readxl)
#hair cortisol
df_hair<-read_xlsx(path="C:/Users/nico/PowerFolders/project_sega/data/hair_cortisol_1439_23082023.xlsx",
                   skip=11) #starting reading in line 12

  #change variable to appropriate length
  df_hair$`externe ID`<-substr(df_hair$`externe ID`,nchar(df_hair$`externe ID`)-7,nchar(df_hair$`externe ID`))
  
  #change names
  names(df_hair)<-c('id','hair_mass','hair_comments','hair_cortisol','hair_cortisone')
  
#saliva cortisol
df_saliva<-read_xlsx(path="C:/Users/nico/PowerFolders/project_sega/data/saliva_cortisol_20230912_Results_SEGA_Salimetrics.xlsx",
                     skip=9) #starting reading in line 10

  names(df_saliva)<-c('saliva_sample_id','id','saliva_date','saliva_time',
                      'saliva_cortisol_1','saliva_cortisol_2','saliva_cortisol_avg','saliva_cv','saliva_comments')

  #categorize sample
  df_saliva$saliva_when<-rep(c('before','after'),nrow(df_saliva)/2)
  
  #get time as numeric variable
  df_saliva$saliva_time_num<-as.numeric(format(df_saliva$saliva_time, "%H"))
  
  
### - create SEGA id variable to match df_hair and df_saliva ####
  
sega_id<-ifelse(nchar(levels(df_trial$id))==2,
                paste0('0',levels(df_trial$id)),
                ifelse(nchar(levels(df_trial$id))==1,
                       paste0('00',levels(df_trial$id)),levels(df_trial$id)))

sega_id<-paste0('SEGA_',sega_id)

levels(df_trial$id)<-sega_id


### - adapt SEGA id variable to mean pd data

df_meanpd$sega_id<-ifelse(nchar(df_meanpd$id)==2,
             paste0('0',df_meanpd$id),
             ifelse(nchar(df_meanpd$id)==1,
                    paste0('00',df_meanpd$id),df_meanpd$id))

df_meanpd$sega_id<-paste0('SEGA_',df_meanpd$sega_id)


### - estimate pupillometry measures from trial data #####
df_agg<-aggregate(df_trial[,names(df_trial) %in% c('rpd_auc','rpd_high','rpd_low','rpd')],
                by=with(df_trial,list(id,oddball,manipulation,reverse)),
                mean,na.rm=T)

names(df_agg)<-c('id','oddball','manipulation','reverse','rpd_auc','rpd_high','rpd_low','rpd')


### - saliva change ####

saliva_cortisol_change<-with(df_saliva,saliva_cortisol_avg[saliva_when=='after']-
                               saliva_cortisol_avg[saliva_when=='before'])
  
id<-unique(df_saliva$id)

df_saliva_change<-data.frame(id,saliva_cortisol_change)


### - merge data ####
df_agg_hair<-merge(df_agg,df_hair,by='id')
df_agg_hair<-merge(df_agg_hair,df_meanpd,by.x='id',by.y='sega_id')
df_agg_hair<-merge(df_agg_hair,df_saliva_change,by='id',all.x=T)


df_agg_saliva<-merge(df_agg,df_saliva,by='id')
df_agg<-merge(df_agg_saliva,df_hair,by='id')
df_agg<-merge(df_agg,df_meanpd,by.x='id',by.y='sega_id')
df_agg<-merge(df_agg,df_saliva_change,by='id')


### - check data distribution ####

hist(df_hair$hair_cortisol)

ggplot(df_hair,aes(hair_cortisone,hair_cortisol,color=hair_mass))+geom_point()+
  labs(x='hair cortisone (pg/mg, log-scaled)', y='hair cortisol (pg/mg, log-scaled)', color='sample mass (mg)')+
  geom_smooth(color='grey',method='lm')+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme_bw()

hist(df_hair$hair_cortisol[df_hair$hair_cortisol<5])

hist(df_saliva$saliva_cortisol_avg)
hist(df_saliva$saliva_cortisol_avg[df_saliva$saliva_when=='before'])
hist(df_saliva$saliva_cortisol_avg[df_saliva$saliva_when=='after'])
hist(df_saliva$saliva_cv)

ggplot(df_saliva,aes(saliva_cortisol_1,saliva_cortisol_2,color=scale(saliva_cv)))+geom_point()+
  labs(x='saliva sample A (nmol/l, log-scaled)', y='saliva sample B (nmol/l, log-scaled)', color='variation index (z)')+
  geom_smooth(color='grey',method='lm')+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme_bw()



### - data analysis ####

length(unique(df_agg$id))
names(df_agg)
require(performance)

#association of hair variables
summary(lm(scale(hair_cortisol)~scale(hair_mass)+scale(hair_cortisone),df_agg[!duplicated(df_agg$id),]))
###--> hair cortisone and cortisol highly correlated (0.98)

#association of saliva variables
lmm<-lmer(scale(saliva_cortisol_avg)~saliva_when*saliva_time_num+(1|id),df_agg)
summary(lmm)
r2_nakagawa(lmm)


confint(contrast(emmeans(lmm,~saliva_when),'pairwise'))
###--> saliva cortisol reduction after experiment
contrast(emtrends(lmm,~saliva_when,var = 'saliva_time_num'),'pairwise')
###--> time of day is more strongly associated with cortisol after experiment

##effect of experiment
ggplot(df_agg,aes(x=factor(saliva_when,level=c('before','after')),log(saliva_cortisol_avg),fill=saliva_when))+
  geom_violin(adjust=2,alpha=0.5)+geom_boxplot(notch=T,width=0.4)+
  labs(x='saliva sample',y='cortisol (nmol/l, log-scaled)')+
  theme_bw()

ggplot(df_agg,aes(x=factor(saliva_when,level=c('before','after')),log(saliva_cortisol_avg),group=id))+
  geom_point()+geom_line()+
  labs(x='saliva sample',y='cortisol (nmol/l), log-scaled')+
  theme_bw()


###asocaition of hair and saliva
lmm<-lmer(scale(saliva_cortisol_avg)~saliva_when*saliva_time_num*scale(hair_cortisol)+(1|id),df_agg)
summary(lmm)
r2_nakagawa(lmm)

emtrends(lmm,~saliva_time_num|saliva_when,var = 'hair_cortisol',at=list(saliva_time_num=c(10,14,18)))
#before the experiment: 
# - higher hair cortisol is associated with higher saliva cortisol in morning measurements
# - higher hair cortisol is associated with lower saliva cortisol in evening measurements

    hist(df_agg$saliva_time_num)
    ggplot(df_agg[df_agg$saliva_time_num<12,],aes(hair_cortisol,saliva_cortisol_avg))+geom_point()+geom_smooth(method='lm')
    ggplot(df_agg[df_agg$saliva_time_num>16,],aes(hair_cortisol,saliva_cortisol_avg))+geom_point()+geom_smooth(method='lm')
    ###--> to few data points to draw strong conclusion

###asocaition of pupillometric and saliva and hair
lm<-lm(scale(mean_pd)~scale(hair_cortisol),
          df_agg_hair[!(duplicated(df_agg_hair$id)),])
summary(lm)

lm<-lm(scale(mean_pd)~scale(saliva_cortisol_avg)*saliva_when,
       df_agg[!(duplicated(interaction(df_agg$id,df_agg$saliva_when))),])
summary(lm)
###--> pd_baseline not associated with saliva and hair cortisol

### - SEPR and hair and saliva cortisol ####

df_agg$hair_cortisol_z<-scale(df_agg$hair_cortisol) #z standardize before to use emtrends later
lmm<-lmer(scale(rpd)~(oddball*manipulation*reverse)*hair_cortisol_z+(1|id),df_agg)
anova(lmm)
r2_nakagawa(lmm)

contrast(emmeans(lmm,~oddball),'pairwise')
###--> substantial oddball effect
emtrends(lmm,~manipulation,var = 'hair_cortisol_z')
confint(contrast(emtrends(lmm,~manipulation,var = 'hair_cortisol_z'),'pairwise'))
###--> higher hair cortisol is associated with higher rpd responses to all trials before manipulation

ggplot(df_agg[df_agg$manipulation=='before' &
                df_agg$hair_cortisol<15,],aes(hair_cortisol,rpd))+geom_point()+geom_smooth(method='lm')+
  labs(x='hair corisol (pg/mg)',y='pupillary response (z)')+
  theme_bw()
###--> rather driven by outliers

df_agg$saliva_cortisol_z<-scale(df_agg$saliva_cortisol_avg) #z standardize before to use emtrends later
lmm<-lmer(scale(rpd)~(oddball*manipulation*reverse)*saliva_cortisol_z+saliva_cv+(1|id),
          df_agg)
anova(lmm)
r2_nakagawa(lmm)

plot(emtrends(lmm,~oddball|manipulation+reverse,var='saliva_cortisol_z'))
##--> higher saliva cortisol is associated with higher oddball response before manipulation
##--> higher saliva cortisol is associated with lower oddball response after manipulation

  ggplot(df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward' &
                  df_agg$saliva_cortisol_avg<10,],
         aes(x=saliva_cortisol_avg,y=rpd,color=oddball))+geom_point()+geom_smooth(method='lm')+
    labs(x='salivary cortisol (nmol/l)',y='pupillary repsonse (mm)')+
    theme_bw()
  #--> higher saliva cortisol differentiates pupillary responses

  ggplot(df_agg[df_agg$manipulation=='after' & df_agg$reverse=='forward' &
                  df_agg$saliva_cortisol_avg<10,],
         aes(x=saliva_cortisol_avg,y=rpd,color=oddball))+geom_point()+geom_smooth(method='lm')+
    labs(x='saliva cortisol (nmol/l)',y='pupillary repsonse (mm)')+
    theme_bw()
  ###--> higher saliva cortisol might attenuate the effect of the manipulation

  
### - BPS and hair and saliva cortisol ####  
  
lmm<-lmer(scale(rpd_low)~(oddball*manipulation*reverse)*hair_cortisol_z+(1|id),df_agg)
anova(lmm)
r2_nakagawa(lmm)

contrast(emmeans(lmm,~reverse),'pairwise')
emtrends(lmm,~reverse,var='hair_cortisol_z')
confint(contrast(emtrends(lmm,~reverse,var='hair_cortisol_z'),'pairwise'))
###--> higher hair cortisol leads to higher BPS in reverse condition


lmm<-lmer(scale(rpd_low)~(oddball*manipulation*reverse)*saliva_cortisol_z+(1|id),
          df_agg[df_agg$saliva_when=='before',])
anova(lmm)
r2_nakagawa(lmm)

cbind(round(fixef(lmm)['saliva_cortisol_z'],2),
      round(confint(lmm,parm = 'saliva_cortisol_z'),2))
#--> higher saliva cortisol before is associated with lower BPS

lmm<-lmer(scale(rpd_low)~(oddball*manipulation*reverse)*saliva_cortisol_z+(1|id),
          df_agg[df_agg$saliva_when=='after',])
anova(lmm)
r2_nakagawa(lmm)

contrast(emtrends(lmm,~manipulation,var='saliva_cortisol_z'),'pairwise')
#--> the effect of saliva cortisol after is attenuated with the manipulation

### cortisol change ####

lm<-lm(scale(saliva_cortisol_change)~scale(rpd_low),df_agg[!duplicated(df_agg$id),])
summary(lm)


summary(lm(scale(saliva_cortisol_change)~scale(rpd_low)+scale(hair_cortisol),
           df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward',]))
summary(lm(scale(saliva_cortisol_change)~scale(rpd_low)+scale(hair_cortisol),
           df_agg[df_agg$manipulation=='after' & df_agg$reverse=='forward',]))
summary(lm(scale(saliva_cortisol_change)~scale(rpd_low)+scale(hair_cortisol),
           df_agg[df_agg$manipulation=='after' & df_agg$reverse=='reverse',]))
##--> higher BPS associated with attenuated cortisol decrease

summary(lmer(scale(hair_cortisol)~scale(rpd)*oddball*manipulation*reverse+(1|id),df_agg))


###--> analysis in ERC GRant 2024 ####

##assocaition of hair cortisol and pupillary response
df_agg$hair_cortisol_z<-scale(df_agg$hair_cortisol) #z standardize before to use emtrends later
lmm<-lmer(scale(rpd)~(oddball*manipulation*reverse)*hair_cortisol_z+(1|id),df_agg)

anova(lmm)
r2_nakagawa(lmm)

contrast(emmeans(lmm,~oddball),'pairwise')
###--> substantial oddball effect
emtrends(lmm,~manipulation,var = 'hair_cortisol_z')
confint(contrast(emtrends(lmm,~manipulation,var = 'hair_cortisol_z'),'pairwise'))
###--> higher hair cortisol is associated with higher rpd responses to all trials before manipulation

#association of baseline pupil size and saliva cortisol change 
lm<-lm(scale(saliva_cortisol_change)~scale(rpd_low),df_agg[!duplicated(df_agg$id),])
summary(lm)

#association of saliva cortisol before the task and baseline pupil size
df_agg$saliva_cortisol_z<-scale(df_agg$saliva_cortisol_avg) #z standardize before to use emtrends later
lmm<-lmer(scale(rpd_low)~(oddball*manipulation*reverse)*saliva_cortisol_z+(1|id),
          df_agg[df_agg$saliva_when=='before',])
anova(lmm)
r2_nakagawa(lmm)

cbind(round(fixef(lmm)['saliva_cortisol_z'],2),
      round(confint(lmm,parm = 'saliva_cortisol_z'),2))
#--> higher saliva cortisol before is associated with lower BPS

#association of saliva cortisol and pupillary response
df_agg$saliva_cortisol_z<-scale(df_agg$saliva_cortisol_avg) #z standardize before to use emtrends later
lmm<-lmer(scale(rpd)~(oddball*manipulation*reverse)*saliva_cortisol_z+(1|id),
          df_agg)
anova(lmm)
r2_nakagawa(lmm)

plot(emtrends(lmm,~oddball|manipulation+reverse,var='saliva_cortisol_z'))
##--> higher saliva cortisol is associated with higher oddball response before manipulation
##--> higher saliva cortisol is associated with lower oddball response after manipulation

require(wesanderson) #custom color palettes
custom_condition_colors <- wes_palette('Darjeeling1',2,type='discrete') #reverse custom colors to match color coding in other figures

  
ggplot(df_agg[df_agg$manipulation=='before' & df_agg$reverse=='forward' &
                df_agg$saliva_cortisol_avg<10,],
       aes(x=saliva_cortisol_avg,y=rpd,color=oddball))+geom_point()+geom_smooth(method='lm')+
  scale_color_manual(values = custom_condition_colors)+
  labs(x='salivary cortisol (nmol/l)',y='pupillary repsonse (mm)')+
  theme_bw()
#--> higher saliva cortisol differentiates pupillary responses

### ADD demographic data (database export) ####

#data is generated with: data_preprocessing_database.Rmd
df_dem<-readRDS("C:/Users/nico/PowerFolders/project_sega/data/demographics_14092023.rds")

df_agg<-merge(df_agg,df_dem,by='id',all.x=T) # only data with saliva and hair data
df_trial<-merge(df_trial,df_dem,by='id',all.x=T) #all data 
df_trial<-merge(df_trial,df_meanpd,by.x='id',by.y='sega_id',all.x=T) #ad mean PD data

### - association of cortisol with BMI ####
df_agg$BMI<-with(df_agg,Gewicht/((Groeße/100)*(Groeße/100)))
hist(df_agg$BMI)


lm<-lm(scale(saliva_cortisol_avg)~scale(BMI),
        df_agg[df_agg$saliva_when=='before' & df_agg$manipulation=='before' & df_agg$oddball=='oddball' & df_agg$reverse=='forward',])
summary(lm)
###--> higher BMI is associated with higher cortisol before assessments

lm<-lm(scale(saliva_cortisol_change)~scale(BMI),
        df_agg[!duplicated(df_agg$id),])
summary(lm)
###--> higher BMI also associated with stronger cortisol dropoff

lm<-lm(scale(hair_cortisol)~scale(BMI),
       df_agg[!duplicated(df_agg$id),])
summary(lm)
###--> but not associated with hair cortisol

lm<-lm(scale(mean_pd)~scale(BMI),
       df_agg[!duplicated(df_agg$id),])
summary(lm)
###--> but not associated with baseline pupil size


### - analysis Internalizing and Externalizing with core measures ####

### explore in df_trial
lmm<-lmer(scale(rpd)~oddball*(SDQ_int_s+SDQ_ext_s)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$reverse=='forward',])
summary(lmm)
anova(lmm)
r2_nakagawa(lmm)

fixef(lmm)
emtrends(lmm,~oddball,var="SDQ_ext_s")
###--> contrary to hypothesis - SDQ ext self reports are associated with higher pupillary responses to oddballs

lmm<-lmer(scale(rpd_low)~oddball*(SDQ_int_s+SDQ_ext_s)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$reverse=='forward',])
anova(lmm)
###--> no association

lmm<-lmer(scale(rpd)~oddball*(CBCL_T_INT+CBCL_T_EXT+CBCL_T_GES)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$reverse=='forward',])
summary(lmm)
anova(lmm)
###--> no association

lmm<-lmer(scale(rpd_low)~oddball*(CBCL_T_INT+CBCL_T_EXT+CBCL_T_GES)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$reverse=='forward',])
summary(lmm)
anova(lmm)
###--> no association

df_trial$YSR_T_INT_z<-scale(df_trial$YSR_T_INT)
df_trial$YSR_T_EXT_z<-scale(df_trial$YSR_T_EXT)
df_trial$YSR_T_GES_z<-scale(df_trial$YSR_T_GES)
lmm<-lmer(scale(rpd)~oddball*(YSR_T_INT_z+YSR_T_EXT_z+YSR_T_GES_z)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$reverse=='forward',])
summary(lmm)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm,~oddball,var="YSR_T_INT_z")
emtrends(lmm,~oddball,var="YSR_T_GES_z")
## higher self-rated Internalizing is associated with lower pupillary response to oddballs
## higher self-rated p-factor assocaited with higher pupilalry response to oddballs
        
    ###--> explore for an per-participant aggregated measure
    df_agg_firstblock<-aggregate(df_trial[df_trial$manipulation=='before' & df_trial$reverse=='forward',
                            names(df_trial) %in% c('rpd_auc','rpd_high','rpd_low','rpd')],
                  by=with(df_trial[df_trial$manipulation=='before' & df_trial$reverse=='forward',],
                          list(id,oddball)),
                  mean,na.rm=T)
    names(df_agg_firstblock)<-c('id','oddball','rpd_auc','rpd_high','rpd_low','rpd')
    

    df_agg_firstblock<-merge(df_agg_firstblock,df_dem,by='id')
    
    df_agg_firstblock$YSR_T_INT_z<-scale(df_agg_firstblock$YSR_T_INT)
    df_agg_firstblock$YSR_T_EXT_z<-scale(df_agg_firstblock$YSR_T_EXT)
    df_agg_firstblock$YSR_T_GES_z<-scale(df_agg_firstblock$YSR_T_GES)
    df_agg_firstblock$YSR_T_split<-ifelse(df_agg_firstblock$YSR_T_GES>69,'above clinical cutoff','below clinical cutoff')
    df_agg_firstblock<-df_agg_firstblock[!is.na(df_agg_firstblock$YSR_T_INT),]
        
    lmm<-lmer(scale(rpd)~oddball*(YSR_T_INT_z+YSR_T_EXT_z+YSR_T_GES_z)+(1|id),
              df_agg_firstblock)
    summary(lmm)
    anova(lmm)
    r2_nakagawa(lmm)
    
    emtrends(lmm,~oddball,var='YSR_T_INT_z')    
    emtrends(lmm,~oddball,var='YSR_T_GES_z')    
    
    lmm<-lmer(rpd~oddball*(YSR_T_INT+YSR_T_EXT+YSR_T_GES)+(1|id),
              df_agg_firstblock)
    df_agg_firstblock$predicted_rpd<-predict(lmm)
    ggplot(df_agg_firstblock,aes(YSR_T_GES,predicted_rpd,color=oddball))+
      geom_point()+
      labs(x='p-factor (YSR total t-score)',y='pupillary response (mm)')+
      geom_smooth(method='lm')
    
    ggplot(df_agg_firstblock,aes(YSR_T_GES,predicted_rpd,group=interaction(YSR_T_split,oddball),color=oddball))+
      geom_point()+
      geom_vline(xintercept=65,linetype=2)+
      labs(x='p-factor (YSR total t-score)',y='pupillary response (mm)')+
      geom_smooth(method='lm')+theme_bw()
    
    
    require(wesanderson) #custom color palettes
    df_agg_firstblock$custom_group<-with(df_agg_firstblock,
                                               ifelse(YSR_T_GES<65 & YSR_T_INT<65,'non-clinical',
                                                      ifelse(YSR_T_GES>YSR_T_INT,'high p-factor',
                                                             ifelse(YSR_T_GES<=YSR_T_INT,'high int','else'))))
    
    
    custom_condition_colors <- wes_palette('Cavalcanti1',3,type='discrete') #reverse custom colors to match color coding in other figures
    
    ggplot(df_agg_firstblock[df_agg_firstblock$oddball=='oddball',],
           aes(YSR_T_INT,YSR_T_GES,color=custom_group,
               size=predicted_rpd,alpha=predicted_rpd))+
      geom_point()+
      scale_color_manual(values = custom_condition_colors)+
      labs(x='internalizing (t-score)',y='p-factor (t-score)',
           color = "cluster group",size= 'pupillary resp. (mm)',alpha= 'pupillary resp. (mm)')+ theme_bw()
      
    
    
    
    ggplot(df_agg_firstblock,aes(YSR_T_INT_z,predicted_rpd,color=oddball))+
      geom_point()+
      geom_smooth(method='lm')
    
    
    

lmm<-lmer(scale(rpd_low)~oddball*(YSR_T_INT+YSR_T_EXT+YSR_T_GES)+(1|id),
          df_trial[df_trial$manipulation=='before' & df_trial$reverse=='forward',])
summary(lmm)
anova(lmm)


lm<-lm(scale(mean_pd)~YSR_T_INT+YSR_T_EXT+YSR_T_GES,
          df_trial[!(duplicated(df_trial$id)),])
summary(lm)
##--> no associations

