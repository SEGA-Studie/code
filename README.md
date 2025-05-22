# Purpose
Analysis of the Auditory Oddball Task conducted for the SEGA Study (2022–2025): Pupil diameter as a proxy for Locus Coeruleus-Norepinephrine (LC-NE) activity and EEG ERPs (MMN, P3a) as measures of snesory processing on cortical level in an ASD-MHC-CON group comparison.

# Requirements
Installing Biomanager: https://www.bioconductor.org/install/  
Installing libraries (see below)
R version 4.3.0

# The Eye Tracker
model: Tobii Pro Spectrum  
sampling rate: 300 Hz -> 300 samples per sec. -> one sample every 3,3 ms
# Input data
## Eye tracking data
Eye tracking data is stored in a .hdf5 format. 
Event types suported by the eye tracker:  
https://www.psychopy.org/api/iohub/device/eyetracker_interface/Tobii_Implementation_Notes.html
### Variables:
*experiment_id* = 1 (constant)  
*session_ID* = 1 (constant)  
*device_id* = 0 (constant)  
*event_id* = running number  
*type* = 52 (constant)  
*device_time* = timestamp in millisecond format  
*logged_time* = no timestamp but time format, since start of experiment, used for matching with trial data  
*timestamp_tracker*: unix epoch, psychopy timestamp: seconds since 01.01.1070  
*time* = another timestamp in seconds format  
*confidence interval* = constant  
*delay* = looks like time difference value  
*filter* id = 0

## Trial data
trialdata is stored in a .csv format.  
*timestamp_exp* (trial data) + *logged_time* (ET data) represent time since start of experiment, used for matching with logged_time variabale from eye tracker.

## EEG data
Preprocessed EEG data is imported as .txt files.

# Required libraries
- rhdf5
- data.table
- zoo
- pbapply
- lme4
- lmerTest
- emmeans
- hexbin
- gridExtra
- ggpubr
- tidyverse

# Paths
Code is independent of OS although paths have to be adjusted.

# Code structure
## NOTES ON MISSING DATA
- List of all data not included in the analysis for the Auditory Oddball-Paper

## SETUP
- sessionInfo
- loading libraries
- paths

## DATA IMPORT AND RESHAPING
- Load eye-tracking data from HDF5 files and corresponding trial data from CSV files. Extract relevant variables and store them in dedicated lists.
- Check for .hdf5 + .csv file matching
- Introduce a new variable: *grip_strength*.
- Organize trial data into a structured, trial-wise DataFrame (*df_trial*).

## MERGE PUPIL AND TRIAL DATA
- Merge pupil and trial data using the user-defined function *fun_merge_all_ids*, and store the result in a new list of DataFrames (*df_list*).

## NEW VARIABLES, REMOVE 3 STANDARDS + RESHAPE INTO LIST OF TRIALS
- Add two new variables: *trial_index* and _trial_number
- Remove the first 3 trials (standards) of each block which were were included by design to establish standards in the oddball paradigm
- Reshape data into *list_split_trial*, a a list of all trials of all subjects. Each list element contains a data frame of trial data.
- Add a new variable, *ts_trial*, representing the time stamp for each pupil data point within a trial. 

## DATA PREPROCESSING
- *fun_blink_cor* is a UDF that defines blinks between 75–250 ms of NAs and replaces data during blinks + 8 seconds around blinks with NA
- *func_pd_preprocess* is a UDF that detects + removes pupil diameter outliers (<2 mm or > 8 mm), calculates + removes dilation speed outliers (3 x median change values), performs robust scatter plot smoothing, interpolation and averages pupil diameter of both eyes (right, left) to one data point, *pd*.

## DEFINE BPS, RESHAPING INTO DATA FRAME + CALCULATE BASELINES 
- Define Baseline Pupil Size (BPS) as the average pupil diameter during the first 250 ms of each trial. BPS serves as a primary outcome variable in subsequent analyses.
- Reshape *list_split_trial* into a large DataFrame (*df*)
- Pupil diameter of each of the 2 baseline phases before (block counters 3, 5) + 2 baseline phases after the manipulation (block counters 10, 12) are averaged. Visualize the results using a boxplot to assess potential habituation effects across task blocks.

## ADD NEW VARIABLES IN TASK BLOCK STRUCTURE + RESHAPE INTO LIST OF TRIALS
- Split the DataFrame *df* into a list of DataFrames (*list_split_block*), grouped by block_counter.
- Add 2 new variables: *trial_index_in_block* + *trial_number_in_block*
- Recombine *list_split_block* into a single DataFrame (*df)
- Split *df* into *list_split_trial*, a list where each element corresponds to an individual trial.
- Remove unused variables from the global environment to optimize memory usage.

## AREA UNDER THE CURVE (AUC)
- AUC as another approach to calculate the pupil response

## DEFINE SEPR + ADD 2 NEW VARIABLES FOR TRIAL INDEXING
- Define Stimulus-Evoked Pupillary Response (*SEPR*) as the mean pupil diameter within the 500–1500 ms window of each trial.
- Compute Baseline Pupil Size (*BPS*) for the newly created DataFrame *et_df_trial*.
- Add two indexing variables (*trial_number* + *trial_number_in_block*) to each element in list_split_trial for trial identification.

## COMPLEMENT + SAVE *DF_TRIAL*
- Add variables to index specific trials in the data frame *df_trial*, AUC variables, SEPR, baseline means, block baseline corrected SEPR (*rpd_block*), variables for manipulation (before, after), block (forward, reverse), pitch, trial (standard, oddball) + trial-corrected pupil reponse
- save (read) *df_trial* as .Rds file for later analysis or reproducibility

## VISUALIZATION
- some plots to visualize current state of data

## READ SAMPLE DESCRIPTIVES (group, z_grip_strength + IQ)
- *exp_group.csv* is a manually created file containing group assignment, date of data collection, grip strength measurements, and IQ scores.
- The *group* variable was manually added to account for a third group (MHC) introduced during data collection, which was not included in the original task implementation that only distinguished between ASD and TD.
- *z_grip_strength* represents grip strength values standardized against normative data provided in the JAMAR hand dynamometer manual.
- *verbal/non-verbal IQ* are computed as the mean of two respective subtests, respectively.


## NEW: MERGE PUPIL AND EEG data ON SINGLE-TRIAL + AGGREGATED LEVEL
- The preprocessing and merging procedures for single-trial and aggregated EEG data follow analogous steps to ensure consistency across data levels.
- EEG data were separately preprocessed and imported here as .txt files
- read EEG files (MMN, P3a, P3b)
- remove 3 standards at the beginning of each task block analogously to pupil data
- Manual renaming of selected EEG files due to incorrect file naming at recording
- merge pupil + EEG data (trial level: *ET_ERP_TRIAL*; aggregated level: *ET_ERP_subject)*
- Qquestionnaire data (SCQ, SRS, CBCL, YSR, SP2) were added
- dependent variables were z-standardized for subsequent analyses.
- Final datasets (*ET_ERP_TRIAL* and *ET_ERP_subject*) were saved for use in the analysis script.

## SAVE DATA FRAME "DF" AS RDS FILE
- The data frame *df* is saved as .Rds to read it in the analysis script.

## MMN DIFFERENCE WACE
- Read as MMN difference wave preprocessed EEG data
- add questionnaire data: SCQ, SRS, CBCL, YSR, SP2
- z-standardization of dependent variables for subsequent analysis
- save *MMN_diff_df* to read in analysis script