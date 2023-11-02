# Purpose
Analysis of the Auditory Oddball Task conducted for the SEGA Study (2022–2025): Pupil diameter as a proxy for Locus Coeruleus-Norepinephrine (LC-NE) activity in an ASD–NonASD group comparison.

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

# Code structure
## Required libraries
- rhdf5
- data.table
- zoo
- pbapply
- ggplot2
- lme4
- lmerTest
- emmeans
- hexbin
- gridExtra
- dplyr

## Paths
Code is independent of OS although paths have to be adjusted.
## Data import and reshaping
## Data preprocessing
1. invalid pupil sizes (< 2 mm or > 8 mm)
2. dilation speed outliers (median dilation speed = 3 * MAD)
3. blink correction (75 ms > Na <= 250 ms)
4. scatter plot smoothing
5. two-pass approach 
6. interpolation
7. average across both eyes

## Data analysis
### Baseline phase
### Oddball phase
### Manipulation phase
## Visualization
## EEG analysis