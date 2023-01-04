# read datasets in HDF5 files

#INSTALL
# # install package from the Biocmanager repository, not CRAN
 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
# 
 BiocManager::install(version = "3.15")
# 
 BiocManager::install("rhdf5")

require(rhdf5)
 
#browseVignettes("rhdf5")

 
 ###EYE TRACKING DATA -->
#read dataset of HDF5 file - location in container file (name) needs to be provided
df <- h5read(file = 'C:/Users/Nico/PowerFolders/project_sega/data/eyetracking/555555_2022-07-08-1010.hdf5',
            name = "data_collection/events/eyetracker/BinocularEyeSampleEvent")

df <- h5read(file = "C:/Users/Nico/PowerFolders/project_sega/data/AuditoryOddball/3_2022-08-16-1614.hdf5",
             name = "data_collection/events/eyetracker/BinocularEyeSampleEvent")

h5closeAll()

names(df)
table(is.na(df))
apply(df,2,function(x){table(is.na(x))})

#infos on data structure
  #experiment_id = 1 (constant)
  #session_ID = 1 (constant)
  #device_id = 0 (constant)
  #event_id = running sequence 
  #type = 52 (constant)
  #--> device_time = timestamp in microsecond format
  #--> logged_time = no timestamp but time format (since start of exp?)
  #--> time = another timestamp in seconds format
  # confidence interval = constant
  #delay = looks like time difference value
  #filter id = 0
  
attach(df)
names(df)

time[2]-time[1]

#raw pupil data
hist(df$left_pupil_measure1)
hist(df$right_pupil_measure1)

#raw gaze data
hist(df$left_gaze_x[df$left_gaze_x<200 & df$left_gaze_x>-200])
hist(df$right_gaze_x[df$right_gaze_x<200 & df$right_gaze_x>-200])
hist(df$left_gaze_y[df$left_gaze_y<200 & df$left_gaze_y>-200])
hist(df$right_gaze_y[df$right_gaze_y<200 & df$right_gaze_y>-200])

require(ggplot2)
ggplot(df[left_gaze_x<200 & left_gaze_x>-200 & left_gaze_y<200 & left_gaze_y>-200,],aes(x=left_gaze_x,y=left_gaze_y))+geom_hex()

#df_trial <- read.csv('C:/Users/Nico/PowerFolders/project_sega/data/trialdata/555555_2022-07-08-1010.csv')
df_trial <- read.csv("C:/Users/Nico/PowerFolders/project_sega/data/AuditoryOddball/3_2022-08-16-1614.csv")

names(df_trial)
table(df_trial$.thisRepN) #repetitions? range = 0-19
table(df_trial$.thisN) # tirals in block?
table(df_trial$.thisTrialN) #trial in sequence?
table(df_trial$phase) #
table(df_trial$block_counter)

hist(df_trial$stimulus_duration[df_trial$stimulus_duration<1],30)
table(df_trial$gaze_offset_duration)
table(df_trial$trial_nodata_duration)

#how to match
  ##--> df_trial$timestamp_exp ~~ df$logged_time

summary(df_trial$timestamp)
format(df_trial$timestamp[[100]],scientific=F)
format(df_trial$timestamp_tracker[[100]],scientific=F)
#--> unix epoch (seconds since 01.01.1070) - psychopy timestamp

summary(df_trial$timestamp_exp) #with df$logged_time or df$time
df_trial$timestamp_exp[100]

##-->

summary(df$device_time)
format(df$device_time[[100]],scientific=F)


summary(df$logged_time)
df$logged_time[100]
summary(df$time)

### ---------- MERGE DATA  ----------- ####
  #not: only returns et data that can be matched to trial data
  #if participants and trial data are in lists - put below into function

#required transformation
start_ts<-df_trial$timestamp_exp
end_ts<-c(df_trial$timestamp_exp[-1],NA)
et_ts<-df$logged_time
split_trial_data<-split(df_trial,seq(nrow(df_trial)))

#merge function
fun_merge_data<-function(ts_1,ts_2,trial_data){
  
  matched_time<-which(et_ts>=ts_1 & et_ts<ts_2)
  selected_et_data<-df[matched_time,] #select et data for a trial
  repeated_trial_data<-data.frame(sapply(trial_data,function(x){rep(x,length(matched_time))},simplify = F))
  merged_data<-data.frame(repeated_trial_data,selected_et_data)

}

df<-mapply(fun_merge_data,ts_1=start_ts,ts_2=end_ts,trial_data=split_trial_data,SIMPLIFY=F)
df<-df[sapply(df,nrow)!=0] #remove frames without entries
#test<-plyr::rbind.fill(test)
df<-dplyr::bind_rows(df) #faster than rbind.fill

### create ts_trial variable

table(df$phase)
df_oddball<-df[df$phase %in% c('oddball_block','oddball_block_rev'),]
table(df_oddball$trial)

df_oddball$ts_trial<-df_oddball$logged_time-df_oddball$timestamp_exp
#eye tracker time - trial onset time


##### -- visualization

df_oddball$trial_type<-substr(df_oddball$trial,1,7)

require(ggplot2)
ggplot(df_oddball[df_oddball$ts_trial<2.1,],aes(x=ts_trial,y=left_pupil_measure1,group=trial,color=trial))+
  geom_smooth()+theme_bw()+labs(x='time (in s)',y='pupil size (in mm)')+facet_wrap(~block_counter)

ggplot(df_oddball[df_oddball$ts_trial<2.1,],aes(x=ts_trial,y=left_pupil_measure1,group=trial_type,color=trial_type))+
  geom_smooth()+theme_bw()+labs(x='time (in s)',y='pupil size (in mm)')

#by block
ggplot(df_oddball[df_oddball$ts_trial<2.1 & df_oddball$block_counter==3,],aes(x=ts_trial,y=left_pupil_measure1,group=trial_type,color=trial_type))+
  geom_smooth()+theme_bw()+labs(x='time (in s)',y='pupil size (in mm)')
ggplot(df_oddball[df_oddball$ts_trial<2.1 & df_oddball$block_counter==5,],aes(x=ts_trial,y=left_pupil_measure1,group=trial_type,color=trial_type))+
  geom_smooth()+theme_bw()+labs(x='time (in s)',y='pupil size (in mm)')
ggplot(df_oddball[df_oddball$ts_trial<2.1 & df_oddball$block_counter==10,],aes(x=ts_trial,y=left_pupil_measure1,group=trial_type,color=trial_type))+
  geom_smooth()+theme_bw()+labs(x='time (in s)',y='pupil size (in mm)')
ggplot(df_oddball[df_oddball$ts_trial<2.1 & df_oddball$block_counter==12,],aes(x=ts_trial,y=left_pupil_measure1,group=trial_type,color=trial_type))+
  geom_smooth()+theme_bw()+labs(x='time (in s)',y='pupil size (in mm)')

####--> get this for all participants
  #correct for baseline data

ggplot(df_oddball[df_oddball$ts_trial<2.1,],aes(x=.thisN,y=left_pupil_measure1))+
  stat_summary(geom='errorbar')+facet_wrap(~block_counter)



hist(df_oddball$.thisN)
