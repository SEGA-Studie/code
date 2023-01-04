## SCRIPT HEADER --------------------------- 
##
## Script purpose: READ ANALYZE TEST DATA - SEGA PROJECT - AUDITORY ODDBALL
##
##
## Author: Nico Bast
##
## Date Created: `r paste(Sys.Date())`
##
## Copyright (c) Nico Bast, `r paste(format(Sys.Date(), "%Y"))`
## Email: nico.bast@kgu.de
##
## NOTES --------------------------- 

        #relevant data frames
        # df == unaggregated data
        # df_trial == trial aggregated data


        # - INFO ON HDF5 file <-- ET raw data#
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

        # - INFO ON TRIAL DATA: -#
        #timestamp_tracker --> unix epoch (seconds since 01.01.1070) - psychopy timestamp
        #timestamp_exp (trial data) and logged_time (ET data) represent time since start of experiment
          ###--> is used for matching

## ---------------------------

sessionInfo()

# SETUP ####

#- options
  #options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

# #install package rhdf5 from the Biocmanager repository, not CRAN
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# # 
#  BiocManager::install(version = "3.15")
#  BiocManager::install("rhdf5")


 #- required packages

require(rhdf5) #read eye-tracking raw data - from psychopy recordings
require(data.table) #fread uses parallelization and thus much faster than read.csv
require(zoo) #na.approx
require(pbapply) # progress bar for apply functions
require(ggplot2) #visualization

#analysis
require(lme4)
require(lmerTest)
require(emmeans)

#install.packages('plyr')
#install.packages('dplyr')
#install.packages('psych')
#install.packages('gridExtra')
#install.packages('hexbin')


 
 #PATHS###
 
 #check for OS --> define home path (script independent of OS)
 ifelse(Sys.info()['sysname']=='Linux',
        home_path<-'~',
        home_path<-'C:/Users/Nico')
 
 
 #project path
 project_path<-'/PowerFolders/project_sega'
 data_path<-'/PowerFolders/project_sega/data/AuditoryOddball' 
 data_path_eeg<-'/PowerFolders/project_sega/data/AuditoryOddball_EEG'
 
 datapath<-paste0(home_path,data_path)
 datapath_eeg<-paste0(home_path,data_path_eeg)
 
 
 #FUNCTIONS###
 
 func_pd_preprocess<-function(x){
   
   #define variables
   Left_Diameter<-x$left_pupil_measure1
   Right_Diameter<-x$right_pupil_measure1
   #RemoteTime<-x$timestamp
   
   #constant for MAD caluclation
   constant<-3 ##--> if change speed is higher than constant * median change --> values are excluded
   #constant<-3 #default value
   
   # STEP 1 - exclude invalid data ####
   pl <- ifelse((Left_Diameter<2|Left_Diameter>8), NA, Left_Diameter)
   pr <- ifelse((Right_Diameter<2|Right_Diameter>8), NA, Right_Diameter)
   #table(is.na(pl))
   #table(is.na(pr))
   
   # STEP 2 - filtering ####
   ## A) normalized dilation speed, take into account time jumps with Remotetimestamps: ####
   #maximum change in pd compared to last and next pd measurement
   #Left
   # pl.speed1<-diff(pl)/diff(RemoteTime) #compared to last
   # pl.speed2<-diff(rev(pl))/diff(rev(RemoteTime)) #compared to next
   # pl.speed1<-c(NA,pl.speed1)
   # pl.speed2<-c(rev(pl.speed2),NA)
   # pl.speed<-pmax(pl.speed1,pl.speed2,na.rm=T)
   # rm(pl.speed1,pl.speed2)
   # #Right
   # pr.speed1<-diff(pr)/diff(RemoteTime)
   # pr.speed2<-diff(rev(pr))/diff(rev(RemoteTime))
   # pr.speed1<-c(NA,pr.speed1)
   # pr.speed2<-c(rev(pr.speed2),NA)
   # pr.speed<-pmax(pr.speed1,pr.speed2,na.rm=T)
   # rm(pr.speed1,pr.speed2)
   # #median absolute deviation -SPEED
   # #constant<-3
   # pl.speed.med<-median(pl.speed,na.rm=T)
   # pl.mad<-median(abs(pl.speed-pl.speed.med),na.rm = T)
   # pl.treshold.speed<-pl.speed.med+constant*pl.mad #treshold.speed units are mm/microsecond
   # #plot(abs(pl.speed))+abline(h=pl.treshold.speed)
   # pr.speed.med<-median(pr.speed,na.rm=T)
   # pr.mad<-median(abs(pr.speed-pr.speed.med),na.rm = T)
   # pr.treshold.speed<-pr.speed.med+constant*pr.mad #treshold.speed units are mm/microsecond
   # #plot(abs(pr.speed))+abline(h=pr.treshold.speed)
   # #correct pupil dilation for speed outliers
   # pl<-ifelse(abs(pl.speed)>pl.treshold.speed,NA,pl)
   # pr<-ifelse(abs(pr.speed)>pr.treshold.speed,NA,pr)
   # 
   # ## B) delete data around blinks - not applied ####
   # #gaps=missing data sections > 75ms; Leonie: also <=250ms, otherwise not likely to be a blink
   # #to be excluded: samples within 50 ms of gaps -> +-25 (8 data points) oder 50?
   # pl<-fun_blink_cor(pl)
   # pr<-fun_blink_cor(pr)
   # 
   # ## C) normalized dilation size - median absolute deviation -SIZE ####
   # #applies a two pass approach
   # #first pass: exclude deviation from trend line derived from all samples
   # #second pass: exclude deviation from trend line derived from samples passing first pass
   # #-_> reintroduction of sample that might have been falsely excluded due to outliers
   # #estimate smooth size based on sampling rate
   # smooth.length<-150 #measured in ms
   # #take sampling rate into account (300 vs. 120):
   # #smooth.size<-round(smooth.length/mean(diff(RemoteTime)/1000)) #timestamp resolution in microseconds
   # smooth.size<-round(smooth.length/median(diff(RemoteTime),na.rm=T)) #timestamp resolution in milliseconds
   # is.even<-function(x){x%%2==0}
   # smooth.size<-ifelse(is.even(smooth.size)==T,smooth.size+1,smooth.size) #make sure to be odd value (see runmed)
   # #Left
   # pl.smooth<-na.approx(pl,na.rm=F,rule=2) #impute missing values with interpolation
   # #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
   # if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
   # pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
   # #Right
   # pr.smooth<-na.approx(pr,na.rm=F,rule=2) #impute missing values with interpolation
   # #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
   # if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
   # pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
   # #correct pupil dilation for size outliers - FIRST pass
   # pl.pass1<-ifelse((pl>pl.smooth+constant*pl.mad)|(pl<pl.smooth-constant*pl.mad),NA,pl)
   # pr.pass1<-ifelse((pr>pr.smooth+constant*pr.mad)|(pr<pr.smooth-constant*pr.mad),NA,pr)
   # #Left
   # pl.smooth<-na.approx(pl.pass1,na.rm=F,rule=2) #impute missing values with interpolation
   # #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
   # if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
   # pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
   # #Right
   # pr.smooth<-na.approx(pr.pass1,na.rm=F,rule=2) #impute missing values with interpolation
   # #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
   # if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
   # pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
   # #correct pupil dilation for size outliers - SECOND pass
   # pl.pass2<-ifelse((pl>pl.smooth+constant*pl.mad)|(pl<pl.smooth-constant*pl.mad),NA,pl)
   # pr.pass2<-ifelse((pr>pr.smooth+constant*pr.mad)|(pr<pr.smooth-constant*pr.mad),NA,pr)
   # pl<-pl.pass2
   # pr<-pr.pass2
   # 
   # ## D) sparsity filter - not applied ####
   # # STEP 3 - processing valid samples  ####
   # #take offset between left and right into account
   pd.offset<-pl-pr
   pd.offset<-na.approx(pd.offset,rule=2)
   #mean pupil dilation across both eyes
   pl <- ifelse(is.na(pl)==FALSE, pl, pr+pd.offset)
   pr <- ifelse(is.na(pr)==FALSE, pr, pl-pd.offset)
   
   # #interpolation of NA (for <=300ms)
   # pl<-na.approx(pl, na.rm=F, maxgap=90, rule=2)
   # pr<-na.approx(pr, na.rm=F, maxgap=90, rule=2)
   
   pd <- (pl+pr)/2
   # end of function --> return ####
   #detach(x)
   
   x[,'pd']<-pd
   return(x)
 }
 
## ----------- READ ET AND TRIAL DATA AND MERGE ---------####
# READ ET DATA ####
 
 #read data from datapath and store in according objects
 data.files<-list.files(path=datapath,full.names=T)
 data_files_et<-data.files[grepl('.hdf5',data.files)]
 
 #remove files without data
 data_files_et<-data_files_et[!(data_files_et %in% c("C:/Users/Nico/PowerFolders/project_sega/data/AuditoryOddball/t01_2022-08-15-1125.hdf5"))]

 #read all files to list (hdf5 is a container format and et data is stored in BinocularEyeSampleEvent folder within the container) 
 list_et_data<-list(0)
 for(i in 1:length(data_files_et)){
   list_et_data[[i]]<-h5read(file = data_files_et[i], name = "data_collection/events/eyetracker/BinocularEyeSampleEvent")
   print(paste0('read ET data file: ',i))
 }
 
  #close open handles
  h5closeAll()

  #add names to list
  #id.names<-substr(data_files_et,nchar(datapath)+2,nchar(data_files_et)-21) 
  id.names<-substr(data_files_et,nchar(datapath)+2,nchar(data_files_et)-10) 
  id.names
  
  names(list_et_data)<-id.names
  
# DROP NON-INFORMATION (constant) ET VARIABLES    ####
  
  constant_variables<-c('experiment_id','session_id','device_id','type','confidence_interval','filter_id',
                        'left_gaze_z','left_angle_x','left_angle_y','left_raw_x','left_raw_y',
                        'left_pupil_measure1_type','left_pupil_measure2_type','left_ppd_x','left_ppd_y',
                        'left_velocity_x','left_velocity_y','left_velocity_xy','left_pupil_measure2',
                        'right_gaze_z','right_angle_x','right_angle_y','right_raw_x','right_raw_y',
                        'right_pupil_measure1_type','right_pupil_measure2_type','right_ppd_x','right_ppd_y',
                        'right_velocity_x','right_velocity_y','right_velocity_xy','right_pupil_measure2')
  
  list_et_data<-lapply(list_et_data,function(x){x[!(names(x) %in% constant_variables)]})
  
# INITIAL ET DATA QUALITY ####
  
  # #missing gaze data (~10%)
  # sapply(list_et_data,function(x){round(table(is.na(x$left_gaze_x))[2]/sum(table(is.na(x$left_gaze_x))),2)})
  # #missing pupil data (~10%) 
  # sapply(list_et_data,function(x){round(table(is.na(x$left_pupil_measure1))[2]/sum(table(is.na(x$left_pupil_measure1))),2)})
  # 
  # #time gap --> should be 0.0033
  # sapply(list_et_data,function(x){x$time[2]-x$time[1]})
  # sapply(list_et_data,function(x){hist(diff(x$time),50)})
  # 
  # #normally distributed pupil size  
  # sapply(list_et_data,function(x){hist(x$left_pupil_measure1,50)})
  # sapply(list_et_data,function(x){hist(x$right_pupil_measure1,50)})
  # 
  
# READ TRIAL DATA ####
  
  #read data from datapath and store in according objects
  data.files<-list.files(path=datapath,full.names=T)
  data_files_trial<-data.files[grepl('.csv',data.files)]
  
  #remove files without data
  data_files_trial<-data_files_trial[!(data_files_trial %in% c("C:/Users/Nico/PowerFolders/project_sega/data/AuditoryOddball/t01_2022-08-15-1125.csv"))]
  
  #read trial data into list
  list_trial_data<-list(0)
  for(i in 1:length(data_files_trial)){
    list_trial_data[[i]]<-fread(data_files_trial[i])
    print(paste0('read TRIAL data file: ',i))
  }
  
  #convert to data.frame 
  list_trial_data<-lapply(list_trial_data,data.frame)
  
  #remove empty variable
  list_trial_data<-lapply(list_trial_data,function(x){x[!(names(x) %in% c('V26'))]})
  
  #ID as character
  list_trial_data<-lapply(list_trial_data,function(x){x$id<-as.character(x$id);return(x)})
  
  #add names to list
  #id.names<-substr(data_files_trial,nchar(datapath)+2,nchar(data_files_trial)-20) 
  id.names<-substr(data_files_trial,nchar(datapath)+2,nchar(data_files_trial)-9) 
  names(list_trial_data)<-id.names
  
# INITIAL TRIAL DATA QUALITY ####  
  
  # #INFO on data structure
  # sapply(list_trial_data,function(x){sum(table(x$.thisTrialN)[-1])}) # total number of trials (without practice trials [-1])
  # sapply(list_trial_data,function(x){table(x$.thisTrialN)}) 
  #   
  # sapply(list_trial_data,function(x){table(x$.thisRepN)}) #repetitions of trial sequence (4 standard + 1)
  # sapply(list_trial_data,function(x){table(x$.thisN)}) # trials in a sequence
  # sapply(list_trial_data,function(x){table(x$phase)}) 
  # sapply(list_trial_data,function(x){table(x$block_counter)}) 

# CREATE df-trial ####
df_trial<-plyr::rbind.fill(list_trial_data)   
  
hist(df_trial$stimulus_duration[df_trial$stimulus_duration<1],30) #stimulus duration
hist(df_trial$ISI_duration[df_trial$ISI_duration<3],30) #sometimes ISI is too long
table(df_trial$gaze_offset_duration)
table(df_trial$trial_nodata_duration)

hist(df_trial$stimulus_duration[df_trial$stimulus_duration<0.2],col='grey',breaks=30,main='auditory oddball - actual stimulus duration',xlab='time (s)')
psych::describe(df_trial$stimulus_duration[df_trial$stimulus_duration<0.2])

## ---> MERGE ET AND TRIAL DATA ####
  
    #how to match - as there are different timestamp formats
      ##--> df_trial$timestamp_exp ~~ df$logged_time
    
        # summary(df_trial$timestamp)
        # format(df_trial$timestamp[[100]],scientific=F)
        # format(df_trial$timestamp_tracker[[100]],scientific=F)
        # #--> unix epoch (seconds since 01.01.1070) - psychopy timestamp
        # 
        # summary(df_trial$timestamp_exp) #with df$logged_time or df$time
        # hist(df_trial$timestamp_exp)
        
        ##-->
        
## --> MERGE DATA  ##
  #not: only returns et data that can be matched to trial data

#- merge function (per id per trial)
fun_merge_all_ids<-function(et_data,trial_data){

  ##required transformation
  start_ts<-trial_data$timestamp_exp
  end_ts<-c(trial_data$timestamp_exp[-1],NA)
  et_ts<-et_data$logged_time
  split_trial_data<-split(trial_data,seq(nrow(trial_data)))
  
  #merge function per ID and TRIAL
  fun_merge_data<-function(ts_1,ts_2,trial_data_splitted){
    
    matched_time<-which(et_ts>=ts_1 & et_ts<ts_2)
    selected_et_data<-et_data[matched_time,] #select et data for a trial
    repeated_trial_data<-data.frame(sapply(trial_data_splitted,function(x){rep(x,length(matched_time))},simplify = F))
    merged_data<-data.frame(repeated_trial_data,selected_et_data)
    
  }
  
  #merge data of one id and convert to data frame
  df_one_id<-mapply(fun_merge_data,ts_1=start_ts,ts_2=end_ts,trial_data_splitted=split_trial_data,SIMPLIFY=F)
  #df_one_id<-df_one_id[sapply(df_one_id,nrow)!=0] #remove frames without entries
  df_one_id<-dplyr::bind_rows(df_one_id) #faster than rbind.fill
  
}
    
df_list<-pbmapply(fun_merge_all_ids,et_data=list_et_data,trial_data=list_trial_data,SIMPLIFY=F)
#may take a while

## ------------ DATA PREPROCESSING (Creating variables, PD preprocessing) --------#####
##CREATE a trial_index, trial_number, and trial in oddball phase variable####
#INFO: Trial is defined --> .thisTrialN (trials within sequence) x .thisRepn (repetition of sequences in block) x block_counter (blocks)
df_list<-sapply(df_list,function(x){
  #all trials
  x$trial_index<-with(x,droplevels(interaction(.thisTrialN,.thisRepN,block_counter))) #name to indentify individual trials
  x$trial_number<-with(x,rep(seq_along(table(trial_index)),times=table(trial_index))) #all trials including empty trials
  #x$trial_number<-with(x,rep(seq_along(table(trial_number)),times=table(trial_number))) #rerun remove empty trials from trial number
  return(x)})

# #testing
# test<-df_list[[1]]
# trial_index<-with(test,droplevels(interaction(.thisTrialN,.thisRepN,block_counter)))
# trial_number<-rep(seq_along(table(trial_index)),times=table(trial_index))
# table(trial_index)
# 
# list_split_trial<-lapply(test,function(x){split(x,x$trial_number)}) #split by individual trials
# list_split_trial<-unlist(list_split_trial,recursive=F) #remove top level list - every trial is now a list element
# testdf<-list_split_trial[[4]]
# sapply(list_split_trial[[1]],nrow)


##SPLIT data per trial####
list_split_trial<-lapply(df_list,function(x){split(x,x$trial_number)})
list_split_trial<-unlist(list_split_trial,recursive=F) #remove top level list - every trial is now a list element

##CREATE ts_trial variable####
list_split_trial<-lapply(list_split_trial,function(x){
    x$ts_trial<-x$logged_time-x$timestamp_exp
    return(x)
  })


## PD preprocessing ####

#define functions
fun_blink_cor <- function(signal,lower_threshold=23,upper_threshold=75,samples_before=8,samples_after=8) {
  #change NA to 999 for rle()-function
  findna <- ifelse(is.na(signal),999,signal)
  #find blinks:
  #output of rle(): how many times values (NA) are repeated
  repets <- rle(findna)
  #stretch to length of PD vector for indexing
  repets <- rep(repets[["lengths"]], times=repets[["lengths"]])
  #difference between two timestamps~3.33ms -> 75/3.333=22.5 -> wenn 23 Reihen PD=NA, dann blink gap
  #if more than 150ms (45 rows) of NA, missing data due to blink unlikely
  #dummy coding of variables (1=at least 23 consecutive repetitions, 0=less than 23 repetitions)
  repets <- ifelse(repets>=lower_threshold & repets<=upper_threshold, 1, 0)
  #exclude cases where other values than NA (999) are repeated >=23 times by changing dummy value to 0:
  repets[findna!=999 & repets==1] <- 0
  #gives out where changes from 0 to 1 (no NA at least 23xNA) or 1 to 0 (23x NA to no NA) appear
  changes <- c(diff(repets),0)
  #define start (interval before blink/missing data)
  changes.start<-which(changes==1) #where NA-sequence starts
  #gives out row numbers of NA (blink) and previous 8 frames
  start.seq<-unlist(lapply(changes.start, function(x) {seq(max(x-(samples_before-1),1), x)}))
  repets[start.seq]<-1
  #define end (interval after blink/missing data)
  changes.end<-which(changes==-1)+1 #where NA.sequence ends
  #gives out row numbers of NA (blink) and subsequent 8 frames
  end.seq<-unlist(lapply(changes.end, function(x) {seq(x, min(x+(samples_before-1),length(repets)))}))
  repets[end.seq]<-1
  #replace PD data in blink interval (start to end) with NA
  signal[repets==1]<-NA
  return(signal)
}

func_pd_preprocess<-function(x){
  
  #define variables
  Left_Diameter<-x$left_pupil_measure1
  Right_Diameter<-x$right_pupil_measure1
  RemoteTime<-x$ts_trial*1000 #*1000 to retain ms format of timestamp (was in s format)
  
  #constant for MAD caluclation
  constant<-3 ##--> if change speed is higher than constant * median change --> values are excluded
  
  # STEP 1 - exclude invalid data ####
  pl <- ifelse((Left_Diameter<2|Left_Diameter>8), NA, Left_Diameter)
  pr <- ifelse((Right_Diameter<2|Right_Diameter>8), NA, Right_Diameter)
  
  # STEP 2 - filtering ####
  ## A) normalized dilation speed, take into account time jumps with Remotetimestamps: ####
  #maximum change in pd compared to last and next pd measurement
  #Left
  pl.speed1<-diff(pl)/diff(RemoteTime) #compared to last
  pl.speed2<-diff(rev(pl))/diff(rev(RemoteTime)) #compared to next
  pl.speed1<-c(NA,pl.speed1)
  pl.speed2<-c(rev(pl.speed2),NA)
  pl.speed<-pmax(pl.speed1,pl.speed2,na.rm=T)
  rm(pl.speed1,pl.speed2)
  #Right
  pr.speed1<-diff(pr)/diff(RemoteTime)
  pr.speed2<-diff(rev(pr))/diff(rev(RemoteTime))
  pr.speed1<-c(NA,pr.speed1)
  pr.speed2<-c(rev(pr.speed2),NA)
  pr.speed<-pmax(pr.speed1,pr.speed2,na.rm=T)
  rm(pr.speed1,pr.speed2)
  #median absolute deviation -SPEED
  #constant<-3
  pl.speed.med<-median(pl.speed,na.rm=T)
  pl.mad<-median(abs(pl.speed-pl.speed.med),na.rm = T)
  pl.treshold.speed<-pl.speed.med+constant*pl.mad #treshold.speed units are mm/microsecond
  #plot(abs(pl.speed))+abline(h=pl.treshold.speed)
  pr.speed.med<-median(pr.speed,na.rm=T)
  pr.mad<-median(abs(pr.speed-pr.speed.med),na.rm = T)
  pr.treshold.speed<-pr.speed.med+constant*pr.mad #treshold.speed units are mm/microsecond
  #plot(abs(pr.speed))+abline(h=pr.treshold.speed)
  #correct pupil dilation for speed outliers
  pl<-ifelse(abs(pl.speed)>pl.treshold.speed,NA,pl)
  pr<-ifelse(abs(pr.speed)>pr.treshold.speed,NA,pr)
  
  ## B) delete data around blinks - not applied ####
  #gaps=missing data sections > 75ms; Leonie: also <=250ms, otherwise not likely to be a blink
  #to be excluded: samples within 50 ms of gaps -> +-25 (8 data points) oder 50?
  pl<-fun_blink_cor(pl)
  pr<-fun_blink_cor(pr)
  
  ## C) normalized dilation size - median absolute deviation -SIZE ####
  #applies a two pass approach
  #first pass: exclude deviation from trend line derived from all samples
  #second pass: exclude deviation from trend line derived from samples passing first pass
  #-_> reintroduction of sample that might have been falsely excluded due to outliers
  #estimate smooth size based on sampling rate
  smooth.length<-150 #measured in ms
  #take sampling rate into account (300 vs. 120):
  #smooth.size<-round(smooth.length/mean(diff(RemoteTime)/1000)) #timestamp resolution in microseconds
  smooth.size<-round(smooth.length/median(diff(RemoteTime),na.rm=T)) #timestamp resolution in milliseconds
  is.even<-function(x){x%%2==0}
  smooth.size<-ifelse(is.even(smooth.size)==T,smooth.size+1,smooth.size) #make sure to be odd value (see runmed)
  #Left
  pl.smooth<-na.approx(pl,na.rm=F,rule=2) #impute missing values with interpolation
  #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
  #Right
  pr.smooth<-na.approx(pr,na.rm=F,rule=2) #impute missing values with interpolation
  #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
  #correct pupil dilation for size outliers - FIRST pass
  pl.pass1<-ifelse((pl>pl.smooth+constant*pl.mad)|(pl<pl.smooth-constant*pl.mad),NA,pl)
  pr.pass1<-ifelse((pr>pr.smooth+constant*pr.mad)|(pr<pr.smooth-constant*pr.mad),NA,pr)
  #Left
  pl.smooth<-na.approx(pl.pass1,na.rm=F,rule=2) #impute missing values with interpolation
  #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
  #Right
  pr.smooth<-na.approx(pr.pass1,na.rm=F,rule=2) #impute missing values with interpolation
  #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
  #correct pupil dilation for size outliers - SECOND pass
  pl.pass2<-ifelse((pl>pl.smooth+constant*pl.mad)|(pl<pl.smooth-constant*pl.mad),NA,pl)
  pr.pass2<-ifelse((pr>pr.smooth+constant*pr.mad)|(pr<pr.smooth-constant*pr.mad),NA,pr)
  pl<-pl.pass2
  pr<-pr.pass2
  
  ## D) sparsity filter - not applied ####
  # STEP 3 - processing valid samples  ####
  #take offset between left and right into account
  pd.offset<-pl-pr
  pd.offset<-na.approx(pd.offset,rule=2)
  #mean pupil dilation across both eyes
  pl <- ifelse(is.na(pl)==FALSE, pl, pr+pd.offset)
  pr <- ifelse(is.na(pr)==FALSE, pr, pl-pd.offset)
  
  #interpolation of NA (for <=300ms)
  pl<-na.approx(pl, na.rm=F, maxgap=90, rule=2)
  pr<-na.approx(pr, na.rm=F, maxgap=90, rule=2)
  
  pd <- (pl+pr)/2
  # end of function --> return ####
  #detach(x)
  
  x[,'pd']<-pd
  return(x)
}

list_split_trial<-pblapply(list_split_trial, func_pd_preprocess) #apply with progress bar (pb)

##--> bind list of merged data to data.frame####
#df<-dplyr::bind_rows(df_list) #faster than plyr::rbind.fill
df<-dplyr::bind_rows(list_split_trial) #faster than plyr::rbind.fill

    # what does this variable encode - tracking status?
    #convert raw to numeric
    #test$status<-as.numeric(as.character(test$status))
    #table(as.character(df$status))


##split by block and id
list_split_blocks<-split(df,droplevels(interaction(df$block_counter,df$id)))

list_split_blocks<-lapply(list_split_blocks,function(x){
  #all trials
  x$trial_index_in_block<-with(x,droplevels(interaction(.thisTrialN,.thisRepN))) #name to indentify individual trials within BLOCK
  x$trial_number_in_block<-with(x,rep(seq_along(table(trial_index_in_block)),times=table(trial_index_in_block))) #all trials including empty trials
  return(x)})


#melt to data.frame
df<-dplyr::bind_rows(list_split_blocks) #faster than plyr::rbind.fill
#melt to split by trial (for further processing)
list_split_trial<-split(df,droplevels(interaction(df$id,df$trial_number)))

    
table(df$phase)
table(df$trial)

hist(df$ts_trial[df$ts_trial<5])
hist(df$pd)
#preprocessed pd for every participant
ggplot(df[is.finite(df$pd),],aes(x=pd))+geom_histogram(bins=100)+facet_wrap(~id)

##--> reduce ET data to per trial data (and merge with df_trial) ####

#define pupillary response varaibles
rpd_high<-sapply(list_split_trial,function(x){mean(x$pd[x$ts_trial>0.875 & x$ts_trial<1.125])})
rpd_low<-sapply(list_split_trial,function(x){mean(x$pd[x$ts_trial<0.250])})

#additional variables
trial_number<-sapply(list_split_trial,function(x){unique(x$trial_number)})
trial_number_in_block<-sapply(list_split_trial,function(x){unique(x$trial_number_in_block)})

#merge data
merger_id<-sapply(list_split_trial,function(x){interaction(x$id,x$trial_index)[1]})
df_et_trial<-data.frame(merger_id,rpd_high,rpd_low,trial_number,trial_number_in_block)

df_trial$trial_index<-with(df_trial,interaction(.thisTrialN,.thisRepN,block_counter)) #name to identify individual trials
df_trial$merger_id<-interaction(df_trial$id,df_trial$trial_index)

df_trial<-merge(df_trial,df_et_trial,id='merger_id')

#define pupillary response
df_trial$rpd<-df_trial$rpd_high-df_trial$rpd_low

#define additional variables
df_trial$manipulation<-factor(ifelse(df_trial$block_counter<8,'before',
                                     ifelse(df_trial$block_counter>8,'after','manipulation')),levels=c('before','after'))
df_trial$pitch<-ifelse(df_trial$trial=='oddball',df_trial$oddball_frequency,df_trial$standard_frequency)
df_trial$reverse<-ifelse(df_trial$phase=='oddball_block_rev','reverse','forward')
df_trial$oddball<-as.factor(ifelse(grepl('oddball',df_trial$trial),'oddball','standard'))
df_trial$id_type<-ifelse(df_trial$id %in% c('3','555555','t002','t01'),'test','study')


table(df_trial$block_counter,df_trial$phase)
table(df_trial$trial)
table(df_trial$reverse,df_trial$oddball)

##--> CHECK trial-based data ####

#overall
hist(df_trial$rpd,50)
with(df_trial,by(rpd,trial,mean,na.rm=T))
lmm<-lmer(rpd~oddball*manipulation*reverse+trial_number_in_block+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev') & df_trial$id_type=='study',])
summary(lmm)
contrast(emmeans(lmm,~oddball),'pairwise')

#
ggplot(df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev') & is.finite(df_trial$rpd) & df_trial$id_type=='study',],
       aes(x=trial_number,y=rpd,group=oddball,color=oddball))+geom_smooth()+facet_grid(rows = vars(reverse), cols = vars(manipulation))

ggplot(df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev') & df_trial$trial_number_in_block>=10 & is.finite(df_trial$rpd) & df_trial$id_type=='study',],
       aes(x=trial_number_in_block,y=rpd,group=oddball,color=oddball))+geom_smooth()+facet_grid(rows = vars(reverse), cols = vars(manipulation))+theme_bw()

ggplot(df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev') & is.finite(df_trial$rpd) & df_trial$id_type=='study',],
       aes(x=as.factor(trial_number_in_block),y=rpd))+geom_boxplot()+facet_grid(rows = vars(reverse), cols = vars(manipulation,oddball))+theme_bw()


####### --------- PHASE DATA ANALYSIS ------------####
##--> DATA ANALYSIS: BASELINE PHASE ####    

df_baseline<-df[df$phase %in% c('baseline','baseline_calibration'),]

  with(df_baseline[df_baseline$phase=='baseline_calibration',],table(trial))
  #baseline, bloackslide, whiteslide
  with(df_baseline,table(trial,block_counter))

  with(df_baseline,by(timestamp_exp,interaction(baseline_trial_counter,phase,id),mean,na.rm=T))
  #--> first comes baseline_calibration
  
  df_baseline<-df_baseline[is.finite(df_baseline$pd),]
  
        #BASELINE VISUALIZATION
        # #baseline pupil size (start of experiment)
        # require(ggplot2)
        # ggplot(df_baseline[df_baseline$phase=='baseline_calibration',],aes(x=trial,y=pd,fill=id))+geom_boxplot()
        # ggplot(df_baseline[df_baseline$phase=='baseline',],aes(x=block_counter,y=pd,group=block_counter))+geom_boxplot()
        # 
        # ggplot(df_baseline[df_baseline$phase=='baseline_calibration',],aes(x=ts_trial,y=pd,group=trial,color=trial))+geom_smooth()+facet_wrap(~id)
        # 
        # 
        # #baseline pupil size in experiment
        # tiff(file=paste0(home_path,project_path,"/output/figure_audio_effectofmanipulation_testdata.tiff"), # create a file in tiff format in current working directory
        #      width=8, height=4, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure
        # 
        # ggplot(df_baseline,aes(x=as.factor(baseline_trial_counter),y=pd,fill=id))+
        #   geom_boxplot()+facet_wrap(~id)+geom_vline(xintercept=2.5,lty=2)+theme_bw()+
        #   xlab('baseline PD measurements')+ylab('pupil size (mm)')+labs(title='effect of manipulation')
        # 
        # dev.off()
        # 
        # #pd change in baseline calibration
        # tiff(file=paste0(home_path,project_path,"/output/figure_audio_changeofpd_testdata.tiff"), # create a file in tiff format in current working directory
        #      width=8, height=4, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure
        # 
        # ggplot(df_baseline[df_baseline$phase=='baseline_calibration' & df_baseline$ts_trial<6,],aes(x=ts_trial,y=pd,group=trial,color=trial))+
        #   geom_smooth()+facet_wrap(~id)+theme_bw()+xlab('time (s)')+ylab('pupil size (mm)')+labs(title='change of PD')
        # 
        # dev.off()
        # 
        # ggplot(df_baseline[df_baseline$phase=='baseline_calibration' & df_baseline$ts_trial<6,],aes(x=ts_trial,y=pd,group=trial,color=trial))+
        #   geom_smooth()+theme_bw()+xlab('time (s)')+ylab('pupil size (mm)')+labs(title='change of PD')
        # 
        
##--> estimate mean PD per ID --####
    df_meanpd<-with(df_baseline[df_baseline$phase=='baseline_calibration' & df_baseline$trial=='baseline',],by(as.numeric(pd),id,mean))

    df_meanpd<-data.frame(names(df_meanpd),as.numeric(df_meanpd))
    names(df_meanpd)<-c('id','mean_pd')

    
##--> DATA ANALYSIS: ODDBALL PHASE ####    
table(df$phase)
    
##select oddball data
df_oddball<-df[df$phase %in% c('oddball_block','oddball_block_rev'),]
table(df_oddball$trial)

#merge with baseline PD
df_oddball<-merge(df_oddball,df_meanpd,by='id')

#baselinecorrected pd
df_oddball$rpd<-df_oddball$pd-df_oddball$mean_pd
hist(df_oddball$rpd,50)

#correct for non-finite values
df_oddball<-df_oddball[is.finite(df_oddball$rpd),]

#trial type (standard versus oddball)
df_oddball$trial_type<-substr(df_oddball$trial,1,7)

#
table(df_oddball$block_counter,df_oddball$trial)
df_oddball$manipulation<-factor(ifelse(df_oddball$block_counter<8,'before',
                                ifelse(df_oddball$block_counter>8,'after','manipulation')),levels=c('before','after'))

df_oddball$order<-ifelse(grepl('rev',df_oddball$trial),'reverse','normal')
table(df_oddball$order,df_oddball$block_counter)

table(df_oddball$id)
df_oddball$id_type<-ifelse(df_oddball$id %in% c('3','555555','t002','t01'),'test','study')


#define trials in oddball phase
df_oddball$trial_index_oddballphase<-with(df_oddball,interaction(block_counter,.thisRepN,.thisTrialN)) #name to indentify individual trials
df_oddball$trial_number_oddballphase<-with(df_oddball,rep(seq_along(table(trial_index_oddballphase)),times=table(trial_index_oddballphase))) #all trials including empty trials
df_oddball$trial_number_oddballphase<-with(df_oddball,rep(seq_along(table(trial_number_oddballphase)),times=table(trial_number_oddballphase))) #rerun remove empty trials from trial number

table(df_oddball$trial_number)
hist(df_oddball$trial_number_oddballphase)

###pupillary response by trial
ggplot(df_oddball,aes(x=trial_number_oddballphase,y=rpd))+geom_smooth(method='lm')+facet_wrap(~trial+manipulation)

ggplot(df_oddball,aes(rpd, fill=interaction(manipulation)))+geom_density(alpha=0.2)+facet_wrap(~trial)

#define 


##### -- visualization

  tiff(file=paste0(home_path,project_path,"/output/figure_audio_effectofmanipulation_pupilresponse_testdata.tiff"), # create a file in tiff format in current working directory
       width=8, height=4, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure

# ggplot(df_oddball[df_oddball$ts_trial<2 & (df_oddball$id %in% c('5','22','31','44','81')),],
#        aes(x=ts_trial,y=scale(rpd),group=trial_type,color=trial_type))+
#   geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupil response (z)',title='effect of manipulation of pupil response')+facet_wrap(~manipulation+order)

  #restricted to normal trials
ggplot(df_oddball[df_oddball$ts_trial<1.8 & df_oddball$order=='normal',],
       aes(x=ts_trial,y=scale(rpd),group=trial_type,color=trial_type))+
  geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupil response (z)',title='effect of manipulation of pupil response')+facet_wrap(~manipulation)

  dev.off()

  
  # ###--> split for normal versus reverse
  # ggplot(df_oddball[df_oddball$ts_trial<2 & df_oddball$order=='normal',],
  #        aes(x=ts_trial,y=scale(rpd),group=trial_type,color=trial_type))+
  #   geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupil response (z)',title='effect of manipulation of pupil response - normal condition')+facet_wrap(~manipulation)
  # 
  # ggplot(df_oddball[df_oddball$ts_trial<2 & df_oddball$order=='reverse',],
  #        aes(x=ts_trial,y=scale(rpd),group=trial_type,color=trial_type))+
  #   geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupil response (z)',title='effect of manipulation of pupil response - reverse condition')+facet_wrap(~manipulation)
  # 
  
  ## --> NEUROPHYSIOLOGICAL HABITUATION to standard trials within block ####
  ggplot(df_oddball[df_oddball$ts_trial<2 & is.finite(df_oddball$rpd) & df_oddball$trial_type=='standar',],
         aes(x=.thisRepN,y=scale(rpd),group=interaction(manipulation,order),color=interaction(manipulation,order)))+geom_smooth()+theme_bw()
  
  
  tiff(file=paste0(home_path,project_path,"/output/figure_audio_neurophysiological_habituation_standards_Jan2023.tiff"), # create a file in tiff format in current working directory
       width=8, height=4, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure
  
  
  ggplot(df_oddball[df_oddball$ts_trial<2 & is.finite(df_oddball$rpd),],
         aes(x=trial_number_in_block,y=scale(rpd),group=interaction(manipulation,order,trial_type),color=interaction(manipulation,order),linetype=trial_type))+geom_smooth()+theme_bw()+
          xlab('trial within block')+ylab('relative pupil size (z)')
  
  dev.off()
  
### --> gaze data ####
  
#raw pupil data
hist(df_oddball$left_pupil_measure1)
hist(df_oddball$right_pupil_measure1)

#raw gaze data
hist(df_oddball$left_gaze_x[df_oddball$left_gaze_x<200 & df_oddball$left_gaze_x>-200])
hist(df_oddball$right_gaze_x[df_oddball$right_gaze_x<200 & df_oddball$right_gaze_x>-200])
hist(df_oddball$left_gaze_y[df_oddball$left_gaze_y<200 & df_oddball$left_gaze_y>-200])
hist(df_oddball$right_gaze_y[df_oddball$right_gaze_y<200 & df_oddball$right_gaze_y>-200])

require(ggplot2)
require(hexbin)

## OLD STYLE
# df_oddball<-df_oddball[is.finite(df_oddball$left_gaze_x) & is.finite(df_oddball$left_gaze_x),]
# ggplot(df_oddball[df_oddball$left_gaze_x<100 & 
#                     df_oddball$left_gaze_x>-100 & 
#                     df_oddball$left_gaze_y<100 & 
#                     df_oddball$left_gaze_y>-100,],aes(x=left_gaze_x,y=left_gaze_y))+geom_hex(bins=30)+
#   scale_fill_gradientn(colours=rev(rainbow(3)))

#-->refined heatmap
ggplot(df_oddball[df_oddball$ts_trial<2 & is.finite(df_oddball$rpd) & is.finite(df_oddball$left_gaze_x) & is.finite(df_oddball$left_gaze_y),],aes(x=left_gaze_x,y=left_gaze_y))+
  stat_density_2d(aes(fill=after_stat(density)), n=50, geom = "raster", contour = FALSE)+
  xlim(-200,200)+ylim(-200,200)+
  scale_fill_gradientn(colours=rev(rainbow(3)))


## --> DATA ANALYSIS: MANIPULATION ####

##select oddball data
df_manip<-df[df$phase %in% c('manipulation_block'),]

#relevant variables
table(df_manip$trial)
table(df_manip$manipulation_trial_counter)
with(df_manip,table(trial,manipulation_trial_counter))

#merge with baseline PD
df_manip<-merge(df_manip,df_meanpd,by='id')

#baselinecorrected pd
df_manip$rpd<-df_manip$pd-df_manip$mean_pd
hist(df_manip$rpd)

#correct for non-finite values
df_manip<-df_manip[is.finite(df_manip$rpd),]

#define different baseline scenarios
df_manip$trial_type<-with(df_manip,ifelse(manipulation_trial_counter %in% c(3,7,11,15,19),'baseline_after_squeeze',
                           ifelse(manipulation_trial_counter %in% c(5,9,13,17),'baseline_after_relax',
                                  ifelse(manipulation_trial_counter %in% c(2,6,10,14,18),'squeeze',
                                         ifelse(manipulation_trial_counter %in% c(4,8,12,16,20),'relax','baseline_start')))))

#visualization 
#--> 18 is duration of squeeze phase
ggplot(df_manip[df_manip$ts_trial<18 & df_manip$trial!='baseline',],aes(x=ts_trial,y=rpd,group=trial,color=trial))+
  geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupillary response (z)',title='pupil response between manipulation phases')

ggplot(df_manip[df_manip$trial_type!='baseline_start',],aes(x=trial_type,y=rpd,fill=trial_type))+geom_boxplot()+theme_bw()+
  labs(x='manipulation phase (s)',y='standardized pupillary response (z)',title='pupil response between manipulation phases')

ggplot(df_manip,aes(x=manipulation_trial_counter,y=rpd,group=manipulation_trial_counter,fill=trial_type))+geom_boxplot()+theme_bw()

##model
require(lme4)
require(lmerTest)
require(emmeans)

lmm<-lmer(rpd~trial_type+(1|manipulation_trial_counter),data=df_manip[df_manip$trial_type!='baseline_start',])
anova(lmm)
plot(contrast(emmeans(lmm,~trial_type),'pairwise'))

####### ---------- PRELIMINARY EEG ANALYSIS ------####
### READ EEG data ####

data_files_eeg<-list.files(path=datapath_eeg,full.names=T)

list_eeg_data<-lapply(data_files_eeg,read.table, header = TRUE, fill=TRUE, skip=2)

df_MMN_500oddball<-list_eeg_data[[1]]
df_MMN_750oddball<-list_eeg_data[[2]]
df_P3_500oddball<-list_eeg_data[[3]]
df_P3_750oddball<-list_eeg_data[[4]]


df_MMN_500oddball<-reshape2::melt(df_MMN_500oddball,id.vars='File')
df_MMN_750oddball<-reshape2::melt(df_MMN_750oddball,id.vars='File')
df_P3_500oddball<-reshape2::melt(df_P3_500oddball,id.vars='File')
df_P3_750oddball<-reshape2::melt(df_P3_750oddball,id.vars='File')

df_erp<-rbind(df_MMN_500oddball,df_MMN_750oddball,df_P3_500oddball,df_P3_750oddball)

#define variables
df_erp$electrode<-substr(df_erp$variable,1,7)
df_erp$manipulation<-substr(df_erp$variable,nchar(as.character(df_erp$variable))-5,nchar(as.character(df_erp$variable)))
df_erp$pitch<-substr(df_erp$variable,nchar(as.character(df_erp$variable))-12,nchar(as.character(df_erp$variable))-9)
df_erp$reverse<-as.logical(ifelse(nchar(as.character(df_erp$variable))>42,'T','F'))
df_erp$erp<-substr(df_erp$variable,18,20)
df_erp$erp<-ifelse(df_erp$erp=='_MM','MMN',df_erp$erp)
df_erp$trial<-substr(df_erp$variable,22,26)
df_erp$trial<-ifelse(grepl('db',df_erp$trial),'oddball','standard')

#remove characters 
df_erp$electrode<-gsub('_','',df_erp$electrode)
df_erp$manipulation<-gsub('_','',df_erp$manipulation)
df_erp$pitch<-gsub('_','',df_erp$pitch)
df_erp$erp<-gsub('_','',df_erp$erp)

#change decimal - ERP amplitude variable
df_erp$value<-gsub(',','.',df_erp$value)
df_erp$value<-ifelse(df_erp$value=='???',NA,df_erp$value)
df_erp$value<-as.numeric(df_erp$value)
df_erp$erp_amplitude<-df_erp$value
df_erp<-df_erp[,!(names(df_erp)=='value')]

#extract id
df_erp$id<-ifelse(grepl('SEGA_AuditoryOddball',df_erp$File),substr(df_erp$File,21,24),substr(df_erp$File,6,8))
df_erp$id<-as.numeric(gsub('_','',df_erp$id))

#which id has EEG and also ET data?
names(table(df_erp$id)) %in% names(table(df$id)) #3 #4
names(table(df$id)) %in% names(table(df_erp$id)) #32 #82

  
#
ggplot(df_erp[df_erp$electrode=="L.Peak",],aes(x=value))+geom_histogram()+facet_wrap(~manipulation+pitch)

table(df$id)

### converge EEG with eye-tracking oddball data on condition level ####

#prepare ET data - allign variables
df_oddball$manipulation<-ifelse(df_oddball$block_counter<8,'before','after')
df_oddball$pitch<-ifelse(df_oddball$trial=='oddball',df_oddball$oddball_frequency,df_oddball$standard_frequency)
df_oddball$reverse<-as.logical(ifelse(df_oddball$block_counter %in% c(5,12),'T','F'))

#create aggregated varaibles per condition (id,manipulation,reverse,phase)
df_oddball$merger_id<-with(df_oddball,interaction(id,manipulation,reverse,trial_type))

oddball_per_condition_high<-with(df_oddball[df_oddball$ts_trial>0.75 & df_oddball$ts_trial<1.25,],by(as.numeric(rpd),merger_id,mean,na.rm=T))
oddball_per_condition_low<-with(df_oddball[df_oddball$ts_trial>0 & df_oddball$ts_trial<0.25,],by(as.numeric(rpd),merger_id,mean,na.rm=T))
oddball_per_condition<-as.numeric(oddball_per_condition_high-oddball_per_condition_low)

id_per_condition<-as.character(with(df_oddball,by(id,merger_id,head,n=1)))
manipulation_per_condition<-as.character(with(df_oddball,by(manipulation,merger_id,head,n=1)))
reverse_per_condition<-as.logical(with(df_oddball,by(reverse,merger_id,head,n=1)))
trial_type_per_condition<-as.character(with(df_oddball,by(trial_type,merger_id,head,n=1)))


df_rpd_agg<-data.frame(id_per_condition,manipulation_per_condition,reverse_per_condition,oddball_per_condition,trial_type_per_condition)
names(df_rpd_agg)<-c('id','manipulation','reverse','rpd','trial_type')
df_rpd_agg$trial_type<-ifelse(df_rpd_agg$trial_type=='oddball','oddball','standard')

#rather plausible

##--> merge
df_erp$merger_id<-with(df_erp,interaction(id,manipulation,trial,reverse))
df_rpd_agg$merger_id<-with(df_rpd_agg,interaction(id,manipulation,trial_type,reverse))
df_rpd_agg<-df_rpd_agg[,!(names(df_rpd_agg) %in% c("id","manipulation","trial_type","reverse"))]

df_erp_rpd_agg<-merge(df_rpd_agg,df_erp,by='merger_id')

##compare
require(lme4)
summary(lmer(scale(rpd)~trial+manipulation+oddball_pitch+(1|id),df_erp_rpd_agg[df_erp_rpd_agg$reverse==F,]))
summary(lmer(scale(rpd)~trial+manipulation+(1|id),df_erp_rpd_agg[df_erp_rpd_agg$reverse==T,]))
##--> looks fine for normal, not working for reverse

#visualize
#pupillary response
ggplot(df_erp_rpd_agg,aes(x=interaction(manipulation,trial),y=rpd,fill=interaction(manipulation,trial)))+
  geom_violin()+geom_boxplot(alpha=0.4)+facet_wrap(~reverse)

ggplot(df_erp_rpd_agg,aes(rpd,fill=trial))+
  geom_density(alpha=0.2,adjust=2)+facet_wrap(~reverse+manipulation)


#erp
ggplot(df_erp_rpd_agg[df_erp_rpd_agg$electrode=='L.Peak',],aes(x=interaction(manipulation,trial),y=erp_amplitude,fill=interaction(manipulation,trial)))+
  geom_violin()+geom_boxplot(alpha=0.4)+facet_wrap(~reverse)

####display
ggplot(df_erp_rpd_agg[df_erp_rpd_agg$electrode=='L.Peak',],aes(x=rpd,y=erp_amplitude))+geom_point()+geom_smooth(method='lm')+facet_wrap(~trial+reverse+manipulation)
ggplot(df_erp_rpd_agg[df_erp_rpd_agg$electrode=='Pz.Peak',],aes(x=rpd,y=erp_amplitude))+geom_point()+geom_smooth(method='lm')+facet_wrap(~trial+reverse+manipulation)


ggplot(df_erp_rpd_agg[df_erp_rpd_agg$electrode=='L.Peak' & df_erp_rpd_agg$trial=='standard' & df_erp_rpd_agg$reverse==FALSE,],aes(x=scale(rpd),y=scale(erp_amplitude)))+geom_point()+geom_smooth(method='lm',color='chocolate2',fill='chocolate1')+
xlab('pupillary response (z)')+ylab('P3 amplitude (z)')+theme_bw()

ggplot(df_erp_rpd_agg[df_erp_rpd_agg$electrode=='L.Peak',],aes(x=rpd,y=erp_amplitude))+geom_point()+geom_smooth(method='lm')
  

### --> VISUALIZATION PANEL #####
require(gridExtra)

g1<-ggplot(df_baseline[df_baseline$phase=='baseline_calibration' & df_baseline$ts_trial<6,],aes(x=ts_trial,y=pd,group=trial,color=trial))+
  geom_smooth()+theme_bw()+xlab('time (s)')+ylab('pupil size (mm)')+labs(title='change of PD during initial baseline between IDs')

g2<-ggplot(df_manip[df_manip$ts_trial<18 & df_manip$trial!='baseline',],aes(x=ts_trial,y=rpd,group=trial,color=trial))+
  geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupillary response (z)',title='pupil response during manipulation')

g3<-ggplot(df_manip[df_manip$trial_type!='baseline_start',],aes(x=trial_type,y=rpd,fill=trial_type))+geom_boxplot()+theme_bw()+
  labs(x='manipulation phase (s)',y='standardized pupillary response (z)',title='pupil response in manipulation phases')

g4<-ggplot(df_oddball[df_oddball$ts_trial<1.8 & df_oddball$order=='normal',],aes(x=ts_trial,y=scale(rpd),group=trial_type,color=trial_type))+
  geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupil response (z)',title='effect of manipulation on stimulus-evoked pupillary response (before vs. after manipulation)')+
  facet_wrap(~manipulation)+scale_color_brewer(palette="Dark2")


tiff(file=paste0(home_path,project_path,"/output/figure_audio_results_testdata.tiff"), # create a file in tiff format in current working directory
     width=12, height=16, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure

grid.arrange(g1,g2,g3,g4, layout_matrix = rbind(c(1, 1),
                                                c(2, 3),
                                                c(4, 4)))

dev.off()


