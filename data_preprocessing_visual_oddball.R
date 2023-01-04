## SCRIPT HEADER --------------------------- 
##
## Script purpose: ANALYSE TEST DATA - SEGA PROJECT - VISUAL ODDBALL
##     - reads hdf5 data (et)
##     - matches to trial data from psychopy
##     - preliminary data preprocessing
##     - initial data analysis / data visualization
##
## Author: Nico Bast
##
## Date Created: `r paste(Sys.Date())`
##
## Copyright (c) Nico Bast, `r paste(format(Sys.Date(), "%Y"))`
## Email: nico.bast@kgu.de
##
## NOTES --------------------------- 

        # - INFO ON HDF5 file -#
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
require(rhdf5)
require(data.table) #fread uses parallelization and thus much faster than read.csv
require(zoo) #na.approx
require(ggplot2)
 
 #PATHS###
 
 #check for OS --> define home path (script independent of OS)
 ifelse(Sys.info()['sysname']=='Linux',
        home_path<-'~',
        home_path<-'C:/Users/Nico')
 
 
 #project path
 project_path<-'/PowerFolders/project_sega'
 data_path<-'/PowerFolders/project_sega/data/VisualOddball' 
 datapath<-paste0(home_path,data_path)
 
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
   
   # # STEP 2 - filtering ####
   # ## A) normalized dilation speed, take into account time jumps with Remotetimestamps: ####
   # #maximum change in pd compared to last and next pd measurement
   # #Left
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
 
 
# READ ET DATA ####
 
 #read data from datapath and store in according objects
 data.files<-list.files(path=datapath,full.names=T)
 data_files_et<-data.files[grepl('.hdf5',data.files)]
 
 #read all files to list (hdf5 is a container format and et data is stored in BinocularEyeSampleEvent folder within the container) 
 list_et_data<-list(0)
 for(i in 1:length(data_files_et)){
   list_et_data[[i]]<-h5read(file = data_files_et[i], name = "data_collection/events/eyetracker/BinocularEyeSampleEvent")
   print(paste0('read ET data file: ',i))
 }
 
  #close open handles
  h5closeAll()

  #add names to list
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
  
  #missing gaze data (~10%)
  sapply(list_et_data,function(x){round(table(is.na(x$left_gaze_x))[2]/sum(table(is.na(x$left_gaze_x))),2)})
  #missing pupil data (~10%) 
  sapply(list_et_data,function(x){round(table(is.na(x$left_pupil_measure1))[2]/sum(table(is.na(x$left_pupil_measure1))),2)})
  #time gap --> should be 0.0033
  sapply(list_et_data,function(x){hist(diff(x$time[-1]))}) #delete first entry

# READ TRIAL DATA ####
  
  #read data from datapath and store in according objects
  data.files<-list.files(path=datapath,full.names=T)
  data_files_trial<-data.files[grepl('.csv',data.files)]
  
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
  
  #remove empty variable
  list_trial_data<-lapply(list_trial_data,function(x){x$id<-as.character(x$id);return(x)})
  
  #add names to list
  id.names<-substr(data_files_trial,nchar(datapath)+2,nchar(data_files_trial)-9) 
  names(list_trial_data)<-id.names
  
# INITIAL TRIAL DATA QUALITY ####  
  
  sapply(list_trial_data,function(x){sum(table(x$.thisTrialN)[-1])}) # total number of trials (without practice trials [-1])
    sapply(list_trial_data,function(x){table(x$.thisRepN)}) #repetitions of trial sequence (4 standard + 1)
    sapply(list_trial_data,function(x){table(x$.thisN)}) # trials in a sequence
    sapply(list_trial_data,function(x){table(x$phase)}) 
    sapply(list_trial_data,function(x){table(x$block_counter)}) 
  
      ## --> CREATE df-trial ####
      df_trial<-plyr::rbind.fill(list_trial_data)   
        
      hist(df_trial$stimulus_duration[df_trial$stimulus_duration<1],30) #stimulus duration
      hist(df_trial$ISI_duration[df_trial$ISI_duration<3],30) #sometimes ISI is too long
      table(df_trial$gaze_offset_duration)
      table(df_trial$trial_nodata_duration)
      
      ## --> TRIAL DATA DISTRIBUTION ####
      
      hist(df_trial$stimulus_duration[df_trial$stimulus_duration<0.2],col=df_trial$.thisIndex,breaks=20,main='visual oddball - actual stimulus duration',xlab='time (s)')
      
      require(ggplot2)
      ggplot(df_trial[df_trial$stimulus_duration<0.2,],
             aes(x=stimulus_duration,fill=as.factor(.thisTrialN)))+geom_histogram()+theme_bw()+facet_wrap(~as.factor(.thisTrialN))
      
      ggplot(df_trial[df_trial$stimulus_duration<0.2,],
             aes(x=stimulus_duration,fill=as.factor(.thisIndex)))+geom_histogram()+theme_bw()+facet_wrap(~as.factor(trial))
      
      table(df_trial$trial,df_trial$.thisTrialN)
      
      hist(df_trial$ISI_duration[df_trial$ISI_duration<3],col='grey',breaks=20,main='visual oddball - ISI duration',xlab='time (s)')
      
# MERGE ET AND TRIAL DATA ####
  
    #how to match - as there are different timestamp formats
      ##--> df_trial$timestamp_exp ~~ df$logged_time
    
        summary(df_trial$timestamp)
        format(df_trial$timestamp[[100]],scientific=F)
        format(df_trial$timestamp_tracker[[100]],scientific=F)
        #--> unix epoch (seconds since 01.01.1070) - psychopy timestamp
        
        summary(df_trial$timestamp_exp) #with df$logged_time or df$time
        hist(df_trial$timestamp_exp)
        
        ##-->
        
## --> MERGE DATA  ##
  #note: only returns et data that can be matched to trial data

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
    
df_list<-mapply(fun_merge_all_ids,et_data=list_et_data,trial_data=list_trial_data,SIMPLIFY=F)

#define a trial_index and trial_number variable
#INFO: Trial is defined --> .thisTrialN (trials within sequence) x .thisRepn (repetition of sequences in block) x block_counter (blocks)
df_list<-lapply(df_list,function(x){
  x$trial_index<-with(x,interaction(block_counter,.thisRepN,.thisTrialN)) #name to indentify individual trials
  x$trial_number<-with(x,rep(seq_along(table(trial_index)),times=table(trial_index))) #all trials including empty trials
  x$trial_number<-with(x,rep(seq_along(table(trial_number)),times=table(trial_number))) #rerun remove empty trials from trial number
  return(x)})

df<-dplyr::bind_rows(df_list) #faster than plyr::rbind.fill
#df<-plyr::rbind.fill(df_list)

#df<-df[sapply(df,nrow)!=0] #remove frames without entries

    # what does this variable encode - tracking status?
    #convert raw to numeric
    #test$status<-as.numeric(as.character(test$status))
    table(as.character(df$status))


table(df$phase)

###--> create ts_trial variable ####
df$ts_trial<-df$logged_time-df$timestamp_exp
hist(df$ts_trial[df$ts_trial<5])

##--> DATA ANALYSIS: BASELINE PHASE ####    

df_baseline<-df[df$phase %in% c('baseline','baseline_calibration'),]

  with(df_baseline[df_baseline$phase=='baseline_calibration',],table(trial))
  #baseline, bloackslide, whiteslide
  
  with(df_baseline,by(timestamp_exp,interaction(baseline_trial_counter,phase,id),mean,na.rm=T))
  #--> first comes baseline_calibration
  
  df_baseline<-func_pd_preprocess(df_baseline)
  df_baseline<-df_baseline[is.finite(df_baseline$pd),]
  
        # #baseline pupil size (start of experiment)
        # ggplot(df_baseline[df_baseline$phase=='baseline_calibration',],aes(x=trial,y=pd,fill=id))+geom_boxplot()
        # ggplot(df_baseline[df_baseline$phase=='baseline_calibration',],aes(x=ts_trial,y=pd,group=trial,color=trial))+geom_smooth()+facet_wrap(~id)
        # 
        # #baseline pupil size in experiment
        # tiff(file=paste0(home_path,project_path,"/output/figure_visual_effectofmanipulation_testdata.tiff"), # create a file in tiff format in current working directory
        #      width=8, height=4, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure
        # 
        # ggplot(df_baseline,aes(x=as.factor(baseline_trial_counter),y=pd,fill=id))+
        #   geom_boxplot()+facet_wrap(~id)+geom_vline(xintercept=2.5,lty=2)+theme_bw()+
        #   xlab('baseline PD measurements')+ylab('pupil size (mm)')+labs(title='effect of manipulation')
        # 
        # dev.off()
        # 
        # #pd change in baseline calibration
        # tiff(file=paste0(home_path,project_path,"/output/figure_visual_changeofpd_testdata.tiff"), # create a file in tiff format in current working directory
        #      width=8, height=4, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure
        # 
        # ggplot(df_baseline[df_baseline$phase=='baseline_calibration' & df_baseline$ts_trial<6,],aes(x=ts_trial,y=pd,group=trial,color=trial))+
        #   geom_smooth()+facet_wrap(~id)+xlim(c(0,5))+theme_bw()+xlab('time (s)')+ylab('pupil size (mm)')+labs(title='change of PD')
        # 
        # dev.off()
        
        
    #--> estimate mean PD per ID ####
    df_meanpd<-with(df_baseline[df_baseline$phase=='baseline_calibration' & df_baseline$trial=='baseline',],by(as.numeric(pd),id,mean,na.rm=T))
  
    df_meanpd<-data.frame(names(df_meanpd),as.numeric(df_meanpd))
    names(df_meanpd)<-c('id','mean_pd')
    
    
##--> DATA ANALYSIS: ODDBALL PHASE ####    

##select oddball data
df_oddball<-df[df$phase %in% c('oddball_--','oddball_-+','oddball_+-','oddball_++'),]
table(df_oddball$trial)

#basic preprocessing (has to be elaborated in final data analysis)
df_oddball<-func_pd_preprocess(df_oddball)

##trial number oddball phase
df_oddball$trial_index_oddballphase<-with(df_oddball,interaction(block_counter,.thisRepN,.thisTrialN)) #name to indentify individual trials
df_oddball$trial_number_oddballphase<-with(df_oddball,rep(seq_along(table(trial_index_oddballphase)),times=table(trial_index_oddballphase))) #all trials including empty trials
df_oddball$trial_number_oddballphase<-with(df_oddball,rep(seq_along(table(trial_number_oddballphase)),times=table(trial_number_oddballphase))) #rerun remove empty trials from trial number
#hist(df_oddball$trial_number_oddballphase,50)

#estimate mean pd of every trial first 250 ms
id_by_trial<-with(df_oddball,interaction(id,trial_number_oddballphase))
baseline_pd_trial<-with(df_oddball,by(as.numeric(pd),id_by_trial,mean,na.rm=T))
df_meanpd_trial<-data.frame(names(baseline_pd_trial),as.numeric(baseline_pd_trial))
names(df_meanpd_trial)<-c('id_by_trial','mean_pd_trial')

#merge with pd trial (first 250ms)
df_oddball$id_by_trial<-with(df_oddball,interaction(id,trial_number_oddballphase))
df_oddball<-merge(df_oddball,df_meanpd_trial,by='id_by_trial')

#merge with baseline PD
df_oddball<-merge(df_oddball,df_meanpd,by='id')

#baselinecorrected pd
df_oddball$rpd<-df_oddball$pd-df_oddball$mean_pd
hist(df_oddball$rpd)

#baseline per trial corrected pd
df_oddball$rpdt<-df_oddball$pd-df_oddball$mean_pd_trial
hist(df_oddball$rpdt,50)

#correct for non-finite values
df_oddball<-df_oddball[is.finite(df_oddball$rpd),]
df_oddball<-df_oddball[is.finite(df_oddball$rpdt),]

#trial type (standard versus oddball)

table(df_oddball$trial)
table(df_oddball$phase)

##### -- visualization

### ODDBALL RESPONSE PER CONDITION
  tiff(file=paste0(home_path,project_path,"/output/figure_visual_effectofmanipulation_pupilresponse.tiff"), # create a file in tiff format in current working directory
       width=8, height=4, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure

# #rpd - mean pd corrected
# ggplot(df_oddball[df_oddball$ts_trial<2.2,],aes(x=ts_trial,y=scale(rpd),group=interaction(trial,phase),color=phase,lty=trial))+
#   geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupil response (z)',title='effect of manipulation on pupil response')

#rpd per trial 250ms corrected
ggplot(df_oddball[df_oddball$ts_trial<2.2,],aes(x=ts_trial,y=scale(rpdt),group=interaction(trial,phase),color=phase,lty=trial))+
  geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupil response (z)',title='effect of manipulation on pupil response')+
  scale_color_brewer(palette="Dark2")

  
  dev.off()

  
  

#table(df$block_counter,df$phase) #--> manipulation is always block 8
  
### performance data
  
  #convert to numeric
  df_trial$rt<-as.numeric(substr(df_trial$responses_rt,2,nchar(df_trial$responses_rt)-1))
  
  hist(df_trial$rt)
  
  with(df_trial[df_trial$phase %in% c('oddball_--','oddball_-+','oddball_+-','oddball_++'),],by(rt,phase,psych::describe))
  
  ggplot(df_trial[df_trial$phase %in% c('oddball_--','oddball_-+','oddball_+-','oddball_++'),], aes(x = rt, colour = phase)) +
    geom_density()

  summary(lm(scale(rt)~phase,data=df_trial[df_trial$phase %in% c('oddball_--','oddball_-+','oddball_+-','oddball_++'),]))  
  ###--> manipulation check - faster performance in oddball_++ condition
  
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


df_oddball<-df_oddball[is.finite(df_oddball$left_gaze_x) & is.finite(df_oddball$left_gaze_x),]
ggplot(df_oddball[df_oddball$left_gaze_x<200 & 
                    df_oddball$left_gaze_x>-200 & 
                    df_oddball$left_gaze_y<200 & 
                    df_oddball$left_gaze_y>-200,],aes(x=left_gaze_x,y=left_gaze_y))+geom_hex()


##model
require(lme4)
require(lmerTest)
require(emmeans)

lmm<-lmer(rpd~trial*phase+(1|id),data=df_oddball)
anova(lmm)
contrast(emmeans(lmm,~trial|phase),'pairwise')
contrast(emmeans(lmm,~trial|phase),'pairwise')


plot(contrast(emmeans(lmm,~trial|phase),'pairwise'))+theme_bw()
plot(contrast(emmeans(lmm,~phase+trial),'eff'))+theme_bw()+geom_vline(xintercept=0,lty=2)



### --> VISUALIZATION PANEL #####
require(gridExtra)

ggplot(df_oddball[df_oddball$ts_trial<2.2,],aes(x=ts_trial,y=scale(rpd),group=interaction(trial,phase),color=phase,lty=trial))+
  geom_smooth()+theme_bw()+labs(x='trial duration (s)',y='standardized pupil response (z)',title='effect of manipulation on pupil response')+
  scale_color_brewer(palette="Dark2")


#lmm needs to be calculated ahead
g2<-plot(contrast(emmeans(lmm,~phase+trial),'eff'))+theme_bw()+geom_vline(xintercept=0,lty=2)+
  labs(x='marginal effect',title='condition and stimulus on PD')


tiff(file=paste0(home_path,project_path,"/output/figure_visual_results_testdata.tiff"), # create a file in tiff format in current working directory
     width=8, height=16, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure

grid.arrange(g1,g2,layout_matrix=rbind(1,1,2))

dev.off()

