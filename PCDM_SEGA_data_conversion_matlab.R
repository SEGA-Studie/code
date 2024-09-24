#convert data for DCPM

### taken from dataAnalysis.m of DCPM package

# Inputs:
# %
# %              input struct "in" with fields containing cell array with one
# %              cell per run of data. Field names:
# %
# %              - yPos: horizontal gaze position in dva (column vector)
# %              - yPos: vertical gaze position in dva (column vector)
# %                      (the coodinates (0,0) should be at fixation)
# %              - pupilArea: pupil area (column vector)
# %              - startInds: trial start indexes in samples (n x 2 matrix
#                                                             %                with trial start and end times, n is number of trials)
# %              - sampleRate: sampling rate of eye tracker (Hz)
# %              - trialTypes: trial types (e.g., easy / hard, or corr\error)
# %                integers for each trial, e.g. 1,2,3,4,5 (1 x n vector)
# %              - predictionWindow: time window beyond trial onset to make a
# %                prediction within (e.g. 4 sec in a jittered ISI expt)

## packages ####

require(R.matlab)
require(ggplot2)
require(data.table) #rbindlist - fast row binding of lists

require(pbapply) # progress bar for apply functions
require(data.table) #rbindlist - fast row binding of lists


## load data ####

### --> in epoched data ####

df_sega<-readRDS("C:/Users/nico/Desktop/preprocessed_auditory_ETdata.rds")
length(df_sega)
names(df_sega)
###--> list of single trials

  #extract only oddball trial data (exclude baseline data)
  oddball_data<-pblapply(df_sega,function(x){
    trial_type<-head(x$phase,1)
    ifelse(trial_type %in% c('oddball_block','oddball_block_rev'),
           return(TRUE),
           return(FALSE))})
  df_oddball_validate<-df_sega[unlist(oddball_data)]

  #split into list by participant
  ids<-sapply(df_oddball_validate,function(x){head(x$id,1)})
  df_oddball_validate<-split(df_oddball_validate,ids)

  #sample a fraction
  unique_ids<-unique(ids)
  sample_ids<-sample(unique_ids,10)

#--> LOOP across all particiapnts ####
for(participant in seq(sample_ids)){


  # select one partiicpant in epoch data
test_data<-df_oddball_validate[names(df_oddball_validate)==sample_ids[participant]]

#bind to data frame
test_data<-rbindlist(test_data[[1]])


####------- gaze data not included currently  !!! ####


## convert gaze coord from relative space to degress in visual angle
###formula: visual angle = 2 x atan (0.5* gaze_coord/screen_distance*0.1)
screen_width<-345 #mm fixed width of presentation screen EU  AIMS LEAP
screen_height<-259 #mm fixed height of presentation screen EU  AIMS LEAP
degrees_by_radian<-360/(2*pi) #fixed conversion facor
#screen_dist<-test_data$screen_dist
screen_dist<-610
x_norm <- test_data$gazepos.x-0.5
y_norm <- test_data$gazepos.y-0.5
x_dva<-degrees_by_radian*2*atan(x_norm*screen_width/(2*screen_dist))
y_dva<-degrees_by_radian*2*atan(y_norm*screen_height/(2*screen_dist))
# a = 2 * arctan(size / (2 * distance))


###convert pupil size to Arbitary units as recorded by Eye-Link
# % convert pupil size in mm to AREA
# % see https://doi.org/10.3758/s13428-015-0588-x, Experiment 1 end
# % AU = pupil_size / (rescaling_in_radians_of_AU * distance_to_screen )
#radians_per_AU<-1.7*10^(-4) #see https://doi.org/10.3758/s13428-015-0588-x, Experiment 1 end
#pd_AU <- with(test_data,pd/(radians_per_AU*screen_dist))
#pd_AU <- with(test_data,pd_res/(radians_per_AU*screen_dist))
pd_AU <- test_data$pd_res


###recalculations for correct data format
trial_start<-which(test_data$ts_epoch==1)
trial_end<-as.numeric(cumsum(with(test_data,by(ts_epoch,EventCounter_epoch,max))))
trial_number<-as.numeric(with(test_data,by(EventCounter_epoch,EventCounter_epoch,unique)))

event_types<-as.character(with(test_data,by(EventData_epoch,EventCounter_epoch,head,n=1)))
event_types<-as.numeric(sapply(event_types,function(x){switch(x,
                                                   "201" = 1,
                                                   "202" = 2,
                                                   "203" = 2,
                                                   "204" = 2)}))

# SINGLE RUN
matlab_input<-list()
matlab_input[[1]]<-x_dva #0,0 should be center
matlab_input[[2]]<-y_dva  #0,0 should be center
matlab_input[[3]]<-pd_AU
matlab_input[[4]]<-as.matrix(cbind(trial_start,trial_end))
matlab_input[[5]]<-with(test_data,round(1/frequency_rate)[1])
matlab_input[[6]]<-event_types
names(matlab_input)<-c('xPos','yPos','pupilArea','startInds','sampleRate','trialTypes')

# ###SPOLIT FOR DIFFERENT RUNS

#split_point<-700 #median(df$EventCounter)

# #split into unnested lists as this can be handled by matlab import
# split_runs_samples<-with(test_data,ifelse(EventCounter_epoch<split_point,1,2))
# split_runs_trials<-ifelse(trial_number<split_point,1,2)
#
# matlab_input<-list()
# matlab_input[1:2]<-split(x_dva,split_runs_samples) #0,0 should be center
# matlab_input[3:4]<-split(y_dva,split_runs_samples)  #0,0 should be center
# matlab_input[5:6]<-split(pd_AU,split_runs_samples)
# matlab_input[7:8]<-split(trial_start,split_runs_trials)
# matlab_input[9:10]<-split(trial_end,split_runs_trials)
# matlab_input[11:12]<-list(with(test_data,round(1/frequency_rate)[1]),with(test_data,round(1/frequency_rate)[1]))
# matlab_input[13:14]<-split(event_types,split_runs_trials)
# ##reset startinds of run 2 (relative to absolute)
# first_sample_of_run<-matlab_input[[8]][1]-1
# matlab_input[[8]]<-matlab_input[[8]]-first_sample_of_run
# matlab_input[[10]]<-matlab_input[[10]]-first_sample_of_run
#
# names(matlab_input)<-c('xPos_1','xPos_2','yPos_1','yPos_2','pupilArea_1','pupilArea_2','startInds_1s','startInds_2s','startInds_1e','startInds_2e','sampleRate_1','sampleRate_2','trialTypes_1','trialTypes_2')

#matlab_input<-lapply(matlab_input,as.numeric)

writeMat(paste0('C:/Users/nico/PowerFolders/project_matlab_pcdm/PCDM/input/',sample_ids[participant],'.mat'),input=matlab_input) #input - is the name how it appears in matlab workspace
print(paste0('saved: ',sample_ids[participant]))

}


## READ MATLAB OUTPUT AFTER PCDM FITTING ####


mat_output<-readMat("C:/Users/nico/PowerFolders/project_matlab_pcdm/PCDM/input/378195533007_wave2.mat")


# #read and analyze matlab output from DCPM ####
# mat_output<-readMat("C:/Users/nico/PowerFolders/project_matlab_pcdm/PCDM/sample_ouput.mat")
# ?readMat

##output is saved as cell types
require(R.matlab)
mat_output<-readMat("C:/Users/nico/PowerFolders/project_oddball_LEAP/data/pcdm_estimates_2runs_310323.mat")
mat_output<-readMat("C:/Users/nico/PowerFolders/project_oddball_LEAP/data/pcdm_estimates_1run_290323.mat")


#check parameters - taken from one participant
pcdm_parameters<-mat_output[[2]]
names(pcdm_parameters)<-c('interpolateBlinks','downsampleRate','fitTimeseries','change_cutoff_by_samplingRate','low_bandpass_filter_range','blinkint_vel_thres_on','blinkint_window_from_onset','blinkint_min_dur','rel_vel_thres','min_sac_dur')

list_pcdm_estimates<-list()

for (i in seq(length(mat_output[[1]]))){

pcdm_estimate<-mat_output[[1]][i][[1]][[1]] #extract relevant data of single participant
names(pcdm_estimate)<-dimnames(pcdm_estimate)[[1]] #assign names of dimensions
pcdm_estimate<-lapply(pcdm_estimate,unlist) # remove internal list structure

list_pcdm_estimates[[i]]<-pcdm_estimate

}

names(list_pcdm_estimates)<-sapply(list_pcdm_estimates,function(x){x['name']})

#extract pcdm estimates - 1run
gain_pcdm<-unlist(sapply(list_pcdm_estimates,function(x){x['gain']}))
Rsq_pcdm<-unlist(sapply(list_pcdm_estimates,function(x){x['Rsq']}))
offset_pcdm<-unlist(sapply(list_pcdm_estimates,function(x){x['offset']}))
id<-substr(names(list_pcdm_estimates),1,nchar(names(list_pcdm_estimates))-4)
id<-id[id!=""] #remove empty entries

#remove outlier
hist(gain_pcdm[gain_pcdm<0.2 & gain_pcdm>-0.2],30)
hist(Rsq_pcdm)
hist(offset_pcdm)

table(Rsq_pcdm==0) #54 zero fits, and 3 NA -290323
table(is.na(Rsq_pcdm)) #58 zero fits, and 4 NA -270323

gain_pcdm[gain_pcdm>0.2 | gain_pcdm<(-0.2)]<-NA
hist(gain_pcdm)


summary(gain_pcdm)

df_pcdm<-data.frame(id,gain_pcdm,Rsq_pcdm,offset_pcdm)


    #extract pcdm estimates - 2 runs
    gain_pcdm<-data.frame(t(matrix(unlist(sapply(list_pcdm_estimates,function(x){x['gain']})),nrow=2))) ###adapt to others
    Rsq_pcdm<-data.frame(t(matrix(unlist(sapply(list_pcdm_estimates,function(x){x['Rsq']})),nrow=2)))
    offset_pcdm<-data.frame(t(matrix(unlist(sapply(list_pcdm_estimates,function(x){x['offset']})),nrow=2)))
    id<-substr(names(list_pcdm_estimates),1,nchar(names(list_pcdm_estimates))-4)
    id<-id[id!=""] #remove empty entries

    #remove outlier
    hist(gain_pcdm[gain_pcdm<0.2 & gain_pcdm>-0.2],30)
    hist(Rsq_pcdm)
    hist(offset_pcdm)

    gain_pcdm[gain_pcdm>0.2 | gain_pcdm<(-0.2)]<-NA
    hist(gain_pcdm)

    df_pcdm<-data.frame(id,gain_pcdm,Rsq_pcdm,offset_pcdm)
    names(df_pcdm)<-c('id','gain_pcdm_1st_half','gain_pcdm_2nd_half','Rsq_pcdm_1st_half','Rsq_pcdm_2nd_half','offset_pcdm_1st_half','offset_pcdm_2nd_half')


df_timepoint<-merge(df_timepoint,df_pcdm,by='id',all.x=T)

#df_timepoint<-df_timepoint[,1:138]
names(df_timepoint)

with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],cor.test(rpd_auc.201,gain_pcdm))
with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],cor.test(rpd_auc.201,Rsq_pcdm))
with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],cor.test(rpd_auc.201,offset_pcdm))

with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],cor.test(rpd_response.201,gain_pcdm))


with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],cor.test(pd,gain_pcdm))
with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],cor.test(pd_baseline,gain_pcdm))


# sample with Rsq>0.2 ~ n =150; similar in age sex and IQ
with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],table(t1_diagnosis))
with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],t.test(ageyrs~t1_diagnosis))
with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],chisq.test(sex,t1_diagnosis))
with(df_timepoint[df_timepoint$Rsq_pcdm>0.2,],t.test(t1_piq~t1_diagnosis))

#higher gain for standards in ASD
require(effectsize)
#1run
summary(lm(scale(gain_pcdm)~t1_diagnosis,df_timepoint[df_timepoint$Rsq_pcdm>0.2,])) #higher gain in ASD for standards
cohens_d(x=scale(gain_pcdm)~t1_diagnosis,data=df_timepoint[df_timepoint$Rsq_pcdm>0.2,])
#--> d=0.4 in ASD

#2runs --> driven by second run
summary(lm(scale(gain_pcdm_1st_half)~t1_diagnosis,df_timepoint[df_timepoint$Rsq_pcdm_1st_half>0.2,])) #tendency for higher gain in ASD for standards
summary(lm(scale(gain_pcdm_2nd_half)~t1_diagnosis,df_timepoint[df_timepoint$Rsq_pcdm_2nd_half>0.2,])) #tendency for higher gain in ASD for standards
#--> gains between halfs essentially the same
hist(df_pcdm$gain_pcdm_1st_half)

df_timepoint$gain_change<-df_timepoint$gain_pcdm_2nd_half-df_timepoint$gain_pcdm_1st_half
summary(lm(scale(gain_change)~t1_diagnosis,df_timepoint[df_timepoint$Rsq_pcdm_2nd_half>0.2,])) #tendency for higher gain in ASD for standards



    require(ggplot2)
    ggplot(df_timepoint[df_timepoint$Rsq_pcdm>0.2, ],aes(x=gain_pcdm,fill=t1_diagnosis))+geom_histogram(position='jitter',alpha=0.5,bins=30)+theme_bw()


#correlations make sense
summary(lm(scale(gain_pcdm)~scale(sdq_externalising_p),df_timepoint[df_timepoint$Rsq_pcdm>0.2,]))
summary(lm(scale(gain_pcdm)~scale(ssp_audfilt),df_timepoint[df_timepoint$Rsq_pcdm>0.2,]))



summary(lm(Rsq_pcdm~t1_diagnosis,df_timepoint)) #unrelated to Rsq_pcdm

summary(lm(gain_pcdm~t1_diagnosis*schedule_enrol+Rsq_pcdm,df_timepoint)) #tendency for higher gain in ASD for standards

with(df_timepoint,plot(Rsq_pcdm,gain_pcdm))


#gain associated with externalising
summary(lm(scale(gain_pcdm)~scale(sdq_externalising_p),df_timepoint)) #tendency for higher gain in ASD for standards

with(df_timepoint,plot(scale(sdq_externalising_p),scale(gain_pcdm)))

### convert lists to cells

test<-readMat("C:/Users/nico/PowerFolders/project_matlab_pcdm/PCDM/test_dcpm_input5.mat")


