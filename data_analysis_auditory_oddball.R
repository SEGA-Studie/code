### --- LOADS preprocessed auditory oddball data and analysis pupillometry variables

# SETUP ####

# Can be used to skip preprocessing and directly read proprocessed data from .rds file:
df_trial <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ETdata.rds'))

#changed random intercept to a factor
df_trial$id<-as.factor(df_trial$id) #change ID to factor
#df_trial$trial<-as.factor(df_trial$trial)

require(performance)


### DATA ANALYSIS  ####

lmm <- lmer(
  scale(rpd) ~ oddball * manipulation * reverse + trial_number_in_block + (1 | id),
  data = df_trial[
    df_trial$phase %in% c("oddball_block", "oddball_block_rev"), ])
anova(lmm)
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

# WTAS 2024 (analysis + data: 14.11.2023)
# LMM: SEPR
lmm_sepr_wtas <- lmer(scale(rpd)~oddball*group*manipulation + (1|id),
                      data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm_sepr_wtas) # -> 2 sig. main effects (group + oddball)

# contrast oddball
contrast(emmeans(lmm_sepr_wtas, ~ oddball), "pairwise")
# -> stronger SEPR to oddballs vs. standards across groups

# contrast group
contrast(emmeans(lmm_sepr_wtas, ~ group), "pairwise")
# -> stronger SEPR in TD than in ASD group

# LMM: BPS
lmm_bps_wtas <- lmer(scale(rpd_low)~oddball*group*manipulation + (1|id),
                     data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm_bps_wtas) # -> interaction effect of group x manipulation

# contrast manipulation x group interaction
contrast(emmeans(lmm_bps_wtas, ~ manipulation|group), "pairwise")
emmip(lmm_bps_wtas, ~ manipulation|group) # plot: manipulation effect on BPS for both groups
# -> after manipulation higher BPS in both groups

# LMM (SEPR) with only forward blocks
lmm_sepr_forward <- lmer(scale(rpd)~oddball*group*manipulation + (1|id),
                         data=df_trial[df_trial$phase %in% c('oddball_block'),])
anova(lmm_sepr_forward) # -> sig. main effects of group + oddball and trend of manipulation
contrast(emmeans(lmm_sepr_forward, ~ manipulation), "pairwise")
emmip(lmm_sepr_forward, ~ manipulation|group)

# LMM (BPS) with only forward blocks
lmm_bps_forward <- lmer(scale(rpd_low)~oddball*group*manipulation + (1|id),
                        data=df_trial[df_trial$phase %in% c('oddball_block'),])
anova(lmm_bps_forward) # -> sig. main effects of group + oddball and trend of manipulation
contrast(emmeans(lmm_bps_forward, ~ manipulation|group), "pairwise")
emmip(lmm_bps_forward, ~ manipulation|group)

# #these effects are not significant - as of Sept 2023
# contrast(emmeans(lmm,~group),'pairwise') ##--> TD with higher response to all trials
# contrast(emmeans(lmm,~group|manipulation),'pairwise') ##--> higher response in TD before manipulation
# contrast(emmeans(lmm,~manipulation|group),'pairwise') ##--> manipulation has an  effect in ASD
# contrast(emmeans(lmm,~group|oddball+manipulation+reverse),'pairwise')
# ###--> TD compared to ASD with higher response to standards in foward trials before manipulation
# ###--> TD comapred to ASD with higehr repsonse to oddballs in reverse trials before manipulation

# rpd in block progression
lmm1 <- lmer(
  scale(rpd) ~ oddball * reverse * trial_number_in_block + (1|id),
  data = df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm1)
# -> no main effect of trial_number_in_block (p = 0.13)

emtrends(lmm1,~oddball,var = 'trial_number_in_block')
# -> standard + oddball responses decrease with trial_number_in_block

emtrends(lmm1,~reverse,var = 'trial_number_in_block')
# -> pupil responses towards all stimuli in forward + reverse blocks decrease with trial_number_in_block

emtrends(lmm1,~oddball|reverse,var='trial_number_in_block')
# -> In forward blocks, slight trend towards decreasing response to oddballs
# and increasing to standards, which is the other way round in reverse blocks.

hist(df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_reverse"), ]$trial_number_in_block)
table(df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_reverse"), ]$trial_number_in_block)
# -> Each block contains 100 x trial_number_in_block, thus 100 oddball trials per block.

bin_size <- 20
df_trial <- df_trial %>%
  mutate(trial_number_in_block_binned = factor(trial_number_in_block%/%bin_size*20))

ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev"), ],
  aes(x = trial_number_in_block_binned,
      y = scale(rpd))) + 
  geom_boxplot() + 
  facet_grid(
    rows = vars(reverse),
    cols = vars(manipulation)) +
  theme_bw() +
  ggtitle("rpd during blocks, split by forward + reverse \nbefore + after manipulation") 


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

# DATA ANALYSIS: grip strength ####

hist(df_trial$mean_grip_strength)

df_trial$mean_grip_strength_z<-scale(df_trial$mean_grip_strength)
#pupillary response - SEPR 
lmm<-lmer(scale(rpd)~mean_grip_strength_z*manipulation*oddball*group*reverse+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm) #--> no effect of hand grip strength

#pupillary response - BPS
lmm<-lmer(scale(rpd_low)~mean_grip_strength_z*manipulation*oddball*group*reverse+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm)

fixef(lmm) #higher hand grip is assoicated with lower bps
emtrends(lmm,~manipulation|group+reverse,var = 'mean_grip_strength_z') ##--> higher pupil size after manipulation
#the hand grip strength has an effect on rpd only in the TD group


# DATA ANALYSIS: PRELIMINARY EEG analysis ####
# List and read MMN data
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
df_P3a_500oddball <- list_eeg_data[[3]]
df_P3a_750oddball <- list_eeg_data[[4]]

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
  "L-Peak_P3a_P3a_Standard_rev_500_Hz_before",
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

# Correct column names for 750Hz-P3a data frame
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

# Remove empty columns + rows
df_MMN_500oddball <- Filter(function(x)!all(is.na(x)), df_MMN_500oddball)
df_MMN_500oddball <- df_MMN_500oddball[!(df_MMN_500oddball$`L-Peak_MMN_Oddball_500_Hz_before` == "???"), ]
df_MMN_750oddball <- Filter(function(x)!all(is.na(x)), df_MMN_750oddball)
df_MMN_750oddball <- df_MMN_750oddball[!(df_MMN_750oddball$`L-Peak_MMN_Oddball_750_Hz_before` == "???"), ]

df_P3a_500oddball <- Filter(function(x)!all(is.na(x)), df_P3a_500oddball)
df_P3a_500oddball <- df_P3a_500oddball[!(df_P3a_500oddball$`L-Peak_P3a_Oddbal_500_Hz_before` == "???"), ]
df_P3a_750oddball <- Filter(function(x)!all(is.na(x)), df_P3a_750oddball)
df_P3a_750oddball <- df_P3a_750oddball[!(df_P3a_750oddball$`L-Peak_P3a_Oddball_750_Hz_before` == "???"), ]

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

# Diff columns contain MMN-amplitude differences between oddball + standards
df_MMN_500oddball$diff_O500_S750_for_before <- (
  df_MMN_500oddball$Peak_MMN_Oddball_500_Hz_before)-(df_MMN_500oddball$Peak_MMN_Standard_750_Hz_before) # 1st block (for) before
df_MMN_500oddball$diff_O750_S500_rev_before <- (
  df_MMN_500oddball$Peak_MMN_Oddball_rev_750_Hz_before)-(df_MMN_500oddball$Peak_MMN_Standard_rev_500_Hz_before) # 2nd block (rev) before
df_MMN_500oddball$diff_O500_S750_for_after <- (
  df_MMN_500oddball$Peak_MMN_Oddball_500_Hz_after)-(df_MMN_500oddball$Peak_MMN_Standard_750_Hz_after) # 3rd block (for) after
df_MMN_500oddball$diff_O750_S500_rev_after <- (
  df_MMN_500oddball$Peak_MMN_Oddball_rev_750_Hz_after)-(df_MMN_500oddball$Peak_MMN_Standard_rev_500_Hz_after) # 4th block (rev) after

df_MMN_750oddball$diff_O750_S500_for_before <- (
  df_MMN_750oddball$Peak_MMN_Oddball_750_Hz_before)-(df_MMN_750oddball$Peak_MMN_Standard_500_Hz_before) # 1st block (for) before
df_MMN_750oddball$diff_O500_S750_rev_before <- (
  df_MMN_750oddball$Peak_MMN_Oddball_rev_500_Hz_before)-(df_MMN_750oddball$Peak_MMN_Standard_rev_750_Hz_before) # 2nd block (rev) before
df_MMN_750oddball$diff_O750_S500_for_after <- (
  df_MMN_750oddball$Peak_MMN_Oddball_750_Hz_after)-(df_MMN_750oddball$Peak_MMN_Standard_500_Hz_after) # 3rd block (for), after
df_MMN_750oddball$diff_O500_S750_rev_after <- (
  df_MMN_750oddball$Peak_MMN_Oddball_rev_500_Hz_after)-(df_MMN_750oddball$Peak_MMN_Standard_rev_750_Hz_after) # 4th block (rev) after

# Diff columns contain P3a-amplitude differences between oddball + standards
df_P3a_500oddball$diff_O500_S750_for_before <- (
  df_P3a_500oddball$Peak_P3a_Oddbal_500_Hz_before)-(df_P3a_500oddball$Peak_P3a_Standard_750_Hz_before) # 1st block (for) before
df_P3a_500oddball$diff_O750_S500_rev_before <- (
  df_P3a_500oddball$Peak_P3a_Oddball_rev_750_Hz_before)-(df_P3a_500oddball$Peak_P3a_Standard_rev_500_Hz_before) #2nd block (rev) before
df_P3a_500oddball$diff_O500_S750_for_after <- (
  df_P3a_500oddball$Peak_P3a_Oddball_500_Hz_after)-(df_P3a_500oddball$Peak_P3a_Standard_750_Hz_after) # 3rd block (for) after
df_P3a_500oddball$diff_O750_s500_rev_after <- (
  df_P3a_500oddball$Peak_P3a_Oddball_rev_750_after)-(df_P3a_500oddball$Peak_P3a_Standard_rev_500_Hz_after) # 4th block (rev) after

df_P3a_750oddball$diff_O750_S500_for_before <- c(
  df_P3a_750oddball$Peak_P3a_Oddball_750_Hz_before)-(df_P3a_750oddball$Peak_P3a_Standard_500_Hz_before) # 1st block (for) before
df_P3a_750oddball$diff_O500_S750_rev_before <- c(
  df_P3a_750oddball$Peak_P3a_Oddball_rev_500_Hz_before)-(df_P3a_750oddball$Peak_P3a_Standard_rev_750_Hz_before) # 2nd block (rev) before
df_P3a_750oddball$diff_O750_S500_for_after <- (
  df_P3a_750oddball$Peak_P3a_Oddball_750_Hz_after)-(df_P3a_750oddball$Peak_P3a_Standard_500_Hz_after) # 3rd block (for) after
df_P3a_750oddball$diff_O500_S750_rev_after <- (
  df_P3a_750oddball$Peak_P3a_Oddball_rev_500_Hz_after)-(df_P3a_750oddball$Peak_P3a_Standard_rev_750_Hz_after) # 4th block (rev) after

# Diff_df is a subset of df_MMN_500oddball/df_MMN_750oddball
diff_500_MMN <- df_MMN_500oddball[c("id", "diff_O500_S750_for_before", "diff_O750_S500_rev_before", "diff_O500_S750_for_after", "diff_O750_S500_rev_after")]
diff_500_MMN <- reshape2::melt(diff_500_MMN, id = "id")
diff_750_MMN <- df_MMN_750oddball[c("id", "diff_O750_S500_for_before", "diff_O500_S750_rev_before", "diff_O750_S500_for_after", "diff_O500_S750_rev_after")]
diff_750_MMN <- reshape2::melt(diff_750_MMN, id = "id")

diff_500_750_MMN <- rbind(diff_500_MMN, diff_750_MMN)
diff_500_750_MMN$value <- scale(as.numeric(diff_500_750_MMN$value))

# Diff_df is a subset of df_P3a_500oddball/df_P3a_750oddball
diff_500_P3a <- df_P3a_500oddball[c("id", "diff_O500_S750_for_before", "diff_O750_S500_rev_before", "diff_O500_S750_for_after", "diff_O750_s500_rev_after")]
diff_500_P3a <- reshape::melt(diff_500_P3a, id = "id")
diff_750_P3a <- df_P3a_750oddball[c("id", "diff_O750_S500_for_before", "diff_O500_S750_rev_before", "diff_O750_S500_for_after", "diff_O500_S750_rev_after")]
diff_750_P3a <- reshape::melt(diff_750_P3a, id = "id")

diff_500_750_P3a <- rbind(diff_500_P3a, diff_750_P3a)
diff_500_750_P3a$value <- scale(as.numeric(diff_500_750_P3a$value))

# additional columns for block, pitch + manipulation
diff_500_750_MMN$pitch <- as.factor(substr(diff_500_750_MMN$variable, 7, 9)) # MMN 
diff_500_750_MMN$block <- as.factor(substr(diff_500_750_MMN$variable, 16, 18))
diff_500_750_MMN$manipulation <- as.factor(substr(diff_500_750_MMN$variable, 20, 25))

diff_500_750_P3a$pitch <- as.factor(substr(diff_500_750_P3a$variable, 7,9)) # P3a 
diff_500_750_P3a$block <- as.factor(substr(diff_500_750_P3a$variable, 16, 18))
diff_500_750_P3a$manipulation <- as.factor(substr(diff_500_750_P3a$variable, 20, 25))

# LMM: MMN
lmm <- lmer(value ~ pitch * block * manipulation + (1|id), data = diff_500_750_MMN)
anova(lmm)
em <- emmeans(lmm, list(pairwise ~ pitch|manipulation), adjust = "tukey")
plot(em)

bxp_750_500 <- boxplot(
  value ~ variable,
  data = diff_500_750_MMN,
  ylab = "amplitude_diff",
  xlab = "task block",
  cex.axis = 0.5)

# LMM: P3a
lmm <- lmer(value ~ pitch * block * manipulation + (1|id), data = diff_500_750_P3a)
anova(lmm)

bxp_750_500 <- boxplot(
  value ~ variable,
  data = diff_500_750_P3a,
  ylab = "amplitude_diff",
  xlab = "task block",
  cex.axis = 0.5)

# the end_____ 

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

