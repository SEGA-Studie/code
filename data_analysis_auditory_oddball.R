### --- LOADS preprocessed auditory oddball data and analysis ERPs and pupillometry variables

# SETUP ####
# REQUIRED PACKAGES
require(performance) # marginal + conditional R2
require(lme4) # linear-mixed-effects models
require(lmerTest, warn.conflicts = FALSE) # linear-mixed-effects models
require(emmeans, warn.conflicts = FALSE) # estimated marginal means (EMMs)
library(ggpubr, warn.conflicts = FALSE) # ggscatter()-function
library(knitr) # dynamic report generation
library(kableExtra) # table formatting
library(tidyverse)
library(cowplot, warn.conflicts = FALSE) # get only legend from plot
library(simr) # power analysis
library(gridExtra) # arrange ggplots

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
  home_path <- "/Volumes/common/"
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

# Can be used to skip preprocessing and directly read proprocessed data from .rds file:
df_trial <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ETdata.rds')) # trial pupil data including ALL trials
df <- readRDS(paste0(home_path, project_path,'/data/preprocessed_auditory_df.rds')) # pupil data (all data points for visualizing)
et_erp_subject <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ET_ERP_subject.rds')) # eeg + pupil (subject level, only oddball(rev) blocks)
et_erp_trial <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ET_ERP_trial.rds')) # eeg + pupil(single-trial, only oddball(rev) blocks)
MMN_diff <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_MMN_diff.rds'))

# Renaming the variable "trial" into "stimulus" to match terminology in the paper
names(et_erp_subject)[names(et_erp_subject) == "trial"] <- "stimulus"
names(et_erp_trial)[names(et_erp_trial) == "trial"] <- "stimulus"
names(df_trial)[names(df_trial) == "trial"] <- "stimulus"
names(df)[names(df) == "trial"] <- "stimulus"

# DATA ANALYSIS ON SUBJECT LEVEL ####
# Distributions of dependent variables
par(mfrow = c(3,2), mar = c(4, 4, 2, 1))
hist(et_erp_subject$z_rpd,
     main = "Distribution of rpd (500-1500 ms)",
     xlab = "rpd",
     xlim = c(-4, 6),
     breaks = 200)
hist(et_erp_subject$z_rpd_low,
     main = "Distribution of rpd_low (0-250 ms)",
     xlab = "rpd_low",
     xlim = c(-4, 6),
     breaks = 200)
hist(et_erp_subject$z_rpd_block,
     main = "Distribution of rpd_block",
     xlab = "rpd_block",
     xlim = c(-4, 4),
     breaks = 200)
hist(et_erp_subject$z_MMN_amplitude,
     main = "Distribution of MMN amplitude (100-150 ms)",
     xlab = "MMN amplitude (FC1, FC2, FCz, Fz)",
     breaks = 200)
hist(et_erp_subject$z_P3a_amplitude,
     main = "Distribution of P3a amplitude (150-250 ms)",
     xlab = "P3a amplitude (Cz, FCz)",
     ylim = c(0, 30),
     xlim = c(-4, 6),
     breaks = 200)
hist(et_erp_subject$z_P3b_amplitude,
     main = "Distribution of P3b amplitude (250-500 ms)",
     xlab = "P3b amplitude (Pz)",
     xlim = c(-6, 8),
     breaks = 200)
par(mfrow = c(1, 1))

# New data frame: 1 row per subject containing sample description data (consistent across conditions)
sample_description_df <- et_erp_subject %>%
  group_by(SEGA_ID) %>%
  summarise(across(everything(), ~ {
    vals <- na.omit(.x) # removes NAs so we only check actual values
    # keeps the column if all non-NA values are identical + returns that consistent value (or NA if all were NA)
    if (length(unique(vals)) <= 1) unique(vals)[1] else NA 
  }), .groups = "drop") %>% 
  # Keep columns that have at least some non-NA values + removes columns that are completely NA across all subjects
  select(where(~ any(!is.na(.))))
sample_description_df$block_baseline_mean <- NULL

# Sample size
nrow(sample_description_df)
table(sample_description_df$group)

# Function returning group mean + sd
fun_return_descriptives <- function(group){
  group_df <- sample_description_df[
    sample_description_df$group == group, ]
  # gender
  n <- length(unique(group_df$SEGA_ID))
  male <- (length(which(group_df$gender == "männlich")))
  female <- (length(which(group_df$gender == "weiblich")))
  gender_f_m <- paste(female,"/",male)
  # age
  age_mean <- round(mean(group_df$age, na.rm = TRUE), digits = 1)
  age_sd <- round(sd(group_df$age, na.rm = TRUE), digits = 1)
  age <- paste(age_mean, "(",age_sd,")" )
  # SRS
  SRS_mean <- round(mean(group_df$SRS, na.rm = TRUE), digits = 1)
  SRS_sd <- round(sd(group_df$SRS, na.rm = TRUE), digits = 1)
  SRS <- paste(SRS_mean, "(",SRS_sd,")")
  # CBCL
  CBCL_mean <- round(mean(group_df$CBCL, na.rm = TRUE), digits = 1)
  CBCL_sd <- round(sd(group_df$CBCL, na.rm = TRUE), digits = 1)
  CBCL <- paste(CBCL_mean, "(",CBCL_sd,")")
  # SCQ
  SCQ_mean <- round(mean(group_df$SCQ, na.rm = TRUE), digits = 1)
  SCQ_sd <- round(sd(group_df$SCQ, na.rm = TRUE), digits = 1)
  SCQ <- paste(SCQ_mean, "(", SCQ_sd, ")")
  # YSR
  YSR_mean <- round(mean(group_df$YSR, na.rm = TRUE), digits = 1)
  YSR_sd <- round(sd(group_df$YSR, na.rm = TRUE), digits = 1)
  YSR <- paste(YSR_mean, "(", YSR_sd, ")")
  # SP2
  SP2_mean <- round(mean(group_df$SP2, na.rm = TRUE), digits = 1)
  SP2_sd <- round(sd(group_df$SP2, na.rm = TRUE), digits = 1)
  SP2 <- paste(SP2_mean, "(", SP2_sd, ")")
  # grip strength
  grip_strength_mean <- round(mean(group_df$z_grip_strength, na.rm = TRUE), digits = 1)
  grip_strength_sd <- round(sd(group_df$z_grip_strength, na.rm = TRUE), digits = 1)
  grip_strength <- paste(grip_strength_mean, "(", grip_strength_sd, ")")
  # verbal IQ
  verbal_IQ_mean <- round(mean(group_df$verbal_IQ, na.rm = TRUE), digits = 1)
  verbal_IQ_sd <- round(sd(group_df$verbal_IQ, na.rm = TRUE), digits = 1)
  verbal_IQ <- paste(verbal_IQ_mean, "(", verbal_IQ_sd, ")")
  # non-verbal IQ
  non_verbal_IQ_mean <- round(mean(group_df$non_verbal_IQ, na.rm = TRUE), digits = 1)
  non_verbal_IQ_sd <- round(sd(group_df$non_verbal_IQ, na.rm = TRUE), digits = 1)
  non_verbal_IQ <- paste(non_verbal_IQ_mean, "(", non_verbal_IQ_sd, ")")
  # included trials eeg
  included_trials_eeg_mean <- round((mean(group_df$included_trials_eeg))*100, digits = 1)
  included_trials_eeg_sd <- round((sd(group_df$included_trials_eeg))*100, digits = 1 )
  included_trials_eeg <- paste(included_trials_eeg_mean, "(", included_trials_eeg_sd, ")")
  # included trials et
  included_trials_et_mean <- round((mean(group_df$included_trials_et))*100, digits = 1)
  included_trials_et_sd <- round((sd(group_df$included_trials_et))*100, digits = 1)
  included_trials_et <- paste(included_trials_et_mean, "(", included_trials_et_sd, ")")
  
  group_description <- data.frame(
    n,
    gender_f_m,
    age,
    verbal_IQ,
    non_verbal_IQ,
    CBCL,
    YSR,
    SRS,
    SCQ,
    SP2,
    grip_strength,
    included_trials_et,
    included_trials_eeg)
  t(group_description)
}

asd_description <- fun_return_descriptives(group = "ASD")
colnames(asd_description) <- "ASD"
con_description <- fun_return_descriptives(group = "CON")
colnames(con_description) <- "CON"
mhc_description <- fun_return_descriptives(group = "MHC")
colnames(mhc_description) <- "MHC"

## p-values for descriptive statistics
SCQ_anova <- aov(SCQ ~ group, data = sample_description_df)
SCQ_anova_p <- summary(SCQ_anova)[[1]][["Pr(>F)"]][[1]]
age_anova <- aov(age ~ group, data = sample_description_df)
age_anova_p <- summary(age_anova)[[1]][["Pr(>F)"]][[1]]
CBCL_anova <- aov(CBCL ~ group, data = sample_description_df)
CBCL_anova_p <- summary(CBCL_anova)[[1]][["Pr(>F)"]][[1]]
SRS_anova <- aov(SRS ~ group, data = sample_description_df)
SRS_anova_p <- summary(SRS_anova)[[1]][["Pr(>F)"]][[1]]
YSR_anova <- aov(YSR ~ group, data = sample_description_df)
YSR_anova_p <- summary(YSR_anova)[[1]][["Pr(>F)"]][[1]]
SP2_anova <- aov(SP2 ~ group, data = sample_description_df)
SP2_anova_p <- summary(SP2_anova)[[1]][["Pr(>F)"]][[1]]
grip_strength_anova <- aov(z_grip_strength ~ group, data = sample_description_df)
grip_strength_anova_p <- summary(grip_strength_anova)[[1]][["Pr(>F)"]][[1]]
verbal_IQ_anova <- aov(verbal_IQ ~ group, data = sample_description_df)
verbal_IQ_anova_p <- summary(verbal_IQ_anova)[[1]][["Pr(>F)"]][[1]]
non_verbal_IQ_anova <- aov(non_verbal_IQ ~ group, data = sample_description_df)
non_verbal_IQ_anova_p <- summary(non_verbal_IQ_anova)[[1]][["Pr(>F)"]][[1]]
gender_chi2 <- chisq.test(sample_description_df$gender, sample_description_df$group)
gender_chi2_p <- gender_chi2$p.value
included_trials_et_anova <- aov(included_trials_et ~ group, data = sample_description_df)
included_trials_et_p <- summary(included_trials_et_anova)[[1]][["Pr(>F)"]][[1]]
included_trials_eeg_anova <- aov(included_trials_eeg ~ group, data = sample_description_df)
included_trials_eeg_p <- summary(included_trials_eeg_anova)[[1]][["Pr(>F)"]][[1]]

p_values <- c(
  NA,
  gender_chi2_p,
  age_anova_p,
  SRS_anova_p,
  CBCL_anova_p,
  SCQ_anova_p,
  SP2_anova_p,
  YSR_anova_p,
  grip_strength_anova_p,
  verbal_IQ_anova_p,
  non_verbal_IQ_anova_p,
  included_trials_et_p,
  included_trials_eeg_p)

p_value <- sapply(p_values, function(p) {
  if (is.na(p)) {
    NA
  } else if (p < 0.001) {
    "< 0.001"
  } else {
    format(p, scientific = F)
    round(p, 3)
  }
})

sample_table <- cbind(asd_description, con_description, mhc_description, p_value)
# Redo table row and column names
rownames(sample_table)[rownames(sample_table) == "grip_strength"] <- "grip strength [z]"
rownames(sample_table)[rownames(sample_table) == "gender_f_m"] <- "gender (f/m)"
rownames(sample_table)[rownames(sample_table) == "verbal_IQ"] <- "verbal IQ"
rownames(sample_table)[rownames(sample_table) == "non_verbal_IQ"] <- "non verbal IQ"
rownames(sample_table)[rownames(sample_table) == "included_trials_eeg"] <- "included eeg trials [%]"
rownames(sample_table)[rownames(sample_table) == "included_trials_et"] <- "included et trials [%]"
colnames(sample_table)[colnames(sample_table) == "p_value"] <- "p value"

sample_description_table <- kable(sample_table, caption = "Sample description", digits = 4) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% kable_styling
sample_description_table

# Supplements: No of trials per condition
fun_trials_per_cond <- function(group, variable){
  group_df <- droplevels(et_erp_trial[et_erp_trial$group == group, ])
  # subjects starting with 500 Hz oddball
  ## block 1
  df_before_forward_O500_S750 <- group_df %>% filter(
    (manipulation == "before" & phase == "oddball_block" & stimulus == "oddball" & pitch == "500") |
      (manipulation == "before" & phase == "oddball_block" & stimulus == "standard" & pitch == "750"))
  sample_size <- length(unique(df_before_forward_O500_S750$SEGA_ID))
  max_trials <- sample_size *100
  before_forward_O500_S750 <- paste(round(((sum(!is.na(df_before_forward_O500_S750[variable]), na.rm = T))/max_trials)*100, 2), "%")
  ## block 2
  df_before_reverse_O750_S500 <- group_df %>% filter(
    (manipulation == "before" & phase == "oddball_block_rev" & stimulus == "oddball" & pitch == "750") |
      (manipulation == "before" & phase == "oddball_block_rev" & stimulus == "standard" & pitch == "500"))
  before_reverse_O750_S500 <- paste(round(((sum(!is.na(df_before_reverse_O750_S500[variable]), na.rm = T))/max_trials)*100, 2), "%")
  ## block 3
  df_after_forward_O500S750 <- group_df %>% filter(
    (manipulation == "after" & phase == "oddball_block" & stimulus == "oddball" & pitch == "500") |
      (manipulation == "after" & phase == "oddball_block" & stimulus == "standard" & pitch == "750"))
  after_forward_O500S750 <- paste(round(((sum(!is.na(df_after_forward_O500S750[variable]), na.rm = T))/max_trials)*100, 2), "%")
  # block 4
  df_after_reverse_O750S500 <- group_df %>% filter(
    (manipulation == "after" & phase == "oddball_block_rev" & stimulus == "oddball" & pitch == "750") |
      (manipulation == "after" & phase == "oddball_block_rev" & stimulus == "standard" & pitch == "500"))
  after_reverse_O750S500 <- paste(round(((sum(!is.na(df_after_reverse_O750S500[variable]), na.rm = T))/max_trials)*100, 2), "%")
  # subjects starting with 750 Hz oddball
  ## block 1
  df_before_forward_O750_S500 <- group_df %>% filter(
    (manipulation == "before" & phase == "oddball_block" & stimulus == "oddball" & pitch == "750") |
      (manipulation == "before" & phase == "oddball_block" & stimulus == "standard" & pitch == "500"))
  sample_size <- length(unique(df_before_forward_O750_S500$SEGA_ID))
  max_trials <- sample_size *100
  before_forward_O750_S500 <- paste(round(((sum(!is.na(df_before_forward_O750_S500[variable]), na.rm = T))/max_trials)*100, 2), "%")
  ## block 2
  df_before_reverse_O500S750 <- group_df %>% filter(
    (manipulation == "before" & phase == "oddball_block_rev" & stimulus == "oddball" & pitch == "500") |
      (manipulation == "before" & phase == "oddball_block_rev" & stimulus == "standard" & pitch == "750"))
  before_reverse_O500S750 <- paste(round(((sum(!is.na(df_before_reverse_O500S750[variable]), na.rm = T))/max_trials)*100, 2), "%")
  ## block 3
  df_after_forward_O750_S500 <- group_df %>% filter(
    (manipulation == "after" & phase == "oddball_block" & stimulus == "oddball" & pitch == "750") |
      (manipulation == "after" & phase == "oddball_block" & stimulus == "standard" & pitch == "500"))
  after_forward_O750_S500 <- paste(round(((sum(!is.na(df_after_forward_O750_S500[variable]), na.rm = T))/max_trials)*100, 2), "%")
  ## block 4
  df_after_reverse_O500S750 <- group_df %>% filter(
    (manipulation == "after" & phase == "oddball_block_rev" & stimulus == "oddball" & pitch == "500") |
      (manipulation == "after" & phase == "oddball_block_rev" & stimulus == "standard" & pitch == "750"))
  after_reverse_O500S750 <- paste(round(((sum(!is.na(df_after_reverse_O500S750[variable]), na.rm = T))/max_trials)*100, 2), "%")
  
   trial_per_condition <- data.frame(
     before_forward_O500_S750,
     before_reverse_O750_S500,
     after_forward_O500S750,
     after_reverse_O750S500,
     before_forward_O750_S500,
     before_reverse_O500S750,
     after_forward_O750_S500,
     after_reverse_O500S750)
   t(trial_per_condition)
   }

# Call function for each group for EEG
asd_trials_per_cond_eeg <- fun_trials_per_cond(group = "ASD", variable = "z_MMN_amplitude")
colnames(asd_trials_per_cond_eeg) <- "ASD"
con_trials_per_cond_eeg <- fun_trials_per_cond(group = "CON", variable = "z_MMN_amplitude")
colnames(con_trials_per_cond_eeg) <- "CON"
mhc_trials_per_cond_eeg <- fun_trials_per_cond(group = "MHC", variable = "z_MMN_amplitude")
colnames(mhc_trials_per_cond_eeg) <- "MHC"
# Call function for each group for ET
asd_trials_per_cond_et <- fun_trials_per_cond(group = "ASD", variable = "z_rpd")
colnames(asd_trials_per_cond_et) <- "ASD"
con_trials_per_cond_et <- fun_trials_per_cond(group = "CON", variable = "z_rpd")
colnames(con_trials_per_cond_et) <- "CON"
mhc_trials_per_cond_et <- fun_trials_per_cond(group = "MHC", variable = "z_rpd")
colnames(mhc_trials_per_cond_et) <- "MHC"

condition_table <- cbind(
  asd_trials_per_cond_eeg, con_trials_per_cond_eeg, mhc_trials_per_cond_eeg,
  asd_trials_per_cond_et, con_trials_per_cond_et, mhc_trials_per_cond_et)

condition_trials_table <- kable(
  condition_table, align = "c", caption = "Number of included trials per condition") %>% 
  kable_classic(full_width = F, html_font = "Cambria") %>%
  kable_styling() %>%
  add_header_above(c(" " = 1, "Eye Tracking" = 3, "EEG" = 3))
condition_trials_table

## Plot SCQ
ggplot(sample_description_df, aes(x = group, y = SCQ), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(24, 25.5, 21)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("SCQ") +
  theme(plot.title = element_text(face = "bold"))

## Plot SRS
ggplot(sample_description_df, aes(x = group, y = SRS), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(42, 45, 40)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("SRS") +
  theme(plot.title = element_text(face = "bold"))

## Plot CBCL
ggplot(sample_description_df, aes(x = group, y = CBCL), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(87, 90, 84)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("CBCL") +
  theme(plot.title = element_text(face = "bold"))

## Plot: YSR
ggplot(sample_description_df, aes(x = group, y = YSR), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(87, 90, 84)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("YSR") +
  theme(plot.title = element_text(face="bold"))

## Plot: SP2-Auditiv
ggplot(sample_description_df, aes(x = group, y = SP2), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(43, 46, 40)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("SP2-Auditiv") +
  theme(plot.title = element_text(face="bold"))

## Plot: Age
ggplot(sample_description_df, aes(x = group, y = age), col = group) +
  geom_boxplot(fill = "grey") +
  geom_jitter() +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(19.5, 20, 19)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("Age") +
  theme(plot.title = element_text(face="bold"))

## Plot: IQ-verbal
plot_verbal_iq <- ggplot(sample_description_df, aes(x = group, y = verbal_IQ), col = group) + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(135, 140, 130)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("IQ-verbal") +
  theme(plot.title = element_text(face="bold")) +
  theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))

## Plot: IQ-non-verbal
plot_non_verbal_iq <- ggplot(sample_description_df, aes(x = group, y = non_verbal_IQ), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(135, 140, 130)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("IQ-non-verbal") +
  theme(plot.title = element_text(face="bold")) +
  theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
  
grid.arrange(plot_verbal_iq, plot_non_verbal_iq, ncol = 2)

# Trial duration by group (REVIEWER 2025-10)
## on TRIAL level
et_erp_trial$stimulus_isi_duration <- et_erp_trial$stimulus_duration + et_erp_trial$ISI_duration
ggplot(et_erp_trial, aes(x = stimulus_isi_duration)) +
  geom_histogram(binwidth = 0.020, fill = "grey", color = "white") +
  labs(x = "Duration (s)", y = "Count") +
  coord_cartesian(xlim = c(1.85, 3.4)) +
  facet_wrap(~ group) +   # separate histogram per group
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"))

et_erp_trial %>%
  group_by(group) %>%
  summarise(
    mean_stimulus_isi = round(mean(stimulus_isi_duration, na.rm = TRUE), 1),
    sd_stimulus_isi = round(sd(stimulus_isi_duration, na.rm = TRUE), 1))

stimulus_isi_duration_anova <- aov(stimulus_isi_duration ~ group, data = et_erp_trial)
stimulus_isi_duration_p <- summary(stimulus_isi_duration_anova)[[1]][["Pr(>F)"]][[1]]

## on PARTICIPANT level
et_erp_trial <- et_erp_trial %>% # same aggregated value repeated in each trial within participant
  group_by(SEGA_ID) %>%
  mutate(stimulus_isi_duration_aggregated = mean(stimulus_isi_duration, na.rm = TRUE)) %>%
  ungroup()

df_stimulus_isi_duration_aggregated <- et_erp_trial %>% # new data frame with 1 participant per row for anova
  group_by(SEGA_ID, group) %>%
  summarise(stimulus_isi_duration_aggregated = first(stimulus_isi_duration_aggregated),
            .groups = "drop")

ggplot(df_stimulus_isi_duration_aggregated, aes(x = stimulus_isi_duration_aggregated)) +
  geom_histogram(binwidth = 0.2, fill = "grey", color = "white") +
  labs(x = "Duration (s)", y = "Count") +
  coord_cartesian(xlim = c(1.85, 3.4)) +
  facet_wrap(~ group) +   # separate histogram per group
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"))

stats_stimulus_isi_duration_aggregated <- df_stimulus_isi_duration_aggregated %>%
  group_by(group) %>%
  summarise(
    mean_stimulus_isi_duration = mean(stimulus_isi_duration_aggregated, na.rm = TRUE),
    sd_stimulus_isi_duration = sd(stimulus_isi_duration_aggregated, na.rm = TRUE))
stats_stimulus_isi_duration_aggregated

anova_stimulus_isi_duration_aggregated <- aov(
  stimulus_isi_duration_aggregated ~ group, data = df_stimulus_isi_duration_aggregated)
summary(anova_stimulus_isi_duration_aggregated)

ggplot(df_stimulus_isi_duration_aggregated, aes(x = group, y = stimulus_isi_duration_aggregated, fill = group)) +
  geom_boxplot(outlier.shape = 16, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "black") +
  labs(
    x = NULL, # remioves x-axs label "group"
    y = "Stimulus_ISI_Duration (aggregated on particiopant level)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = "14"))

# Power Analysis
n_size<- 150 # sample size
k_size<-4*100 # number of trials --> auditory oddball: 
subj <- factor(1:n_size) 
trial_id <- 1:k_size # remains integer variable
condition <- c("before", "after")
group <- c("ASD", "MHC", "CON")

subj_full <- rep(subj, k_size)
trial_full <- rep(trial_id, each=n_size)
condition_full <- rep(condition, each=k_size*n_size/2)
group_full <- rep(group, k_size*n_size/3)

covars <- data.frame(id=subj_full, trial=trial_full, condition=condition_full, group=group_full)
fixed <- c(0.1, 0.2, 0.2, 0.2, 0.1, 0.1) # effect sizes
rand <- list(0.05) # random intercept
res <- 1.2 # residual standard deviation

model <- makeLmer(
  y ~ group*condition + (1|id), fixef=fixed, VarCorr=rand, sigma=res, data=covars) # create model
summary(model) # model

power_interaction <- powerSim(model, nsim=1000, test = fcompare(y ~ group + condition))
print(power_interaction)
power_curve_interaction <- powerCurve(model, test = fcompare(y ~ group + condition), along = "id")
plot(power_curve_interaction)

# EXPLORATORY (Reviewer's comment): Correlation SP2-P3a amplitude
## new df with mean P3a amplitude across conditions per participant
df_SP2_P3a <- et_erp_subject %>%
  group_by(SEGA_ID) %>%                             
  summarise(
    mean_z_P3a_amplitude = mean(z_P3a_amplitude),  
    SP2 = first(SP2)) %>% # just take the first value, is the same
  ungroup()


corr_Sp2_P3a <- df_SP2_P3a %>%
  summarise(
    cor = cor(mean_z_P3a_amplitude, SP2, use = "complete.obs"),
    p_value = cor.test(mean_z_P3a_amplitude, SP2)$p.value,
    conf_low = cor.test(mean_z_P3a_amplitude, SP2)$conf.int[1],
    conf_high = cor.test(mean_z_P3a_amplitude, SP2)$conf.int[2]
  )
print(corr_Sp2_P3a)

ggplot(df_SP2_P3a, aes(x = SP2, y = mean_z_P3a_amplitude)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  labs(
    x = "SP2",
    y = "Mean z_P3a amplitude",
    title = "Correlation between SP2 and Mean z_P3a amplitude"
  ) +
  annotate(
    "text", 
    x = min(df_SP2_P3a$SP2), 
    y = max(df_SP2_P3a$mean_z_P3a_amplitude), 
    label = paste0("r = ", round(corr_test$estimate, 2),
                   "\n95% CI [", round(corr_test$conf.int[1],2),
                   ", ", round(corr_test$conf.int[2],2), "]\n",
                   "p = ", round(corr_test$p.value, 3)),
    hjust = 0, vjust = 1, size = 4
  ) +
  theme_minimal(base_size = 14)


# CORRELATION: BLOCK-BASELINE vs. TRIALBASELINE
corr_block_trial_bl <- et_erp_subject %>%
  summarise(
    cor = cor(block_baseline_mean, rpd_low, use = "complete.obs"),
    p_value = cor.test(block_baseline_mean, rpd_low)$p.value,
    conf_low = cor.test(block_baseline_mean, rpd_low)$conf.int[1],
    conf_high = cor.test(block_baseline_mean, rpd_low)$conf.int[2]
  )
print(corr_block_trial_bl)

# Scatter plot with regression line and annotated correlation
ggplot(et_erp_subject, aes(x = block_baseline_mean, y = rpd_low)) +
  geom_point(alpha = 0.6, color = "#0072B2", size = 3) +          # points
  geom_smooth(method = "lm", se = TRUE, color = "#D55E00") +      # regression line + CI
  annotate("text", x = Inf, y = -Inf, 
           hjust = 1.1, vjust = -0.5,
           label = sprintf("r = %.2f, p < .001", corr_block_trial_bl$cor),  # add correlation
           size = 5) +
  labs(
    x = "separate baseline phase",
    y = "BPS",
    title = "Correlation between separate baseline phase and BPS"
  ) +
  theme_minimal(base_size = 14)

# RESULT 1: SEPR ON SUBJECT LEVEL
lmm <- lmer(
  z_rpd ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

## Post-hoc: stimulus
contrast(emmeans(lmm, ~ stimulus), "pairwise")
confint(contrast(emmeans(lmm, ~ stimulus), "pairwise"))

# RESULT 2: BPS ON SUBJECT LEVEL
lmm <- lmer(
  z_rpd_low ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm) 
r2_nakagawa(lmm) 

## post-hoc: manipulation * group
contrast(emmeans(lmm, ~ manipulation|group), method = "revpairwise")
confint(contrast(emmeans(lmm, ~ manipulation|group), method = "revpairwise"))
emmip(lmm, ~ manipulation | group, linearg = list(linetype = "blank"), CIs = T)

## Plot WTAS: manipulation * group
emm_data <- emmip(lmm, ~ manipulation | group, CIs = TRUE, plotit = FALSE)
ggplot(emm_data, aes(x = manipulation, y = yvar, color = group)) +
  geom_point(size = 3) +  # Increased point size
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2, linewidth = 1) +
  facet_wrap(~ group) +
  labs(x = "Manipulation", y = "Estimated Marginal Mean") +
  theme_bw() +
  scale_color_manual(values = c("ASD" = "red", "CON" = "#045D5D", "MHC" = "#A05000")) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none") +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(text = element_text(size = 14))

### exploratory: custom contrast: group comparisons of BPS before manipulation
emm <- emmeans(lmm, ~ group * manipulation)
contrast(emm, method = list(
  "ASD.before vs CON.before" = c(-1, 1, 0, 0, 0, 0)), infer = T)
contrast(emm, method = list(
  "ASD.before vs MHC.before" = c(-1, 0, 1, 0, 0, 0)), infer = T)
contrast(emm, method = list(
  "CON.before vs MHC.before" = c(0, 0, 1, 0, 0, 0)), infer = T)

## Post-hoc: manipulation * block
contrast(emmeans(lmm, ~ manipulation|block), "revpairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "revpairwise"))
emmip(lmm, ~ manipulation | block, linearg = list(linetype = "blank"), CIs = T)

## Plot: manipulation x block interaction  
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)

interaction_plot_BPS <- ggplot(plot_data) +
  geom_crossbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
    alpha = 0.8, 
    position = position_dodge(width = 0.7),  
    width = 0.3) +  
  geom_errorbar(aes(
    x = manipulation, y = emmean, color = block, fill = block,
    ymin = lower.CL, ymax = upper.CL), 
    position = position_dodge(width = 0.7), 
    width = 0.2) +  
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) +   
  theme_bw() +
  labs(title = "Manipulation x block interaction on BPS",
       x = "Manipulation",
       y = "BPS [z]") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none") 
print(interaction_plot_BPS)

## Plot BPS for WTAS 2025: manipulation x block interaction
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)

interaction_plot_BPS <- ggplot(plot_data) +
  geom_crossbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
    alpha = 0.8, 
    position = position_dodge(width = 0.9),  
    width = 0.3,
    linewidth = 0.8) +  
  geom_errorbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = lower.CL, ymax = upper.CL), 
    position = position_dodge(width = 0.9), 
    width = 0.3,
    linewidth = 0.8) +  
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) +   
  theme_bw() +
  labs(title = "Manipulation x block interaction on BPS",
       x = "Manipulation",
       y = "BPS [z]") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none") +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(text = element_text(size = 14))
print(interaction_plot_BPS)

## Plot BPS: manipulation x block interaction (SEGA_Paper)
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)

interaction_plot_BPS <- ggplot(plot_data) +
  geom_crossbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
    alpha = 0.8, 
    position = position_dodge(width = 0.9),  
    width = 0.3,
    linewidth = 0.8) +  
  geom_errorbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = lower.CL, ymax = upper.CL), 
    position = position_dodge(width = 0.9), 
    width = 0.3,
    linewidth = 0.8) +  
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) +   
  theme_bw() +
  labs(x = "Manipulation", y = "BPS [z]") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(text = element_text(size = 18))
print(interaction_plot_BPS)

## custom contrast: Manipulation effect in block 3? (before.reverse vs. after.forward)
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs forward.after" = c(0, -1, 1, 0)), infer = T)
emmip(emm, ~  block * manipulation, linearg = list(linetype = "blank"), CI = TRUE) +
  theme_bw() + labs(x = "task block", y = "Estimated Marginal Means")

### custom contrast: Manipulation effect in block 4? (before.reverse vs. after.reverse).
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs reverse.after" = c(0, -1, 0, 1)), infer = T)
emmip(emm, ~  block * manipulation, linearg = list(linetype = "blank"), CI = TRUE) +
  theme_bw() + labs(x = "task block", y = "Estimated Marginal Means")

### custom contrast: Habituation from block 1 to block 2?
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs forward.before" = c(-1, 1, 0, 0)), infer = T)

### custom contrast: Habituation from block 3 to 4
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.after vs forward.after" = c(0, 0, -1, 1)), infer = T)

# RESULT 3: MMN AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  z_MMN_amplitude ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

## Post-hoc: stimulus
contrast(emmeans(lmm, ~ stimulus), "pairwise")
confint(contrast(emmeans(lmm, ~ stimulus), "pairwise"))

## Post-hoc: manipulation * block
contrast(emmeans(lmm, ~ manipulation|block), "revpairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "revpairwise"))
emm <- emmeans(lmm, ~ block * manipulation)
emmip(emm, ~  block * manipulation, linearg = list(linetype = "blank"), CI = TRUE) +
  theme_bw() + labs(x = "task block", y = "Estimated Marginal Means")

## Plot 1: manipulation x block interaction  
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)

interaction_plot_MMN <- ggplot(plot_data) +
  geom_crossbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
                alpha = 0.8, 
                position = position_dodge(width = 0.7),  
                width = 0.3) +  
  geom_errorbar(aes(
    x = manipulation, y = emmean, color = block, fill = block,
    ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(width = 0.7), 
                width = 0.2) +  
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) +   
  theme_bw() +
  labs(title = "Manipulation x block interaction on MMN amplitude",
       x = "Manipulation",
       y = "MMN_amplitude [z]") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none") 
print(interaction_plot_MMN)

legend <- cowplot::get_legend(interaction_plot_MMN + theme(legend.position="right"))
legend_plot <- ggdraw() + draw_plot(legend)
print(legend_plot)

## Slightly modified plot for WTAS 2025
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)

interaction_plot_MMN <- ggplot(plot_data) +
  geom_crossbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
    alpha = 0.8, 
    position = position_dodge(width = 0.9),  
    width = 0.3,
    linewidth = 0.8) +  
  geom_errorbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = lower.CL, ymax = upper.CL), 
    position = position_dodge(width = 0.9), 
    width = 0.3,
    linewidth = 0.8) +  
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) +   
  theme_bw() +
  labs(title = "Manipulation x block interaction on MMN amplitude",
       x = "Manipulation",
       y = "MMN_amplitude [z]") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none") +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(text = element_text(size = 14))
print(interaction_plot_MMN)

## Plot MMN amplitude (SEGA-Paper)
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)

interaction_plot_MMN <- ggplot(plot_data) +
  geom_crossbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
    alpha = 0.8, 
    position = position_dodge(width = 0.9),  
    width = 0.3,
    linewidth = 0.8) +  
  geom_errorbar(aes(
    x = manipulation, y = emmean, color = block,
    ymin = lower.CL, ymax = upper.CL), 
    position = position_dodge(width = 0.9), 
    width = 0.3,
    linewidth = 0.8) +  
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) +   
  theme_bw() +
  labs(x = "Manipulation", y = "MMN_amplitude [z]") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(text = element_text(size = 18))
print(interaction_plot_MMN)

## Custom contrast: Manipulation effect in block 3?  (block 2 vs. 3).
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs forward.after" = c(0, -1, 1, 0)), infer = T)

## Custom contrast: Manipulation effect in block 4? (block 2 vs 4).
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs reverse.after" = c(0, -1, 0, 1)), infer = T)

## Correlation with age
corr_by_group <- et_erp_subject %>%
  group_by(group) %>%
  summarise(cor = cor(age, z_MMN_latency, use = "complete.obs"),
            p_value = cor.test(age, z_MMN_latency)$p.value)
print(corr_by_group)

# RESULT 4: MMN LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  z_MMN_latency ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

## Correlation with age
corr_by_group <- et_erp_subject %>%
  group_by(group) %>%
  summarise(cor = cor(age, z_MMN_amplitude, use = "complete.obs"),
            p_value = cor.test(age, z_MMN_latency)$p.value)
print(corr_by_group)

# RESULT 5: MMN DIFF AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  z_MMN_diff_amplitude ~  manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = MMN_diff)
anova(lmm)
r2_nakagawa(lmm)

# RESULT 6: MMN DIFF LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  z_MMN_diff_latency ~  manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = MMN_diff)
anova(lmm)
r2_nakagawa(lmm)

# RESULT 7: P3A AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  z_P3a_amplitude ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

## Post-hoc: stimulus * group * block
emm <- emmeans(lmm, ~ group * stimulus * block)
stimulus_diff <- contrast(emm, method = "pairwise", by = c("group", "block"))
stimulus_diff_emm <- as.emm_list(stimulus_diff)
group_diff <- contrast(stimulus_diff_emm, method = "pairwise", by = "block")
summary(group_diff)
confint(group_diff)

## WTAS: Plot with stimulus * group
emm_data <- emmip(lmm, ~ stimulus | group, CIs = TRUE, plotit = FALSE)
ggplot(emm_data, aes(x = xvar, y = yvar, color = stimulus)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2, size = 1) +
  facet_wrap(~ group) +
  labs(x = "Stimulus", y = "Estimated Marginal Mean") +
  theme_bw() +
  scale_color_manual(values = c("Oddball" = "red", "Standard" = "black")) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none") +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(text = element_text(size = 14))

## Post-hoc: manipulation * block
emmeans(lmm, ~ manipulation|block)
contrast(emmeans(lmm, ~ manipulation|block), "revpairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "revpairwise"))
emm <- emmeans(lmm, ~ block * manipulation)
emmip(emm, ~  block * manipulation, linearg = list(linetype = "blank"), CI = TRUE) +
  theme_bw() + labs(x = "task block", y = "Estimated Marginal Means")

## Plot: manipulation * block
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)

interaction_plot_P3a <- ggplot(plot_data, aes(x = manipulation, y = emmean, group = block, color = block)) +
  geom_crossbar(aes(ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
                alpha = 0.8, 
                position = position_dodge(width = 0.9),  
                width = 0.3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(width = 0.9), 
                width = 0.2) +
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) + 
  scale_x_discrete(expand = expansion(mult = c(0.5, 0.5))) +
  theme_bw() +
  labs(title = "Manipulation x block interaction on P3a amplitude",
       x = "Manipulation",
       y = "P3a_amplitude [z]") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")
print(interaction_plot_P3a)

## Slightly modified for WTAS 2025
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)
interaction_plot_P3a <- ggplot(plot_data, aes(x = manipulation, y = emmean, group = block, color = block)) +
  geom_crossbar(aes(ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
                alpha = 0.8, 
                position = position_dodge(width = 0.9),
                width = 0.3,
                size = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(width = 0.9), 
                width = 0.2,
                size = 0.8) +
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) + 
  scale_x_discrete(expand = expansion(mult = c(0.5, 0.5))) +
  theme_bw() +
  labs(title = "Manipulation x block interaction on P3a amplitude",
       x = "Manipulation",
       y = "P3a_amplitude [z]") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none") +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(text = element_text(size = 14))
print(interaction_plot_P3a)

legend <- cowplot::get_legend(interaction_plot_MMN + theme(legend.position="right"))
legend_plot <- ggdraw() + draw_plot(legend)
print(legend_plot)

## Plot P3a: Manipulation effect (SEGA Paper)
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)
interaction_plot_P3a <- ggplot(plot_data, aes(x = manipulation, y = emmean, group = block, color = block)) +
  geom_crossbar(aes(ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
                alpha = 0.8, 
                position = position_dodge(width = 0.9),
                width = 0.3,
                size = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(width = 0.9), 
                width = 0.2,
                size = 0.8) +
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) + 
  scale_x_discrete(expand = expansion(mult = c(0.5, 0.5))) +
  theme_bw() +
  labs(x = "Manipulation", y = "P3a_amplitude [z]") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(text = element_text(size = 18))
print(interaction_plot_P3a)

## custom contrast: Manipulation effect in block 3? (before.reverse vs. after.forward)
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs forward.after" = c(0, -1, 1, 0)), infer = T)

## custom contrast: Manipulation effect in block 4? (before.reverse vs. after.reverse)
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs reverse.after" = c(0, 1, 0, -1)), infer = T)

# RESULT 8: P3A LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  z_P3a_latency ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

## Post-hoc: stimulus
contrast(emmeans(lmm, ~ stimulus), "pairwise")
confint(contrast(emmeans(lmm, ~ stimulus), "pairwise"))

# RESULT 9: P3B AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  z_P3b_amplitude ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

## Post-hoc: stimulus
contrast(emmeans(lmm, ~ stimulus), "pairwise")
confint(contrast(emmeans(lmm, ~ stimulus), "pairwise"))

## Post-hoc: manipulation
contrast(emmeans(lmm, ~ manipulation), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation), "pairwise"))
emmip(lmm, ~ manipulation, linearg = list(linetype = "blank"), CI = T)


# RESULT 10: P3B LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  z_P3b_latency ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

# Post-hoc: stimulus * manipulation * group
contrast(emmeans(lmm, ~ group|manipulation*stimulus), "pairwise")
confint(contrast(emmeans(lmm, ~ group|manipulation*stimulus), "pairwise"))
emmip(lmm, ~ group|manipulation*stimulus, CI = T)
contrast(emmeans(lmm, ~ manipulation|group*stimulus), "pairwise")

# RESULT 11: ASSOCIATIONS OF PUPILLOMETRIC MEASURES
lmm <- lmer(
  z_rpd ~ z_rpd_low * stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

## Post-hoc: BPS
emtrends(lmm, ~ z_rpd_low, var = "z_rpd_low")
summary(emtrends(lmm, ~ z_rpd_low, var = "z_rpd_low"), infer = T)

# RESULT 12: ASSOCIATIONS: PUPIL DATA- ERPs ON SUBJECT LEVEL
## MMN amplitude
lmm <- lmer(
  z_MMN_amplitude ~ z_rpd * z_rpd_low * stimulus *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

## Post-hoc: SEPR * BPS * stimulus
emt <- emtrends(lmm, ~ z_rpd_low * z_rpd | stimulus, var = 'z_rpd', at = list(z_rpd_low = c(-2,-1,0,1,2)))
summary(emt, infer = T)

## P3a amplitude
lmm <- lmer(
  z_P3a_amplitude ~ z_rpd * z_rpd_low * stimulus *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

### post-hoc: SEPR x manipulation x group x stimulus
emt <- emtrends(lmm, ~ z_rpd  * manipulation | group * stimulus, var = 'z_rpd')
summary(emt, infer = T)

### post-hoc: BPS x stimulus x group
emt <- emtrends(lmm, ~ z_rpd_low * group | stimulus, var = 'z_rpd_low', at = list(z_rpd_low = c(-2,-1,0,1,2)))
summary(emt, infer = T)

lmm <- lmer(
  z_P3b_amplitude ~ z_rpd * z_rpd_low * stimulus *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

crl <- cor(et_erp_subject[et_erp_subject$stimulus == "Oddball", c(
  "z_rpd_low",
  "z_rpd",
  "z_MMN_amplitude",
  "z_P3a_amplitude",
  "z_P3b_amplitude"
)],
use="complete.obs")
crl
corrplot::corrplot(
  crl,
  method = "circle",
  title = "Association between pupil data and ERPs in oddball trials",
  mar=c(0,0,1,0))

crl <- cor(et_erp_subject[et_erp_subject$stimulus == "Standard", c(
  "z_rpd_low",
  "z_rpd",
  "z_MMN_amplitude",
  "z_P3a_amplitude",
  "z_P3b_amplitude"
)],
use="complete.obs")
crl
corrplot::corrplot(
  crl,
  method = "circle",
  title = "Association between pupil data and ERPs in standard trials",
  mar=c(0,0,1,0))

# DATA ANALYSIS ON TRIAL LEVEL ####
# Distributions of dependent variables
par(mfrow = c(2,2), mar = c(4, 4, 2, 1))
hist(et_erp_trial$z_rpd,
     main = "Distribution of rpd (500-1500 ms)",
     xlab = "z_rpd",
     xlim = c(-6, 6),
     breaks = 200)
hist(et_erp_trial$z_rpd_low,
     main = "Distribution of rpd_low (0-250 ms)",
     xlab = "rpd_low",
     xlim = c(-4, 6),
     breaks = 200)
hist(et_erp_trial$z_MMN_amplitude,
     main = "Distribution of MMN amplitude (100-150 ms)",
     xlab = "MMN amplitude (FC1, FC2, FCz, Fz)",
     xlim = c(-6, 6),
     breaks = 200)
hist(et_erp_trial$z_P3a_amplitude,
     main = "Distribution of P3a amplitude (150-250 ms)",
     xlab = "P3a amplitude (Cz, FCz)",
     xlim = c(-4, 6),
     breaks = 200)
par(mfrow = c(1, 1))

# RESULT 12: SEPR ON TRIAL LEVEL
lmm <- lmer(
  z_rpd ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 


## Post-hoc: stimulus
contrast(emmeans(lmm, ~ stimulus), "pairwise")
confint(contrast(emmeans(lmm, ~ stimulus), "pairwise"))

# EXPLORATORY: SEPR OVER TASK
lmm <- lmer(
  z_rpd ~ stimulus * manipulation * group * block * oddball_trial_counter + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emt <- emtrends(lmm, pairwise ~ manipulation, var = "oddball_trial_counter")
summary(emt, infer = T)

emt <- emtrends(lmm,~ stimulus * group * block, var = 'oddball_trial_counter')
summary(emt, infer = T)

et_erp_trial$binned_trial_number <- cut(
  et_erp_trial$oddball_trial_counter, c(
    1,  103, 206, 309, 412))
ggplot(et_erp_trial[et_erp_trial$stimulus == "oddball", ]) + geom_boxplot(aes(binned_trial_number, z_rpd)) + 
  facet_grid(cols = vars(group)) + 
  ggtitle("SEPR in oddball trials during the task, \nsplit by group") +
  theme(plot.title = element_text(face = "bold", size = 14))

# EXPLORATORY: SEPR OVER BLOCK
lmm <- lmer(
  z_rpd ~ stimulus * manipulation * group * block * trial_number_in_block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emt <- emtrends(lmm, pairwise ~ manipulation, var = "trial_number_in_block")
summary(emt, infer = T)

emt <- emtrends(lmm,~ stimulus * group * block, var = 'trial_number_in_block')
summary(emt, infer = T)

bin_size <- 20
df_trial <- df_trial %>%
  mutate(trial_number_in_block_binned = factor(trial_number_in_block%/%bin_size*20))
ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev") & df_trial$stimulus == "oddball", ],
  aes(x = trial_number_in_block_binned,
      y = scale(rpd))) + 
  geom_boxplot() + 
  facet_grid(
    rows = vars(block),
    cols = vars(manipulation)) +
  theme_bw() +
  ggtitle("SEPR in oddball trials during blocks, \nsplit by forward + reverse and before + after manipulation") +
  theme(plot.title = element_text(face = "bold", size = 14))

et_erp_trial$binned_trial_number <- cut(
  et_erp_trial$trial_number_in_block, c(
    0,  20, 40, 60, 80, 100))
ggplot(et_erp_trial[et_erp_trial$stimulus == "oddball", ]) + geom_boxplot(aes(binned_trial_number, z_rpd)) + 
  facet_grid(cols = vars(group)) + 
  ggtitle("SEPR in oddball trials during blocks, split by group") +
  theme(plot.title = element_text(face = "bold", size = 14))

# RESULT 13: BPS ON TRIAL LEVEL
lmm <- lmer(
  z_rpd_low ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

## Post-hoc: manipulation * group
contrast(emmeans(lmm, ~ manipulation|group), "revpairwise")
confint(contrast(emmeans(lmm, ~ manipulation|group), "revpairwise"))
emm <- emmeans(lmm, ~ manipulation * group)
emmip(emm, ~  manipulation|group, linearg = list(linetype = "blank"), CI = TRUE) +
  theme_bw() + labs(x = "task block", y = "Estimated Marginal Means")

## Post-hoc: manipulation * block
contrast(emmeans(lmm, ~ manipulation|block), "revpairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "revpairwise"))
emm <- emmeans(lmm, ~ block * manipulation)
emmip(emm, ~  block * manipulation, linearg = list(linetype = "blank"), CI = TRUE) +
  theme_bw() + labs(x = "task block", y = "Estimated Marginal Means")

## custom contrast: manipulation effect in block 3?  (block 2 vs. 3).
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs forward.after" = c(0, -1, 1, 0)), infer = T)

## custom contrast: manipulation effect in block 4? (block 2 vs 4).
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs reverse.after" = c(0, -1, 0, 1)), infer = T)

## Post-hoc: group * block
contrast(emmeans(lmm, ~ block|group), "pairwise")
confint(contrast(emmeans(lmm, ~ block|group), "pairwise"))
emm <- emmeans(lmm, ~ block | group)
emmip(emm, ~  block | group, linearg = list(linetype = "blank"), CI = TRUE) +
  theme_bw() + labs(x = "task block", y = "Estimated Marginal Means")

# EXPLORATORY: BPS OVER TASK
lmm <- lmer(
  z_rpd_low ~ stimulus * manipulation * group * block * oddball_trial_counter + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emt <- emtrends(lmm,~ manipulation * group * block, var = 'oddball_trial_counter')
summary(emt, infer = T)

et_erp_trial$binned_trial_number <- cut(
  et_erp_trial$oddball_trial_counter, c(
    1,  103, 206, 309, 412))
ggplot(et_erp_trial[et_erp_trial$stimulus == "standard", ]) + geom_boxplot(aes(binned_trial_number, z_rpd_low)) + 
  facet_grid(cols = vars(group)) + 
  ggtitle("BPS in standard trials during the task, split by group") +
  theme(plot.title = element_text(face = "bold", size = 14))

# EXPLORATORY: BPS OVER BLOCK
lmm <- lmer(
  z_rpd_low ~ stimulus * manipulation * group * block * trial_number_in_block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emt <- emtrends(lmm,~ manipulation * group * block, var = 'trial_number_in_block')
summary(emt, infer = T)

bin_size <- 20
df_trial <- df_trial %>%
  mutate(trial_number_in_block_binned = factor(trial_number_in_block%/%bin_size*20))
ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev") & df_trial$stimulus == "standard", ],
  aes(x = trial_number_in_block_binned,
      y = scale(rpd_low))) + 
  geom_boxplot() + 
  facet_grid(
    rows = vars(block),
    cols = vars(manipulation)) +
  theme_bw() +
  ggtitle("BPS in standard trials during blocks, \nsplit by forward + reverse and before + after manipulation") +
  theme(plot.title = element_text(face = "bold", size = 14))

# RESULT 8: MMN AMPLITUDE ON TRIAL LEVEL
lmm <- lmer(
  z_MMN_amplitude ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

## Post-hoc: stimulus
contrast(emmeans(lmm, ~ stimulus), "pairwise")
confint(contrast(emmeans(lmm, ~ stimulus), "pairwise"))

# RESULT 9: MMN LATENCY ON TRIAL LEVEL
lmm <- lmer(
  z_MMN_latency ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

## Post-hoc: stimulus
contrast(emmeans(lmm, ~ stimulus), "pairwise")
confint(contrast(emmeans(lmm, ~ stimulus), "pairwise"))

## Post-hoc: manipulation * block
contrast(emmeans(lmm, ~ manipulation|block), "revpairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "revpairwise"))
emm <- emmeans(lmm, ~ block * manipulation)
emmip(emm, ~  block * manipulation, linearg = list(linetype = "blank"), CI = TRUE) +
  theme_bw() + labs(x = "task block", y = "Estimated Marginal Means")

## Plot Manipulation effect on MMN Latency 
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)

interaction_plot_MNN_latency <- ggplot(plot_data, aes(x = manipulation, y = emmean, group = block, color = block)) +
  geom_crossbar(aes(ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
                alpha = 0.8, 
                position = position_dodge(width = 0.7),  
                width = 0.3) +  
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                position = position_dodge(width = 0.7), 
                width = 0.2) +  
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) +   
  theme_bw() +
  labs(title = "Manipulation x block interaction on MMN latency",
       x = "Manipulation",
       y = "MMN latency [z]") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")
print(interaction_plot_MNN_latency)

## Slightly modified plot for WTAS 2025
plot_data <- as.data.frame(emmeans(lmm, ~ block * manipulation))
plot_data <- na.omit(plot_data)

interaction_plot_MNN_latency <- ggplot(plot_data, aes(x = manipulation, y = emmean, group = block, color = block)) +
  geom_crossbar(aes(ymin = emmean - 1 * SE, ymax = emmean + 1 * SE), 
                alpha = 0.8, 
                position = position_dodge(width = 0.9),  
                width = 0.3,
                linewidth = 0.8) +  
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                position = position_dodge(width = 0.9), 
                width = 0.2,
                linewidth = 0.8) +  
  scale_color_manual(values = c("blue", "orange")) +  
  scale_fill_manual(values = c("blue", "orange")) +   
  theme_bw() +
  labs(title = "Manipulation x block interaction on MMN latency",
       x = "Manipulation",
       y = "MMN latency [z]") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none") +
  theme(axis.text = element_text(face = "bold")) +
  theme(text = element_text(size = 14))
print(interaction_plot_MNN_latency)

grid.arrange(interaction_plot_BPS,
             interaction_plot_MMN,
             interaction_plot_P3a)

# RESULT 10: P3A AMPLITUDE ON TRIAL LEVEL
lmm <- lmer(
  z_P3a_amplitude ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

## Post-hoc: stimulus * group * block
emm <- emmeans(lmm, ~ group * stimulus * block)
stimulus_diff <- contrast(emm, method = "pairwise", by = c("group", "block"))
stimulus_diff_emm <- as.emm_list(stimulus_diff)
group_diff <- contrast(stimulus_diff_emm, method = "pairwise", by = "block")
summary(group_diff)
confint(group_diff)

# Post-hoc: manipulation * block
contrast(emmeans(lmm, ~ manipulation|block), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "pairwise"))
emm <- emmeans(lmm, ~ block * manipulation)
emmip(emm, ~  block * manipulation, linearg = list(linetype = "blank"), CI = TRUE) +
  theme_bw() + labs(x = "task block", y = "Estimated Marginal Means")

## custom contrast: Manipulation effect in block 3?  (block 2 vs. 3).
emm <- emmeans(lmm, ~ block * manipulation)
contrast(emm, method = list(
  "reverse.before vs forward.after" = c(0, -1, 1, 0)), infer = T)

## custom contrast: Manipulation effect in block 4? (block 2 vs 4).
contrast(emm, method = list(
  "reverse.before vs reverse.after" = c(0, 1, 0, -1)), infer = T)

# RESULT 11: P3A LATENCY ON TRIAL LEVEL
lmm <- lmer(
  z_P3a_latency ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

# Post-hoc: stimulus
contrast(emmeans(lmm, ~ stimulus), "pairwise")
confint(contrast(emmeans(lmm, ~ stimulus), "pairwise"))

# RESULT 14: ASSOCIATION BETWEEN PUPILLOMETRIC MEASURES
ggscatter(et_erp_trial[et_erp_trial$stimulus == "oddball", ], 
          x = "z_rpd_low",
          y = "z_rpd",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "scaled BPS (rpd_low)",
          ylab = "scaled SEPR (rpd)",
          title = "Association between SEPR + BPS in oddball trials")

ggscatter(et_erp_trial[et_erp_trial$stimulus == "standard", ], 
          x = "z_rpd_low",
          y = "z_rpd",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "scaled BPS (rpd_low)",
          ylab = "scaled SEPR (rpd)",
          title = "Association between SEPR + BPS in standard trials")

# Associations of the pupillometric measures
lmm <- lmer(
  z_rpd ~ z_rpd_low * stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

### post-hoc: BPS
summary(emtrends(lmm, ~ z_rpd_low, var = "z_rpd_low"), infer = T)

### post-hoc: BPS x manipulation x block
emt <- emtrends(lmm, ~ manipulation|block, var = "z_rpd_low")
contrast(emt, "revpairwise")
confint(contrast(emt, "revpairwise"))

### EXPLORATORY (Reviewer's comment, 2025-10-27)
# lmer()-, but not poly()-function can handle NAs in the predictors, so
# et_erp_trial_for_model_fit excludes NAs for these exploratory model fit comparisons. 
et_erp_trial_for_model_fit <- na.omit(et_erp_trial[ , c(
  "z_rpd",
  "z_rpd_low",
  "stimulus",
  "manipulation",
  "group",
  "block",
  "SEGA_ID",
  "age",
  "gender")])

# linear
lmm_linear_fit <- lmer(
  z_rpd ~ z_rpd_low * stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial_for_model_fit, REML = F)
anova(lmm_linear_fit)
r2_nakagawa(lmm_linear_fit) 

# quadratic
lmm_quadratic_fit <- lmer(
  z_rpd ~ poly(z_rpd_low, 2) * stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial_for_model_fit, REML = F)
anova(lmm_quadratic_fit)
r2_nakagawa(lmm_quadratic_fit) 
# cubic
lmm_cubic_fit <- lmer(
  z_rpd ~ poly(z_rpd_low, 3) * stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial_for_model_fit, REML = F)
anova(lmm_cubic_fit)
r2_nakagawa(lmm_cubic_fit) 

table_model_compare <- anova(lmm_linear_fit, lmm_quadratic_fit, lmm_cubic_fit)

table_model_compare <- cbind(
  c("linear_fit", "quadratic_fit", "cubic_fit"),
  table_model_compare)

table_formatted_BPS_SEPR <- table_model_compare %>% 
  kbl(caption = "Model comparision of BPS-SEPR association",
      col.names = c("","number of parameters",'AIC','BIC','log likelihood','deviance','Chi-squared','df','p-value'),
      row.names = F) %>%
  kable_classic(full_width = F, html_font = "Cambria")

table_formatted_BPS_SEPR

# RESULT 15: Correlation: Pupil-ERP
crl <- cor(et_erp_trial[et_erp_trial$stimulus == "oddball", c(
  "z_rpd_low",
  "z_rpd",
  "z_MMN_amplitude",
  "z_P3a_amplitude"
)],
use="complete.obs")
crl
corrplot::corrplot(
  crl,
  method = "circle",
  title = "Association between pupil data and ERPs in oddball trials",
  mar=c(0,0,1,0))

crl <- cor(et_erp_trial[et_erp_trial$stimulus == "standard", c(
  "z_rpd_low",
  "z_rpd",
  "z_MMN_amplitude",
  "z_P3a_amplitude"
)],
use="complete.obs")
crl
corrplot::corrplot(
  crl,
  method = "circle",
  title = "Association between pupil data and ERPs in standard trials",
  mar=c(0,0,1,0))

# RESULT 16: Association between pupillometric indices and ERPs
## MMN amplitude
lmm <- lmer(
  z_MMN_amplitude ~ z_rpd * z_rpd_low * stimulus *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

### post-hoc: SEPR
summary(emtrends(lmm, ~ z_rpd, var = "z_rpd"), infer = T)

### post-hoc: BPS
summary(emtrends(lmm, ~ z_rpd_low, var = "z_rpd_low"), infer = T)

### post-hoc: SEPR * BPS *group
summary(emtrends(lmm, ~  z_rpd * z_rpd_low|group,  var = 'z_rpd', at = list(z_rpd_low = c(-2,0, 2))), infer = T)

emtr <- emtrends(lmm, ~ group | z_rpd_low, var = "z_rpd", at = list(z_rpd_low = c(-2, 0, 2)))
contrast_summary <- summary(pairs(emtr), infer = TRUE)
print(contrast_summary)

## P3a amplitude
lmm <- lmer(
  z_P3a_amplitude ~ z_rpd * z_rpd_low * stimulus *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

### Post-hoc: SEPR
summary(emtrends(lmm, ~ z_rpd, var = "z_rpd"), infer = T)

## BPS * group
summary(emtrends(lmm, ~ z_rpd_low|group,  var = 'z_rpd_low'), infer = T)

# RESULT 17: EXPLORATORY ANALYSIS OF MODEL FIT FOR "trial_number_in_block" ON BPS
linear_fit <- lmer(
  z_rpd_low ~ stimulus * manipulation * group * trial_number_in_block * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial, REML = F)
anova(linear_fit)

quadratic_fit <- lmer(
  z_rpd_low ~ stimulus * manipulation * group * poly(trial_number_in_block, 2) * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial, REML = F)
anova(quadratic_fit)
r2_nakagawa(quadratic_fit) 

cubic_fit <- lmer(
  z_rpd_low ~ stimulus * manipulation * group * poly(trial_number_in_block, 3) * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial, REML = F)
anova(cubic_fit)
r2_nakagawa(cubic_fit) 

table_model_compare <- anova(linear_fit, quadratic_fit, cubic_fit)

table_model_compare <- cbind(
  c("linear_fit", "quadratic_fit", "cubic_fit"),
  table_model_compare
)

table_formatted_BPS <- table_model_compare %>% 
  kbl(caption = "Model comparision of trial_number_in_block on BPS",
      col.names = c("","number of parameters",'AIC','BIC','log likelihood','deviance','Chi-squared','df','p-value'),
      row.names = F) %>%
  kable_classic(full_width = F, html_font = "Cambria")

table_formatted_BPS

# RESULT 18: EXPLORATIVE ANALYSIS OF MODEL FIT FOR "trial_number_in_block" on SEPR**
linear_fit <- lmer(
  z_rpd ~ stimulus * manipulation * group * trial_number_in_block * block+ (1|SEGA_ID) + pitch + age + gender,
  data = et_erp_trial, REML = F)
anova(linear_fit)
r2_nakagawa(linear_fit) 

quadratic_fit <- lmer(
  z_rpd ~ stimulus * manipulation * group * poly(trial_number_in_block, 2) * block + (1|SEGA_ID) + pitch + age + gender,
  data = et_erp_trial, REML = F)
anova(quadratic_fit)
r2_nakagawa(quadratic_fit) 

cubic_fit <- lmer(
  z_rpd ~ stimulus * manipulation * group * poly(trial_number_in_block, 3) * block + (1|SEGA_ID) + pitch + age + gender,
  data = et_erp_trial, REML = F)
anova(cubic_fit)
r2_nakagawa(cubic_fit) 

table_model_compare <- anova(linear_fit, quadratic_fit, cubic_fit)

table_model_compare <- cbind(
  c("linear_fit", "quadratic_fit", "cubic_fit"),
  table_model_compare
)

table_formatted_SEPR <- table_model_compare %>% 
  kbl(caption = "Model comparision of trial_number_in_block on SEPR",
      col.names = c("","number of parameters",'AIC','BIC','log likelihood','deviance','Chi-squared','df','p-value'),
      row.names = F) %>%
  kable_classic(full_width = F, html_font = "Cambria")

table_formatted_SEPR

# RESULT 19: ANALYSIS OF GRIP STRENGTH EFFECT
## MNN amplitude
lmm <- lmer(
  z_MMN_amplitude ~ stimulus * group * block * z_grip_strength + (1|SEGA_ID) + age + gender,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm)

## MMN latency
lmm <- lmer(
  z_MMN_latency ~ stimulus * group * block * z_grip_strength+ (1|SEGA_ID) + age + gender,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm)

## P3a amplitude
lmm <- lmer(
  z_P3a_amplitude ~ stimulus * group * block * z_grip_strength + (1|SEGA_ID) + age + gender,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm) 

## P3a latency
lmm <- lmer(
  z_P3a_latency ~ stimulus * group * block * z_grip_strength+ (1|SEGA_ID) + age + gender,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm) 

## SEPR
lmm <- lmer(
  z_rpd ~ stimulus * group * block * z_grip_strength + (1|SEGA_ID) + age + gender,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm)

## BPS
lmm <- lmer(
  z_rpd_low ~ stimulus * group * block * z_grip_strength + (1|SEGA_ID) + age + gender,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm) 
r2_nakagawa(lmm)

# RESULT 20: Covariates Models (on subject + trial level)
## MMN amplitude
lmm <- lmer(
  z_MMN_amplitude ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender + verbal_IQ + non_verbal_IQ,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

lmm <- lmer(
  z_MMN_amplitude ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender + verbal_IQ + non_verbal_IQ,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

## MMN latency
lmm <- lmer(
  z_MMN_latency ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age + verbal_IQ + non_verbal_IQ,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

## Post-hoc: age
emtrends(lmm, ~ age, var = "age")

lmm <- lmer(
  z_MMN_latency ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age + verbal_IQ + non_verbal_IQ,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

## MMN diff amplitude
lmm <- lmer(
  z_MMN_diff_amplitude ~  manipulation * group * block + (1|SEGA_ID) + gender + age + verbal_IQ + non_verbal_IQ,
  data = MMN_diff)
anova(lmm)
r2_nakagawa(lmm)

## MMN diff latency
lmm <- lmer(
  z_MMN_diff_latency ~  manipulation * group * block + (1|SEGA_ID) + gender + age + verbal_IQ + non_verbal_IQ,
  data = MMN_diff)
anova(lmm)
r2_nakagawa(lmm)

## P3a amplitude
lmm <- lmer(
  z_P3a_amplitude ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age + verbal_IQ + non_verbal_IQ,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

## Post-hoc: age
emtrends(lmm, ~ age, var = "age")

lmm <- lmer(
  z_P3a_amplitude ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age + verbal_IQ + non_verbal_IQ,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)


## P3a latency
lmm <- lmer(
  z_P3a_latency ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age + verbal_IQ + non_verbal_IQ,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

lmm <- lmer(
  z_P3a_latency ~ stimulus * manipulation * group * block + (1|SEGA_ID) + gender + age + verbal_IQ + non_verbal_IQ,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

## SEPR
lmm <- lmer(
  z_rpd ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender + verbal_IQ + non_verbal_IQ,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

lmm <- lmer(
  z_rpd ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender + verbal_IQ + non_verbal_IQ,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

## BPS
lmm <- lmer(
  z_rpd_low ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender + verbal_IQ + non_verbal_IQ,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

## Post-hoc: age
emtrends(lmm, ~ age, var = "age")

lmm <- lmer(
  z_rpd ~ stimulus * manipulation * group * block + (1|SEGA_ID) + age + gender + verbal_IQ + non_verbal_IQ,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

# ANALYSIS OF df (ONLY PUPIL DATA-ALL DATA POINTS) ####
# Baseline phase
df_baseline <- df[df$phase %in% c("baseline", "baseline_calibration"), ]
df_baseline <- df_baseline[is.finite(df_baseline$pd), ]

# Number of pupil data per baseline trial (white, black and grey slide)
with(df_baseline[df_baseline$phase == "baseline_calibration", ], table(stimulus))
# Check: Block counter should be in line with aforementioned table
with(df_baseline, table(stimulus, block_counter))
with(df_baseline, by(
  timestamp_exp, interaction(baseline_trial_counter, phase, id),
  mean, na.rm = TRUE))

# mean_pd (subject baseline) is defined as mean of baseline trial during calibration phase
df_meanpd <- with(
  df_baseline[df_baseline$phase == "baseline_calibration" &
                df_baseline$stimulus == "baseline" , ], by(as.numeric(pd), id, mean))
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
      y = pd,
      group = group,
      color = group)) +
  geom_smooth() +
  theme_bw() +
  xlab("trial duration [s]") +
  ylab("pupil dilation [mm]") +
  labs(title = "LAPR for group and trial") +
  facet_wrap(~ factor(stimulus,
                      levels = c("baseline", "baseline_whiteslide", "baseline_blackslide")))

ggsave(
  "output/baseline_calibration_pd.tiff",
  device = "tiff",
  width = 6,
  height = 4,
  units = "in",
  compression = "lzw",
  dpi = 800
)

# Oddball phase
df_oddball <- df[df$phase %in% c("oddball_block", "oddball_block_rev"), ]
table(df_oddball$stimulus)

# merge with baseline PD
df_oddball <- merge(df_oddball, df_meanpd, by = "id")

# subject-baseline corrected pd
df_oddball$subj_corr_rpd <- df_oddball$pd - df_oddball$mean_pd
hist(df_oddball$subj_corr_rpd, 50)

# correct for non-finite values
df_oddball <- df_oddball[is.finite(df_oddball$subj_corr_rpd), ]
df_oddball <- df_oddball[is.finite(df_oddball$trial_corr_rpd), ]

# trial type (standard versus oddball)
df_oddball$trial_type <- substr(df_oddball$stimulus, 1, 7)

table(df_oddball$block_counter, df_oddball$stimulus)
df_oddball$manipulation <- factor(ifelse(df_oddball$block_counter < 8, "before",
                                         ifelse(df_oddball$block_counter > 8, "after", "manipulation")),
                                  levels = c("before", "after"))

df_oddball$order <- ifelse(grepl("rev", df_oddball$stimulus), "reverse", "normal")
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

# Plot: Pupillary response across experiment
ggplot(
  df_oddball,
  aes(x = trial_number_oddballphase, y = rpd)) +
  geom_smooth(method = "lm") + facet_wrap(~stimulus + manipulation)
dev.off()

# Plot: Density of manipulation per condition
ggplot(df_oddball,
       aes(rpd, fill = interaction(manipulation))) +
  geom_density(alpha = 0.2) + facet_wrap(~stimulus)
dev.off()

# Plot: Pupil response over a trial (only forward)
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

# Plot: Pupil response over a trial
ggplot(
  df_oddball[df_oddball$ts_trial < 1.8, ],
  aes(x = ts_trial,
      y = scale(rpd),
      group = interaction(order, manipulation),
      color = order,
      linetype = trial_type)) +
  geom_smooth() +
  theme_bw() +
  labs(x = "trial duration (s)",
       y = "standardized pupil response (z)",
       title = "effect of manipulation of pupil response by block")
dev.off()


# Plot: Pupil response over a trial across groups (SEGA Paper Fig 2)
plot_dgkjp <- ggplot(
  df_oddball[df_oddball$ts_trial < 1.9, ],
  aes(x = ts_trial,
      y = pd,
      group = trial_type,
      color = trial_type)) + 
  geom_smooth(formula = y ~ x + poly(x, 5)) +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16)) + 
  labs(x = "trial duration [s]",
       y = "standardized pupil response [z]",
       title = "Oddball Effect on pupillary response") +
    scale_color_manual(values = c("red", "black"))
plot_dgkjp
dev.off()


# get legend only and save it
legend_plot_dgkjp <- get_legend(plot_dgkjp)
as_ggplot(legend_plot_dgkjp)
ggsave(
  "output/plot_legend.tiff",
  device = "tiff",
  width = 3,
  height = 2,
  units = "in",
  compression = "lzw",
  dpi = 800
)

# Plot: Effect of manipulation on pupillary response
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

lmm_interaction <- lmerTest::lmer(
  rpd ~ trial_type * manipulation * group + (1|id),
  data = df_oddball)
anova(lmm_interaction)
confint(contrast(emmeans::emmeans(
  lmm_interaction, ~ group + trial_type + manipulation), method = 'pairwise'))
confint(contrast(emmeans::emmeans(
  lmm_interaction, ~ manipulation + trial_type | group), method = 'pairwise'))

# Pupil response to manipulation and conditions within block
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

# Manipulation phase
df_manip <- df[df$phase %in% c("manipulation_block"), ]

table(df_manip$stimulus)
table(df_manip$manipulation_trial_counter)
with(df_manip, table(stimulus, manipulation_trial_counter))

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
  df_manip[df_manip$ts_trial < 18 & df_manip$stimulus != "baseline", ],
  aes(x = ts_trial,
      y = rpd,
      group = stimulus,
      color = stimulus)) +
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

lmm <- lmer(
  rpd ~ trial_type + (1 | manipulation_trial_counter),
  data = df_manip[df_manip$trial_type != "baseline_start", ])
anova(lmm)
plot(contrast(emmeans(lmm, ~ trial_type), "pairwise"))

# Grip strength ####
hist(df_trial$mean_grip_strength)
df_trial$mean_grip_strength_z<-scale(df_trial$mean_grip_strength)
#pupillary response - SEPR 
lmm<-lmer(scale(rpd)~mean_grip_strength_z*manipulation*stimulus*group*block+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm) #--> no effect of hand grip strength

#pupillary response - BPS
lmm<-lmer(scale(rpd_low)~mean_grip_strength_z*manipulation*stimulus*group*block+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm)
fixef(lmm) #higher hand grip is assoicated with lower bps
emtrends(lmm,~manipulation|group+block,var = 'mean_grip_strength_z') ##--> higher pupil size after manipulation
#the hand grip strength has an effect on rpd only in the TD group

# Visualization
g1 <- ggplot(
  df_baseline[df_baseline$phase == "baseline_calibration" &
                df_baseline$ts_trial < 6, ],
  aes(x = ts_trial,
      y = pd,
      group = stimulus,
      color = stimulus)) +
  geom_smooth() +
  theme_bw() +
  xlab("time (s)") +
  ylab("pupil size (mm)") +
  labs(title = "change of PD during initial baseline between IDs")
g2 <- ggplot(
  df_manip[df_manip$ts_trial < 18 & df_manip$stimulus != "baseline", ],
  aes(x = ts_trial,
      y = rpd,
      group = stimulus,
      color = stimulus)) +
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