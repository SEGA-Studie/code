### --- LOADS preprocessed auditory oddball data and analysis ERPs and pupillometry variables

# SETUP ####
# REQUIRED PACKAGES
require(performance) # marginal + conditional R2
require(lme4) # linear-mixed-effects models
require(lmerTest, warn.conflicts = FALSE) # linear-mixed-effects models
require(emmeans, warn.conflicts = FALSE) # estimated marginal means (EMMs)
require(ggplot2, warn.conflicts = FALSE) # creating graphs
require(dplyr, warn.conflicts = FALSE) # for %>% operator
library(ggpubr, warn.conflicts = FALSE) # ggscatter()-function
library(knitr) # dynamic report generation
library(kableExtra) # table formatting
library(stringr) # for string modifications
library(tidyverse)
library(cowplot, warn.conflicts = FALSE) # get only legend from plot
library(ggpubr, warn.conflicts = FALSE) # save legend from plot

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

# Can be used to skip preprocessing and directly read proprocessed data from .rds file:
df_trial <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ETdata.rds')) # trial pupil data including ALL trials
df <- readRDS(paste0(home_path, project_path,'/data/preprocessed_auditory_df.rds')) # pupil data (all data points for visualizing)
et_erp_subject <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ET_ERP_subject.rds')) # eeg + pupil (subject level, only oddball(rev) blocks)
et_erp_trial <- readRDS(paste0(home_path,project_path,'/data/preprocessed_auditory_ET_ERP_trial.rds')) # eeg + pupil(single-trial, only oddball(rev) blocks)

# changes random intercept to a factor
et_erp_subject$SEGA_ID <- as.factor(et_erp_subject$SEGA_ID)
df_trial$id<-as.factor(df_trial$id) 
et_erp_trial$SEGA_ID <- as.factor(et_erp_trial$SEGA_ID)

# changes covariates to a factor
et_erp_subject$gender <- as.factor(et_erp_subject$gender)
et_erp_trial$gender <- as.factor(et_erp_trial$gender)
df_trial$group <- as.factor(df_trial$group)
df_trial$pitch <- as.factor(df_trial$pitch)
df_trial$block <- as.factor(df_trial$block)

# DATA ANALYSIS ON SUBJECT LEVEL ####
# Sample size
length(unique(et_erp_subject$SEGA_ID))
table(et_erp_subject$group)/8

# Distributions of dependent variables
hist(et_erp_subject$rpd,
     main = "Distribution of rpd (500-1500 ms)",
     xlab = "rpd",
     xlim = c(-0.1, 0.2),
     breaks = 200)

hist(et_erp_subject$rpd_low,
     main = "Distribution of rpd_low (0-250 ms)",
     xlab = "rpd_low",
     xlim = c(2, 5.5),
     breaks = 200)

hist(et_erp_subject$rpd_block,
     main = "Distribution of rpd_block",
     xlab = "rpd_block",
     breaks = 200)

hist(et_erp_subject$MMN_amplitude,
     main = "Distribution of MMN amplitude (100-150 ms)",
     xlab = "MMN amplitude (FC1, FC2, FCz, Fz)",
     ylim = c(0, 30),
     xlim = c(-15, 4.5),
     breaks = 200)

hist(et_erp_subject$P3a_amplitude,
     main = "Distribution of P3a amplitude (150-250 ms)",
     xlab = "P3a amplitude (Cz, FCz)",
     ylim = c(0, 30),
     xlim = c(-5, 15),
     breaks = 200)

hist(et_erp_subject$P3b_amplitude,
     main = "Distribution of P3a amplitude (250-500 ms)",
     xlab = "P3b amplitude (Pz)",
     ylim = c(0, 30),
     xlim = c(-10, 20),
     breaks = 200)

# QUESTIONNAIRE DATA AND DEMOGRAPHICS
exp_groups <- read.csv("exp_groups.csv", header = TRUE)
exp_groups$SEGA_ID <- as.factor(exp_groups$SEGA_ID)
exp_groups$group <- as.factor(exp_groups$group)

## FSK
fsk_import <- read.csv("FSK.csv", header = T, sep = ";", dec = ",", fill = T)
fsk <- fsk_import[c("ID_Studie", "FSK_ges")]
colnames(fsk)[colnames(fsk) == "ID_Studie"] <- "SEGA_ID" # same column name as in exp_groups for matching
fsk$SEGA_ID <- sub(".*SEGA_", "", fsk$SEGA_ID)
fsk$SEGA_ID <- str_replace(fsk$SEGA_ID, "^0+", "")
fsk$SEGA_ID <- as.factor(fsk$SEGA_ID)

## CBCL
cbcl_import <- read.csv("CBCL.csv", header = T, sep = ";", dec = ",", fill = T)
cbcl <- cbcl_import[c("ID_Studie", "CBCL_T_INT", "CBCL_T_EXT", "CBCL_T_GES")]
colnames(cbcl)[colnames(cbcl) == "ID_Studie"] <- "SEGA_ID" # same column name as in exp_groups for matching
cbcl$SEGA_ID <- sub(".*SEGA_", "", cbcl$SEGA_ID)
cbcl$SEGA_ID <- str_replace(cbcl$SEGA_ID, "^0+", "")
cbcl$SEGA_ID <- as.factor(cbcl$SEGA_ID)

## SRS
srs_import <- read.csv("SRS.csv", header = T, sep = ";", dec = ",", fill = T)
srs <- srs_import[c("ID_Studie", "Gesamtwert_N_k_RW")]
colnames(srs)[colnames(srs) == "ID_Studie"] <- "SEGA_ID" # same column name as in exp_groups for matching
### add questionnaire name to column names 
colnames(srs)[colnames(srs) == "Gesamtwert_N_k_RW"] <- "SRS_Gesamtwert_N_k_RW"
colnames(srs)[colnames(srs) == "Gesamt_TW_ASS"] <- "SRS_Gesamt_TW_ASS"
colnames(srs)[colnames(srs) == "Gesamt_TW_AB"] <- "SRS_Gesamt_TW_AB"
srs$SEGA_ID <- sub(".*SEGA_", "", srs$SEGA_ID)
srs$SEGA_ID <- str_replace(srs$SEGA_ID, "^0+", "")
srs$SEGA_ID <- as.factor(srs$SEGA_ID)

### Create df for sample characteristics
age_gender <- et_erp_subject[, c("gender", "age", "SEGA_ID")]
# Combine in one df: sample_characteristics
list_sample <- list(age_gender, exp_groups, fsk, cbcl, srs)
sample_characteristics <- list_sample %>% reduce(full_join, by = "SEGA_ID")
# et_erp_subjects brings 6 rows with identical age + gender per subject -> delete duplicates
sample_characteristics <- sample_characteristics[!duplicated(sample_characteristics), ]
# Only keep sample characteristics for subjects in ERP_ET data
sample_characteristics <- sample_characteristics[sample_characteristics$SEGA_ID %in% et_erp_subject$SEGA_ID, ]


fun_return_descriptives <- function(group){
  group_df <- sample_characteristics[sample_characteristics$group == group, ]
  n <- count(group_df)
  male <- length(which(group_df$gender == "männlich"))
  female <- length(which(group_df$gender == "weiblich"))
  gender_f_m <- paste(female,"/",male)
  age_mean <- round(mean(group_df$age, na.rm = TRUE), digits = 1)
  age_sd <- round(sd(group_df$age, na.rm = TRUE), digits = 1)
  age <- paste(age_mean, "(",age_sd,")" )
  srs_mean <- round(mean(group_df$SRS_Gesamtwert_N_k_RW, na.rm = TRUE), digits = 1)
  srs_sd <- round(sd(group_df$SRS_Gesamtwert_N_k_RW, na.rm = TRUE), digits = 1)
  srs <- paste(srs_mean, "(",srs_sd,")")
  cbcl_mean <- round(mean(group_df$CBCL_T_GES, na.rm = TRUE), digits = 1)
  cbcl_sd <- round(sd(group_df$CBCL_T_GES, na.rm = TRUE), digits = 1)
  cbcl <- paste(cbcl_mean, "(",cbcl_sd,")")
  fsk_mean <- round(mean(group_df$FSK_ges, na.rm = TRUE), digits = 1)
  fsk_sd <- round(sd(group_df$FSK_ges, na.rm = TRUE), digits = 1)
  fsk <- paste(fsk_mean, "(", fsk_sd, ")")
  group_characteristics <- data.frame(n, gender_f_m, age, srs, cbcl, fsk)
  t(group_characteristics)
}
## Descriptive statistics
asd_characteristics <- fun_return_descriptives(group = "ASD")
colnames(asd_characteristics) <- "ASD"
con_characteristics <- fun_return_descriptives(group = "CON")
colnames(con_characteristics) <- "CON"
mhc_characteristics <- fun_return_descriptives(group = "MHC")
colnames(mhc_characteristics) <- "MHC"

fsk_anova <- aov(FSK_ges ~ group, data = sample_characteristics)
fsk_anova_p <- summary(fsk_anova)[[1]][["Pr(>F)"]][[1]]
age_anova <- aov(age ~ group, data = sample_characteristics)
age_anova_p <- summary(age_anova)[[1]][["Pr(>F)"]][[1]]
age_anova_p <- round(age_anova_p, digits = 2)
cbcl_anova <- aov(CBCL_T_GES ~ group, data = sample_characteristics)
cbcl_anova_p <- summary(cbcl_anova)[[1]][["Pr(>F)"]][[1]]
srs_anova <- aov(SRS_Gesamtwert_N_k_RW ~ group, data = sample_characteristics)
srs_anova_p <- summary(srs_anova)[[1]][["Pr(>F)"]][[1]]
# don't know how to extract p-value from chi2-test -> inserted manually
gender_chi2 <- chisq.test(sample_characteristics$gender, sample_characteristics$group)
gender_chi2
p_values <- c(NA, 0.005, age_anova_p, srs_anova_p, cbcl_anova_p, fsk_anova_p)

sample_table <- cbind(asd_characteristics, con_characteristics, mhc_characteristics, p_values)

sample_characteristics_table <- kable(sample_table, caption = "Sample characteristics", digits = 2) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% kable_styling
sample_characteristics_table

## Plot FSK
boxplot_fsk <- boxplot(
  FSK_ges ~ group, data = sample_characteristics, main = "FSK")
stripchart(sample_characteristics$FSK_ges ~ sample_characteristics$group, vertical = TRUE, method = "jitter",
           pch = 19, add = TRUE, col = 1:length(levels(sample_characteristics$group)))
ggplot(sample_characteristics, aes(x = group, y = FSK_ges), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(24, 25.5, 21)) +
  theme_bw() + 
  theme_classic()

## Plot SRS
boxplot_srs <- boxplot(
  SRS_Gesamtwert_N_k_RW ~ group, data = sample_characteristics, main = "SRS")
stripchart(sample_characteristics$SRS_Gesamtwert_N_k_RW ~ sample_characteristics$group, vertical = TRUE, method = "jitter",
           pch = 19, add = TRUE, col = 1:length(levels(sample_characteristics$group)))
ggplot(sample_characteristics, aes(x = group, y = SRS_Gesamtwert_N_k_RW), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(42, 45, 40)) +
  theme_bw() + 
  theme_classic()

## Plot CBCL
boxplot_cbcl <- boxplot(
  CBCL_T_GES ~ group, data = sample_characteristics, main = "CBCL")
stripchart(sample_characteristics$CBCL_T_GES ~ sample_characteristics$group, vertical = TRUE, method = "jitter",
           pch = 19, add = TRUE, col = 1:length(levels(sample_characteristics$group)))
ggplot(sample_characteristics, aes(x = group, y = CBCL_T_GES), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(87, 90, 84)) +
  theme_bw() + 
  theme_classic()

# RESULT 1: MMN AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  scale(MMN_amplitude) ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

emmeans(
  lmm, list(pairwise ~ trial), adjust = "tukey")

emmeans(
  lmm, list(pairwise ~ manipulation|block), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ block|manipulation), adjust = "tukey")
plot(emmeans(
  lmm, list(pairwise ~ manipulation|block), adjust = "tukey"))


# RESULT 2: MMN LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  scale(MMN_latency) ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

contrast(emmeans(lmm, ~ trial * manipulation|group), "pairwise")
emmip(lmm, ~ manipulation |group * trial, CI = T)

# RESULT 3: P3A AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  scale(P3a_amplitude) ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

emm <- emmeans(
  lmm, list(pairwise ~ trial ), adjust = "tukey")
emmip(lmm, ~ trial, CIs = T)
contrast(emmeans(lmm, ~ trial|group), "pairwise")
contrast(emmeans(lmm, ~ group|trial), "pairwise")
emmip(lmm, ~ trial|group, CI = T)
contrast(emmeans(lmm, ~ manipulation|block), "pairwise")
emmip(lmm, ~ manipulation|block, CI = T)
contrast(emmeans(lmm, ~ group * trial|block), "pairwise")
emmip(lmm, ~ trial|group * block, CI = T)

# RESULT 4: P3A LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  scale(P3a_latency) ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

emmeans(
  lmm, list(pairwise ~ trial), adjust = "tukey")
plot(emmeans(
  lmm, list(pairwise ~ manipulation), adjust = "tukey"))
emmip(lmm, ~ manipulation, CIs = T)

# RESULT 5: P3B AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  scale(P3b_amplitude) ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

emmeans(
  lmm, list(pairwise ~ trial), adjust = "tukey")

# RESULT 6: P3B LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  scale(P3b_latency) ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

emmeans(
  lmm, list(pairwise ~ trial), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ group|trial), adjust = "tukey")
emmip(lmm, ~ group | trial, CI = T)
emmeans(
  lmm, list(pairwise ~ group|trial * manipulation), adjust = "tukey")
emmip(lmm, ~ group | trial * manipulation, CI = T)


emmeans(
  lmm, list(pairwise ~ manipulation|trial * block), adjust = "tukey")
emmip(lmm, ~ manipulation| trial * block, CI = T)

# RESULT 7: SEPR ON SUBJECT LEVEL
lmm <- lmer(
  scale(rpd) ~ trial * manipulation * group * block + (1|SEGA_ID),
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

emmeans(
  lmm, list(pairwise ~ trial), adjust = "tukey")

# RESULT 6: BPS ON SUBJECT LEVEL
lmm <- lmer(
  scale(rpd_low) ~ trial * manipulation * group * block + (1|SEGA_ID),
  data = et_erp_subject)
anova(lmm) 
r2_nakagawa(lmm) 

emmeans(
  lmm, list(pairwise ~ manipulation), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ manipulation|block), adjust = "tukey")
emmip(lmm, ~ manipulation | block, CIs = T)
emmip(lmm, ~ manipulation, CIs = T)

# RESULT 7: DOES PUPUIL DATA PREDICT ERPs?
lmm <- lmer(
  scale(MMN_amplitude) ~ rpd * rpd_low * trial *  manipulation * group + (1|SEGA_ID),
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

lmm <- lmer(
  scale(P3a_amplitude) ~ rpd * rpd_low * trial *  manipulation * group + (1|SEGA_ID),
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

lmm <- lmer(
  scale(P3b_amplitude) ~ rpd * rpd_low * trial *  manipulation * group + (1|SEGA_ID),
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

# RESULT 8: ASSOCIATION BETWEEN PUPIL DATA AND ERPs ON SUBJECT LEVEL
crl <- cor(et_erp_subject[et_erp_subject$trial == "Oddball", c(
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

crl <- cor(et_erp_subject[et_erp_subject$trial == "Standard", c(
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
# RESULT 8: MMN AMPLITUDE ON TRIAL LEVEL
lmm <- lmer(
  scale(MMN_amplitude) ~ trial * manipulation * group * block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

emmeans(
  lmm, list(pairwise ~ trial), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ manipulation|block), adjust = "tukey")

# RESULT 9: MMN LATENCY ON TRIAL LEVEL
lmm <- lmer(
  scale(MMN_latency) ~ trial * manipulation * group * block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

emmeans(
  lmm, list(pairwise ~ trial), adjust = "tukey")

# RESULT 10: P3A AMPLITUDE ON TRIAL LEVEL
lmm <- lmer(
  scale(P3a_amplitude) ~ trial * manipulation * group * block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ trial | group), "pairwise")
contrast(emmeans(lmm, ~ group | trial), "pairwise")
plot(emmeans(lmm, ~ trial|group))
contrast(emmeans(lmm, ~ manipulation | block), "pairwise")
plot(emmeans(lmm, ~ manipulation|block))

# RESULT 11: P3A LATENCY ON TRIAL LEVEL
lmm <- lmer(
  scale(P3a_latency) ~ trial * manipulation * group * block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emmeans(
  lmm, list(pairwise ~ trial), adjust = "tukey")
emmip(lmm, ~block|group * manipulation|trial, CIs = T)

# RESULT 12: SEPR ON TRIAL LEVEL
lmm <- lmer(
  scale(rpd) ~ trial * manipulation * group * block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

emmeans(
  lmm, list(pairwise ~ trial), adjust = "tukey")

lmm <- lmer(
  scale(rpd) ~ trial * manipulation * group * block * oddball_trial_counter + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

plot(emtrends(lmm,~ group, var = 'oddball_trial_counter'))
emtrends(lmm,~ group, var = 'oddball_trial_counter')

ggplot(
  et_erp_trial,
  aes(
    x = trial_number_in_block,
    y = scale(rpd),
    group = trial,
    color = trial)) +
  geom_smooth() +
  theme_bw() + 
  ggtitle("SEPR (rpd) over block: A group comparison.") +
  facet_grid(cols = vars(group))

ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev") &
             df_trial$trial_number_in_block >= 10 & is.finite(
               df_trial$rpd), ],
  aes(x = trial_number_in_block, y = rpd, group = trial, color = trial)) +
  geom_smooth() +
  facet_grid(rows = vars(block), cols = vars(manipulation)) + theme_bw()

ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev") &
             is.finite(df_trial$rpd), ],
  aes(x = as.factor(trial_number_in_block), y = rpd)) +
  geom_boxplot() +
  facet_grid(
    rows = vars(block),
    cols = vars(manipulation, trial)) +
  theme_bw()

bin_size <- 20
df_trial <- df_trial %>%
  mutate(trial_number_in_block_binned = factor(trial_number_in_block%/%bin_size*20))
ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev"), ],
  aes(x = trial_number_in_block_binned,
      y = scale(rpd))) + 
  geom_boxplot() + 
  facet_grid(
    rows = vars(block),
    cols = vars(manipulation)) +
  theme_bw() +
  ggtitle("rpd during blocks, split by forward + reverse \nbefore + after manipulation")

et_erp_trial$binned_trial_number <- cut(
  et_erp_trial$oddball_trial_counter, c(
    1,  25,  50, 75,  100, 125,  150, 175, 200,  225, 250, 275, 300, 325, 350, 375, 400))
ggplot(et_erp_trial) + geom_boxplot(aes(binned_trial_number, rpd)) + 
  facet_grid(cols = vars(group))

# RESULT 13: BPS ON TRIAL LEVEL
lmm <- lmer(
  scale(rpd_low) ~ trial * manipulation * group * block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

emmeans(
  lmm, list(pairwise ~ manipulation), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ group|manipulation), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ group|block), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ block|group), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ manipulation|block), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ group|manipulation), adjust = "tukey")
emmeans(
  lmm, list(pairwise ~ manipulation|group), adjust = "tukey")
emmip(lmm, ~ manipulation|group, CIs = T)

ggplot(
  et_erp_trial,
  aes(
    x = trial_number_in_block,
    y = scale(rpd_low),
    group = trial,
    color = trial)) +
  geom_smooth() +
  theme_bw() + 
  ggtitle("BPS (rpd_low) over block: A group comparison.") +
  facet_grid(cols = vars(group))

et_erp_trial$binned_trial_number <- cut(
  et_erp_trial$oddball_trial_counter, c(
    1,  25,  50, 75,  100, 125,  150, 175, 200,  225, 250, 275, 300, 325, 350, 375, 400))
ggplot(et_erp_trial[et_erp_trial$trial == "standard", ]) + geom_boxplot(aes(binned_trial_number, rpd_low)) + 
  facet_grid(cols = vars(group))

lmm <- lmer(
  scale(rpd_low) ~ trial * manipulation * group * block * oddball_trial_counter + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

# RESULT 14: ASSOCIATION BETWEEN BPS AND SEPR ON TRIAL LEVEL
ggscatter(et_erp_trial[et_erp_trial$trial == "oddball", ], 
          x = "z_rpd_low",
          y = "z_rpd",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "scaled BPS (rpd_low)",
          ylab = "scaled SEPR (rpd)",
          title = "Association between SEPR + BPS in oddball trials")

ggscatter(et_erp_trial[et_erp_trial$trial == "standard", ], 
          x = "z_rpd_low",
          y = "z_rpd",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "scaled BPS (rpd_low)",
          ylab = "scaled SEPR (rpd)",
          title = "Association between SEPR + BPS in standard trials")

lmm <- lmer(
  scale(rpd) ~ scale(rpd_low) * trial * manipulation * group * block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

emtrends(lmm, ~ group, var = "rpd_low")
plot(emtrends(lmm, ~ group, var = "rpd_low"))

# RESULT 15: ASSOCIATION BETWEEN PUPIL DATA AND ERPs ON TRIAL LEVEL
crl <- cor(et_erp_trial[et_erp_trial$trial == "oddball", c(
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

crl <- cor(et_erp_trial[et_erp_trial$trial == "standard", c(
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

# RESULT 16: DOES PUPIL DATA PREDICT ERPs?
## MMN amplitude
lmm <- lmer(
  scale(MMN_amplitude) ~ rpd * rpd_low * trial *  manipulation * group + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

summary(lmm)
emtrends(lmm, ~ group, var = c("rpd"))

## P3a amplitude
lmm <- lmer(
  scale(P3a_amplitude) ~ rpd * rpd_low * trial *  manipulation * group + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm, ~ manipulation | group, var = c("rpd_low"))
plot(emtrends(lmm, ~ manipulation | group, var = c("rpd_low")))

emtrends(lmm, ~ manipulation, var = c("rpd_low"))

# RESULT 17: EXPLORATORY ANALYSIS OF MODEL FIT FOR "trial_number_in_block" ON BPS
linear_fit <- lmer(
  scale(rpd_low) ~ trial * manipulation * group * trial_number_in_block * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial, REML = F)
anova(linear_fit)

quadratic_fit <- lmer(
  scale(rpd_low) ~ trial * manipulation * group * poly(trial_number_in_block, 2) * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial, REML = F)
anova(quadratic_fit)
r2_nakagawa(quadratic_fit) 

cubic_fit <- lmer(
  scale(rpd_low) ~ trial * manipulation * group * poly(trial_number_in_block, 3) * block + (1|SEGA_ID) + age + gender,
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
  scale(rpd) ~ trial * manipulation * group * trial_number_in_block * block+ (1|SEGA_ID) + pitch + age + gender,
  data = et_erp_trial, REML = F)
anova(linear_fit)
r2_nakagawa(linear_fit) 

quadratic_fit <- lmer(
  scale(rpd) ~ trial * manipulation * group * poly(trial_number_in_block, 2) * block + (1|SEGA_ID) + pitch + age + gender,
  data = et_erp_trial, REML = F)
anova(quadratic_fit)
r2_nakagawa(quadratic_fit) 

cubic_fit <- lmer(
  scale(rpd) ~ trial * manipulation * group * poly(trial_number_in_block, 3) * block + (1|SEGA_ID) + pitch + age + gender,
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

# RESULT 19: EXPLORATORY ANALYSIS OF GRIP STRENGTH EFFECT (Z-VALUES)
lmm <- lmer(
  scale(MMN_amplitude) ~ trial * group * block * z_handdynamometer + (1|SEGA_ID) + age + gender,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm)

lmm <- lmer(
  scale(MMN_latency) ~ trial * group * block * z_handdynamometer+ (1|SEGA_ID) + gender + age,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm)

lmm <- lmer(
  scale(P3a_amplitude) ~ trial * group * block * z_handdynamometer + (1|SEGA_ID) + gender + age,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm) 

lmm <- lmer(
  scale(P3a_latency) ~ trial * group * block * z_handdynamometer+ (1|SEGA_ID) + gender + age,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm) 

lmm <- lmer(
  scale(rpd) ~ trial * group * block * z_handdynamometer + (1|SEGA_ID),
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm) 

lmm <- lmer(
  scale(rpd_low) ~ trial * group * block * z_handdynamometer + (1|SEGA_ID),
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm) 
r2_nakagawa(lmm)

# New df et_erp_subject shows same result as df_trial, is valide ####
lmm <- lmer(
  scale(rpd) ~ trial * manipulation * block + trial_number_in_block + (1 | id),
  data = df_trial[
    df_trial$phase %in% c("oddball_block", "oddball_block_rev"), ])
anova(lmm)

lmm <- lmer(
  scale(rpd) ~ trial * manipulation * block + trial_number_in_block + (1 | SEGA_ID),
  data = et_erp_trial)
anova(lmm)

# ANALYSIS OF df (ONLY PUPIL DATA-ALL DATA POINTS) ####
# Baseline phase
df_baseline <- df[df$phase %in% c("baseline", "baseline_calibration"), ]
df_baseline <- df_baseline[is.finite(df_baseline$pd), ]

# Number of pupil data per baseline trial (white, black and grey slide)
with(df_baseline[df_baseline$phase == "baseline_calibration", ], table(trial))
# Check: Block counter should be in line with aforementioned table
with(df_baseline, table(trial, block_counter))
with(df_baseline, by(
  timestamp_exp, interaction(baseline_trial_counter, phase, id),
  mean, na.rm = TRUE))

# mean_pd (subject baseline) is defined as mean of baseline trial during calibration phase
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
      y = pd,
      group = group,
      color = group)) +
  geom_smooth() +
  theme_bw() +
  xlab("trial duration [s]") +
  ylab("pupil dilation [mm]") +
  labs(title = "LAPR for group and trial") +
  facet_wrap(~ factor(trial,
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

# Plot: Pupillary response across experiment
ggplot(
  df_oddball,
  aes(x = trial_number_oddballphase, y = rpd)) +
  geom_smooth(method = "lm") + facet_wrap(~trial + manipulation)
dev.off()

# Plot: Density of manipulation per condition
ggplot(df_oddball,
       aes(rpd, fill = interaction(manipulation))) +
  geom_density(alpha = 0.2) + facet_wrap(~trial)
dev.off()

# Plot: Pupil response over a trial
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

# Plot: Pupil response over a trial across groups
plot_dgkjp <- ggplot(
  df_oddball[df_oddball$ts_trial < 1.8, ],
  aes(x = ts_trial,
      y = scale(rpd),
      group = trial_type,
      color = trial_type,
      linetype = trial_type,)) +
  geom_smooth() +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14))
  labs(x = "trial duration (s)",
       y = "standardized pupil response [z]",
       title = "Oddball Effect on pupillary response")
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

lmm <- lmer(
  rpd ~ trial_type + (1 | manipulation_trial_counter),
  data = df_manip[df_manip$trial_type != "baseline_start", ])
anova(lmm)
plot(contrast(emmeans(lmm, ~ trial_type), "pairwise"))

# Grip strength ####
hist(df_trial$mean_grip_strength)
df_trial$mean_grip_strength_z<-scale(df_trial$mean_grip_strength)
#pupillary response - SEPR 
lmm<-lmer(scale(rpd)~mean_grip_strength_z*manipulation*trial*group*block+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
anova(lmm) #--> no effect of hand grip strength

#pupillary response - BPS
lmm<-lmer(scale(rpd_low)~mean_grip_strength_z*manipulation*trial*group*block+(1|id),data=df_trial[df_trial$phase %in% c('oddball_block','oddball_block_rev'),])
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