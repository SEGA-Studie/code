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

# DATA ANALYSIS ON SUBJECT LEVEL ####
# Sample size
length(unique(et_erp_subject$SEGA_ID))
table(et_erp_subject$group)/8

excluded_ids <- c("117", "152", "36") # to balance age between groups
et_erp_subject <- subset(et_erp_subject, ! (SEGA_ID %in% excluded_ids))
et_erp_trial <- subset(et_erp_trial, ! (SEGA_ID %in% excluded_ids))

# Sample Size after removing subjects (age balance)
length(unique(et_erp_subject$SEGA_ID))
table(et_erp_subject$group)/8

# Distributions of dependent variables
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

fun_return_descriptives <- function(group){
  group_df <- et_erp_subject[et_erp_subject$group == group, c("gender", "age", "SRS", "SCQ", "CBCL", "YSR", "SP2", "z_handdynamometer", "SEGA_ID", "verbal_IQ", "non_verbal_IQ")]
  n <- length(unique(group_df$SEGA_ID))
  male <- (length(which(group_df$gender == "männlich"))/8)
  female <- (length(which(group_df$gender == "weiblich"))/8)
  gender_f_m <- paste(female,"/",male)
  age_mean <- round(mean(group_df$age, na.rm = TRUE), digits = 1)
  age_sd <- round(sd(group_df$age, na.rm = TRUE), digits = 1)
  age <- paste(age_mean, "(",age_sd,")" )
  SRS_mean <- round(mean(group_df$SRS, na.rm = TRUE), digits = 1)
  SRS_sd <- round(sd(group_df$SRS, na.rm = TRUE), digits = 1)
  SRS <- paste(SRS_mean, "(",SRS_sd,")")
  CBCL_mean <- round(mean(group_df$CBCL, na.rm = TRUE), digits = 1)
  CBCL_sd <- round(sd(group_df$CBCL, na.rm = TRUE), digits = 1)
  CBCL <- paste(CBCL_mean, "(",CBCL_sd,")")
  SCQ_mean <- round(mean(group_df$SCQ, na.rm = TRUE), digits = 1)
  SCQ_sd <- round(sd(group_df$SCQ, na.rm = TRUE), digits = 1)
  SCQ <- paste(SCQ_mean, "(", SCQ_sd, ")")
  YSR_mean <- round(mean(group_df$YSR, na.rm = TRUE), digits = 1)
  YSR_sd <- round(sd(group_df$YSR, na.rm = TRUE), digits = 1)
  YSR <- paste(YSR_mean, "(", YSR_sd, ")")
  SP2_mean <- round(mean(group_df$SP2, na.rm = TRUE), digits = 1)
  SP2_sd <- round(sd(group_df$SP2, na.rm = TRUE), digits = 1)
  SP2 <- paste(SCQ_mean, "(", SP2_sd, ")")
  grip_strength_mean <- round(mean(group_df$z_handdynamometer, na.rm = TRUE), digits = 1)
  grip_strength_sd <- round(sd(group_df$z_handdynamometer, na.rm = TRUE), digits = 1)
  grip_strength <- paste(grip_strength_mean, "(", grip_strength_sd, ")")
  verbal_IQ_mean <- round(mean(group_df$verbal_IQ, na.rm = TRUE), digits = 1)
  verbal_IQ_sd <- round(sd(group_df$verbal_IQ, na.rm = TRUE), digits = 1)
  verbal_IQ <- paste(verbal_IQ_mean, "(", verbal_IQ_sd, ")")
  non_verbal_IQ_mean <- round(mean(group_df$non_verbal_IQ, na.rm = TRUE), digits = 1)
  non_verbal_IQ_sd <- round(sd(group_df$non_verbal_IQ, na.rm = TRUE), digits = 1)
  non_verbal_IQ <- paste(non_verbal_IQ_mean, "(", non_verbal_IQ_sd, ")")
  group_description <- data.frame(
    n,
    gender_f_m,
    age,
    SRS,
    CBCL,
    SCQ,
    YSR,
    SP2,
    grip_strength,
    verbal_IQ,
    non_verbal_IQ)
  t(group_description)
}

## Descriptive statistics
asd_description <- fun_return_descriptives(group = "ASD")
colnames(asd_description) <- "ASD"
con_description <- fun_return_descriptives(group = "CON")
colnames(con_description) <- "CON"
mhc_description <- fun_return_descriptives(group = "MHC")
colnames(mhc_description) <- "MHC"

SCQ_anova <- aov(SCQ ~ group, data = et_erp_subject)
SCQ_anova_p <- summary(SCQ_anova)[[1]][["Pr(>F)"]][[1]]
age_anova <- aov(age ~ group, data = et_erp_subject)
age_anova_p <- summary(age_anova)[[1]][["Pr(>F)"]][[1]]
age_anova_p <- age_anova_p
CBCL_anova <- aov(CBCL ~ group, data = et_erp_subject)
CBCL_anova_p <- summary(CBCL_anova)[[1]][["Pr(>F)"]][[1]]
SRS_anova <- aov(SRS ~ group, data = et_erp_subject)
SRS_anova_p <- summary(SRS_anova)[[1]][["Pr(>F)"]][[1]]
YSR_anova <- aov(YSR ~ group, data = et_erp_subject)
YSR_anova_p <- summary(YSR_anova)[[1]][["Pr(>F)"]][[1]]
SP2_anova <- aov(SP2 ~ group, data = et_erp_subject)
SP2_anova_p <- summary(SP2_anova)[[1]][["Pr(>F)"]][[1]]
grip_strength_anova <- aov(z_handdynamometer ~ group, data = et_erp_subject)
grip_strength_anova_p <- summary(grip_strength_anova)[[1]][["Pr(>F)"]][[1]]
verbal_IQ_anova <- aov(verbal_IQ ~ group, data = et_erp_subject)
verbal_IQ_anova_p <- summary(verbal_IQ_anova)[[1]][["Pr(>F)"]][[1]]
non_verbal_IQ_anova <- aov(non_verbal_IQ ~ group, data = et_erp_subject)
non_verbal_IQ_anova_p <- summary(non_verbal_IQ_anova)[[1]][["Pr(>F)"]][[1]]
gender_chi2 <- chisq.test(et_erp_subject$gender, et_erp_subject$group)
gender_chi2_p <- gender_chi2$p.value
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
  non_verbal_IQ_anova_p)

p_value <- c()
for (p in p_values){
  ifelse (p < 0.001, p <- "< 0.001",
          ifelse(p < 0.01, p <- "< 0.01",
                 ifelse(p < 0.05, p <- "< 0.05", p <- p )))
  p_value <- c(p_value, p)
}

sample_table <- cbind(asd_description, con_description, mhc_description, p_value)
sample_description_table <- kable(sample_table, caption = "Sample description", digits = 4) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% kable_styling
sample_description_table

## Plot SCQ
ggplot(et_erp_subject, aes(x = group, y = SCQ), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(24, 25.5, 21)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("SCQ") +
  theme(plot.title = element_text(face = "bold"))

## Plot SRS
ggplot(et_erp_subject, aes(x = group, y = SRS), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(42, 45, 40)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("SRS") +
  theme(plot.title = element_text(face = "bold"))

## Plot CBCL
ggplot(et_erp_subject, aes(x = group, y = CBCL), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(87, 90, 84)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("CBCL") +
  theme(plot.title = element_text(face = "bold"))

## Plot: YSR
ggplot(et_erp_subject, aes(x = group, y = YSR), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(87, 90, 84)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("YSR") +
  theme(plot.title = element_text(face="bold"))

## Plot: SP2-Auditiv
ggplot(et_erp_subject, aes(x = group, y = SP2), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(43, 46, 40)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("SP2-Auditiv") +
  theme(plot.title = element_text(face="bold"))

## Plot: Age
age_df <- et_erp_subject[!duplicated(et_erp_subject$SEGA_ID), c("SEGA_ID", "group", "age")]
ggplot(age_df, aes(x = group, y = age), col = "group") + 
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
IQ_df <- et_erp_subject[!duplicated(et_erp_subject$SEGA_ID), c("SEGA_ID", "group", "verbal_IQ", "non_verbal_IQ")]
ggplot(IQ_df, aes(x = group, y = verbal_IQ), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(135, 140, 130)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("IQ-verbal") +
  theme(plot.title = element_text(face="bold"))

## Plot: IQ-non-verbal
ggplot(IQ_df, aes(x = group, y = non_verbal_IQ), col = "group") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(135, 140, 130)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("IQ-non-verbal") +
  theme(plot.title = element_text(face="bold"))

# RESULT 1: MMN AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  z_MMN_amplitude ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

contrast(emmeans(lmm, ~ trial), "pairwise")
confint(contrast(emmeans(lmm, ~ trial), "pairwise"))

contrast(emmeans(lmm, ~ manipulation|block), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "pairwise"))
emmip(lmm, ~ manipulation | block, CIs = T)

# RESULT 2: MMN LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  z_MMN_latency ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

# RESULT 3: P3A AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  z_P3a_amplitude ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ trial|group), "pairwise")
confint(contrast(emmeans(lmm, ~ trial|group), "pairwise"))
emmip(lmm, ~ trial|group, CI = T)

emmeans(lmm, ~ manipulation|block)
contrast(emmeans(lmm, ~ manipulation|block), "pairwise")
emmip(lmm, ~ manipulation|block, CI = T)

# RESULT 4: P3A LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  z_P3a_latency ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ trial), "pairwise")
confint(contrast(emmeans(lmm, ~ trial), "pairwise"))

# RESULT 5: P3B AMPLITUDE ON SUBJECT LEVEL
lmm <- lmer(
  z_P3b_amplitude ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ trial), "pairwise")
confint(contrast(emmeans(lmm, ~ trial), "pairwise"))
emmip(lmm, ~ trial, CI = T)

contrast(emmeans(lmm, ~ manipulation), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation), "pairwise"))
emmip(lmm, ~ manipulation, CI = T)

# RESULT 6: P3B LATENCY ON SUBJECT LEVEL
lmm <- lmer(
  z_P3b_latency ~ trial * manipulation * group * block + (1|SEGA_ID) + gender + age,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ group|manipulation*trial), "pairwise")
confint(contrast(emmeans(lmm, ~ group|manipulation*trial), "pairwise"))
emmip(lmm, ~ manipulation|group*trial, CI = T)

# RESULT 7: SEPR ON SUBJECT LEVEL
lmm <- lmer(
  z_rpd ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ trial), "pairwise")
confint(contrast(emmeans(lmm, ~ trial), "pairwise"))

# RESULT 6: BPS ON SUBJECT LEVEL
lmm <- lmer(
  z_rpd_low ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm) 
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ manipulation|group), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation|group), "pairwise"))
emmip(lmm, ~ manipulation | group, CIs = T)

contrast(emmeans(lmm, ~ manipulation|block), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "pairwise"))
emmip(lmm, ~ manipulation | block, CIs = T)

# Does BPS predict SEPR?
lmm <- lmer(
  z_rpd ~ z_rpd_low * trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm, ~ group, var = "z_rpd_low")
contrast(emtrends(lmm, ~ group, var = "z_rpd_low"))
confint(contrast(emtrends(lmm, ~ group, var = "z_rpd_low")))

# RESULT 8: CORRELATION: PUPIL DATA AND ERPs ON SUBJECT LEVEL
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

# RESULT 7: DOES PUPUIL DATA PREDICT ERPs?
lmm <- lmer(
  z_MMN_amplitude ~ z_rpd * z_rpd_low * trial *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm, ~ z_rpd_low * z_rpd | trial, var = 'z_rpd', at = list(z_rpd_low = c(-2,-1,0,1,2)))
contrast(emtrends(lmm, ~ z_rpd_low * z_rpd | trial, var = 'z_rpd', at = list(z_rpd_low = c(-2,-1,0,1,2))))

lmm <- lmer(
  z_P3a_amplitude ~ z_rpd * z_rpd_low * trial *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

### post-hoc: SEPR x manipulation x group
emtrends(lmm, ~ z_rpd  * manipulation | group * trial, var = 'z_rpd')
contrast(emtrends(lmm, ~ z_rpd  * manipulation | group * trial, var = 'z_rpd'))

### post-hoc: BPS x trial x group
emtrends(lmm, ~ z_rpd_low * trial * group, var = 'z_rpd_low')
contrast(emtrends(lmm, ~ z_rpd_low * trial * group, var = 'z_rpd_low'))
confint(contrast(emtrends(lmm, ~ z_rpd_low * trial * group, var = 'z_rpd_low')))

lmm <- lmer(
  z_P3b_amplitude ~ z_rpd * z_rpd_low * trial *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_subject)
anova(lmm)
r2_nakagawa(lmm)

# DATA ANALYSIS ON TRIAL LEVEL ####
# Distributions of dependent variables
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

# RESULT 8: MMN AMPLITUDE ON TRIAL LEVEL
lmm <- lmer(
  z_MMN_amplitude ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ trial), "pairwise")
confint(contrast(emmeans(lmm, ~ trial), "pairwise"))
emmip(lmm, ~ trial, CIs = T)

contrast(emmeans(lmm, ~ manipulation|block), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "pairwise"))
emmip(lmm, ~ manipulation | block, CIs = T)

# RESULT 9: MMN LATENCY ON TRIAL LEVEL
lmm <- lmer(
  z_MMN_latency ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ trial), "pairwise")
confint(contrast(emmeans(lmm, ~ trial), "pairwise"))

# RESULT 10: P3A AMPLITUDE ON TRIAL LEVEL
lmm <- lmer(
  z_P3a_amplitude ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ trial | group))
emmeans(lmm, ~ trial | group)
contrast(emmeans(lmm, ~ manipulation|block), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "pairwise"))
emmip(lmm, ~ manipulation | block, CIs = T)

# RESULT 11: P3A LATENCY ON TRIAL LEVEL
lmm <- lmer(
  z_P3a_latency ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

contrast(emmeans(lmm, ~ trial), "pairwise")
confint(contrast(emmeans(lmm, ~ trial), "pairwise"))
emmip(lmm, ~block|group * manipulation|trial, CIs = T)

# RESULT 12: SEPR ON TRIAL LEVEL
lmm <- lmer(
  z_rpd ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

emmeans(lmm, list(pairwise ~ trial))

# EXPLORATORY: SEPR OVER TASK
lmm <- lmer(
  z_rpd ~ trial * manipulation * group * block * oddball_trial_counter + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm,~ manipulation, var = 'oddball_trial_counter')
emtrends(lmm,~ group * group * block, var = 'oddball_trial_counter')

et_erp_trial$binned_trial_number <- cut(
  et_erp_trial$oddball_trial_counter, c(
    1,  103, 206, 309, 412))
ggplot(et_erp_trial[et_erp_trial$trial == "oddball", ]) + geom_boxplot(aes(binned_trial_number, z_rpd)) + 
  facet_grid(cols = vars(group)) + 
  ggtitle("SEPR in oddball trials during the task, \nsplit by group") +
  theme(plot.title = element_text(face = "bold", size = 14))

# EXPLORATORY: SEPR OVER BLOCK
lmm <- lmer(
  z_rpd ~ trial * manipulation * group * block * trial_number_in_block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm,~ manipulation, var = 'trial_number_in_block')

bin_size <- 20
df_trial <- df_trial %>%
  mutate(trial_number_in_block_binned = factor(trial_number_in_block%/%bin_size*20))
ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev") & df_trial$trial == "oddball", ],
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
ggplot(et_erp_trial[et_erp_trial$trial == "oddball", ]) + geom_boxplot(aes(binned_trial_number, z_rpd)) + 
  facet_grid(cols = vars(group)) + 
  ggtitle("SEPR in oddball trials during blocks, split by group") +
  theme(plot.title = element_text(face = "bold", size = 14))

# RESULT 13: BPS ON TRIAL LEVEL
lmm <- lmer(
  z_rpd_low ~ trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

contrast(emmeans(lmm, ~ manipulation|block), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation|block), "pairwise"))
emmip(lmm, ~ manipulation | block, CIs = T)

contrast(emmeans(lmm, ~ manipulation|group), "pairwise")
confint(contrast(emmeans(lmm, ~ manipulation|group), "pairwise"))
emmip(lmm, ~ manipulation | group, CIs = T)

contrast(emmeans(lmm, ~ trial), "pairwise")
confint(contrast(emmeans(lmm, ~ trial), "pairwise"))

# EXPLORATORY: BPS OVER TASK
lmm <- lmer(
  z_rpd_low ~ trial * manipulation * group * block * oddball_trial_counter + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm,~ manipulation * group|block, var = 'oddball_trial_counter')

et_erp_trial$binned_trial_number <- cut(
  et_erp_trial$oddball_trial_counter, c(
    1,  103, 206, 309, 412))
ggplot(et_erp_trial[et_erp_trial$trial == "standard", ]) + geom_boxplot(aes(binned_trial_number, z_rpd_low)) + 
  facet_grid(cols = vars(group)) + 
  ggtitle("BPS in standard trials during the task, split by group") +
  theme(plot.title = element_text(face = "bold", size = 14))

# EXPLORATORY: BPS OVER BLOCK
lmm <- lmer(
  z_rpd_low ~ trial * manipulation * group * block * trial_number_in_block + (1|SEGA_ID),
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm,~ manipulation * group|block, var = 'trial_number_in_block')

bin_size <- 20
df_trial <- df_trial %>%
  mutate(trial_number_in_block_binned = factor(trial_number_in_block%/%bin_size*20))
ggplot(
  df_trial[df_trial$phase %in% c("oddball_block", "oddball_block_rev") & df_trial$trial == "standard", ],
  aes(x = trial_number_in_block_binned,
      y = scale(rpd_low))) + 
  geom_boxplot() + 
  facet_grid(
    rows = vars(block),
    cols = vars(manipulation)) +
  theme_bw() +
  ggtitle("BPS in standard trials during blocks, \nsplit by forward + reverse and before + after manipulation") +
  theme(plot.title = element_text(face = "bold", size = 14))

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

## Does BPS predict SEPR?
lmm <- lmer(
  z_rpd ~ z_rpd_low * trial * manipulation * group * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm) 

### post-hoc: BPS x manipulation x block
emtrends(lmm, ~ manipulation * block, var = "z_rpd_low")
contrast(emtrends(lmm, ~ manipulation * block, var = "z_rpd_low"))
confint(contrast(emtrends(lmm, ~ manipulation * block, var = "z_rpd_low")))
plot(emtrends(lmm, ~ manipulation * block, var = "z_rpd_low"))

emtrends(lmm, ~ block|manipulation, var = "z_rpd_low")
contrast(emtrends(lmm, ~ block|manipulation, var = "z_rpd_low"))
confint(contrast(emtrends(lmm, ~ block|manipulation, var = "z_rpd_low")))

### post-hoc: BPS x manipulation x group
emtrends(lmm, ~ group|manipulation, var = "z_rpd_low")
contrast(emtrends(lmm, ~ group|manipulation, var = "z_rpd_low"))
confint(contrast(emtrends(lmm, ~ group|manipulation, var = "z_rpd_low")))
plot(emtrends(lmm, ~ group|manipulation, var = "z_rpd_low"))

# RESULT 15: Correlation: Pupil-ERP
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
  z_MMN_amplitude ~ z_rpd * z_rpd_low * trial *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm, ~ z_rpd_low * z_rpd | group, var = 'z_rpd', at = list(z_rpd_low = c(-2,-1,0,1,2)))
contrast(emtrends(lmm, ~ z_rpd_low * z_rpd | group, var = 'z_rpd', at = list(z_rpd_low = c(-2,-1,0,1,2))))
confint(contrast(emtrends(lmm, ~ z_rpd_low * z_rpd | group, var = 'z_rpd', at = list(z_rpd_low = c(-2,-1,0,1,2)))))

## P3a amplitude
lmm <- lmer(
  z_P3a_amplitude ~ z_rpd * z_rpd_low * trial *  manipulation * group + (1|SEGA_ID) + age + gender,
  data = et_erp_trial)
anova(lmm)
r2_nakagawa(lmm)

emtrends(lmm, ~ manipulation | group, var = c("z_rpd_low"))
contrast(emtrends(lmm, ~ manipulation | group, var = c("z_rpd_low")))
confint(contrast(emtrends(lmm, ~ manipulation | group, var = c("z_rpd_low"))))

emtrends(lmm, ~ z_rpd_low * z_rpd | manipulation, var = 'z_rpd', at = list(z_rpd_low = c(-2,-1,0,1,2)))
contrast(emtrends(lmm, ~ z_rpd_low * z_rpd | manipulation, var = 'z_rpd', at = list(z_rpd_low = c(-2,-1,0,1,2))))

# RESULT 17: EXPLORATORY ANALYSIS OF MODEL FIT FOR "trial_number_in_block" ON BPS
linear_fit <- lmer(
  z_rpd_low ~ trial * manipulation * group * trial_number_in_block * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial, REML = F)
anova(linear_fit)

quadratic_fit <- lmer(
  z_rpd_low ~ trial * manipulation * group * poly(trial_number_in_block, 2) * block + (1|SEGA_ID) + age + gender,
  data = et_erp_trial, REML = F)
anova(quadratic_fit)
r2_nakagawa(quadratic_fit) 

cubic_fit <- lmer(
  z_rpd_low ~ trial * manipulation * group * poly(trial_number_in_block, 3) * block + (1|SEGA_ID) + age + gender,
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
  z_rpd ~ trial * manipulation * group * trial_number_in_block * block+ (1|SEGA_ID) + pitch + age + gender,
  data = et_erp_trial, REML = F)
anova(linear_fit)
r2_nakagawa(linear_fit) 

quadratic_fit <- lmer(
  z_rpd ~ trial * manipulation * group * poly(trial_number_in_block, 2) * block + (1|SEGA_ID) + pitch + age + gender,
  data = et_erp_trial, REML = F)
anova(quadratic_fit)
r2_nakagawa(quadratic_fit) 

cubic_fit <- lmer(
  z_rpd ~ trial * manipulation * group * poly(trial_number_in_block, 3) * block + (1|SEGA_ID) + pitch + age + gender,
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
  z_MMN_amplitude ~ trial * group * block * z_handdynamometer + (1|SEGA_ID) + age + gender,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm)

lmm <- lmer(
  z_MMN_latency ~ trial * group * block * z_handdynamometer+ (1|SEGA_ID) + gender + age,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm)

lmm <- lmer(
  z_P3a_amplitude ~ trial * group * block * z_handdynamometer + (1|SEGA_ID) + gender + age,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm) 

emtrends(lmm,~ trial * block, var = 'z_handdynamometer')


lmm <- lmer(
  z_P3a_latency ~ trial * group * block * z_handdynamometer+ (1|SEGA_ID) + gender + age,
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm) 

lmm <- lmer(
  z_rpd ~ trial * group * block * z_handdynamometer + (1|SEGA_ID),
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm)
r2_nakagawa(lmm) 

lmm <- lmer(
  z_rpd_low ~ trial * group * block * z_handdynamometer + (1|SEGA_ID),
  data = et_erp_subject[et_erp_subject$manipulation == "after", ])
anova(lmm) 
r2_nakagawa(lmm)

# New df et_erp_subject shows same result as df_trial, is valide ####
lmm <- lmer(
z_rpd ~ trial * manipulation * block + trial_number_in_block + (1 | id),
  data = df_trial[
    df_trial$phase %in% c("oddball_block", "oddball_block_rev"), ])
anova(lmm)

lmm <- lmer(
  z_rpd ~ trial * manipulation * block + trial_number_in_block + (1 | SEGA_ID),
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