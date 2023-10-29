# STAT_MODEL_TYPE1.R
# This is the main code that runs the first type of model, namely to compare the 
# cortical activity of the PD group to the HC group for the various gait events 
# when no freezing occurred.
# INPUT:
# - data_nrlgait.csv (created by create_dataframe_nrlgait.m)
# - participants.xlsx
# OUTPUT:
# - data_nrl.Rdata: the reorganized tabulated data in R format
# - folder 'models': containing all the created (raw) models for each condition and each ROI
# - folder 'figures' (optional): containing figures of the conditional effects of each condition and each ROI
# - folder 'summary': containing a summary for each condition of all fixed effect of all ROIs with the probability intervals
# - folder 'priors_v2': containing similar models as the folder 'models' but with normal distributed priors for the sensitivity analysis
# - folder 'models_LvR': containing all the created (raw) models of check_LR_diff for each condition and each ROI
# DEPENDENCIES:
# - run_model.R
# - check_LR_diff.R

# load libraries & functions
# .libPaths("/home/hcockx/R/x86_64-pc-linux-gnu-library/4.1", FALSE)
library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(xlsx)
library(tidybayes)
library(modelr)
library(emmeans)
library(tibble)
# library(tidyverse)
setwd("/vol/dcn/biophysics/prompt/freezing_fnirs/scripts/finalfinal")
source("run_model.R")
source("stat_table.R")
source("check_LR_diff.R")

# set options
options(contrasts = c("contr.sum","contr.poly"))
theme_set(theme_classic(base_size = 12, base_family = ""))
options(mc.cores = parallel::detectCores())

## Read & prepare the data
participants <- read.xlsx('participants.xlsx', sheetName = 'Summary', endRow = 47)
participants[is.na(participants)]<-0
data <- read.csv('data_nrlgait.csv')
data <- merge(data, participants, by = "ID")

# read as factors
data$ID <-as.factor(data$ID)
# data$group <- factor(data$group, levels = c("PD", "HC"))
data$group <- as.factor(data$group)
data$trial <- as.factor(data$trial)
data$trigger <- as.factor(data$trigger)
data$walkingvsstanding <- as.factor(data$walkingvsstanding)
data$direction <- as.factor(data$direction)
data$channel <- as.factor(data$channel)
data$ROI <- factor(data$ROI, levels=c("M1", "PMC", "SMA", "dlPFC", "PPC"))
data$hemisphere <- as.factor(data$hemisphere)
data$time <- factor(data$time, levels=c("baseline", "pre", "post", "post2"))# change order of time
data$X.TF.turn <- data$X.TF.turn.left + data$X.TF.turn.right

# reshape data to contain a column for each time point
data <- pivot_wider(data, names_from = "time", values_from = c("HbO", "HbR"))

# baseline correct the data & multiply by 100
data$HbO_pre <- (data$HbO_pre - data$HbO_baseline)
data$HbO_post <- (data$HbO_post - data$HbO_baseline)
data$HbO_post2 <- (data$HbO_post2 - data$HbO_baseline)
# (do the same for HbR)

# exclude PD10 (MSA) & PD89 (too many bad channels)
data <- droplevels(subset(data, ID!="PD10" & ID!="PD89"))

# center & standardize FOG severity & MDS-UPDRS for the PD group
data_PD <- droplevels(subset(data, group == "PD"))
data_PD$MDS.UPDRS_sd <- scale(data_PD$MDS.UPDRS)
data_PD$X.TF.total_sd <- scale(data_PD$X.TF.total)
data_PD$X.TF.turn_sd <- scale(data_PD$X.TF.turn)
data_PD$X.TF.door_sd <- scale(data_PD$X.TF.door)
cov_PD <- data_PD %>% select(c("trial", "channel", "MDS.UPDRS_sd", "X.TF.total_sd", "X.TF.turn_sd", "X.TF.door_sd"))
data <- merge(data, cov_PD, by = c("trial", "channel"), all.y = TRUE, all.x = TRUE, no.dups = TRUE)

save(data, file = "data_nrl.Rdata")

### Run the analysis
# start to walk
data_start <- droplevels(subset(data, data$trigger=="start")) # select the data
summary_start <- run_model(data_start, "HbO_post ~ 1 + group + (1|ID)", c(), c(), "ROI", "start", 0)
save(summary_start, file = "summary/nrl/start_ROI.Rdata")
write.table(summary_start, file = "summary/nrl/start_ROI.csv", sep =",", row.names = FALSE)
summary_start_sens <- run_model(data_start, "HbO_post ~ 1 + group + (1|ID)", c(), c(), "ROI", "start", 1) # run the models once more with normal distributed priors for the sensitivity analysis and save these models in priors_v2
summary_start_chan <- run_model(data_start, "HbO_post ~ 1 + group + (1|ID)", c(), c(), "channels", "start", 0) # run model for each channel (for plotting purposes)
save(summary_start_chan, file = "summary/nrl/start_channel.Rdata")
write.table(summary_start_chan, file = "summary/nrl/start_channel.csv", sep =",", row.names = FALSE)
# test correlations with FOG sev. & MDS-UPDRS
data_start_PD <- droplevels(subset(data_PD, trigger=="start"))
summary_start_PDcov <- run_model(data_start_PD, "HbO_post ~ 1 + X.TF.total_sd + MDS.UPDRS_sd + (1|ID)", c(), c(), "ROI", "start_PDcov", 0)
save(summary_start_PDcov, file = "summary/PDcov/start_PDcov_ROI.Rdata")
write.table(summary_start_PDcov, file = "summary/PDcov/start_PDcov_ROI.csv", sep =",", row.names = FALSE)

# during walking (= [7 10]s)
summary_walking <- run_model(data_start, "HbO_post2 ~ 1 + group + (1|ID)", c(), c(), "ROI", "walking", 0)
save(summary_walking, file = "summary/nrl/walking_ROI.Rdata")
write.table(summary_walking, file = "summary/nrl/walking_ROI.csv", sep =",", row.names = FALSE)
summary_walking_sens <- run_model(data_start, "HbO_post2 ~ 1 + group + (1|ID)", c(), c(), "ROI", "walking", 1)
summary_walking_chan <- run_model(data_start, "HbO_post2 ~ 1 + group + (1|ID)", c(), c(), "channels", "walking", 0)
save(summary_walking_chan, file = "summary/nrl/walking_channel.Rdata")
write.table(summary_walking_chan, file = "summary/nrl/walking_channel.csv", sep =",", row.names = FALSE)
# test correlations with FOG sev. & MDS-UPDRS
summary_walking_PDcov <- run_model(data_start_PD, "HbO_post2 ~ 1 + X.TF.total_sd + MDS.UPDRS_sd + (1|ID)", c(), c(), "ROI", "walking_PDcov", 0)
save(summary_walking_PDcov, file = "summary/PDcov/walking_PDcov_ROI.Rdata")
write.table(summary_walking_PDcov, file = "summary/PDcov/walking_PDcov_ROI.csv", sep =",", row.names = FALSE)

# stopping
data_stop <- droplevels(subset(data, trigger == "stop"))
summary_stop <- run_model(data_stop, "HbO_post ~ 1 + group + (1|ID)", c(), c(), "ROI", "stop", 0)
save(summary_stop, file = "summary/nrl/stop_ROI.Rdata")
write.table(summary_stop, file = "summary/nrl/stop_ROI.csv", sep =",", row.names = FALSE)
summary_stop_sens <- run_model(data_stop, "HbO_post ~ 1 + group + (1|ID)", c(), c(), "ROI", "stop", 1)
summary_stop_chan <- run_model(data_stop, "HbO_post ~ 1 + group + (1|ID)", c(), c(), "channels", "stop", 0)
save(summary_stop_chan, file = "summary/nrl/stop_channel.Rdata")
write.table(summary_stop_chan, file = "summary/nrl/stop_channel.csv", sep =",", row.names = FALSE)
# test correlations with FOG sev. & MDS-UPDRS
data_stop_PD <- droplevels(subset(data_PD, trigger == "stop"))
summary_stop_PDcov <- run_model(data_stop_PD, "HbO_post ~ 1 + X.TF.total_sd + MDS.UPDRS_sd + (1|ID)", c(), c(), "ROI", "stop_PDcov", 0)
save(summary_stop_PDcov, file = "summary/PDcov/stop_PDcov_ROI.Rdata")
write.table(summary_stop_PDcov, file = "summary/PDcov/stop_PDcov_ROI.csv", sep =",", row.names = FALSE)

# standing
summary_standing <- run_model(data_stop, "HbO_post2 ~ 1 + group + (1|ID)", c(), c(), "ROI", "standing", 0) # no LR differences
save(summary_standing, file = "summary/nrl/standing_ROI.Rdata")
write.table(summary_standing, file = "summary/nrl/standing_ROI.csv", sep =",", row.names = FALSE)
summary_standing_sens <- run_model(data_stop, "HbO_post2 ~ 1 + group + (1|ID)", c(), c(), "ROI", "standing", 1) 
summary_standing_chan <- run_model(data_stop, "HbO_post2 ~ 1 + group + (1|ID)", c(), c(), "channels", "standing", 0)
save(summary_standing_chan, file = "summary/nrl/standing_channel.Rdata")
write.table(summary_standing_chan, file = "summary/nrl/standing_channel.csv", sep =",", row.names = FALSE)
# test correlations with FOG sev. & MDS-UPDRS
summary_standing_PDcov <- run_model(data_stop_PD, "HbO_post2 ~ 1 + X.TF.total_sd + MDS.UPDRS_sd + (1|ID)", c(), c(), "ROI", "standing_PDcov", 0)
save(summary_standing_PDcov, file = "summary/PDcov/standing_PDcov_ROI.Rdata")
write.table(summary_standing_PDcov, file = "summary/PDcov/standing_PDcov_ROI.csv", sep =",", row.names = FALSE)

## door
data_door <- droplevels(subset(data, trigger == "door"))
data_door <- droplevels(subset(data_door, walkingvsstanding == "walking"))
summary_door <- run_model(data_door, 'HbO_post ~ 1 + group + (1|ID)', c(), c(), "ROI", "door", 0)
save(summary_door, file = "summary/nrl/door_ROI.Radata")
write.table(summary_door, file = "summary/nrl/door_ROI.csv", sep =",", row.names = FALSE)
summary_door_sens <- run_model(data_door, 'HbO_post ~ 1 + group + (1|ID)', c(), c(), "ROI", "door", 1)
summary_door_chan <- run_model(data_door, 'HbO_post ~ 1 + group + (1|ID)', c(), c(), "channels", "door",0)
save(summary_door_chan, file = "summary/nrl/door_channel.Rdata")
write.table(summary_door_chan, file = "summary/nrl/door_channel.csv", sep =",", row.names = FALSE)
# make new subgroup of patients that froze at the door
levels(data_door$group) <- c("HC", "FOG-", "FOG+")
data_door$group[which(data_door$group == "PD")] <- "FOG-"
data_door$group[which(data_door$X.TF.door>0)] <- "FOG+"
summary_door_sub <- run_model(data_door, 'HbO_post ~ 1 + group + (1|ID)', c(), c(), "ROI", "door_sub",0)
save(summary_door_sub, file = "summary/door_sub_ROI.Rdata")
write.table(summary_door_sub, file = "summary/door_sub_ROI.csv", sep =",")
summary_door_sub_sens <- run_model(data_door, 'HbO_post ~ 1 + group + (1|ID)', c(), c(), "ROI", "door_sub",1)
# test correlations with FOG sev. & MDS-UPDRS
data_door_PD <- droplevels(subset(data_PD, trigger == "door"))
data_door_PD <- droplevels(subset(data_door_PD, walkingvsstanding == "walking"))
summary_door_PDcov <- run_model(data_door_PD, "HbO_post ~ 1 + X.TF.door_sd + MDS.UPDRS_sd + (1|ID)", c(), c(), "ROI", "door_PDcov",0)
save(summary_door_PDcov, file = "summary/PDcov/door_PDcov_ROI.Rdata")
write.table(summary_door_PDcov, file = "summary/PDcov/door_PDcov_ROI.csv", sep =",", row.names = FALSE)
# correlation of doorway passing with stopping
stop_door<- merge(summary_stop_chan, summary_door_chan, by = c("channel", "factor"), remove = FALSE)
stop_door_I <- subset(stop_door, factor == "Intercept")
stop_door_HC <- subset(stop_door, factor == "HC")
stop_door_PD <- subset(stop_door, factor == "PD")
cor.test(stop_door_I$Estimate.x, stop_door_I$Estimate.y, method = "pearson")
cor.test(stop_door_HC$Estimate.x, stop_door_HC$Estimate.y, method = "pearson")
cor.test(stop_door_PD$Estimate.x, stop_door_PD$Estimate.y, method = "pearson")
stop_door$factor <- factor(stop_door$factor, levels=c('PD', 'HC'))
ggplot(data = subset(stop_door, factor %in% c("HC", "PD")), aes(x = Estimate.x, y = Estimate.y)) + 
  geom_point(aes(color = factor), show.legend = TRUE) + 
  geom_smooth(method = "lm", se = FALSE, aes(color = factor), show.legend = FALSE) +
  facet_wrap(~factor) + labs(x = "stop", y = "doorway") + 
  theme(strip.text.x = element_blank(), text = element_text(size = 20, family = 'Arial', color = "black")) +
  ylim(-0.5, 0.3) + coord_fixed(ratio = 1) +
  theme(legend.position = c(0.95, 0.2), legend.title = element_blank()) +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
  scale_color_manual(values = c("#007172", "#C13700"))
ggsave('figures/correlations/stop_vs_door.eps',device = cairo_ps, width = 18, height = 8, units = 'cm')

## turn
data_turn <- droplevels(subset(data, trigger == "turn"))
data_turn <- droplevels(subset(data_turn, walkingvsstanding == "walking"))
ROI_LR_turn <- check_LR_diff(data_turn, 'HbO_post ~ 1 + hemisphere*direction + (1 |ID)', TRUE, "turn") # 95% CI of direction does not exclude 0
summary_turn <- run_model(data_turn, 'HbO_post ~ 1 + group + (1|ID)', c(), c(), "ROI", "turn", 0)
save(summary_turn, file = "summary/nrl/turn_ROI.Rdata")
write.table(summary_turn, file = "summary/nrl/turn_ROI.csv", sep =",", row.names = FALSE)
summary_turn_sens <- run_model(data_turn, 'HbO_post ~ 1 + group + (1|ID)', c(), c(), "ROI", "turn", 1)
summary_turn_chan <- run_model(data_turn, 'HbO_post ~ 1 + group + (1|ID)', c(), c(), "channels", "turn", 0)
save(summary_turn_chan, file = "summary/nrl/turn_channel.Rdata")
write.table(summary_turn_chan, file = "summary/nrl/turn_channel.csv", sep =",", row.names = FALSE)
# test correlations with FOG sev. & MDS-UPDRS
data_turn_PD <- droplevels(subset(data_PD, trigger == "turn"))
data_turn_PD <- droplevels(subset(data_turn_PD, walkingvsstanding == "walking"))
summary_turn_PDcov <- run_model(data_turn_PD, "HbO_post ~ 1 + X.TF.turn_sd + MDS.UPDRS_sd + (1|ID)", c(), c(), "ROI", "turn_PDcov", 0)
save(summary_turn_PDcov, file = "summary/PDcov/turn_PDcov_ROI.Rdata")
write.table(summary_turn_PDcov, file = "summary/PDcov/turn_PDcov_ROI.csv", sep =",", row.names = FALSE)
# correlation of turning with stopping
stop_turn<- merge(summary_stop_chan, summary_turn_chan, by = c("channel", "factor"), remove = FALSE)
stop_turn_I <- subset(stop_turn, factor == "Intercept")
stop_turn_HC <- subset(stop_turn, factor == "HC")
stop_turn_PD <- subset(stop_turn, factor == "PD")
cor.test(stop_turn_I$Estimate.x, stop_turn_I$Estimate.y, method = "pearson")
cor.test(stop_turn_HC$Estimate.x, stop_turn_HC$Estimate.y, method = "pearson")
cor.test(stop_turn_PD$Estimate.x, stop_turn_PD$Estimate.y, method = "pearson")
stop_turn$factor <- factor(stop_turn$factor, levels=c('PD', 'HC'))
ggplot(data = subset(stop_turn, factor %in% c("HC", "PD")), aes(x = Estimate.x, y = Estimate.y)) + 
  geom_point(aes(color = factor), show.legend = FALSE) + 
  geom_smooth(method = "lm", se = FALSE, aes(color = factor), show.legend = FALSE) +
  facet_wrap(~factor) + labs(x = "stop", y = "turn") + 
  theme(strip.text.x = element_blank(), text = element_text(size = 20, family = 'Arial', color = "black")) +
  ylim(-0.5, 0.3) + coord_fixed(ratio = 1) +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
  scale_color_manual(values = c("#007172", "#C13700"))
ggsave('figures/correlations/stop_vs_turn.eps',device = cairo_ps, width = 18, height = 8, units = 'cm')

