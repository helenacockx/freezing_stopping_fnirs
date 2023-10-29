# STAT_MODEL_TYPE2.R
# This is the main code that runs the second type of model, namely to compare the 
# cortical activity during a freezing event to a voluntary stop and to a successful
# event of the same type as the freezing event (e.g., turning freeze vs stop vs successful turn). 
# INPUT:
# - data_FOG.csv (created by create_dataframe_FOG.m)
# - data_nrl.Rdata (created by stat_model_type1.R)
# - participants.xlsx
# OUTPUT:
# - folder 'models': containing all the created (raw) models for each condition and each ROI
# - folder 'figures' (optional): containing figures of the conditional effects of each condition and each ROI
# - folder 'summary': containing a summary for each condition of all fixed effect of all ROIs with the probability intervals
# - folder 'priors_v2': containing similar models as the folder 'models' but with normal distributed priors for the sensitivity analysis
# DEPENDENCIES:
# - run_model.R

# load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(xlsx)
library(tidybayes)
library(modelr)
library(emmeans)
library(tibble)
library(coda)
# library(lme4)
# library(tidyverse)
setwd("/vol/dcn/biophysics/prompt/freezing_fnirs/scripts/finalfinal")
source("run_model.R")
source("stat_table.R")

# set options
options(contrasts = c("contr.sum","contr.poly"))
theme_set(theme_classic(base_size = 12, base_family = ""))
options(mc.cores = parallel::detectCores())

## Load the data
# load the FOG data
data_FOG<- read.csv('data_FOG.csv')
participants <- read.xlsx('participants.xlsx', sheetName = 'Summary', endRow = 47)
participants[is.na(participants)]<-0
data_FOG <- merge(data_FOG, participants, by = "ID")

# read as factors
data_FOG$ID <-as.factor(data_FOG$ID)
data_FOG$trial <- as.factor(data_FOG$trial)
data_FOG$trigger <- as.factor(data_FOG$trigger)
data_FOG$walkingvsstanding <- as.factor(data_FOG$walkingvsstanding)
data_FOG$direction <- as.factor(data_FOG$direction)
data_FOG$channel <- as.factor(data_FOG$channel)
data_FOG$ROI <- factor(data_FOG$ROI, levels=c("M1", "PMC", "SMA", "dlPFC", "PPC"))
data_FOG$hemisphere <- as.factor(data_FOG$hemisphere)
data_FOG$time <- factor(data_FOG$time, levels=c("baseline", "pre", "post", "post2"))
data_FOG$X.TF.turn <- data_FOG$X.TF.turn.left + data_FOG$X.TF.turn.right

# reshape data to contain a column for each time point
data_FOG <- pivot_wider(data_FOG, names_from = "time", values_from = c("HbO", "HbR"))

# baseline correct the data 
data_FOG$HbO_pre <- (data_FOG$HbO_pre - data_FOG$HbO_baseline)
data_FOG$HbO_post <- (data_FOG$HbO_post - data_FOG$HbO_baseline)
data_FOG$HbO_post2 <- (data_FOG$HbO_post2 - data_FOG$HbO_baseline)
# (do the same for HbR)

# exclude PD10 & PD89
data_FOG <- droplevels(subset(data_FOG, ID!="PD10" & ID!="PD89"))

# combine with the successful events
data_FOG$condition <- as.factor("FOG")

# center & standardize FOG severity & MDS-UPDRS 
# (just to create the same columns as the normal data, we will not use this here)
data_FOG$MDS.UPDRS_sd <- scale(data_FOG$MDS.UPDRS)
data_FOG$X.TF.total_sd <- scale(data_FOG$X.TF.total)
data_FOG$X.TF.turn_sd <- scale(data_FOG$X.TF.turn)
data_FOG$X.TF.door_sd <- scale(data_FOG$X.TF.door)

# load the stops & successful events
load("data_nrl.Rdata")
data_nrl <- droplevels(subset(data, group == "PD"))
data_nrl <- subset(data_nrl, select = -c(group))

## run the models
# turn
data_FOG_turn <- subset(data_FOG, (trigger == "turn" & walkingvsstanding =="walking"))
data_stop_turn <- subset(data_nrl, (trigger == "stop"))
data_stop_turn$condition <- as.factor("stop")
data_nrl_turn <- subset(data_nrl, (trigger == "turn" & walkingvsstanding == "walking"))
data_nrl_turn$condition <- as.factor("nrl")
data_all_turn <- rbind(data_FOG_turn, data_stop_turn, data_nrl_turn)
data_all_turn$group <- data_all_turn$condition # rename condition to group, so we can still use run_model
summary_FOGturn <- run_model(data_all_turn, "HbO_post ~ 1 + group + (1 + group|ID)", c(), c(), "ROI", "FOGturn", 0)
save(summary_FOGturn, file = "summary/FOG/FOGturn_ROI.Rdata")
write.table(summary_FOGturn, file = "summary/FOG/FOGturn_ROI.csv", sep =",", row.names = FALSE)
summary_FOGturn_sens <- run_model(data_all_turn, "HbO_post ~ 1 + group + (1 + group|ID)", c(), c(), "ROI", "FOGturn", 1)  # run the models once more with normal distributed priors for the sensitivity analysis and save these models in priors_v2

# door
data_FOG_door <- subset(data_FOG, (trigger == "door" & walkingvsstanding =="walking"))
data_stop_door <- subset(data_nrl, (trigger == "stop"))
data_stop_door$condition <- as.factor("stop")
data_nrl_door <- subset(data_nrl, (trigger == "door" & walkingvsstanding == "walking"))
data_nrl_door$condition <- as.factor("nrl")
data_all_door <- rbind(data_FOG_door, data_stop_door, data_nrl_door)
data_all_door$group <- data_all_door$condition # rename condition to group, so we can still use run_model
summary_FOGdoor <- run_model(data_all_door, "HbO_post ~ 1 + group + (1 + group|ID)", c(), c(), "ROI", "FOGdoor", 0)
save(summary_FOGdoor, file = "summary/FOG/FOGdoor_ROI.Rdata")
write.table(summary_FOGdoor, file = "summary/FOG/FOGdoor_ROI.csv", sep =",", row.names = FALSE)
summary_FOGdoor_sens <- run_model(data_all_door, "HbO_post ~ 1 + group + (1 + group|ID)", c(), c(), "ROI", "FOGdoor", 1)
