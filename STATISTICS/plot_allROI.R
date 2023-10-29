# PLOT_allROI
# This script creates one plot with the estimated effects for all ROIs and all groups
# for the condition you fill in (in conditions <- c(...))
# There are different parts in the script for type 1 vs type 2 models
# INPUT:
# - folder 'models': containing all the created (raw) models for each condition and each ROI (by stat_model_type1/2.R)
# OUTPUT:
# - figures are saved in folder 'figures' in a folder of their condition

setwd("/vol/dcn/biophysics/prompt/freezing_fnirs/scripts/finalfinal")
library(ggplot2)
library(dplyr)
library(brms)
library(tidybayes)
theme_set(theme_classic(base_size = 12, base_family = ""))

## Type 1 models
conditions <- c('turn')
# collect model info
ROI <- factor(c('M1', 'PMC', 'SMA', 'dlPFC', 'PPC'), levels = c('M1', 'PMC', 'SMA', 'dlPFC', 'PPC'))
for (c in conditions) {
  draws_all <- data.frame()
  for (r in ROI){
    model <- readRDS(file.path('models', c, sprintf('model_%s.rds', r)))
    draws1 <- model %>% spread_draws(b_Intercept, b_group1) 
    draws1$group <- "HC"
    draws2 <- draws1
    draws2$b_group1 <- -draws1$b_group1
    draws2$group <- "PD"
    draws <- rbind(draws1, draws2)  
    draws <- draws %>% mutate(group_mean = b_Intercept + b_group1)
    draws$group <- factor(draws$group, levels=c("PD", "HC"))
    draws$ROI <- r
    draws_all <- rbind(draws_all, draws)
    }
}
draws_all$ROI[draws_all$ROI == 'dlPFC'] <- 'PFC'
draws_all$ROI <- factor(draws_all$ROI, levels = c('M1', 'PMC', 'SMA', 'PFC', 'PPC'))
draws_all$group <- factor(draws_all$group, levels = c("PD", "HC"))
# plot
draws_all %>%
  ggplot(aes(y = group_mean, x = ROI, fill = group)) +  
  # stat_eye(show.legend = FALSE, position = position_dodge(width = 0.7), scale = 0.8, alpha = 0.8, interval_size_range = c(0.2, 0.8)) +
  stat_eye(position = position_dodge(width = 0.7), scale = 0.8,  alpha = 0.8, interval_size_domain = c(1, 6)) +
  labs(title = NULL, x = NULL, y = NULL) +
  geom_hline(yintercept =0) +
  theme(text = element_text(size = 20, family = 'Arial', color = "black")) +
  theme(legend.position = c(0.94, 0.95), legend.title = element_blank()) +
  # scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  scale_fill_manual(values = c("#007172", "#C13700")) + 
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.4, -0.2, 0, 0.2, 0.4))
  # ylim(-0.7, 0.7)
  # labs(y = "estimated \u0394HbO (z-score)")
if (c %in% c('turn', 'door')){
  ggsave(file.path('figures', c, 'allROI.eps'),device = cairo_ps, width = 31, height = 8, units = 'cm')
} else if (c %in% c('stop', 'standing', 'start', 'walking')) {
  ggsave(file.path('figures', c, 'allROI.eps'),device = cairo_ps, width = 22, height = 8, units = 'cm')
}

## Type 1 model, subgroup analysis
conditions <- c('door_sub')
# collect model info
ROI <- factor(c('M1', 'PMC', 'SMA', 'dlPFC', 'PPC'), levels = c('M1', 'PMC', 'SMA', 'dlPFC', 'PPC'))
for (c in conditions) {
  draws_all <- data.frame()
  for (r in ROI){
    model <- readRDS(file.path('models', c, sprintf('model_%s.rds', r)))
    draws1 <- model %>%
      spread_draws(b_Intercept, b_group1, b_group2) 
    draws1$group <- 'HC'
    draws1$b_group <- draws1$b_group1
    draws2 <- draws1
    draws2$b_group <- draws1$b_group2
    draws2$group <- 'FOG-'
    draws3 <- draws1
    draws3$b_group <- -draws1$b_group1 -draws1$b_group2
    draws3$group <- 'FOG+'
    draws <- rbind(draws1, draws2, draws3)
    draws <- draws %>% mutate(group_mean = b_Intercept + b_group)
    draws$group <- factor(draws$group, levels = c("FOG+", "FOG-", "HC"))
    draws$ROI <- r
    draws_all <- rbind(draws_all, draws)
  }
}
draws_all$ROI[draws_all$ROI == 'dlPFC'] <- 'PFC'
draws_all$ROI <- factor(draws_all$ROI, levels = c('M1', 'PMC', 'SMA', 'PFC', 'PPC'))
# draws_all$group <- factor(draws_all$group, levels = c("PD", "HC"))
# plot
draws_all %>%
  ggplot(aes(y = group_mean, x = ROI, fill = group)) +  
  stat_eye(position = position_dodge(width = 0.7), scale = 1, alpha = 0.8, interval_size_domain = c(1, 6)) +
  labs(title = NULL, x = NULL, y = NULL) +
  geom_hline(yintercept =0) + 
  theme(text = element_text(size = 20,family = 'Arial', color = "black")) +
  theme(legend.position = c(0.95, 0.93), legend.title = element_blank()) +
  # scale_fill_manual(values= c("#BD5051", "#F8766D", "#00BFC4")) +
  scale_fill_manual(values= c("#02494F", "#37A6A9", "#C13700")) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.4, -0.2, 0, 0.2, 0.4))
  # ylim(-0.6, 0.6)
ggsave(file.path('figures', c, 'allROI.eps'),device = cairo_ps, width = 31, height = 8, units = 'cm')

## type 2 models
conditions <- c('FOGturn')
# collect model info
ROI <- factor(c('M1', 'PMC', 'SMA', 'dlPFC', 'PPC'), levels = c('M1', 'PMC', 'SMA', 'dlPFC', 'PPC'))
for (c in conditions) {
  draws_all <- data.frame()
  for (r in ROI){
    model <- readRDS(file.path('models', c, sprintf('model_%s.rds', r)))
    draws1 <- model %>%
      spread_draws(b_Intercept, b_group1, b_group2) 
    draws1$group <- 'FOG'
    draws1$b_group <- draws1$b_group1
    draws2 <- draws1
    draws2$b_group <- draws1$b_group2
    draws2$group <- 'stop'
    draws3 <- draws1
    draws3$b_group <- -draws1$b_group1 -draws1$b_group2
    draws3$group <- 'nrl'
    draws <- rbind(draws1, draws2, draws3)
    draws <- draws %>% mutate(group_mean = b_Intercept + b_group)
    draws$ROI <- r
    draws_all <- rbind(draws_all, draws)
  }
}
draws_all$ROI[draws_all$ROI == 'dlPFC'] <- 'PFC'
draws_all$ROI <- factor(draws_all$ROI, levels = c('M1', 'PMC', 'SMA', 'PFC', 'PPC'))
draws_all$group[draws_all$group == "FOG"] <- "freezing"
draws_all$group[draws_all$group == "nrl"] <- "successful"
draws_all$group <- factor(draws_all$group, c("freezing", "stop", "successful"))
# plot
draws_all %>%
  ggplot(aes(y = group_mean, x = ROI, fill = group)) +  
  stat_eye(position = position_dodge(width = 0.7), scale = 0.8, alpha = 0.8, interval_size_domain = c(1, 6)) +
  labs(title = NULL, x = NULL, y = NULL) +
  geom_hline(yintercept =0) +  
  theme(text = element_text(size = 20,family = 'Arial', color = "black")) +
  theme(legend.position = c(0.94, 0.93), legend.title = element_blank()) +
  scale_fill_manual(values= c("#D55E00", "#2271B2", "#359B73")) + 
  ylim(-0.7, 0.7)
ggsave(file.path('figures', c, 'allROI.eps'),device = cairo_ps, width = 31, height = 8, units = 'cm')
