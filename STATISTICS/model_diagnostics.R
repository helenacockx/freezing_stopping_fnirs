# MODEL_DIAGNOSTICS
# This script runs an analysis on the created models to check them for assumptions. 
# It includes caterpillar plots, Gelman-Rubin statistics, checking for divergent transitions,
# effective sample size and a sensitivity analysis.
# INPUT:
# - folder 'models': containing all the created (raw) models for each condition and each ROI (by stat_model_type1/2.R)
# - folder 'priors_v2': ontaining similar models as the folder 'models' but with normal distributed priors for the sensitivity analysis
# OUTPUT:
# - folder 'diagnostics': containing:
#   . summary_diagn_***.Rdata and ***.csv: containing numeric values for the checks in a table
#   . ***.pdf: containing a visual representation of the checks

library(coda)
library(ggplot2)
library(bayesplot)
library(brms)
library(dplyr)
library(ggdist)
setwd("/vol/dcn/biophysics/prompt/freezing_fnirs/scripts/finalfinal")

# create empty table
diagn <- data.frame(matrix(ncol = 7, nrow = 60))
colnames(diagn) <- c('condition', 'ROI', 'parameter', 'Rhat', 'ESS_ratio', 'ESS', 'n_diverg')

# loop over models & collect diagnostic info
conditions <- c('start', 'walking', 'stop', 'standing', 'turn', 'door')
ROI <- c('M1', 'PMC', 'SMA', 'dlPFC', 'PPC')
idx <- 1
for (c in conditions) {
  destination = file.path('diagnostics', sprintf('%s.pdf', c)) # save all plots to pdf
  pdf(file = destination)
  for (r in ROI) {
    model <- readRDS(file.path('models', c, sprintf('model_%s.rds', r)))
    p <- ggplot() + annotate('text', x=4, y = 25, size = 10, label = r) + theme_void()
    print(p)
    plot(model)
    print(pp_check(model))
    model.temp <- as.mcmc(model)
    gelman.plot(model.temp[, "b_Intercept", ])
    gelman.plot(model.temp[, "b_group1", ])
    print(mcmc_acf(model, pars = variables(model)[c(1,2)]))
    # sensitivity analysis
    draws_v1 <- as_draws_df(model) 
    draws_v1$prior <- 'v1'
    model_v2 <- readRDS(file.path('priors_v2', 'models', c, sprintf('model_%s.rds', r)))
    draws_v2 <- as_draws_df(model_v2)
    draws_v2$prior <- 'v2'
    draws <- rbind(draws_v1, draws_v2)
    p <- draws %>%
      ggplot(aes(x = b_Intercept, fill = prior, color = prior)) + 
      stat_slab(alpha = 0.3) + stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single")) +
      labs(title = sprintf("Sensitivity Analysis %s", r))
    print(p)
    p <- draws %>%
      ggplot(aes(x = b_group1, fill = prior, color = prior)) + 
      stat_slab(alpha = 0.3) + stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single")) +
      labs(title = sprintf("Sensitivity Analysis %s", r))
    print(p)
    # compute diagnostics 
    rhat <- rhat(model)
    ESS_ratio <- neff_ratio(model)
    ESS <- round(ESS_ratio * (10000 - 500)*4, digits= 0) # ratio * number of post-warmup samples (https://bookdown.org/arjkurz/Statistical_Rethinking_recoded/multilevel_models.html)
    np <- nuts_params(model)
    n_diverg <- sum(subset(np, Parameter == "divergent__")$Value)
    # fill in table
      diagn$condition[c(idx,idx+1)] <- c
      diagn$ROI[c(idx,idx+1)] <- r
      diagn$parameter[c(idx,idx+1)] <- c('Intercept', 'Group')
      diagn$Rhat[c(idx,idx+1)] <- rhat[1:2]
      diagn$ESS_ratio[c(idx,idx+1)] <- ESS_ratio[1:2]
      diagn$ESS[c(idx,idx+1)] <- ESS[1:2]
      diagn$n_diverg[c(idx, idx+1)] <- n_diverg
      idx <- idx + 2
  }
  dev.off()
}
save(diagn, file = "diagnostics/summary_diagn_model1.Rdata")
write.table(diagn, file = "diagnostics/summary_diagn_model1.csv", sep =",", row.names = FALSE)

## Model 1 subgroup
# create empty table
diagn <- data.frame(matrix(ncol = 7, nrow = 15))
colnames(diagn) <- c('condition', 'ROI', 'parameter', 'Rhat', 'ESS_ratio', 'ESS', 'n_diverg')
# loop over models & collect diagnostic info
conditions <- c('door_sub')
ROI <- c('M1', 'PMC', 'SMA', 'dlPFC', 'PPC')
idx <- 1
for (c in conditions) {
  destination = file.path('diagnostics', sprintf('%s.pdf', c)) # save all plots to pdf
  pdf(file = destination)
  for (r in ROI) {
    model <- readRDS(file.path('models', c, sprintf('model_%s.rds', r)))
    p <- ggplot() + annotate('text', x=4, y = 25, size = 10, label = r) + theme_void()
    print(p)
    plot(model)
    print(pp_check(model))
    model.temp <- as.mcmc(model)
    gelman.plot(model.temp[, "b_Intercept", ])
    gelman.plot(model.temp[, "b_group1", ])
    gelman.plot(model.temp[, "b_group2",])
    print(mcmc_acf(model, pars = variables(model)[c(1,2)]))
    # sensitivity analysis
    draws_v1 <- as_draws_df(model) 
    draws_v1$prior <- 'v1'
    model_v2 <- readRDS(file.path('priors_v2', 'models', c, sprintf('model_%s.rds', r)))
    draws_v2 <- as_draws_df(model_v2)
    draws_v2$prior <- 'v2'
    draws <- rbind(draws_v1, draws_v2)
    p <- draws %>%
      ggplot(aes(x = b_Intercept, fill = prior, color = prior)) + 
      stat_slab(alpha = 0.3) + stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single")) +
      labs(title = sprintf("Sensitivity Analysis %s", r))
    print(p)
    p <- draws %>%
      ggplot(aes(x = b_group1, fill = prior, color = prior)) + 
      stat_slab(alpha = 0.3) + stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single")) +
      labs(title = sprintf("Sensitivity Analysis %s", r))
    print(p)
    p <- draws %>%
      ggplot(aes(x = b_group2, fill = prior, color = prior)) + 
      stat_slab(alpha = 0.3) + stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single")) +
      labs(title = sprintf("Sensitivity Analysis %s", r))
    print(p)
    # compute diagnostics 
    rhat <- rhat(model)
    ESS_ratio <- neff_ratio(model)
    ESS <- round(ESS_ratio * (10000 - 500)*4, digits= 0) # ratio * number of post-warmup samples (https://bookdown.org/arjkurz/Statistical_Rethinking_recoded/multilevel_models.html)
    np <- nuts_params(model)
    n_diverg <- sum(subset(np, Parameter == "divergent__")$Value)
    # fill in table
    diagn$condition[seq(idx,idx+2)] <- c
    diagn$ROI[seq(idx,idx+2)] <- r
    diagn$parameter[seq(idx,idx+2)] <- c('Intercept', 'Group1', 'Group2')
    diagn$Rhat[seq(idx,idx+2)] <- rhat[1:3]
    diagn$ESS_ratio[seq(idx,idx+2)] <- ESS_ratio[1:3]
    diagn$ESS[seq(idx,idx+2)] <- ESS[1:3]
    diagn$n_diverg[seq(idx,idx+2)] <- n_diverg
    idx <- idx + 3
  }
  dev.off()
}
save(diagn, file = "diagnostics/summary_diagn_model1_sub.Rdata")
write.table(diagn, file = "diagnostics/summary_diagn_model1_sub.csv", sep =",", row.names = FALSE)


## FOG table
# create empty table
diagn <- data.frame(matrix(ncol = 7, nrow = 30))
colnames(diagn) <- c('condition', 'ROI', 'parameter', 'Rhat', 'ESS_ratio', 'ESS', 'n_diverg')

# loop over models & collect diagnostic info
conditions <- c('FOGturn', 'FOGdoor')
ROI <- c('M1', 'PMC', 'SMA', 'dlPFC', 'PPC')
idx <- 1
for (c in conditions) {
  destination = file.path('diagnostics', sprintf('%s.pdf', c)) # save all plots to pdf
  pdf(file = destination)
  for (r in ROI) {
    model <- readRDS(file.path('models', c, sprintf('model_%s.rds', r)))
    p <- ggplot() + annotate('text', x=4, y = 25, size = 10, label = r) + theme_void()
    print(p)
    plot(model)
    print(pp_check(model))
    model.temp <- as.mcmc(model)
    gelman.plot(model.temp[, "b_Intercept", ])
    gelman.plot(model.temp[, "b_group1", ])
    gelman.plot(model.temp[, "b_group2",])
    print(mcmc_acf(model, pars = variables(model)[c(1,2)]))
    # sensitivity analysis
    draws_v1 <- as_draws_df(model)
    draws_v1$prior <- 'v1'
    model_v2 <- readRDS(file.path('priors_v2', 'models', c, sprintf('model_%s.rds', r)))
    draws_v2 <- as_draws_df(model_v2)
    draws_v2$prior <- 'v2'
    draws <- rbind(draws_v1, draws_v2)
    p <- draws %>%
      ggplot(aes(x = b_Intercept, fill = prior, color = prior)) +
      stat_slab(alpha = 0.3) + stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single")) +
      labs(title = sprintf("Sensitivity Analysis %s", r))
    print(p)
    p <- draws %>%
      ggplot(aes(x = b_group1, fill = prior, color = prior)) +
      stat_slab(alpha = 0.3) + stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single")) +
      labs(title = sprintf("Sensitivity Analysis %s", r))
    print(p)
      p <- draws %>%
        ggplot(aes(x = b_group2, fill = prior, color = prior)) +
        stat_slab(alpha = 0.3) + stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single")) +
        labs(title = sprintf("Sensitivity Analysis %s", r))
      print(p)
    # compute diagnostics 
    rhat <- rhat(model)
    ESS_ratio <- neff_ratio(model)
    ESS <- round(ESS_ratio * (10000 - 500)*4, digits= 0) # ratio * number of post-warmup samples (https://bookdown.org/arjkurz/Statistical_Rethinking_recoded/multilevel_models.html)
    np <- nuts_params(model)
    n_diverg <- sum(subset(np, Parameter == "divergent__")$Value)
    # fill in table
      diagn$condition[seq(idx,idx+2)] <- c
      diagn$ROI[seq(idx,idx+2)] <- r
      diagn$parameter[seq(idx,idx+2)] <- c('Intercept', 'Condition1', 'Condition2')
      diagn$Rhat[seq(idx,idx+2)] <- rhat[1:3]
      diagn$ESS_ratio[seq(idx,idx+2)] <- ESS_ratio[1:3]
      diagn$ESS[seq(idx,idx+2)] <- ESS[1:3]
      diagn$n_diverg[seq(idx,idx+2)] <- n_diverg
      idx <- idx + 3
  }
  dev.off()
}
save(diagn, file = "diagnostics/summary_diagn_model2.Rdata")
write.table(diagn, file = "diagnostics/summary_diagn_model2.csv", sep =",", row.names = FALSE)
