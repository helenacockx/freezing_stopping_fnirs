# RUN_MODEL.R
# This function runs the actual bayesian model for the given data with or without averaging it over channels belonging to the same ROI
# INPUT:
# - data_input: subset of data for a given condition (see: stat_modeltype1/2.R)
# - brm_form: formula for the Bayesian analysis (see brmsformula)
# - ROI_hemi: ROIs for which the data should be split in a left and right ROI (not used here)
# - ROI_dir: ROIs for which the data containing a condition in different directions (e.g. left and right turn) should be analysed seperately (not used here)
# - test_ROI: "ROI" or "channels". If ROI, then average the data of channels belonging to the same ROI, if channels, do not average the data.
# - cond: name of the condition being tested (for saving purposes)
# - sensitivity_analysis: 1 or 0. If 1, then calculate the model starting from priors with standardized distribution (default = 0)
# OUTPUT (see also stat_model_type1/2.R):
# - folder 'models': containing all the created (raw) models for each condition and each ROI
# - folder 'figures' (optional): containing figures of the conditional effects of each condition and each ROI
# - folder 'priors_v2': containing similar models as the folder 'models' but with normal distributed priors for the sensitivity analysis
# DEPENDENCIES:
# - stat_table.R: creates a summary for all fixed effect of the tested ROI with the probability intervals

run_model <- function(data_input, brm_form, ROI_hemi, ROI_dir, test_ROI, cond, sensitivity_analysis) {
  # average over channels belonging to the same ROI, hemispheres taken together
  data_ROI <- data_input %>%
    group_by(ID, group, trial, trigger, walkingvsstanding, direction, ROI, X.TF.total_sd, X.TF.door_sd, X.TF.turn_sd,  NFOGQ..total.without.door., MDS.UPDRS_sd, MoCA.total, TMT.b.a, HADS, anxiety.protocol) %>%
    summarise(HbO_pre = mean(HbO_pre),
              HbO_post = mean(HbO_post),
              HbO_post2 = mean(HbO_post2)) %>%
    ungroup()
  # average over channels belonging to the same ROI, for each hemisphere separately
  data_ROIlr <- data_input %>%
    group_by(ID, group, trial, trigger, walkingvsstanding, direction, ROI, hemisphere, X.TF.total_sd, X.TF.turn_sd, NFOGQ..total.without.door., MDS.UPDRS_sd, MoCA.total, TMT.b.a, HADS, anxiety.protocol) %>%
    summarise(HbO_pre = mean(HbO_pre),
              HbO_post = mean(HbO_post),
              HbO_post2 = mean(HbO_post2)) %>%
    ungroup()
  # run the model
  if (test_ROI == "channels"){
    save_dir <- paste(getwd(),  "models", cond, "channels", sep="/")
    fig_dir <- paste(getwd(), "figures", cond, "channels", sep = "/")
    ROIs <- unique(data_input$channel)
  } else if (test_ROI == "ROI") {
    save_dir <- paste(getwd(),  "models", cond, sep="/")
    fig_dir <- paste(getwd(), "figures", cond, sep = "/")
    if (sensitivity_analysis){
      save_dir <- paste(getwd(), "priors_v2", "models", cond, sep="/")
      fig_dir <- paste(getwd(), "priors_v2", "figures", cond, sep = "/")
    }
    ROIs <- c("M1", "PMC",  "SMA", "dlPFC", "PPC")
  }
  dir.create(save_dir)
  dir.create(fig_dir)
  first <- TRUE
  ROI_finished <- FALSE
  if (sensitivity_analysis){
    priors <- c(
      prior(normal(0, 1),class = Intercept),
      prior(normal(0,1), class = b, coef = group1),
      prior(cauchy(0,1), class = sigma)
    )
  }
  for (i in ROIs) {
      for (j in c("left", "right")) { # direction of turn
          for (k in c("left", "right")) { # hemisphere
            if (i %in% ROI_dir & i %in% ROI_hemi) { # if hemisphere & direction dependent
            data_selROI <- droplevels(subset(data_ROIlr, (ROI == i & direction == j & hemisphere == k)))
            name <- paste(i, k, "turn", j, sep = "_")
          } else if (i %in% ROI_hemi) { # if hemisphere dependent
          data_selROI <-
            droplevels(subset(data_ROIlr, (ROI == i & hemisphere == k)))
          name <- paste(i, k, sep = "_")
          if (j=="right"){
            ROI_finished <- TRUE
          }
        } else if (!(i %in% ROI_dir) & !(i %in% ROI_hemi)){ # if not hemisphere or direction dependent
          if (test_ROI == "channels"){
            data_selROI <- droplevels(subset(data_input, channel == i))
            } else if (test_ROI == "ROI") {
              data_selROI <- droplevels(subset(data_ROI, ROI == i))
            }
          name <- i
          if (k=="right"){
            ROI_finished <- TRUE
          }
        }
      if (!ROI_finished){
        filename <- paste(save_dir, "/", "model_", name, sep = "")
        mod_selROI <- brm(brm_form, data = data_selROI,family = gaussian,
          iter = 10000, warmup = 500, chains = 4, control = list(adapt_delta = 0.99), 
          file = filename)
        print(name)
        print(summary(mod_selROI, prob = 0.95))
        summary_ROI <- stat_table(mod_selROI)
        if (test_ROI == "channels"){
          summary_ROI <- cbind(channel = c(name), ROI = unique(data_selROI$ROI), summary_ROI)
        } else if (test_ROI == "ROI"){
          summary_ROI <- cbind(ROI = c(name), summary_ROI)
        }
        if (first) {
          summary_table  <- summary_ROI
          } else {
            summary_table <- rbind(summary_table, summary_ROI)
            }
        first <- FALSE
      # plot(mod_selROI)
      # print(pp_check(mod_selROI))
      # png(file = paste(paste(paste(fig_dir, cond, sep = "/"), name, sep = "_"), "png", sep = "."))
      # print(plot(conditional_effects(mod_selROI)))
      # dev.off()
      
      }
          }# loop over direction
      } # loop over hemisphere
    ROI_finished <- FALSE
  } # loop over ROI
  return(summary_table)
}