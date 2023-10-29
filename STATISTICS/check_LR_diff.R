# CHECK_LR_DIFF
# This function tests wether there are statistical differences between the two
# hemispheres or whether there is a significant interaction effect of the hemisphere
# with the direction (of the turn)
# INPUT:
# - data_input: subset of data for a given condition
# - brm_form: formula for the Bayesian analysis (see brmsformula)
# - direction: whether or not also to test for an interaction with direction
# - cond: name of the condition being tested
# OUTPUT:
# - ROI_LR: list of ROIs that show differences between hemispheres ("ROI_hemi")
# and that show differences between hemispheres dependent on the direction ("ROI_direction")
# - folder 'models_LvR': containing all the generated models for this test
# DEPENDENCIES:
# - stat_table.R

check_LR_diff <- function(data_input, brm_form, direction, cond){
  # average over channels belonging to the same ROI, for each hemisphere separately
  data_ROIlr <- data_input %>%
    group_by(ID, group, trial, trigger, walkingvsstanding, direction, ROI, hemisphere) %>%
    summarise(HbO_pre = mean(HbO_pre),
              HbO_post = mean(HbO_post),
              HbO_post2 = mean(HbO_post2)) %>%
    ungroup()

  save_dir <- paste(getwd(), "models_LvR", cond, sep="/")
  dir.create(save_dir)
  ROIs <- c("M1", "PMC",  "SMA", "dlPFC", "PPC")
  ROI_hemi <- c()
  ROI_direction <- c()
  first <- TRUE
  for (i in ROIs) {
    data_selROI <-droplevels(subset(data_ROIlr, ROI == i)) 
    mod_selROI <- brm(brm_form, data = data_selROI, family = gaussian, iter = 10000, warmup = 500, chains = 4,control = list(adapt_delta = 0.99),
                      file = paste(paste(save_dir, "model", sep = "/"),i, sep = "_"))
    print(i)
    print(summary(mod_selROI))
    # print(plot(mod_selROI))
    # # print(pp_check(mod_selROI))
    # print(plot(conditional_effects(mod_selROI)))
    summary_ROI <- stat_table(mod_selROI)
    summary_ROI <- cbind(ROI = c(i, summary_ROI))
    if (first) {
      summary_table  <- summary_ROI
    } else {
      summary_table <- rbind(summary_table, summary_ROI)
    }
    first <- FALSE
    I <- fixef(mod_selROI, pars = "hemisphere1")
    if (!(0>I[1,3] & 0<I[1,4])) { # check if hemisphere is "significant"
      ROI_hemi <- c(ROI_hemi, i)
    }
    if (direction) {
      I <- fixef(mod_selROI, pars = c("direction1", "hemisphere1:direction1"))
      if (!(0>I[2,3] & 0<I[2,4])){
        ROI_hemi <- c(ROI_hemi, i)
        ROI_direction <- c(ROI_direction, i)
      } else if (!(0>I[1,3] & 0<I[1,4])) {
        ROI_direction <- c(ROI_direction, i)
      }
    }
  }
  ROI_LR <- list("ROI_hemi" = ROI_hemi, "ROI_direction" = ROI_direction)
  return(ROI_LR)
}
