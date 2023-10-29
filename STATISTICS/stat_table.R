# STAT_TABLE.R
# function that creates a summary for all fixed effect of the tested ROI with the probability intervals.
# INPUT:
# - model: brms-model (as created by run_model)
# OUTPUT: 
# - ROI_table: summary for all fixed effect of the tested ROI with the probability intervals

stat_table <- function(model) {
  
# fixed effects
fe95 <- rownames_to_column(as.data.frame(brms::fixef(model)))
fe99 <- rownames_to_column(as.data.frame(brms::fixef(model, probs = c(0.005, 0.995))))
fe <- bind_cols(fe95, fe99[,4:5])
fe <- fe[-c(3)]
fe <- fe %>% rename(factor = rowname)
# calculate probabilities
for (i in c(1:length(fe$factor))){
  H <- hypothesis(model, paste(fe$factor[i], ">0"))
  prob <- H[["hypothesis"]][["Post.Prob"]]
  if (prob<0.5){
    prob <- 1-prob
  }
  fe$Prob[i] <- prob
}

if (length(fe$factor)==1){
  em <- data.frame()
} else if (fe$factor[2]=="hemisphere1" | any(grepl("TF", fe$factor)) | any(grepl("NFOGQ", fe$factor))){
  em <- data.frame()
} else if (length(fe$factor)==2){
  # emmeans
  em95 <- as.data.frame(emmeans(model, "group"))
  em99 <- as.data.frame(emmeans(model, "group"), level = .99)
  em <- bind_cols(em95, em99[,3:4])
  # calculate probabilities
  for (i in c(1:length(em$group))){
    if (em$group[i]=="HC"|em$group[i]=="nrl"){
      H <- hypothesis(model, "Intercept + group1 > 0")
    } else if (em$group[i]=="PD"| em$group[i]=="FOG") {
      H <- hypothesis(model, "Intercept - group1 > 0")
    }
      prob <- H[["hypothesis"]][["Post.Prob"]]
    if (prob<0.5){
      prob <- 1-prob
    }
    em$Prob[i] <- prob
  }
  colnames(em) <- colnames(fe)

} else if (length(fe$factor)==3){
  # only keep the intercept
  fe <- fe[1,]
  # test for between group differences
  em95 <- as.data.frame(emmeans(model, specs = pairwise ~ group))
  em99 <- as.data.frame(emmeans(model, specs = pairwise ~ group, level = .99))
  em <- bind_cols(em95, em99[,4:5])
  em[em=="."] <- ""
  em <- em %>%
    unite("factor", group:contrast, sep = "", na.rm = TRUE)
  em$Prob <- NA # just for now
  colnames(em) <- colnames(fe)
} else if (length(fe$factor)==4){
  em95 <- as.data.frame(emmeans(model, specs = ~group|trigger))
  em99 <- as.data.frame(emmeans(model, specs = ~group|trigger, level = .99))
  em <- bind_cols(em95, em99[,4:5])
  em <- em %>%
    unite('factor', group:trigger, remove = TRUE)
  colnames(em) <- colnames(fe)
}

# combine 
ROI_table <- rbind(fe, em)

# check whether it HDI includes 0
sign95 <- !(ROI_table$Q2.5 < 0 & ROI_table$Q97.5 >0) 
sign99 <- !(ROI_table$Q0.5 < 0 & ROI_table$Q99.5 >0)
if (any(is.na(ROI_table$Prob))){
  ROI_table$sign <- ifelse(sign95, ifelse(sign99, '**', '*'), '')
} else {
  ROI_table$sign <- ifelse(ROI_table$Prob>0.95, ifelse(sign95, ifelse(sign99, '**', '*'), '.'), '')
}
# 



return(ROI_table)
}