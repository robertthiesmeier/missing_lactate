############################################################
# Multiple imputation of missing lactate values
# Minimal R version for reproducible analysis
############################################################

# Packages
library(mice)    # multiple imputation
library(pROC)    # AUC / ROC
library(ggplot2) # plots

############################################################
# 1) Prepare data
############################################################

# Use your own data frame name here
# Example: dat <- read.csv("your_data.csv")

# Outcome: in-hospital mortality
dat$mort <- ifelse(dat$cicu_dispo == 3, 1, 0)

# Keep only shock types used in the manuscript
dat_mi <- subset(dat, shock_type %in% c(1, 2, 3))

# Quick check of missingness in variables used for imputation
vars_to_check <- c(
  "lactate_baseline", "sex", "past_medical_cv___3", "past_medical_cv___6",
  "past_medical___1", "hf", "new_acs", "mech_support",
  "sofa_max24", "carrest", "shock_type", "mort"
)

colSums(is.na(dat_mi[vars_to_check]))

# Score component based on observed lactate
dat_mi$iab_lact_with_missing <- ifelse(
  is.na(dat_mi$lactate_baseline), NA,
  ifelse(dat_mi$lactate_baseline > 5, 2, 0)
)

############################################################
# 2) Multiple imputation of lactate using PMM
############################################################

n_imp <- 100   # number of imputations
set.seed(18325)

# Keep only variables needed for the imputation model
imp_data <- dat_mi[, c(
  "lactate_baseline", "sex", "past_medical_cv___3", "past_medical_cv___6",
  "past_medical___1", "hf", "new_acs", "mech_support",
  "sofa_max24", "carrest", "shock_type", "mort",
  "iam_age", "iam_cve", "iam_glu", "iam_creat", "ia_pci",
  "score_3cat", "scai_shock"
)]

# Impute only lactate_baseline
meth <- make.method(imp_data)
meth[] <- ""
meth["lactate_baseline"] <- "pmm"

pred <- make.predictorMatrix(imp_data)
pred[,] <- 0
pred["lactate_baseline", ] <- 1
pred["lactate_baseline", "lactate_baseline"] <- 0

imp <- mice(
  imp_data,
  m = n_imp,
  method = meth,
  predictorMatrix = pred,
  maxit = 10,
  seed = 18325,
  printFlag = FALSE
)

############################################################
# 3) AUC after imputation, by shock type
############################################################

roc_cs1 <- numeric(n_imp)
roc_cs2 <- numeric(n_imp)
roc_cs3 <- numeric(n_imp)

for (i in 1:n_imp) {
  
  d <- complete(imp, i)
  
  # Recalculate lactate component and total score after imputation
  d$iab_lact <- ifelse(d$lactate_baseline > 5, 2, 0)
  
  d$IABP_score_imp <- d$iam_age + d$iam_cve + d$iam_glu +
    d$iam_creat + d$iab_lact + d$ia_pci
  
  d$score_3cat_imp <- ifelse(
    d$IABP_score_imp %in% c(0, 1, 2), 1,
    ifelse(d$IABP_score_imp %in% c(3, 4), 2, 3)
  )
  
  # Shock type 1
  d1 <- subset(d, shock_type == 1)
  fit1 <- glm(mort ~ factor(score_3cat_imp), data = d1, family = binomial())
  p1 <- predict(fit1, type = "response")
  roc_cs1[i] <- as.numeric(roc(d1$mort, p1, quiet = TRUE)$auc)
  
  # Shock type 2
  d2 <- subset(d, shock_type == 2)
  fit2 <- glm(mort ~ factor(score_3cat_imp), data = d2, family = binomial())
  p2 <- predict(fit2, type = "response")
  roc_cs2[i] <- as.numeric(roc(d2$mort, p2, quiet = TRUE)$auc)
  
  # Shock type 3
  d3 <- subset(d, shock_type == 3)
  fit3 <- glm(mort ~ factor(score_3cat_imp), data = d3, family = binomial())
  p3 <- predict(fit3, type = "response")
  roc_cs3[i] <- as.numeric(roc(d3$mort, p3, quiet = TRUE)$auc)
}

mean(roc_cs1)
mean(roc_cs2)
mean(roc_cs3)

# Histograms
par(mfrow = c(1, 3))

hist(roc_cs1, main = "AMI-CS", xlab = "C-statistic", col = "white", border = "black")
abline(v = mean(roc_cs1), lwd = 2)

hist(roc_cs2, main = "Non-AMI-CS", xlab = "C-statistic", col = "white", border = "black")
abline(v = mean(roc_cs2), lwd = 2)

hist(roc_cs3, main = "Mixed", xlab = "C-statistic", col = "white", border = "black")
abline(v = mean(roc_cs3), lwd = 2)

############################################################
# 4) NRI after imputation (optional section)
############################################################

nri_cs1 <- numeric(n_imp)
nri_cs2 <- numeric(n_imp)
nri_cs3 <- numeric(n_imp)

nri_plus_cs1 <- numeric(n_imp)
nri_plus_cs2 <- numeric(n_imp)
nri_plus_cs3 <- numeric(n_imp)

for (i in 1:n_imp) {
  
  d <- complete(imp, i)
  
  d$iab_lact <- ifelse(d$lactate_baseline > 5, 2, 0)
  
  d$IABP_score_imp <- d$iam_age + d$iam_cve + d$iam_glu +
    d$iam_creat + d$iab_lact + d$ia_pci
  
  d$score_3cat_imp <- ifelse(
    d$IABP_score_imp %in% c(0, 1, 2), 1,
    ifelse(d$IABP_score_imp %in% c(3, 4), 2, 3)
  )
  
  for (s in 1:3) {
    
    ds <- subset(d, shock_type == s)
    
    event <- ds$mort == 1
    
    up_event <- event & (ds$score_3cat_imp > ds$score_3cat)
    down_event <- event & (ds$score_3cat_imp < ds$score_3cat)
    
    up_nonevent <- !event & (ds$score_3cat_imp > ds$score_3cat)
    down_nonevent <- !event & (ds$score_3cat_imp < ds$score_3cat)
    
    pr_upcase <- mean(up_event[event], na.rm = TRUE)
    pr_downcase <- mean(down_event[event], na.rm = TRUE)
    pr_downcontrol <- mean(down_nonevent[!event], na.rm = TRUE)
    pr_upcontrol <- mean(up_nonevent[!event], na.rm = TRUE)
    
    nri_plus <- pr_upcase - pr_downcase
    nri_minus <- pr_downcontrol - pr_upcontrol
    nri_total <- nri_plus + nri_minus
    
    if (s == 1) {
      nri_cs1[i] <- 100 * nri_total
      nri_plus_cs1[i] <- 100 * nri_plus
    }
    if (s == 2) {
      nri_cs2[i] <- 100 * nri_total
      nri_plus_cs2[i] <- 100 * nri_plus
    }
    if (s == 3) {
      nri_cs3[i] <- 100 * nri_total
      nri_plus_cs3[i] <- 100 * nri_plus
    }
  }
}

mean(nri_cs1)
mean(nri_cs2)
mean(nri_cs3)

mean(nri_plus_cs1)
mean(nri_plus_cs2)
mean(nri_plus_cs3)

# Simple plots
par(mfrow = c(2, 3))

hist(nri_cs1, main = "Total NRI: AMI-CS", xlab = "NRI (%)", col = "white", border = "black")
abline(v = mean(nri_cs1), lwd = 2)

hist(nri_cs2, main = "Total NRI: Non-AMI-CS", xlab = "NRI (%)", col = "white", border = "black")
abline(v = mean(nri_cs2), lwd = 2)

hist(nri_cs3, main = "Total NRI: Mixed", xlab = "NRI (%)", col = "white", border = "black")
abline(v = mean(nri_cs3), lwd = 2)

hist(nri_plus_cs1, main = "NRI+ cases: AMI-CS", xlab = "NRI+ (%)", col = "white", border = "black")
abline(v = mean(nri_plus_cs1), lwd = 2)

hist(nri_plus_cs2, main = "NRI+ cases: Non-AMI-CS", xlab = "NRI+ (%)", col = "white", border = "black")
abline(v = mean(nri_plus_cs2), lwd = 2)

hist(nri_plus_cs3, main = "NRI+ cases: Mixed", xlab = "NRI+ (%)", col = "white", border = "black")
abline(v = mean(nri_plus_cs3), lwd = 2)

############################################################
# 5) Main regression models
############################################################

# Recalculate passive variables inside each imputed dataset
completed_list <- vector("list", n_imp)

for (i in 1:n_imp) {
  d <- complete(imp, i)
  
  d$iab_lact_with_missing <- ifelse(d$lactate_baseline > 5, 2, 0)
  
  d$IABP_score <- d$iam_age + d$iam_cve + d$iam_glu +
    d$iam_creat + d$iab_lact_with_missing + d$ia_pci
  
  completed_list[[i]] <- d
}

# Put back into mids object for pooled analysis
imp2 <- as.mids(do.call(rbind, completed_list))

# Crude pooled model
fit_mi_crude <- with(
  imp2,
  glm(mort ~ IABP_score, family = binomial(), subset = shock_type %in% c(1, 2, 3))
)
summary(pool(fit_mi_crude), conf.int = TRUE, exponentiate = TRUE)

# Non-imputed analysis: treat missing lactate as 0
dat_cc <- dat_mi
dat_cc$iab_lact_with_missing <- ifelse(
  is.na(dat_cc$lactate_baseline), 0,
  ifelse(dat_cc$lactate_baseline > 5, 2, 0)
)
dat_cc$IABP_score <- dat_cc$iam_age + dat_cc$iam_cve + dat_cc$iam_glu +
  dat_cc$iam_creat + dat_cc$iab_lact_with_missing + dat_cc$ia_pci

fit_cc_crude <- glm(
  mort ~ IABP_score,
  data = dat_cc,
  family = binomial()
)
exp(cbind(OR = coef(fit_cc_crude), confint(fit_cc_crude)))

# Adjusted pooled model
fit_mi_adj <- with(
  imp2,
  glm(mort ~ IABP_score + sofa_max24 + scai_shock,
      family = binomial(),
      subset = shock_type %in% c(1, 2, 3))
)
summary(pool(fit_mi_adj), conf.int = TRUE, exponentiate = TRUE)

# Adjusted non-imputed analysis
fit_cc_adj <- glm(
  mort ~ IABP_score + sofa_max24 + scai_shock,
  data = dat_cc,
  family = binomial()
)
exp(cbind(OR = coef(fit_cc_adj), confint(fit_cc_adj)))