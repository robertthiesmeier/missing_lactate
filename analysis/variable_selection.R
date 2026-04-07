# =============================================================================
# Predictors of Missing Lactate Values in Cardiogenic Shock (CS) Patients
# -----------------------------------------------------------------------------
# Purpose: Identify clinical predictors of missing lactate measurements using
#           stepwise logistic regression with BIC-penalised selection (k = log n
#           approximated as 6.6).  Model performance is evaluated on a held-out
#           test set using the C-statistic (AUC).
# Data: CCCTN registry, waves 1–7
# =============================================================================

rm(list = ls())  

library(haven) 
library(dplyr)  
library(janitor)   
library(caTools)    
library(pROC)       

##### Load data #####
# Primary dataset: all predictors with a low fraction of missing values
data1 <- read_dta("Z:/Robert/CCCTN/data/ccctn1to7_missLactate_allpredictors_low_miss.dta")

# Extended dataset: includes predictors with higher missingness (e.g. SOFA variables)
data2 <- read_dta("Z:/Robert/CCCTN/data/ccctn1to7_missLactate_allpredictors.dta")


##### Descriptives #####
# Count of cardiogenic shock (CS) patients in each dataset
cat("CS distribution in data1:\n"); print(table(data1$cs))
cat("CS distribution in data2:\n"); print(table(data2$cs))

# Missing lactate values stratified by CS status (row percentages with counts)
tab_missing_lact <- data1 %>%
  tabyl(cs, miss_lact) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns()

cat("\nMissing lactate by CS status:\n")
print(tab_missing_lact)

###### Train/test split (80 / 20, stratified by sex) ######

split <- sample.split(data1$sex, SplitRatio = 0.8)
traindata <- subset(data1, split == TRUE)
testdata <- subset(data1, split == FALSE)

cat(sprintf("\nTraining set: n = %d | Test set: n = %d\n",
            nrow(traindata), nrow(testdata)))


###### Variable selection via backward stepwise logistic regression ######
# Penalty k = 6.6 approximates BIC at typical sample sizes, favouring parsimonious models over AIC (k = 2)

stepwise_model <- step(
  glm(
    miss_lact ~
      age +
      as.factor(sex) +
      as.factor(smoking_status) +
      as.factor(past_medical_cv___1) +  
      as.factor(past_medical_cv___2) +
      as.factor(past_medical_cv___3) +
      as.factor(past_medical_cv___4) +
      as.factor(past_medical_cv___6) +
      as.factor(past_medical___1) +
      as.factor(past_medical___2) +
      as.factor(past_medical___3) +
      as.factor(hf) +                   
      as.factor(new_acs) +           
      as.factor(mech_support) + 
      as.factor(carrest) +  
      as.factor(sofa_max24) + 
      eGFR_bl +  
      as.factor(cs),
    data = traindata,
    family = binomial(link = "logit")
  ),
  direction = "backward",
  trace = 0, k = 6.6
)

print(stepwise_model)

###### Model evaluation on test set ######

# Predicted probabilities from the selected model
testdata$pred_prob <- predict(stepwise_model, newdata = testdata, type = "response")

# Binary classification at 0.5 probability threshold
testdata$pred_class <- as.integer(testdata$pred_prob > 0.5)

cat("\nPredicted class distribution (test set):\n")
print(table(Predicted = testdata$pred_class))

cat("\nObserved outcome distribution (test set):\n")
print(table(Observed = testdata$miss_lact))
roc_obj <- roc(testdata$miss_lact, testdata$pred_class)
auc_val <- round(auc(roc_obj), 3)

cat(sprintf("\nC-statistic (AUC) on test set: %s\n", auc_val))
