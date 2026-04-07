# Addressing Informative Missingness of Lactate in Critical Care Registries: A Multiple Imputation Framework for Cardiogenic Shock Research

This repository contains Stata and R code illustrating the workflow of adressing missing lactate values in the Critical Care Cardiology Trials Network (CCCTN). We provide the analysis scripts that were used for the manuscript (Thiesmeier et al.). The workflow follow this structure including the analysis files: 

1. [Identification of key predictors and internal validation: stepwise selection with backward elimination](analysis/_1_variable_selection.R)
2. [Multiple Imputation of lactate values with predictive mean matching](analysis/_2_imputation_lactate.do)
3. IABP-SHOCK II risk score application: Imputation of missing lactate values and reclassification of the risk score

   Discriminitive performance check with

   a) [AUC](analysis/_3_a_AUC_IABP_SHOCK_II_risk_score_application.do)

   b) [the Net Re-classification Index (NRI)](analysis/_3_b_NRI_IABP_SHOCK_II_risk_score_application.do)
5. [Validation using analysis after MI](analysis/_4_IABP_SHOCK_II_risk_score_application_ANALYSIS.do)

Code for different imputation techniques (linear regression with log-transformed values & quantile regression) is provided in an additional [file](analysis/different_imputations_techniques.do). Additional material for the [descriptive analysis](analysis/descriptives.do) is also provided.

Note that the orginal individual-level data cannot be shared.
