# Addressing Informative Missingness of Lactate in Critical Care Registries: A Multiple Imputation Framework for Cardiogenic Shock Research

This repository contains Stata and R code illustrating the workflow of adressing missing lactate values in the Critical Care Cardiology Trials Network (CCCTN). We provide the analysis scripts that were used for the manuscript (Thiesmeier et al.). The workflow follow this structure: 

1. Identification of key predictors and internal validation: stepwise selection with backward elimination (file _1_)
2. Multiple Imputation of lactate values (file _2_)
3. IABP-SHOCK II risk score application: Imputation of missing lactate values and reclassification of the risk score

   Discriminitive performance check with

   a) C-statistic (file _3_)

   b) the net classification index (NRI) (file _4_)
5. Validation using analysis after MI (file _5_)

Code for different imputation techniques (linear regression with log-transformed values & quantile regression) is provided in an additional file.

Note that the orginal individual-level data cannot be shared.
