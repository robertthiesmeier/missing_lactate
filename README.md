# Addressing Informative Missingness of Lactate in Critical Care Registries: A Multiple Imputation Framework for Cardiogenic Shock Research

This repository contains Stata and R code illustrating the workflow of adressing missing lactate values in the Critical Care Cardiology Trials Network (CCCTN). We provide the analysis scripts that were used for the manuscript (Thiesmeier et al.). The workflow follow this structure: 

1. Identification of key predictors and internal validation: stepwise selection with backward elimination
2. Multiple Imputation of lactate values
3. IABP-SHOCK II risk score application: Imputation of missing lactate values and reclassification of the risk score

   Discriminitive performance check with

   a) C-statistic

   b) the net classification index (NRI) 
5. Validation using analysis after MI
