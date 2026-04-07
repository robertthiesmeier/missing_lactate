# Addressing Informative Missingness of Lactate in Critical Care Registries: A Multiple Imputation Framework for Cardiogenic Shock Research

This repository contains Stata and R code illustrating the workflow of adressing missing lactate values in the Critical Care Cardiology Trials Network (CCCTN). We provide the analysis scripts that were used for the manuscript. 

We also provide a generic code script that can be adapted for specific project needs. The code is commented so that each step can be followed.

The outline for the workflor presented in the manuscript is as follows: 
1. Pre-select clinically relevant predictors of missing biomarker
2. Identification of key predictors and internal validation: stepwise selection with backward elimination
3. Multiple Imputation of lactate values
4. IABP-SHOCK II risk score application: Imputation of missing lactate values and reclassification of the risk score
   a. Discriminitive performance check: C-statistics & the net classification index (NRI) 
5. Validation using analysis after MI
