# Addressing Informative Missingness of Lactate in Critical Care Registries: A Multiple Imputation Framework for Cardiogenic Shock Research

This repository contains Stata and R code illustrating the workflow of adressing missing lactate values in the Critical Care Cardiology Trials Network (CCCTN). We provide the analysis scripts that were used for the manuscript. 

We also provide a generic code script that can be adapted for specific project needs. The code is commented so that each step can be followed.

The outline for the workflor presented in the manuscript is as follows: 
1. Pre-select clinically relevant predictors of missing biomarker
2. Select predictors based on variable selection procedure
3. Generate M number of imputations for missing biomarker conditional on selected predictors
4. Validation with the Net-reclassification Index
5. Validation using AUC
6. Validation using analysis after MI
