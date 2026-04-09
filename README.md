# Addressing Informative Missingness of Lactate in Critical Care Registries: A Multiple Imputation Framework for Cardiogenic Shock Research

This repository contains Stata and R code illustrating a workflow for adressing missing lactate values in the Critical Care Cardiology Trials Network (CCCTN). We provide the analysis Stata scripts that were used for the manuscript (Thiesmeier et al. under review), as well as a R simplified template of the key analytica steps. The code can be adapted to users data. 

### Multiple imputtaion workflow for missing lactate values

1. [Identification of key predictors and internal validation](analysis/_1_variable_selection.do)
2. [Multiple Imputation of lactate values with predictive mean matching](analysis/_2_imputation_lactate.do)
3. Assessment of the imputation effect on study-specific settings. In our example, this was conducted assessing the reclassification of IABP-SHOCK II risk score when lactate is imputed through our workflow vs naive imputations (all missing data imputed to 0).

   Study-specific files: 
   a) [AUC](analysis/_3_a_AUC_IABP_SHOCK_II_risk_score_application.do)
   
   b) [Net Re-classification Index (NRI)](analysis/_3_b_NRI_IABP_SHOCK_II_risk_score_application.do)
   
   c) [Validation using analysis after MI](analysis/_4_IABP_SHOCK_II_risk_score_application_ANALYSIS.do)

*Additional files*: Code for different imputation techniques (linear regression with log-transformed values & quantile regression) is provided in an additional [file](analysis/different_imputations_techniques.do). Additional material for the [descriptive analysis](analysis/descriptives.do) is also provided.

 ### [R code to replicate the workflow](analysis/R_imputation_template)
