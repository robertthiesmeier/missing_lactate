/*****************************
Variable selection procedure for mising lactate framework
******************************/ 

* =============================================================================
* Stepwise backward logistic regression + AUC on test set
* =============================================================================

*** Load dataset here (CCCTN (wave 1-7) dataset was used in the manuscript)
use "your_data.dta", clear

*** Train/test split (80/20)
set seed 12345
gen rand = runiform()
gen train = (rand <= 0.8)

*** Backward stepwise logistic regression (BIC-penalised)

sw logit miss_lact age i.sex i.smoking_status ///
    i.past_medical_cv___1 i.past_medical_cv___2 ///
    i.past_medical_cv___3 i.past_medical_cv___4 ///
    i.past_medical_cv___6 i.past_medical___1 ///
    i.past_medical___2 i.past_medical___3 ///
    i.hf i.new_acs i.mech_support i.carrest ///
    i.sofa_max24 eGFR_bl i.cs ///
    if train == 1, pr(0.001) pe(0.001)

*** Predict on test set
predict pred_prob if train == 0, pr
gen pred_class = (pred_prob > 0.5) if train == 0

*** AUC on test set
roctab miss_lact pred_prob if train == 0
