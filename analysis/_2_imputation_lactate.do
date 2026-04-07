****************************************************************************************
*** Multiple imputation of missing lactate values ***
****************************************************************************************

// generate the outcome: mortality (in-hospital) (cicu_dispo==3)
gen mort = . 
replace mort = 0 if (cicu_dispo!=3)
replace mort = 1 if (cicu_dispo==3)
tab mort

// check missigness of predictors
mdesc lactate_baseline sex past_medical_cv___3 past_medical_cv___6 past_medical___1 hf new_acs mech_support sofa_max24 carrest  if inlist(shock_type, 1, 2, 3)

*** Multiple imputation ***
local imp = 100 // change the number of imputations here

// generate iab_lact with missing lactate
cap drop iab_lact_with_missing
gen iab_lact_with_missing = . 
replace iab_lact_with_missing = 0 if lactate_baseline != . 
replace iab_lact_with_missing = 2 if lactate_baseline > 5 & lactate_baseline != .
tab iab_lact_with_missing

*** set up the MI environment ***
mi set wide 
mi register imputed lactate_baseline
mi register passive IABP_score iab_lact_with_missing

mi impute pmm lactate_baseline ///
	sex past_medical_cv___3 past_medical_cv___6 ///
	past_medical___1 hf new_acs mech_support i.sofa_max24 ///
	carrest i.shock_type mort ///
	if inlist(shock_type, 1, 2, 3), knn(10) add(`imp') rseed(18325)
