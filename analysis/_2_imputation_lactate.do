****************************************************************************************
*** Multiple imputation of missing lactate values ***
****************************************************************************************

*** Predictive Mean Matching (PMM) is used in the manuscript ***

// generate the outcome: mortality (in-hospital) 
gen mort = . 
replace mort = 0 if (cicu_dispo!=3)
replace mort = 1 if (cicu_dispo==3)
tab mort

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

*** use PMM to create multiple imputed datasets
mi impute pmm lactate_baseline /// 
		age /// 
		sex /// 
		smoking_status /// 
		past_medical_cv___1 /// 
		past_medical_cv___2 /// 
		past_medical_cv___3 /// 
		past_medical_cv___4 /// 
		past_medical_cv___6 /// 
		past_medical___1 /// 
		past_medical___2 /// 
		past_medical___3 /// 
		hf /// 
		new_acs /// 
		mech_support /// 
		carrest /// 
		sofa_max24 /// 
		eGFR_bl /// 
		i.shock_type /// 
		glu_baseline /// 
		creat_baseline /// 
		ia_pci ///
		mort ///
		if inlist(shock_type, 1, 2, 3), add(`imp') rseed(18325) knn(10)

****************************************************************************************

