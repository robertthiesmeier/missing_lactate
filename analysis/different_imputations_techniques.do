*** Linear regression with log-transformed values ***

mi set wide 
gen log_lac = log(lactate_baseline)
mi register imputed log_lac
mi register passive lactate_baseline IABP_score iab_lact_with_missing

mi impute reg log_lac /// 
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
	mort if inlist(shock_type, 1, 2, 3), add(100) rseed(18325)

mi passive: gen lact_aftermi = exp(log_lac) 
mi passive: replace iab_lact_with_missing = 0 
mi passive: replace iab_lact_with_missing = 2 if lact_aftermi > 5
mi passive: replace IABP_score = (iam_age + iam_cve + iam_glu + iam_creat + iab_lact_with_missing + ia_pci)

// crude 
mi estimate, post : logit mort IABP_score if inlist(shock_type, 1, 2, 3)
mi estimate, or

logit mort IABP_score if inlist(shock_type, 1, 2, 3), or 

// adjusted for sofa score and scai 
mi estimate, post : logit mort IABP_score sofa_max24  scai_shock if inlist(shock_type, 1, 2, 3)
mi estimate, or

logit mort IABP_score sofa_max24  scai_shock  if inlist(shock_type, 1, 2, 3), or 

**************************************************************************************************************
ssc install mi_impute_from, replace 
ssc install qrprocess, replace 
  
*** Quantile imputation using mi impute from ***
gen qlist = mod(_n-1, 99)/100 + 0.01 in 1/99 
qui levelsof qlist in 1/99 , local(levels)

qui qrprocess lactate_baseline /// 
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
	mort , q(`levels')
	
mat ib = e(b)
mat iV = e(V)

mi set wide 
mi register imputed lactate_baseline
mi register passive IABP_score iab_lact_with_missing

mi impute from lactate_baseline , add(100) b(ib) v(iV) imodel(qreg) seed(32125) force

mi passive: replace iab_lact_with_missing = 0 
mi passive: replace iab_lact_with_missing = 2 if lactate_baseline > 5
mi passive: replace IABP_score = (iam_age + iam_cve + iam_glu + iam_creat + iab_lact_with_missing + ia_pci)

// crude 
mi estimate, post : logit mort IABP_score if inlist(shock_type, 1, 2, 3)
mi estimate, or

logit mort IABP_score if inlist(shock_type, 1, 2, 3), or 

// adjusted for sofa score and scai 
mi estimate, post : logit mort IABP_score sofa_max24  scai_shock if inlist(shock_type, 1, 2, 3)
mi estimate, or

logit mort IABP_score sofa_max24  scai_shock  if inlist(shock_type, 1, 2, 3), or 

**************************************************************************************************************
*** Logistic regression for binary lactate value ***

gen lact_hi = (lactate_baseline > 5) if !missing(lactate_baseline)   // 1 = >5, 0 = ≤5, . = missing

mi set wide 
mi register passive iab_lact_with_missing IABP_score
mi register imputed lact_hi

mi impute logit lact_hi /// 
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
	if inlist(shock_type, 1, 2, 3), add(100) rseed(18325)

mi passive: replace iab_lact_with_missing = 0 
mi passive: replace iab_lact_with_missing = 2 if lact_hi == 1
mi passive: replace IABP_score = (iam_age + iam_cve + iam_glu + iam_creat + iab_lact_with_missing + ia_pci)
		
mi estimate, or: logit mort IABP_score if inlist(shock_type,1,2,3)

