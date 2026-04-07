**** Start after imputing missing lactate values ****

*** Multiple imputation analysis (on the first 100 imputations) ***
qui mi passive: replace iab_lact_with_missing = 0 
qui mi passive: replace iab_lact_with_missing = 2 if lactate_baseline > 5
qui mi passive: replace IABP_score = (iam_age + iam_cve + iam_glu + iam_creat + iab_lact_with_missing + ia_pci)

// crude analysis not adjusting for sofa score and scai 
mi estimate, post imp(1/50): logit mort IABP_score if inlist(shock_type, 1, 2, 3)
mi estimate, or

// analysis model without the imputed values treating missing lactate = 0
logit mort IABP_score if inlist(shock_type, 1, 2, 3), or 

// adjusted for sofa score and scai 
mi estimate, post imp(1/50): logit mort IABP_score sofa_max24  scai_shock if inlist(shock_type, 1, 2, 3)
mi estimate, or

// analysis model without the imputed values treating missing lactate = 0
logit mort IABP_score sofa_max24  scai_shock  if inlist(shock_type, 1, 2, 3), or 
