*** Start after MI of lactate ***

*******************************************************************************
*** Results for all imputations ***
*******************************************************************************

// Loop over ALL imputations and summarise the results
cap prog drop imp_roc 
program define imp_roc, rclass 
	
	syntax [, imp(real 1) stype(real 1 )]
	
	gen iab_lact_`imp' = 0
	replace iab_lact_`imp' = 2 if _`imp'_lactate_baseline > 5

	gen IABP_score_imp_`imp' = (iam_age + iam_cve + iam_glu + iam_creat + iab_lact_`imp' + ia_pci)

	gen score_3cat_imp_`imp' = . 
	replace score_3cat_imp_`imp' = 1 if inlist(IABP_score_imp_`imp', 0, 1, 2)  
	replace score_3cat_imp_`imp' = 2 if inlist(IABP_score_imp_`imp', 3, 4)
	replace score_3cat_imp_`imp' = 3 if inlist(IABP_score_imp_`imp', 5, 6, 7, 8, 9)

	qui logit mort i.score_3cat_imp_`imp' if inlist(shock_type, `stype')
	predict p_mort_cs1_`imp', pr

	*su p_mort_cs1_`imp' if (score_3cat_imp_`imp' == 1) & inlist(shock_type, `stype')
	*su p_mort_cs1_`imp' if (score_3cat_imp_`imp' == 2) & inlist(shock_type, `stype')
	*su p_mort_cs1_`imp' if (score_3cat_imp_`imp' == 3) & inlist(shock_type, `stype')

	roctab mort p_mort_cs1_`imp' if inlist(shock_type, `stype')
	return scalar roc_cs`stype'_imp = r(area) 
	
	drop iab_lact_`imp' IABP_score_imp_`imp' score_3cat_imp_`imp' p_mort_cs1_`imp'
end 

forvalues k = 1/3 {
    cap drop roc_imp_final_`k'
    qui gen roc_imp_final_`k' = . 

    forvalues i = 1/`imp' {
        di "." _cont
        qui imp_roc, imp(`i') stype(`k')
        qui replace roc_imp_final_`k' = r(roc_cs`k'_imp) if _n == `i'
    }
}

local roc_amics = round(roc_cs1, 0.001) 
local roc_nonamics = round(roc_cs2, 0.001)
local roc_mixed = round(roc_cs3, 0.001)

qui su roc_imp_final_1
local mean = round(r(mean), 0.001) 
di in red "`= `mean''"
hist roc_imp_final_1,  ///
        fcolor(%0) lcolor(black) bin(10)  xtitle("C-statistic") /// 
        title("AMI-CS") aspect(1) xlab(#5, labsize(small) nogrid) /// 
        ylab( , labsize(small) nogrid) ytitle("") xline(`mean', lcolor(black) lwidth(thick)) /// 
        name("fig1", replace) ///
		note("C-statistic without imputed lactate: `=`roc_amics'' ", size(small)) 
		
	
qui su roc_imp_final_2
local mean = round(r(mean), 0.001) 
di in red "`= `mean''"
hist roc_imp_final_2,  ///
        fcolor(%0) lcolor(black) bin(10)  xtitle("C-statistic") /// 
        title("Non-AMI-CS") aspect(1) xlab(#5, labsize(small) nogrid) /// 
        ylab( , labsize(small) nogrid) ytitle("") xline(`mean', lcolor(black) lwidth(thick)) /// 
        name("fig2", replace) ///
		note("C-statistic without imputed lactate: `=`roc_nonamics'' ", size(small))

qui su roc_imp_final_3
local mean = round(r(mean), 0.001) 
di in red "`= `mean''"
hist roc_imp_final_3,  ///
        fcolor(%0) lcolor(black) bin(10) xtitle("C-statistic") /// 
        title("Mixed") aspect(1) xlab(#5, labsize(small) nogrid) /// 
        ylab( , labsize(small) nogrid) ytitle("") xline(`mean', lcolor(black) lwidth(thick)) /// 
        name("fig3", replace) ///
		note("C-statistic without imputed lactate: `=`roc_mixed'' ", size(small))
		
graph combine fig1 fig2 fig3, ycommon col(3)
graph export "Z:\Robert\CCCTN\results\fig_roc.jpg", as(jpg) width(4000) replace 

