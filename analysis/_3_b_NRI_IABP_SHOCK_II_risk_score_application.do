*** Start after MI of missing lactate ***

********** NRI program *********** (see comments next to the code indicating the logic)
cap prog drop nri 
program define nri, rclass

syntax [, imp(real 1) stype(real 1)] // for the loop over number of imputations and CS types

	preserve
		qui keep if inlist(shock_type, `stype')

		cap drop event up_event down_event up_nonevent down_nonevent
			 
		*** new calculation of risk score after imputation ***
		// generate the new category with imputed lactate (score is recalculated each time)
		qui gen iab_lact_`imp' = 0
		replace iab_lact_`imp' = 2 if _`imp'_lactate_baseline > 5
		
		qui gen IABP_score_imp_`imp' = ///
			(iam_age + iam_cve + iam_glu + iam_creat + iab_lact_`imp' + ia_pci)

		qui gen score_3cat_imp_`imp' = . 
		qui replace score_3cat_imp_`imp' = 1 if inlist(IABP_score_imp_`imp', 0, 1, 2)  
		qui replace score_3cat_imp_`imp' = 2 if inlist(IABP_score_imp_`imp', 3, 4)
		qui replace score_3cat_imp_`imp' = 3 if inlist(IABP_score_imp_`imp', 5, 6, 7, 8, 9)

		qui logit mort i.score_3cat_imp_`imp'
		predict p_mort_cs1_`imp', pr
		
		*** NRI estimatation ***
		gen event = mort == 1

		* movements for cases
		// check if imputed score is higher than old score
		gen up_event = event & (score_3cat_imp_`imp' > score_3cat)
		
		// check if imputed score lower than old score
		gen down_event = event & (score_3cat_imp_`imp' < score_3cat)

		* movements for controls (same as above but for non-events)
		gen up_nonevent = !event & (score_3cat_imp_`imp' > score_3cat)
		gen down_nonevent = !event & (score_3cat_imp_`imp' < score_3cat)

		* calculate the probabilities
		qui su up_event if event
		scalar pr_upcase = r(mean)

		qui su down_event if event
		scalar pr_downcase = r(mean)

		qui su down_nonevent if !event
		scalar pr_downcontrol = r(mean)

		qui su up_nonevent if !event
		scalar pr_upcontrol = r(mean)

		*** NRI components ***
		scalar nri_plus = pr_upcase - pr_downcase
		scalar nri_minus = pr_downcontrol - pr_upcontrol
		scalar nri = nri_plus + nri_minus

		// table with results
		di "---------------------------------"
		di "Pr(Up|Case) = " %6.5f pr_upcase
		di "Pr(Down|Case) = " %6.5f pr_downcase
		di "Pr(Down|Control) = " %6.5f pr_downcontrol
		di "Pr(Up|Control) = " %6.5f pr_upcontrol
		
		di "---------------------------------"
		di "NRI+ (cases) = " %6.5f nri_plus
		di "NRI- (controls) = " %6.5f nri_minus
		di "Total NRI = " %6.5f nri
		
		*** return final results ***
		ret scalar nri_plus = nri_plus
		ret scalar nri_minus = nri_minus
		ret scalar nri_total = nri
		
		drop iab_lact_`imp' IABP_score_imp_`imp' score_3cat_imp_`imp' p_mort_cs1_`imp'
	restore
end

************************************************************
*** Begin: NRI caclulation ***
************************************************************

* looping over all imputations to get the average NRI (call the nri program)
// here we get total NRI and NRI+ (cases) 
// (cases who move up in the risk categories due to "good" lactate imputation)

cap drop nri_cs*
forv k = 1/3{
	qui gen nri_cs`k' = . // total NRI
	qui gen nri_cs`k'_nriplus = . // NRI+
}

forv k = 1/`imp'{
	qui nri , imp(`k') stype(1)
	qui replace nri_cs1 = r(nri_total)*100 if _n == `k'
	qui replace nri_cs1_nriplus = r(nri_plus)*100 if _n == `k'
	
	qui nri , imp(`k') stype(2)
	qui replace nri_cs2 = r(nri_total)*100 if _n == `k'
	qui replace nri_cs2_nriplus = r(nri_plus)*100 if _n == `k'

	qui nri , imp(`k') stype(3)
	qui replace nri_cs3 = r(nri_total)*100 if _n == `k'
	qui replace nri_cs3_nriplus = r(nri_plus)*100 if _n == `k'
}
* optional: bootstrap to get CIs and p values 

cap drop b_nri_*
forv k = 1/3{
	gen b_nri_total_cs`k' = . // collect estimate
	gen b_nri_pos_cs`k' = . 
	gen b_nri_total_cs`k'_p = . // collect p value for NRI and NRI+
	gen b_nri_pos_cs`k'_p = .
}
 
local reps = 5 // reps for bootstrap (increase for final analysis)
forv k = 1/`imp'{
	
	// CS1 event
	qui bootstrap r(nri_total) r(nri_plus), reps(`reps') seed(1): nri , imp(`k') stype(1)
		qui replace b_nri_total_cs1 = r(table)[1,1]*100 if _n == `k'
		qui replace b_nri_pos_cs1 = r(table)[1,2]*100 if _n == `k'
		qui replace b_nri_total_cs1_p = r(table)[4,1] if _n == `k'
		qui replace b_nri_pos_cs1_p = r(table)[4,2] if _n == `k'
		
	// CS2 event
	qui bootstrap r(nri_total) r(nri_plus), reps(`reps') seed(1): nri , imp(`k') stype(2)
		qui replace b_nri_total_cs2 = r(table)[1,1]*100 if _n == `k'
		qui replace b_nri_pos_cs2 = r(table)[1,2]*100 if _n == `k'
		qui replace b_nri_total_cs2_p = r(table)[4,1] if _n == `k'
		qui replace b_nri_pos_cs2_p = r(table)[4,2] if _n == `k'
	
	// CS3 event
	qui bootstrap r(nri_total) r(nri_plus), reps(`reps') seed(1): nri , imp(`k') stype(3)
		qui replace b_nri_total_cs3 = r(table)[1,1]*100 if _n == `k'
		qui replace b_nri_pos_cs3 = r(table)[1,2]*100 if _n == `k'
		qui replace b_nri_total_cs3_p = r(table)[4,1] if _n == `k'
		qui replace b_nri_pos_cs3_p = r(table)[4,2] if _n == `k'
}

*** Graphs **** get a total of 6 graphs and put them all together

*** total NRI ***
su nri_cs1
local mean : di %4.3f r(mean)
hist nri_cs1,  ///
        fcolor(%0) lcolor(black) bin(10) xtitle("NRI") /// 
        aspect(1) xlab(#5, labsize(small) nogrid) /// 
        ylab( , labsize(small) nogrid) ytitle("") xtitle("") title("AMI-CS") /// 
		xline(`mean', lcolor(black) lwidth(thick)) /// 
		note("Average NRI: `mean'%", size(small)) name(nri_cs1, replace)

su nri_cs2
local mean : di %4.3f r(mean)
hist nri_cs2,  ///
        fcolor(%0) lcolor(black) bin(10) xtitle("NRI") /// 
        aspect(1) xlab(#5, labsize(small) nogrid) /// 
        ylab( , labsize(small) nogrid) ytitle("") xtitle("") title("Non-AMI-CS") /// 
		xline(`mean', lcolor(black) lwidth(thick)) /// 
		note("Average NRI: `mean'%", size(small)) name(nri_cs2, replace)
	
su nri_cs3
local mean: di %4.3f r(mean)
hist nri_cs3,  ///
        fcolor(%0) lcolor(black) bin(10) xtitle("NRI") /// 
        aspect(1) xlab(#5, labsize(small) nogrid) /// 
        ylab( , labsize(small) nogrid) ytitle("") xtitle("") title("Mixed") /// 
		xline(`mean', lcolor(black) lwidth(thick)) /// 
		note("Average NRI: `mean'%", size(small)) name(nri_cs3, replace)
	
*** NRI pos cases ***
su nri_cs1_nriplus
local mean: di %4.3f r(mean)
hist nri_cs1_nriplus,  ///
        fcolor(%0) lcolor(black) bin(10) xtitle("NRI") /// 
        aspect(1) xlab(#5, labsize(small) nogrid) /// 
        ylab( , labsize(small) nogrid) ytitle("") xtitle("") title("AMI-CS") /// 
		xline(`mean', lcolor(black) lwidth(thick)) /// 
		note("Average NRI+ cases: `mean'%", size(small)) name(nri_cs1_nriplus, replace)

su nri_cs2_nriplus
local mean: di %4.3f r(mean)
hist nri_cs2_nriplus,  ///
        fcolor(%0) lcolor(black) bin(10) xtitle("NRI") /// 
        aspect(1) xlab(#5, labsize(small) nogrid) /// 
        ylab( , labsize(small) nogrid) ytitle("") xtitle("") title("Non-AMI-CS") /// 
		xline(`mean', lcolor(black) lwidth(thick)) /// 
		note("Average NRI+ cases: `mean'%", size(small)) name(nri_cs2_nriplus, replace)
	
su nri_cs3_nriplus
local mean: di %4.3f r(mean)
hist nri_cs3_nriplus,  ///
        fcolor(%0) lcolor(black) bin(10) xtitle("NRI") /// 
        aspect(1) xlab(#5, labsize(small) nogrid) /// 
        ylab( , labsize(small) nogrid) ytitle("") xtitle("") title("Mixed") /// 
		xline(`mean', lcolor(black) lwidth(thick)) /// 
		note("Average NRI+ cases: `mean'%", size(small)) name(nri_cs3_nriplus, replace)
		
*** combine all grpahs
graph combine ///
	nri_cs1 nri_cs2 nri_cs3 /// 
	nri_cs1_nriplus nri_cs2_nriplus nri_cs3_nriplus, ///
	row(2) name(combined_nri, replace) /// 
	title("{bf:Total NRI}", size(large) pos(12)) /// 
	note("{bf:NRI+ cases}", size(large) pos(0))

*** table for p-values in all three CS categories ***
// average p-value for NRI value and NRI+ value (results from the bootstrap)
forv k = 1/3{
	qui count if b_nri_total_cs`k'_p < 0.05 
	di as text "NRI TOTAL: Prop. of p-value < 0.05 for CS `k': " as result in red r(N)/100 
	
	qui count if b_nri_pos_cs`k'_p < 0.05 
	di as text "NRI+: Prop. of p-value < 0.05 for CS `k': " as result in red r(N)/100 
	di in text "-------------------------------------"
}
************************************************************
*** END: NRI calculation ***
************************************************************
