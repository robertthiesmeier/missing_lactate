/*==============================================================================
Purpose: 
1. Multiple imputation for missing lactate values using PMM
2. NRI applied to the IABP score to validate

Data: CCCTN, wave 1-4 (as in orginal IABP score analysis)
==============================================================================*/

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

*** ----------------------------------------------------------------------- ***

use "/Users/robert/Library/CloudStorage/OneDrive-KarolinskaInstitutet/PhD/Research/TIMI/TIMI_research/method_comp/CCCTN/data/ccctn1to6_missLactate.dta", clear

******************************************************************
***** Building the IABP score as in the orginal publication ******

drop if Campaign == 5 | Campaign == 6 | Campaign == 7  

*** Additional variables needed for the imputation ***
gen hf = . 
replace hf = 1 if (dx_initial_primary == 9)
replace hf = 0 if (dx_initial_primary != 9)
label var hf "Primary presenting cardiac problem: HF"

recode dx_initial_primary (2 = 1) (1 3/20 = 0), gen(new_acs)
label var new_acs "Primary presenting cardiac problem: ACS"

gen carrest = 0 if arrest == 0 | arrest == 2
replace carrest = 1 if arrest ==1 
label var carrest "Cardiac arrest prior to hospitalization"

// AMI (ACS) subjects
cap drop dx_initial_ACS
gen dx_initial_ACS = .
replace dx_initial_ACS = 1 if dx_initial_primary_org == 2 | dx_initial_secondary == 2
replace dx_initial_ACS = 0 if !missing(dx_initial_primary_org) & !missing(dx_initial_secondary) & dx_initial_ACS != 1

// MS with any CS component
cap drop shock_MSwCScomp
gen shock_MSwCScomp = 0
replace shock_MSwCScomp = 1 if shock == 1 & shock_etiol == 4 & shock_mixed_what___1 == 1
tab shock_MSwCScomp

cap drop shock_pureCS
gen shock_pureCS = 0
replace shock_pureCS = 1 if shock == 1 & shock_etiol == 1
tab shock_pureCS

// pure CS or MS with CS+(DS/HS)
gen shock_allCS = 0
replace shock_allCS = 1 if shock_pureCS == 1 | shock_MSwCScomp == 1
tab shock_allCS

gen shock_allCS_AtAdm = 0
replace shock_allCS_AtAdm = 1 if Shock_AtAdm == 1 & shock_allCS == 1

gen dx_AMICS = .
replace dx_AMICS = 1 if dx_initial_ACS == 1 & shock_allCS_AtAdm == 1
replace dx_AMICS = 0 if !missing(dx_initial_ACS) & dx_AMICS != 1
tab dx_AMICS

// non-ANI CS 
cap drop non_AMICS
gen non_AMICS = . 
replace non_AMICS = 0 
replace non_AMICS = 1 if dx_AMICS != 1 & shock_pureCS == 1 & shock_pureCS != 0
tab non_AMICS 

// Mixed shock
tab shock_etiol if shock_etiol == 4

***** IABP score caclulation  ******

cap drop iam_age
gen iam_age = .
replace iam_age = 1 if age > 73
replace iam_age = 0 if age <= 73 & age != .

cap drop iam_cve
gen iam_cve = .
replace iam_cve = 2 if past_medical_cv___4 == 1
replace iam_cve = 0 if past_medical_cv___4 != 1 & past_medical_cv___4 != .

cap drop iam_glu
gen iam_glu = .
replace iam_glu = 1 if glu_baseline > 191
replace iam_glu = 0 if glu_baseline <= 191 & glu_baseline != .

cap drop iam_creat
gen iam_creat = .
replace iam_creat = 1 if creat_baseline > 1.5
replace iam_creat = 0 if creat_baseline <= 1.5 & creat_baseline != .

cap drop iab_lact
gen iab_lact = 0
replace iab_lact = 2 if lactate_baseline > 5 & lactate_baseline != .

cap drop ia_pci
gen ia_pci = 0
replace ia_pci = 2 if ppci == 1 & pci_success == 0 & pci_success != .
replace ia_pci = . if ppci == 1 & pci_success == .

cap drop IABP_score
gen IABP_score = (iam_age + iam_cve + iam_glu + iam_creat + iab_lact + ia_pci)

label variable IABP_score "IABP Shock II score (Zero values were assigned for missing Lactate values only)"

tab IABP_score

count if IABP_score == .
drop if IABP_score == . 

*******************************************************************************
*** Variables used to stratify the score ****

* Acute coronary syndrome (Acute MI) subjects from both primary and secondary;
cap drop dx_initial_ACS
gen dx_initial_ACS = 0
replace dx_initial_ACS = 1 if dx_initial_primary_org == 2 | dx_initial_secondary == 2

gen shock_type = .
replace shock_type = 1 if shock_etiol == 1 & dx_initial_ACS == 1   // CS with AMI (AMI-CS)
replace shock_type = 2 if shock_etiol == 1 & dx_initial_ACS == 0   // CS without AMI (nonAMI-CS)
replace shock_type = 3 if shock_etiol == 4                         // Mixed shock

label define shock_type 1 "AMI-CS" 2 "Non-AMI CS" 3 "Mixed"
label values shock_type shock_type

*******************************************************************************
*** variable of the score as in orginal publication ***

gen score_3cat = . 
replace score_3cat = 1 if inlist(IABP_score, 0, 1, 2)  
replace score_3cat = 2 if inlist(IABP_score, 3, 4)
replace score_3cat = 3 if inlist(IABP_score, 5, 6, 7, 8, 9)
tab score_3cat

tab  score_3cat shock_type , col

bysort score_3cat shock_type: gen count = _N
bysort shock_type: egen total_count = total(count)
gen percent = (count / total_count) * 100

*******************************************************************************
*** Visualisation of the score ***
graph bar (mean) percent,  over(score_3cat) over(shock_type) stack asyvars percent legend(label(1 "Low (0-2)") label(2 "Intermediate (3-4)") label(3 "High (5-9)") pos(12) row(1) size(small)) ytitle("Distribution of Risk Category by IABP SHOCK-II score", size(small)) ylab(#10, labsize(small)) bar(1, color(ltblue)) bar(2, color(gold)) bar(3, color(red%80)) blabel(bar, format(%9.1f) position(center)) title("IABP-SCHOCK II Risk Score Distribution", size(small)) name(fig1, replace)
graph export "Z:\Robert\CCCTN\results\fig_bargraph_orginal_score.jpg", as(jpg) width(4000) replace

*******************************************************************************
*** missing lactate ***
count if  inlist(shock_type, 1, 2, 3)
local total = r(N)
count if lactate_baseline == . & inlist(shock_type, 1, 2, 3)
local miss = r(N)

di "% missing lactate in CS: " %3.2f `miss'/`total' *100 
su lactate_baseline if  inlist(shock_type, 1, 2, 3), detail

count if inlist(shock_type, 1, 2, 3) & lactate_baseline !=. 
count if lactate_baseline > 5 & inlist(shock_type, 1, 2, 3) & lactate_baseline !=. 

*******************************************************************************
*** Imputation of lactate ***
*******************************************************************************

// gen the outcome: mortality (in-hospital) (cicu_dispo==3)
gen mort = . 
replace mort = 0 if (cicu_dispo!=3)
replace mort = 1 if (cicu_dispo==3)
tab mort

// check missigness of predictors
mdesc lactate_baseline sex past_medical_cv___3 past_medical_cv___6 past_medical___1 hf new_acs mech_support sofa_max24 carrest  if inlist(shock_type, 1, 2, 3)

*******************************************************************************
*** Multiple imputation ***
*******************************************************************************
local imp = 100 // change the number of imputations here

// generate iab_lact with missing lactate
cap drop iab_lact_with_missing
gen iab_lact_with_missing = . 
replace iab_lact_with_missing = 0 if lactate_baseline != . 
replace iab_lact_with_missing = 2 if lactate_baseline > 5 & lactate_baseline != .
tab iab_lact_with_missing

*** set up the environment ***
mi set wide 
mi register imputed lactate_baseline
mi register passive IABP_score iab_lact_with_missing

mi impute pmm lactate_baseline ///
	sex past_medical_cv___3 past_medical_cv___6 ///
	past_medical___1 hf new_acs mech_support i.sofa_max24 ///
	carrest i.shock_type mort ///
	if inlist(shock_type, 1, 2, 3), knn(10) add(`imp') rseed(18325) // orginal seed: 18325

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
 
local reps = 5 // reps for bootstrap
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

graph export "Z:\Robert\CCCTN\results\NRI_fig.jpg", as(jpg) width(4000) replace

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
