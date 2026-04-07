/*==============================================================================
 Purpose: Create a table with descriptives statistics startified by CS status for
 all variables included and Figure 1 (distribution of lactate)

 Data: CCCTN, wave 1-7
==============================================================================*/

**********************************************************************
***** Table 1: Descriptives ***** 
  
dtable ///
	i.cs /// 
	i.shock_type /// 
	age ///
	i.sex ///
	i.smoking_status /// medical history
	i.past_medical_cv___1 ///
	i.past_medical_cv___2 ///
	i.past_medical_cv___3 ///
	i.past_medical_cv___4 ///
	i.past_medical___1 ///
	i.past_medical___2 ///
	i.past_medical___3 ///
	i.hf /// 
	i.mech_support /// 
	i.new_acs /// 
	i.carrest /// 
	sofa_max24 /// 
	eGFR_bl /// 
	, cont(age eGFR_bl) by(miss_lact) /// 
    title("Table1. Baseline characteristics by status of lactate measurement availability (column %)") /// 
    export(table1.docx, replace)  titlestyles(font(,bold))

**********************************************************************
***** Figure 1: Distribution of lactate stratified by CS status *****
su lactate_baseline if cs == 0 , d
local median_nocs = r(p50)
su lactate_baseline if cs == 1 , d
local median_cs = r(p50)
two (hist lactate_baseline if cs == 0 & lactate_baseline < 10, fcolor(%0) lcolor(black) bin(50)) ///
	(hist lactate_baseline if cs == 1 & lactate_baseline < 10, fcolor(black%30) lcolor(%0) bin(50)), legend(label(1 "no CS") label(2 "CS") row(1) pos(6)) name(fig1, replace) xlabel(0.0 `=`median_nocs'' `=`median_cs'' 4.0 5.0 6.0 7.0 8.0 9.0 10, nogrid) ylab(, nogrid) xtitle("Lactate at baseline (mmol/L)") ytitle("Distribution") ///
	xline(`=`median_nocs'', lwidth(medthick) lpattern(dash) lcolor(black)) ///
	xline(`=`median_cs'', lpattern(solid) lcolor(black) lwidth(medthick)) 
