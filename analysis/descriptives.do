/*==============================================================================
 Purpose: Distribution of imputed lactate values in different subtypes of cardiogenic shock
 Data: CCCTN, wave 1-4
==============================================================================*/

// start after imputations have been performed
**********************************************************************

	// AMI-CS
	local plot ""
	forv k = 1/10{ 
		local hist = "(hist _`k'_lactate_baseline if _mi_miss == 1 & shock_type == 1 & _`k'_lactate_baseline < 10, bin(20) color(black%5) lcolor(%0))"
		local plot = "`plot' `hist'"
	}

	tw ///
		`plot' ///
		(hist lactate_baseline if shock_type == 1 & lactate_baseline < 10, bin(20) fcolor(%0) lcolor(black) lwidth(medthick)) ,name(plot1, replace) legend(label(1 "Imputed values") label(11 "Observed values") order(1 11) pos(6) row(1)) ytitle("") xtitle("Lactate at baseline (mmol/L)") title("AMI-CS") xlab(, nogrid) ylab(, nogrid)  nodraw
		
	// Non-AMI CS
	local plot ""
	forv k = 1/10{ 
		local hist = "(hist _`k'_lactate_baseline if _mi_miss == 1 & shock_type == 2 & _`k'_lactate_baseline < 10, bin(20) color(black%5) lcolor(%0))"
		local plot = "`plot' `hist'"
	}

	tw ///
		`plot' ///
		(hist lactate_baseline if shock_type == 2 & lactate_baseline < 10, bin(20) fcolor(%0) lcolor(black) lwidth(medthick)) ,name(plot2, replace) legend(label(1 "Imputed values") label(11 "Observed values") order(1 11) pos(6) row(1)) ytitle("") xtitle("Lactate at baseline (mmol/L)") title("Non-AMI CS") xlab(, nogrid) ylab(, nogrid)  nodraw

	// MIxed
	local plot ""
	forv k = 1/10{ 
		local hist = "(hist _`k'_lactate_baseline if _mi_miss == 1 & shock_type == 3 & _`k'_lactate_baseline < 10, bin(20) color(black%5) lcolor(%0))"
		local plot = "`plot' `hist'"
	}

	tw ///
		`plot' ///
		(hist lactate_baseline if shock_type == 3 & lactate_baseline < 10, bin(20) fcolor(%0) lcolor(black) lwidth(medthick)) ,name(plot3, replace) legend(label(1 "Imputed values") label(11 "Observed values") order(1 11) pos(6) row(1)) ytitle("") xtitle("Lactate at baseline (mmol/L)") title("Mixed") xlab(, nogrid) ylab(, nogrid)  nodraw

	grc1leg2  plot1 plot2 plot3, ycommon col(2)
