****************************************************************************************
*** Multiple imputation of missing lactate values ***
*** This MI procedure is illustartive and used for Figure 2 in the manuscript.
****************************************************************************************
ssc install grc1leg2, replace
use "Z:\Robert\CCCTN\data\ccctn1to7_missLactate_allpredictors_low_miss.dta", clear
  
*** set up the MI environment ***
mi set wide 
mi register imputed lactate_baseline

// imputation with the final model
mi impute pmm lactate_baseline sex past_medical_cv___3 past_medical_cv___6 past_medical___1 hf new_acs mech_support i.sofa_max24 carrest cs, add(1) force rseed(14325) knn(10)

tw ///
	(hist _1_lactate_baseline if _mi_miss == 1 &  cs == 0 & _1_lactate_baseline < 10, bin(30) color(black%30) lcolor(%0)) ///
	(hist lactate_baseline if _mi_miss == 0 & cs == 0 & lactate_baseline < 10, bin(30) fcolor(%0) lcolor(black)) ,name(plot1, replace) legend(label(1 "Single imputation") label(2 "Observed values") pos(6) row(1)) ytitle("") xtitle("Lactate at baseline (mmol/L)") title("No CS") xlab(, nogrid) ylab(, nogrid)

tw ///
	(hist _1_lactate_baseline if _mi_miss == 1 &  cs == 1 & _1_lactate_baseline < 10, bin(30) color(black%30) lcolor(%0)) ///
	(hist lactate_baseline if _mi_miss == 0 & cs == 1 & lactate_baseline < 10, bin(30) fcolor(%0) lcolor(black)) ,name(plot2, replace) legend(label(1 "Single imputation") label(2 "Observed values") pos(6) row(1)) ytitle("") xtitle("Lactate at baseline (mmol/L)") title("CS") xlab(, nogrid) ylab(, nogrid)
	
grc1leg2  plot1 plot2, ycommon col(2)

graph export "Z:\Robert\CCCTN\results\fig2_imputaion_lactate.png", replace width(4000)
