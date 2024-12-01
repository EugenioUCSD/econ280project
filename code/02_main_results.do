* Set your working directory
global root "C:\Users\eugen\OneDrive - UC San Diego\Desktop\UC San Diego\2024\Fall 2024\Computation 280\Replication\ECON_280 - Part 4"

* Load the dataset
use "$root\data\data_reg.dta", clear

*
* This file uses Chodorow-Reich & Wieland's two sample two stage least squares method
* ts2sls.ado is in the directory
* ts2sls.sthlp, the help file, is also in the directory


clear all
set more off
set matsize 800

local install_packages = "no"

if "`install_packages'" == "yes" {
ssc install reghdfe
ssc install ranktest
ssc install moremata
ssc install carryforward
ssc install ivreg2
ssc install ftools
ssc install ivreghdfe
ssc install binscatter
ssc install estout
}
*****************************************************************************
****************Step 0: Set the local name with your name ********************
*****************************************************************************


*The dataset contains these variables year state quarter mean_une statecode date constant infl_reg rp qt_bartik_sa
use "$root\data\data_reg.dta", clear


*****************************************************************************
****************Step 1: Calculate present values ****************************
* The point here is to compute the PV approach for unemployment and prices. We do this for a benchmark value of beta = 0.99, and robustness values beta = 0.9, 0.95. The variables rp_sum_* and u_sum_* contain the outcome.
*****************************************************************************
*****************************************************************************
*****************************************************************************

* Set value of beta
global beta 0.99


* Vary truncation length
foreach truncation_length in 10 20 30 40 50 60 {

  * Calculate present value of unemployment
  quietly capture drop u_sum_`truncation_length'
  quietly generate u_sum_`truncation_length' = mean_une

  forvalues i = 1/`truncation_length' {

  	quietly replace u_sum_`truncation_length' = u_sum_`truncation_length' + $beta^`i'*F`i'.mean_une

  }


  * Calculate present value of relative prices. Similar procedure than for unemployment.
  quietly capture drop rp_sum_`truncation_length'
  quietly generate rp_sum_`truncation_length' = rp

  forvalues i = 1/`truncation_length' {

  	quietly replace rp_sum_`truncation_length' = rp_sum_`truncation_length' + $beta^`i'*F`i'.rp

  }

}

* Calculate present values for beta = 0.90, 0.95
foreach beta_local in 75 80 85 90 95 {

	local beta_local_adj = `beta_local' / 100

	quietly capture drop u_sum_20_`beta_local'
	quietly generate u_sum_20_`beta_local' = mean_une

	forvalues i = 1/20 {

		quietly replace u_sum_20_`beta_local' = u_sum_20_`beta_local' + `beta_local_adj'^`i'*F`i'.mean_une

	}
	quietly capture drop rp_sum_20_`beta_local'
	quietly generate rp_sum_20_`beta_local' = rp

	forvalues i = 1/20 {

		quietly replace rp_sum_20_`beta_local' = rp_sum_20_`beta_local' + `beta_local_adj'^`i'*F`i'.rp

	}
}



***********************************************************************
***********************************************************************
***********************************************************************
**************************** Step 2: Run Regression ******************
***********************************************************************
***********************************************************************
***********************************************************************


* Generate some useful variables. infl_reg_sign just changes the sign of inflation so that the coefficients have the sign as in the paper. infl_reg_time_agg divides by 4 because the model is written in quarterly terms and the data are 12-month inflation rates. d20_qt_bartik_sa takes 20 quarter differences of the seasonally adjusted bartik instrument. L4_* variables are 4 lags of the respective variables.

generate infl_reg_sign = -1 * infl_reg
generate infl_reg_time_agg = -1 * infl_reg / 4
generate d20_qt_bartik_sa = qt_bartik_sa - L20.qt_bartik_sa
generate L4_d20_qt_bartik_sa = L4.d20_qt_bartik_sa
generate L4_rp = L4.rp
generate L4_mean_une = L4.mean_une

save "$root\data\data_reg_R.dta", replace


** Table 1: Full sample estimates of the Regional Phillips Curve
* The first three regressions estimate a regression of inflation on the lag of unemployment with different fixed effects. The fourth one runs an IV specification where unemployment is instrumented with the lag of the bartik instrument.
eststo clear

	quietly eststo: reghdfe infl_reg_sign L4.mean_une L4.rp, absorb(constant) cluster(statecode)
	quietly eststo: reghdfe infl_reg_sign L4.mean_une L4.rp, absorb(i.statecode) cluster(statecode)
	quietly eststo: reghdfe infl_reg_sign L4.mean_une L4.rp, absorb(i.statecode i.date) cluster(statecode)
	eststo: ivreghdfe infl_reg_sign (L4.mean_une = L4.d20_qt_bartik_sa) L4.rp, absorb(i.statecode i.date) cluster(statecode)

esttab, se keep(L4.mean_une)

disp("Number of observations is ")
count if !missing(infl_reg_sign) & !missing(L4.mean_une) & !missing(L4.rp)

estout using "$root/output/psi_full_sample.tex", style(tex) keep(L4.mean_une) varlabels(L4.mean_une "$\psi$") cells(b(star fmt(%9.3f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(L4.mean_une)
eststo clear



* Estimates of Kappa + Lambda
* Here we run regressions using the ts2sls function of Chodorow-Reich and Weiland (2019). The function lets us to run 2sls regressions in samples with gaps.
quietly {

	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), absorb(constant) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), absorb(statecode) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_mean_une L4_rp) i.date, absorb(statecode) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_d20_qt_bartik_sa L4_rp) i.date, absorb(statecode) cluster(statecode)

}
esttab, se keep(u_sum_20 rp_sum_20)

disp("Number of observations is ")
count if !missing(infl_reg_time_agg) & !missing(L4_mean_une) & !missing(L4_rp) & !missing(u_sum_20) & !missing(rp_sum_20)

estout using "$root/output/kappa_full_sample.tex", style(tex) keep(u_sum_20) varlabels(u_sum_20 "$\kappa$") cells(b(star fmt(%9.4f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(u_sum_20)
estout using "$root/output/lambda_full_sample.tex", style(tex) keep(rp_sum_20) varlabels(rp_sum_20 "$\lambda$") cells(b(star fmt(%9.4f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(rp_sum_20)


* Estimates of Kappa + Lambda
* Here we run regressions using the ts2sls function of Chodorow-Reich and Weiland (2019). The function lets us to run 2sls regressions in samples with gaps.
quietly {

	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), absorb(constant) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), absorb(statecode) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_mean_une L4_rp) i.date, absorb(statecode) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_d20_qt_bartik_sa L4_rp) i.date, absorb(statecode) cluster(statecode)

}
esttab, se keep(u_sum_20 rp_sum_20)

disp("Number of observations is ")
count if !missing(infl_reg_time_agg) & !missing(L4_mean_une) & !missing(L4_rp) & !missing(u_sum_20) & !missing(rp_sum_20)

estout using "$root/output/kappa_full_sample.tex", style(tex) keep(u_sum_20) varlabels(u_sum_20 "$\kappa$") cells(b(star fmt(%9.4f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(u_sum_20)
estout using "$root/output/lambda_full_sample.tex", style(tex) keep(rp_sum_20) varlabels(rp_sum_20 "$\lambda$") cells(b(star fmt(%9.4f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(rp_sum_20)


* Kappa as Beta varies. Different estimates of kappa depending on the discount rate. For the benchmark value (0.99), and for 0.80, and 0.75.
eststo clear
quietly {

	capture drop endo_var_1
	capture drop endo_var_2

	generate endo_var_1 = u_sum_20
	generate endo_var_2 = rp_sum_20
	eststo: ts2sls infl_reg_time_agg (endo_var_1 endo_var_2 = L4_d20_qt_bartik_sa L4_rp) i.date, absorb(statecode) cluster(statecode)

	replace endo_var_1 = u_sum_20_80
	replace endo_var_2 = rp_sum_20_80
	eststo: ts2sls infl_reg_time_agg (endo_var_1 endo_var_2 = L4_d20_qt_bartik_sa L4_rp) i.date, absorb(statecode) cluster(statecode)

	replace endo_var_1 = u_sum_20_75
	replace endo_var_2 = rp_sum_20_75
	eststo: ts2sls infl_reg_time_agg (endo_var_1 endo_var_2 = L4_d20_qt_bartik_sa L4_rp) i.date, absorb(statecode) cluster(statecode)

	drop endo_var_1
	drop endo_var_2

}
esttab, se keep(endo_var_1 endo_var_2)

estout using "$root/output/beta_vary.tex", style(tex) keep(endo_var_1) varlabels(endo_var_1 "$\kappa$") cells(b(star fmt(%9.4f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(endo_var_1)
