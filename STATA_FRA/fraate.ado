/*
Estimate Average Treatment Effect after Full Regression Adjustment Pre-processing

Inputs:
  varlist: outcome name of outcome whose ATE is being estimated
  treatlevel: value of W corresponding to the treatment whose ATE is being estimated
  control: value of W corresponding to the control group

Outputs:
   Point estimate and standard error

*/

program fraate, rclass
version 16
syntax varlist [if] [in], treatment(varlist) treatlevel(integer) controllevel(integer)

	qui levels `treatment', local(treat_levels)
	local outcome_cols `varlist'
	local treatlevel = `treatlevel'
	local controllevel = `controllevel'

 	qui gen u = u_y_`treatlevel' - u_y_`controllevel'
	
	qui su u
	qui local mean = round(r(mean),.0001)
	qui local se = round(r(sd)/sqrt(_N),.0001)
	
	di "Point Estimate: `mean'"
	di "Standard Error: `se'"

end 
