/*
Estimate Local Average Treatment Effect when Experiment Assignment W is Instrument for Treatment

Inputs:
   varlist: name of outcome whose LATE is being estimated
   endogcol: treatment, which experiment assignment instruments for
   treatlevel: value of W corresponding to treatment whose LATE is being estimated
   controllevel: value of W corresponding to control group
   
Outputs:
   Point estimate and standard error

*/

program fralate
version 16
syntax varlist [if] [in], endogvar(varlist) treatlevel(integer) controllevel(integer)

	local outcomecol `varlist'
	local treatlevel = `treatlevel'
	local controllevel = `controllevel'

	gen u_num = u_`outcomecol'_`treatlevel' - u_`outcomecol'_`controllevel'
	gen u_denom = u_`endogcol'_`treatlevel' - u_`endogcol'_`controllevel'
	qui su u_num
	local pe_num = r(mean)
	qui su u_denom
	local pe_denom = r(mean)
	local pe = round(`pe_num'/`pe_denom',.0001)
	qui corr u_num u_denom, covariance 
	matrix VCV = r(C)
	local D1 = 1/`pe_denom' 
	local D2 = -(`pe_num'/(`pe_denom'^2))
	matrix D = J(2,1,.)
	matrix D[1,1] = `D1'
	matrix D[2,1] = `D2'
	
	matrix DT = J(1,2,.)
	matrix DT[1,1] = `D1'
	matrix DT[1,2] = `D2'
	
	matrix temp = DT*VCV*D
	local se = round(sqrt(temp[1,1]),.0001)
	
	di "Point Estimate: `pe'"
	di "Standard Error: `se'"

end
