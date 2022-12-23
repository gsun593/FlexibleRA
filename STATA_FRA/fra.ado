/*
Perform Flexible Regression Adjustment Pre-Processing

Inputs:
  varlist: name of outcome(s)
  treatment: name of treatment(s)
  covariates: name of covariate(s)
  n_folds: number of folds for sample splitting
  method: regression method used for regression adjustment ("linear" or "rf") or custom ML model supplied by user. Custom ML model should be a function of the form: [function] `outcome' `covariates'. Output should be compatible with predict function. Please note if using "rf", package "rforest" must be installed prior to use.

Outputs:
  Original data with new variables of the form - 
  m_`outcomename'_`treatmentname': fitted value of conditional expectation, E[outcome | X, treatment] for the outcome and treatment named
  u_`outcomename'_`treatmentname': "influence function" for mean potential outcome, E[outcome(treatment)]. Mean of this column is the regression adjusted estimator for E[outcome(treatment)] and variance-covariance matrix of these columns is asymptotically valid estimator of covariance matrix of the regression adjusted point estimates

*/

program fra, rclass
version 16
syntax varlist [if] [in], treatment(varlist) [covariates(varlist)] nfolds(integer) seed(integer) method(name) [mlfuncoptions(string)] [endog]

	qui set seed `seed'
	qui gen rand = runiformint(0, _N)
	qui sort rand 
	gen ord = _n
	qui drop rand
	egen folds = seq(), from(1) to(`nfolds')
	qui levels `treatment', local(treat_levels)
	local outcome_cols `varlist'

	// Linear Regression
	qui if("`method'"== "linear"){
		foreach y in `outcome_cols'{
		foreach treat in `treat_levels'{
			gen m_`y'_`treat'= .
			forval f = 1/`nfolds'{
				reg `y' `covariates' if folds != `f' & w == `treat'
				predict m_`y'_`treat'_`f' if folds == `f' 
				replace m_`y'_`treat' = m_`y'_`treat'_`f' if folds == `f'
				drop m_`y'_`treat'_`f' 
			}
		}
		}
	
	}
	
	// Random Forest
	qui if("`method'" == "rf"){
		foreach y in `outcome_cols'{
		foreach treat in `treat_levels'{
			gen m_`y'_`treat'=.
			forval f=1/`nfolds'{
				display "`y' `treat'"
				rforest `y' `covariates' if folds != `f' & w == `treat', type(class)
				predict m_`y'_`treat'_`f' if folds == `f' 
				replace m_`y'_`treat' = m_`y'_`treat'_`f' if folds == `f'
				drop m_`y'_`treat'_`f' 
			}
		}
		}
	}
	
	// Custom ML function supplied
	qui if ("`method'" != "linear" & "`method'" != "rf") {
		foreach y in `outcome_cols'{
		foreach treat in `treat_levels'{
			gen m_`y'_`treat'=.
			forval f=1/`nfolds'{
				display "`y' `treat'"
				`method' `y' `covariates' if folds != `f' & w == `treat', `MLfuncoptions'
				predict m_`y'_`treat'_`f' if folds == `f' 
				replace m_`y'_`treat' = m_`y'_`treat'_`f' if folds == `f'
				drop m_`y'_`treat'_`f' 
			}
		}
		}
	}
	
	// For each outcome/treatment pair, create column for influence function of the form:
	// 1 / prob(treatment) * (Y - E[Y|X,treatment]) * 1{treatment} + E[Y|X,treatment]
	qui levels `treatment', local(treat_levels)
	qui local outcome_cols `varlist'
	
	gen total = _N
	qui foreach treat in `treat_levels'{
		bysort w: egen num_`treat' = count(w)
		local prop_`treat' = num_`treat'/total
		foreach y in `outcome_cols'{
			gen u_`y'_`treat' = (1/`prop_`treat'')*(`y'- m_`y'_`treat') if w==`treat'
			replace  u_`y'_`treat' = 0 if w!= `treat'
			replace  u_`y'_`treat' =  u_`y'_`treat'+ m_`y'_`treat' 
		}
	}
	
	qui drop num* total
	
end 
 
 
 
 
 
 
 
 
 
 
 
 