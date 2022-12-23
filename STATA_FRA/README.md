# Installation instructions:
	
In your STATA command box use `sysdir\'92` to find your personal ado folder file path

Copy `fra.ado`, `fraate.ado`, and `fraalate.ado` your personal ado directory

Restart STATA

After reopening the functions associated with each .ado file will be available for use


# Description of Functions

`fra.ado`: Perform Flexible Regression Adjustment Pre-Processing

Inputs:

	varlist: name of outcome(s)
	
	treatment: name of treatment(s)
	
	covariates: name of covariate(s)
	
	n_folds: number of folds for sample splitting
	
	method: regression method used for regression adjustment ("linear" or "rf") or custom ML model 	supplied by user. Custom ML model should be a function of the form: [function] `outcome' 	`covariates'. Output should be compatible with predict function. Please note if using "rf", package 	"rforest" must be installed prior to use.	
	
Outputs:

	Original data with new variables of the form - 
	
	m_`outcomename'_`treatmentname': fitted value of conditional expectation, E[outcome | X, treatment] 	for the outcome and treatment named
	
	u_`outcomename'_`treatmentname': "influence function" for mean potential outcome,
	
	E[outcome(treatment)]. Mean of this column is the regression adjusted estimator for 	E[outcome(treatment)] and variance-covariance matrix of these columns is asymptotically valid 	estimator of covariance matrix of the regression adjusted point estimates
	
\
`fraate.ado`: Estimate Average Treatment Effect after Full Regression Adjustment Pre-processing

Inputs:

	varlist: outcome name of outcome whose ATE is being estimated
	
	treatlevel: value of W corresponding to the treatment whose ATE is being estimated
	
	control: value of W corresponding to the control group
	

Outputs:

	Point estimate and standard error


\
`fraalate.ado`: Estimate Local Average Treatment Effect when Experiment Assignment W is Instrument for Treatment

Inputs:

	varlist: name of outcome whose LATE is being estimated
	
	endogcol: treatment, which experiment assignment instruments for
	
	treatlevel: value of W corresponding to treatment whose LATE is being estimated
	
	controllevel: value of W corresponding to control group
	

Outputs:

 	Point estimate and standard error
