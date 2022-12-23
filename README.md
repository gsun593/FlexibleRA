# FlexibleRA
R and STATA Functions to Perform Flexible Regression Adjustment


R Code contains the following functions and a simulation example

```
Perform Flexible Regression Adjustment Pre-Processing
FRA(dat, outcome_cols, treat_col, covariate_cols, n_folds, method)
Inputs:
  dat: data frame with outcomes, treatments, and covariates
  outcome_cols: column names for outcomes of interest
  treat_col: column name of treatment
  covariate_cols: column names of covariates
  n_folds: number of folds for sample splitting
  method: regression method used for regression adjustment
  ML_func: Custom ML model supplied by user. Should be of the form ML_func(formula, data).
           Output should have a predict function.
  
Output:
  dat_with_FRA: original dataframe with extra columns of the form
      'm_{otcuome name}_{treatment name}': fitted value of conditional expectation,
        E[outcome | X, treatment] for the outcome and treatment named
      'u_{outcome name}_{treatment name}': "influence function" for mean potential outcome,
        E[outcome(treatment)]. Mean of this column is the regression adjusted estimator for
        E[outcome(treatment)] and variance-covariance matrix of these columns is asymptotically
        valid estimator of covariance matrix of the regression adjusted point estimates




Estimate Average Treatment Effect after Full Regression Adjustment Pre-processing
FRA_ATE(dat_with_FRA, outcome_col, treat_lvl, ctrl_lvl)
Inputs:
  dat_with_FRA: dataframe with regression adjusted columns
  outcome_col: name of outcome whose ATE is being estimated
  treat_lvl: value of W corresponding to "treatment"
  ctrl_lvl: value of W corresponding to "control"
  
Output:
  Vector with point estimate and standard error





Estimate function of potential outcome means after regression adjustment
FRA_theta(para_func, dat_with_FRA, outcome_treats)
Inputs:
  param_func: function of potential outcome means being estimated
  dat_with_FRA: dataframe with regression adjusted columns
  outcome_treats: vector of strings of the form '{outcome name}_{treatment name}' which
  are the inputs into param_func
Output:
  Vector with point estimate and standard error
```

STATA Package Found in foler titled STATA_FRA
