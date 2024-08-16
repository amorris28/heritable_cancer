RK_get_regPHS <- function(Age, status , Xvars){
  # requires:
  # Age [vector, nevents x 1] - age at diagnosis for observed events, age at last follow-up for censored events 
  # status [vector, nevents x 1] - status indicator, 0 = censored events, 1 = observed events 
  # Xvars [matrix, nevents x nvars] - Matrix of predictors, rows = events (censored and observed), columns = variables 
  # returns: 
  # beta [vector, nvars x 1] - SNP coefficients of lasso-regularized Cox regression model (selected at minimum deviance)
  # R dependencies: 
  # survival
  # glmnet
  
  require(survival)
  require(glmnet)
  
  # convert Xvars to matrix 
  Xvars = as.matrix(Xvars) 
  
  # Cross-validated regularized Cox-PH model with varying lambda
  cat(sprintf('%s | estimating regularized cox....\n', Sys.time()))
  cv.fit = cv.glmnet(Xvars, Surv(Age, status), family = 'cox', maxit = 1000)
  cat(sprintf('%s | done!\n', Sys.time()))
  
  # select coefficients at minumum deviance 
  beta = coef(cv.fit, s = cv.fit$lambda.min)
  
  return(beta)
}