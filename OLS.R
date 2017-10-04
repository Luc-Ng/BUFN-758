###########################################################################
# This function estimates the regression coefficients in multivariate
# linear regressions. It also computes the respective standard errors. 
#
# INPUTS: - y: vector of dependent variable (Tx1)
#         - X: Matrix of independent variables (T x K) including a constant
# OUTPUTS:- A list containing:
#         - Beta vector (K x 1): vector of regression coefficients
#         - Yhat vector (T x 1): vector of fitted values
#         - res vector (T x 1): vector of residuals
#         - varcovarBeta matrix (K x K): variance-covariance matrix of
#            the regression coefficients
#         - varRes (1 x 1): Estimate of the residual variance
#         - tstat (K x 1): t-statistic for the beta coefficients
#
###########################################################################

# Obtaining the coefficients estimates
OLS <- function(X,y) {

  
  beta=solve(t(X) %*% X) %*% t(X) %*% y 
           
  # Obtaining the fitted values from the regression
  yhat= X%*% beta         
 
  # Obtaining the number of observations and the number of regressors 
  # (including the the constant) employed in the analysis
  nobs=dim(X)[1]
  nreg=dim(X)[2]
  
  # Compute the vector of residuals
  res=y-yhat                           
 
  # Computing estimates of the error variance
  shat=as.numeric((t(res) %*% res)/(nobs-nreg))          
       
  # Computing the variance co-variance matrix of the beta estimates
  varcovar= shat * solve(t(X) %*% X)                
 
  # Computing simple t-statistic for the null of beta equal to zero slide 16 of the multiple linear regression
  tstat=beta/sqrt(diag(varcovar))
 
 
  # Computing the residual sum of squares
  RSS=t(res) %*% res             
 
  # Computing the total sum of squares
  TSS=t((y-mean(y))) %*% (y-mean(y)) 
 
  # Computing the R-squared
  Rq=1-(RSS/TSS)      
 
  # Computing the adjusted R-squared
  Rqad=1- ((nobs-1)/(nobs-nreg))*(1-Rq);                          
 
  S=list( beta = beta, yhat = yhat, res = res, varcovarBeta = varcovar,
         varRes = shat, tstat = tstat, Rq = Rq, Rqadj = Rqad)
} 
         
         
         
         