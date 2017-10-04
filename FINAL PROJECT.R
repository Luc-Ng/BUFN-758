###########################################################################
# Define the OLS function 
# This function estimates the alpha and beta coefficients in univariate 
# linear regressions. It also computes the respective standard errors. 
#
# INPUTS: - y: vector of dependent variable (Tx1)
#         - x: vector of independent variable (Tx1)
#         (the constant is automatically added)
# OUTPUTS: - a_est= estimate of the alpha coefficient
#          - b_est= estimate of the beta coefficient       
#          - se_alpha= standard error of the alpha coefficient
#          - se_beta= standard error of the beta coefficient
# Author: Alberto Rossi                                                   #
###########################################################################

univariate_ols <- function(y,x) {
  numerator=sum(x*y)-length(y)*mean(y)*mean(x);    # numerator of the beta coefficient
  denominator=sum(x^2)-length(y)*mean(x)^2;        # denominator of the beta coefficient
  b_est=numerator/denominator;                     # beta coefficient estimate
  a_est=mean(y)-b_est*mean(x);                     # alpha coefficient estimate
  residuals=y-a_est-b_est*x;                       # residuals
  s_squared=(1/(length(y)-2))*sum(residuals^2);    # unbiased estimate of the error variance 
  sum_x_squared=sum(x^2);                         
  sum_x_minus_xbar=sum((x-mean(x))^2);
  
  
  # Standard error of alpha
  se_alpha=sqrt(s_squared)*(sqrt((1/length(y))*sum_x_squared)/sqrt(sum_x_minus_xbar))
  # Standard error of beta
  se_beta=sqrt(s_squared)*(1/sqrt(sum_x_minus_xbar))
  #this will send all values back into the main program
  output=c(a_est=a_est,
           b_est=b_est,
           se_alpha=se_alpha,
           se_beta=se_beta)
}


