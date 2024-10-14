#' Calculates an estimate of Neyman's bound of the cross-term of the variance of average treatment effect estimator under general regression adjsustment in randomized two-arm 
#' experiments based on the Cauchy-Schwarz Inequality.
#' 
#' This function calculates an estimate of Neyman's bound of the cross-term of the variance of average treatment effect estimator under general regression adjsustment in randomized two-arm 
#' experiments based on the Cauchy-Schwarz Inequality.
#' 
#'
#' @param r_obs The residuals of the outcome model,i.e, $Y_i(q) - \hat f_q(X_i)$.
#' @param Z Vector of the treatment assignment.
#' @return estimate of Neyman's Cauchy-Schwarz bound on the cross-term
#' @export
cor_approx_cs <- function(r_obs,Z){
  r_t <- r_obs[Z>0] 
  r_c <- r_obs[Z==0]
  
  m <- sum(Z)
  n <- length(r_obs)
  
  helper_var <- function(x){
    1/(length(x)-1) * sum((x-mean(x))^2)
  }
  var_rt <- helper_var(r_t)
  var_rc <- helper_var(r_c)
  sqrt(var_rt*var_rc)
}
#' Calculates an estimate of the conventional bound of the cross-term of the variance of average treatment effect estimator under general regression adjsustment in randomized two-arm 
#' experiments based on the Cauchy-Schwarz Inequality and the AM-GM inequality.
#' 
#' This function calculates an estimate of Neyman's bound of the cross-term of the variance of average treatment effect estimator under general regression adjsustment in randomized two-arm 
#' experiments based on the Cauchy-Schwarz Inequality.
#' 
#'
#' @param r_obs The residuals of the outcome model,i.e, $Y_i(q) - \hat f_q(X_i)$.
#' @param Z Vector of the treatment assignment.
#' @return estimate of convential bound on the cross-term
#' @export
cor_approx_conventional <- function(r_obs,Z){
  r_t <- r_obs[Z>0] 
  r_c <- r_obs[Z==0]
  
  m <- sum(Z)
  n <- length(r_obs)
  
  helper_var <- function(x){
    1/(length(x)-1) * sum((x-mean(x))^2)
  }
  var_rt <- helper_var(r_t)
  var_rc <- helper_var(r_c)
  1/2*(var_rt+var_rc)
}

#' Calculates an estimates of the cross-term of the variance of average treatment effect estimator under general regression adjsustment in randomized two-arm 
#' experiments under extremal (comonotonic) dependence of the potential outcomes. 
#' 
#' This function calculates a consistent estimate of the cross-term of the variance of average treatment effect estimator under general regression adjsustment in randomized two-arm 
#' experiments under extremal dependence of the potential outcomes
#' 
#'
#' @param r_obs The residuals of the outcome model,i.e, $Y_i(q) - \hat f_q(X_i)$.
#' @param Z Vector of the treatment assignment.
#' @return estimate of the cross-term of the variance average treatment effect estimator under general regression adjsustment in randomized two-arm 
#' experiments under extremal (comonotonic) dependence of the potential outcomes. 
#' @export
cor_approx_sharp <- function(r_obs,Z,upper=TRUE){
  helper_var <- function(x){
    1/(length(x)-1) * sum((x-mean(x))^2)
  }
  
  r_t <- r_obs[Z>0] 
  r_c <- r_obs[Z==0]
  m <- sum(Z)
  
  n <- length(r_obs)
  r_t <- sort(r_t)
  if(upper==TRUE) r_c <- sort(r_c) 
  else r_c <- sort(r_c,decreasing = TRUE)
  p_i <- unique(sort(c(seq(0,n-m,1)/(n-m),seq(0,m,1)/m))) -.Machine$double.eps^.5
  p_i[1] <- .Machine$double.eps^.5
  r_ti <- r_t[ceiling(p_i*m)]
  r_ci <- r_c[ceiling(p_i*(n-m))]
  p_i_minus <- c(NA,p_i[1:(length(p_i)-1)])
  var_rt <- helper_var(r_t)
  var_rc <- helper_var(r_c)
  

  sum(((p_i-p_i_minus)*r_ti*r_ci)[2:length(p_i)]) - mean(r_t)*mean(r_c)
}


