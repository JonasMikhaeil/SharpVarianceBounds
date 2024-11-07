#' Estimation of Sharp Variance Bounds for Randomized Controlled 
#' Experiments. 
#' 
#' @description
#' This function calculates a consistent estimate of the sharp variance bound
#' for general regression-adjusted estimators of the average treatment
#' effect in two-arm randomized controlled experiments. 
#' This variance estimate can be used to construct the asymptotically narrowest conservative Wald-type 
#' confidence interval with the nominal coverage for general regression adjusted estimators for the 
#' treatment effect. For more details, see [(Mikhaeil and Green, 2024)](https://arxiv.org/abs/2411.00191)
#' 
#' @param r_obs The residuals of the outcome model \eqn{Y_i(q) - \hat f_q(X_i)}.
#' @param Z Vector of the treatment assignment.
#' @param upper Logical argument: Determines if upper (TRUE) or lower (FALSE) bounds will be calculated.
#' @return Estimate of upper or lower sharp variance bound
#' @export
sharpvar <- function(r_obs, Z, upper=TRUE){
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
  
  
  1/n*((n-m)/m *var_rt + m/(n-m)*var_rc + 
         2*sum(((p_i-p_i_minus)*r_ti*r_ci)[2:length(p_i)]) 
       - 2*mean(r_t)*mean(r_c))
}
#' Calculates (Neyman's) conventional variance estimator. 
#' 
#' @param r_obs The residuals of the outcome model \eqn{Y_i(q) - \hat f_q(X_i)}.
#' @param Z Vector of the treatment assignment.
#' @return Conventional variance estimator
#' @export
neyman_conventional <- function(r_obs,Z){
  r_t <- r_obs[Z>0] 
  r_c <- r_obs[Z==0]
  
  m <- sum(Z)
  n <- length(r_obs)
  
  helper_var <- function(x){
    1/(length(x)-1) * sum((x-mean(x))^2)
  }
  var_rt <- helper_var(r_t)
  var_rc <- helper_var(r_c)
  (1/m *var_rt + 1/(n-m)*var_rc)
}
#' Calculates Neyman's variance estimator based on the Cauchy-Schwarz Inequality. 
#' 
#' @param r_obs The residuals of the outcome model \eqn{Y_i(q) - \hat f_q(X_i)}.
#' @param Z Vector of the treatment assignment.
#' @return Neyman's Cauchy-Schwarz variance estimator
#' @export
neyman_cs <- function(r_obs,Z){
  r_t <- r_obs[Z>0] 
  r_c <- r_obs[Z==0]
  
  m <- sum(Z)
  n <- length(r_obs)
  
  helper_var <- function(x){
    1/(length(x)-1) * sum((x-mean(x))^2)
  }
  var_rt <- helper_var(r_t)
  var_rc <- helper_var(r_c)
  1/n*((n-m)/m *var_rt + m/(n-m)*var_rc + 2*sqrt(var_rt*var_rc))
}

