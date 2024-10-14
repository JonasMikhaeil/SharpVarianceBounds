
#' Creates a simulated dataset to showcase sharp variance bounds
#' 
#' This function creates the simulated dataset described in (Mikhaeil, Green, 2024) 
#' to showcase the use of and benefits of sharp variance bounds.
#'
#' @param N The number of simulated units.
#' @param sd Simulation parameter
#' @return Simulated dataset containing potential outcomes (Y0,Y1) and covariates X.
#' @export
create_dataset <- function(N = 1000,
                           sd = 0,
                           eps_var = 1) {
  alpha <- 0
  beta <- 1

  
  x <- runif(N,-1, 1)
  eps <- rnorm(N, 0, eps_var)
  y0 <- alpha + beta * x + eps
  y1 <- ifelse(rbinom(N, 1, 1 - sd), y0, 10 + 0.5 * eps)
  
  data.frame(Y0 = y0, Y1 = y1 ,  X = x)
}

#' Randomizes treatment 
#' 
#' This function simulates a CRE by randomizing treatment over a population of potential outcomes. 
#'
#' @param Y0 Potential outcomes under no treatment
#' @param Y1 Potential outcomes under treatment
#' @param x covariates
#' @param x covariates
#' @param N_r ratio of people receiving treatment
#' @return Randomized dataset containing observed outcomes Y_obs, treatment and covariates X
#' @export
randomize <- function(Y0,Y1,x,N_r=0.5){
  N <- length(Y0)
  index <- 1:N
  ind_t <- sample(index,ceiling(N_r*N),replace = FALSE)
  ind_c <- setdiff(index,ind_t)
  Y_t <- Y1[ind_t]
  x_t <- x[ind_t]
  Y_c <- Y0[ind_c]
  x_c <- x[ind_c]
  Y_obs <- c(Y_t,Y_c)
  X<-c(x_t,x_c)
  treatment <- c(rep(1,length(ind_t)),rep(0,length(ind_c)))
  data.frame(Y_obs = Y_obs, treatment = treatment , X= X)
}
#' Calculates the true variance
#' 
#' This function simulates a CRE by randomizing treatment over a population of potential outcomes. 
#'
#' @param eps0 Potential outcomes under no treatment
#' @param eps1 Potential outcomes under treatment
#' @return Randomized dataset containing observed outcomes Y_obs, treatment and covariates X
#' @export
var_true <- function(eps0,eps1){
  n <- length(eps0)
  m <- n/2
  helper_var <- function(x){
    1/(length(x)-1) * sum((x-mean(x))^2)
  }
  helper_cov <- function(x,y){
    1/(length(x)-1) * sum((x-mean(x))*(y-mean(y)))
  }
 
  var_eps1 <- helper_var(eps1)
  var_eps0 <- helper_var(eps0)
  cov_eps1eps0 <- helper_cov(eps1,eps0)
  
  1/n*((n-m)/m *var_eps1 + m/(n-m)*var_eps0 + 
         2*cov_eps1eps0)
}