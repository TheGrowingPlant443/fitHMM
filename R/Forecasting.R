#' @import expm
 

#' @title Forecasting Function
#' 
#' @description Calculate the probability of seeing a value from each state up to time 
#' T + Steps, where T is the length of the observed sequence.  
#' 
#' 
#' @usage Forecasting(HMM, observations, steps = 1)
#' 
#' @param HMM Your specified hidden Markov model by using the function \code{HMM}
#' @param observations Your Data 
#' @param steps Number of steps 
#'
#' @return Vector of probability of being in each state
#' 
Forecasting <- function(HMM, observations, steps=1) {
  AlphaBetaLikelihood <- ForwardBackwardAlgorithm(HMM, observations)
  Alpha <- AlphaBetaLikelihood$Alpha
  Beta <- AlphaBetaLikelihood$Beta
  Likelihood <- AlphaBetaLikelihood$Likelihood_forw
  P_mat <- AlphaBetaLikelihood$P_mat
  Forecasting_predict = colSums(Alpha[,dim(Alpha)[2]] %*% (expm::`%^%`(HMM$TPM,steps)) %*% P_mat)/Likelihood
  return(Forecasting_predict)
}


